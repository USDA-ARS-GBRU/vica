"""A module to train and evaluate Vica models from tfrecords
   of features. It uses The Tensorflow 1.3+ datasets api and estimator api

"""


import os
import tempfile
import functools
import logging
import shutil

import yaml
import tensorflow as tf
import csv

import vica

with open(vica.CONFIG_PATH) as cf:
    config = yaml.safe_load(cf)

def _featureshape(k=5, codonlength=177, minhashlength=267):
    """Determine the shape of the features for each feature type including
    for kmers of different lengths.

    Args:
        k (int): kmer size 4-8
        codonlegnth (int): length of codon feature vectors
        minhashlength (int): length of minhash feature vector

    """
    kmerdim = len(vica.khmer_features.iterate_kmer(k)) - 1
    kmer = tf.feature_column.numeric_column(key='kmer', shape=(kmerdim))
    codon = tf.feature_column.numeric_column(key='codon', shape=(codonlength))
    minhash = tf.feature_column.numeric_column(key='minhash', shape=(minhashlength))
    return kmerdim, kmer, codon, minhash

def _ids_from_tfrecords(filename):
    """Takes a tfrecord filename and returns a list with the labels in order.

    """
    idlist = []
    for example in tf.python_io.tf_record_iterator(filename):
        result = tf.train.Example.FromString(example)
        idstr = result.features.feature["id"].bytes_list.value[0].decode()
        idlist.append(idstr)
    return idlist


def _ids_from_tfrecords_files(file_names):
    """ input: list of tfrecord file names
        output: list with the labels
    """
    id_list = []
    for filename in file_names:
        id_list = id_list + _ids_from_tfrecords(filename)

    return id_list



def base_input_fn(codonlength, minhashlength, kmerdim, shuffle, shuffle_buffer_size, batch, epochs, filenames):
    """The function for feeding and processing training data

    Tensorflow estimators take a function that processes TFrecord Datasets
    into an example iterator that returnes one processed example each time
    it is needed. because an estimator takes a function without arguments
    functools.partial is used to set the parameter values for each application
    (train, evaluate, classify).

    Args:
        codonlength (int): the number of elements in the codon feature set
        minhashlength (int)  the number of elements in the minhash feature set
        kmerdim (int):  the number of elements in the kmer feature set
        shuffle (bool):  Should values be shuffled? (True for train, false
            for evaluate and classify)
        shuffle_buffer_size (int): How large of a record buffer to load for shuffling
        batch (int): the size of the training batch (~32 for training larger,
             ~512 for classification )
        epochs (int): Epochs (number of complete passes through the data)
        filenames(list): List of input tfrecords filenames

    Returns:
        a input function for a Tensorflow estimator

    """
    dataset = tf.data.TFRecordDataset(filenames)
    def parser(record):
        keys_to_features = {"id": tf.FixedLenFeature((), tf.string),
            "label": tf.FixedLenFeature((), tf.int64),
            "kmer": tf.FixedLenFeature([kmerdim], tf.float32),
            "codon": tf.FixedLenFeature([codonlength], tf.float32),
            "minhash": tf.FixedLenFeature([minhashlength], tf.float32)}
        parsed = tf.parse_single_example(record, keys_to_features)
        return {'kmer': parsed['kmer'], 'codon': parsed['codon'], 'minhash': parsed['minhash']}, parsed['label']
    dataset = dataset.map(parser)
    if shuffle:
        dataset = dataset.shuffle(shuffle_buffer_size)
    dataset = dataset.batch(batch)
    dataset = dataset.repeat(epochs)
    iterator = dataset.make_one_shot_iterator()
    features, labels = iterator.get_next()
    return features, labels

# Model definitions
def mk_dnnlogistic_estimator(modeldir, n_classes, minhash, kmer, codon):
    """Specification of Wide and Deep Neural Network model

    Args:
        modeldir (str): path to directory containing model data
        n_classes (int): number of classes in the model
        minhash (obj): Tensorflow feature column object for minhash data
        kmer (obj): Tensorflow feature column object for kmer data
        codon (obj): Tensorflow feature column object for codon data

    Returns:
        (obj): Tensorflow.estimator.DNNLinearCombinedClassifier object

    """
    dnnlogistic_estimator = tf.estimator.DNNLinearCombinedClassifier(
        model_dir = modeldir,
        n_classes=n_classes,
        weight_column=None,
        linear_feature_columns=[minhash],
        linear_optimizer='Ftrl',
        dnn_feature_columns=[kmer, codon],
        dnn_dropout=0.5,
        dnn_activation_fn=tf.nn.relu,
        dnn_hidden_units=[256, 32],
        dnn_optimizer='Adagrad')
    return dnnlogistic_estimator

def mk_dnn_estimator(modeldir, n_classes, kmer, codon):
    """Specification of Deep Neural Network model

    Args:
        modeldir (str): path to directory containing model data
        n_classes (int): number of classes in the model
        kmer (obj): Tensorflow feature column object for kmer data
        codon (obj): Tensorflow feature column object for codon data

    Returns:
        (obj): Tensorflow.estimator.DNNCLassifier object

    """
    dnn_estimator = tf.estimator.DNNClassifier(
        model_dir = modeldir,
        n_classes=n_classes,
        weight_column=None,
        feature_columns=[kmer, codon],
        dropout=0.5,
        activation_fn=tf.nn.relu,
        hidden_units=[256, 32],
        optimizer='Adagrad')
    return dnn_estimator

def train(infiles, out, modeldir, n_classes, configpath):
    """Train a Tensorflow model with training data returning a modeldir.

    Trains a tensorflow model for the number of epoch specified in config file.
    calling the function on the same modeldir resumes training at the last
    checkpoint. The modeldir can be used.

    Args:
        infiles (list): A list with paths to TFrecords file(s) containing labeled test sequences.
        out (str): file where TF Saved model file will be written (Future)
        modeldir (str): a model directory containing a trained model to use
        for evaluation.
        n_classes (int): the number of classes in the model (default 4)
        configpath (str): path to yaml config files

    Returns:
        None

    Todo:
        Convert the model type to Tensorflow SavedModel format and TF serving
        API once the Tensorflow API supports using SavedModels with Python 3.

    Note: the parameter `out` has no effect currently but been reserved
        for use once the SavedModel API is better supported by Tensorflow.

    """
    try:
        logging.info("Beginning tensorflow model training. To see results in real-time run 'tensorboard --logdir={}'".format(modeldir))
        with open(configpath, "r") as cf:
            config = yaml.safe_load(cf)
        kmerdim, kmer, codon, minhash = _featureshape(config["khmer_features"]["ksize"])
        input_fn = functools.partial(base_input_fn,
            codonlength=config["train_eval"]["codonlength"],
            minhashlength=config["train_eval"]["minhashlength"],
            kmerdim=kmerdim,
            shuffle=True,
            shuffle_buffer_size=10000,
            batch=config["train_eval"]["train_batch_size"],
            epochs=config["train_eval"]["epochs"],
            filenames=infiles)
        if config["train_eval"]["model"] == "DNN":
            estimator = mk_dnn_estimator(modeldir=modeldir,
                n_classes=int(n_classes),
                kmer=kmer,
                codon=codon)
            estimator.train(input_fn=input_fn)
        elif config["train_eval"]["model"] == "DNNLogistic":
            estimator = mk_dnnlogistic_estimator(modeldir=modeldir,
                n_classes=int(n_classes),
                minhash=minhash,
                kmer=kmer,
                codon=codon)
            estimator.train(input_fn=input_fn)
    except:
        logging.exception("During tensorflow model training the following exception occurred:")
        raise SystemExit(1)
    try:
        # Save results if successful
        feature_spec={"id": tf.FixedLenFeature((), tf.string),
            "kmer": tf.FixedLenFeature([kmerdim], tf.float32),
            "codon": tf.FixedLenFeature([config["train_eval"]["codonlength"]], tf.float32),
            "minhash": tf.FixedLenFeature([config["train_eval"]["minhashlength"]], tf.float32)}
        serving_input_receiver_fn = tf.estimator.export.build_parsing_serving_input_receiver_fn(feature_spec)
        estimator.export_savedmodel(out, serving_input_receiver_fn)
    except:
        logging.exception("While exporting the model after tensorflow training the following exception occured:")


def evaluate(infiles, out, modeldir, n_classes, configpath):
    """Evaluate a trained model using test data, returning a file of
        predictions and training metrics.


    Classify fasta files or TFrecords containing contigs to identify viral
    contigs. The classifier takes contigs and extracts features then makes
    predictions. It returns a text file of prediction probabilities for each
    model class. This file can be used to filter out viral contigs for analysis.

    Args:
        infiles (list): A list with paths to TFrecords file(s) containing labeled test sequences.
        out (str): a directory where model predictions will be writen
        modeldir (str): a model directory containing a trained model to use
        for evaluation.
        n_classes (int): the number of classes in the model (default 4)
        configpath (str): path to yaml config files

    Returns:
        None

    Todo:
        Convert the model type to Tensorflow SavedModel format and TF serving
        API once the Tensorflow API supports using SavedModels with Python 3.

        Add additional evaluation files: Precision recall plots, confusion matrix.

    """
    try:
        logging.info("Beginning tensorflow model evaluation. To see results in real-time run 'tensorboard --logdir={}'".format(modeldir))
        with open(configpath, "r") as cf:
            config = yaml.safe_load(cf)
        kmerdim, kmer, codon, minhash = _featureshape(config["khmer_features"]["ksize"])
        input_fn = functools.partial(base_input_fn,
            codonlength=config["train_eval"]["codonlength"],
            minhashlength=config["train_eval"]["minhashlength"],
            kmerdim=kmerdim,
            shuffle=False,
            shuffle_buffer_size=0,
            batch=config["train_eval"]["eval_batch_size"],
            epochs=1,
            filenames=infiles)
        if config["train_eval"]["model"] == "DNN":
            estimator = mk_dnn_estimator(modeldir=modeldir,
                n_classes=int(n_classes),
                kmer=kmer,
                codon=codon)
            results = estimator.evaluate(input_fn=input_fn)
            preds = estimator.predict(input_fn=input_fn)
        elif config["train_eval"]["model"] == "DNNLogistic":
            estimator = mk_dnnlogistic_estimator(modeldir=modeldir,
                n_classes=int(n_classes),
                minhash=minhash,
                kmer=kmer,
                codon=codon)
            results = estimator.evaluate(input_fn=input_fn)
            preds = estimator.predict(input_fn=input_fn)
        logging.info("Tensorflow model performance. See also {}.".format(out))
        logging.info(preds)
        if not os.path.exists(out):
            os.mkdir(out)
        predictions = os.path.join(out, "modelpredictions.txt")
        idlist = _ids_from_tfrecords_files(infiles)
        with open(predictions, "w") as outfile:
            csv_writer_instance = csv.writer(outfile, lineterminator='\n')
            header = ["ID", "Class", "Class_id"] + ["Prob_class_" + str(i) for i in range(int(n_classes))]
            csv_writer_instance.writerow(header)
            for idrec, rec in zip(idlist,preds):
                plist = rec['probabilities']
                pliststr = [str(x) for x in plist]
                ll = [idrec, rec['classes'][0].decode("utf-8"), str(rec['class_ids'][0])]
                ll.extend(pliststr)
                csv_writer_instance.writerow(ll)
        for key in sorted(results):
            logging.info('{}: {}'.format(key, results[key]))
    except:
        logging.exception("During tensorflow model evaluation the following exception occured:")


def classify(infile, out, modeldir, n_classes, configpath):
    """Classify fasta sequences or TRrecords data to identify viral sequences

    Classify fasta files or TFrecords containing contigs to identify viral
    contigs. The classifier takes contigs and extracts features then makes
    predictions. It returns a text file of prediction probabilities for each
    model class. This file can be used to filter out viral contigs for analysis.

    Args:
        infile (str): the path to a fasta file of contigs (2 kb or longer is
            recomended) or TFrecords to be classified.
        out (str): a file with contig modelpredictions
        modeldir (str): a model directory containing a trained model to use
        for classification.
        n_classes (int): the number of classes in the model (default 4)
        configpath (str): path to yaml config files

    Returns:
        None

    Todo:
        Convert the model type to Tensorflow SavedModel format and TF serving
        API once the Tensorflow API supports using SavedModels with Python 3.


    """
    logging.info("Beginning Tensorflow classification")
    with open(configpath, "r") as cf:
        config = yaml.safe_load(cf)
    try:
        if infile.endswith(("tfrecord", "TFrecord")):
            logging.info("Classifing data from TFRecord file.")
            getfeat = False
        elif infile.endswith(("fasta","fa","fna")):
            logging.info("Classifing data from fasta file.")
            getfeat =True
    except:
        logging.exception("Files with that suffix are not supported. Please \
            use .fasta, .fna, or .fa files or .tfrecord files created by \
            `vica get_features`")
        raise SystemExit(1)
    if getfeat:
        dtemp = tempfile.mkdtemp()
        tfrecfile = os.path.join(dtemp, "data.tfrecord")
        logging.info("Extracting Features from the sequence data. For more \
            control of options use `vica get_featues`")
        vica.get_features.run(infile=infile,
            output=tfrecfile,
            label=0,
            minhashlocal=None,
            configpath=configpath)
    else:
        tfrecfile=infile
    kmerdim, kmer, codon, minhash = _featureshape(config["khmer_features"]["ksize"])
    input_fn = functools.partial(base_input_fn,
        codonlength=config["train_eval"]["codonlength"],
        minhashlength=config["train_eval"]["minhashlength"],
        kmerdim=kmerdim,
        shuffle=False,
        shuffle_buffer_size=0,
        batch=config["train_eval"]["eval_batch_size"],
        epochs=1,
        filenames=tfrecfile)
    if config["train_eval"]["model"] == "DNN":
        estimator = mk_dnn_estimator(modeldir=modeldir,
            n_classes=int(n_classes),
            kmer=kmer,
            codon=codon)
        preds = estimator.predict(input_fn=input_fn)
        print(preds)
    elif config["train_eval"]["model"] == "DNNLogistic":
        estimator = mk_dnnlogistic_estimator(modeldir=modeldir,
            n_classes=int(n_classes),
            minhash=minhash,
            kmer=kmer,
            codon=codon)
        preds = estimator.predict(input_fn=input_fn)
    idlist = _ids_from_tfrecords(tfrecfile)
    with open(out, "w") as outfile:
        csv_writer_instance = csv.writer(outfile, lineterminator='\n')
        header = ["ID", "Class", "Class_id"] + ["Prob_class_" + str(i) for i in range(int(n_classes))]
        csv_writer_instance.writerow(header)
        for recid, rec in zip(idlist, preds):

            plist = rec['probabilities']
            pliststr = [str(x) for x in plist]
            ll = [recid, rec['classes'][0].decode("utf-8"), str(rec['class_ids'][0])]
            ll.extend(pliststr)
            csv_writer_instance.writerow(ll)
    if getfeat:
        shutil.rmtree(dtemp)
