"""A module to train and evaluate Vica models from tfrecords
   of features. It uses The Tensorflow 1.3+ datasets api and estimator api

"""


import os
import tempfile
import functools
import logging
import shutil
import random

import yaml
import tensorflow as tf
import csv
import ete3
from typing import Tuple
import vica


with open(vica.CONFIG_PATH) as cf:
    config = yaml.safe_load(cf)

    ## Functions for creating Labels from TF record data automatically
    def _ids_from_tfrecords(filename: str) -> list:
        """Takes a tfrecord filename and returns a list with the ids in order.

        """
        idlist = []
        for example in tf.python_io.tf_record_iterator(filename):
            result = tf.train.Example.FromString(example)
            idstr = result.features.feature["id"].bytes_list.value[0].decode()
            idlist.append(idstr)
        return idlist


    def _create_class_lookup(classdict: dict, tf_rec_filenames: list) -> Tuple[dict, dict]:
        """Read  the ids of the records and create a dictionary for fast lookup of labels

        Args:
            clasddict (dict): the dictionary of training classes and their targeted samples
            tf_rec_filenames (list): a list of the tfracords to be used

        Returns:
            labeldict (dict): Dict containing record IDs and the class they belong to

        """
        ncbi = ete3.NCBITaxa()
        class2labels = {}
        for i, val in enumerate(classdict.keys()):
            taxa = ncbi.get_taxid_translator([val])
            class2labels[val] = {"name": taxa[val], "class": i}
        logging.info("Class labels are %s", class2labels)

        idlist = []
        for filename in tf_rec_filenames:
            idlist.extend(_ids_from_tfrecords(filename))
        labeldict = {}
        for item in idlist:
            taxid = item.split("|")[1]
            if not taxid in labeldict:
                try:
                    lineage = ncbi.get_lineage(taxid)
                    classtaxid = list(set(class2labels.keys()).intersection(lineage))
                    assert len(classtaxid) == 1
                    labeldict[taxid] = class2labels[classtaxid[0]]["class"]
                except:
                    labeldict[taxid] = random.randint(0,len(classdict)-1)
                    logging.info("Could not determine the class for taxid %s, randomly assigned class %s" % (str(taxid), labeldict[taxid]))
        return labeldict




    # The input functions for the estimators

    def _base_input_fn(labeldict: dict, shuffle_buffer_size: int,
                      batch: int, epochs: int, filenames: list):
        """The function for feeding and processing training data

        Tensorflow estimators take a function that processes TFrecord Datasets
        into an example iterator that returns one processed example each time
        it is needed. because an estimator takes a function without arguments
        functools.partial is used to set the parameter values for each application
        (train, evaluate, classify).

        Args:

            shuffle_buffer_size (int): How large of a record buffer to load for shuffling, 0 for no shuffling
            batch (int): the size of the training batch (~32 for training larger,
                 ~512 for classification )
            epochs (int): Epochs (number of complete passes through the data)
            filenames(list): List of input tfrecords filenames

        Returns:
            a input function for a Tensorflow estimator

        """

        dataset = tf.data.TFRecordDataset(filenames)

        # define how to parse the tfrecords
        def parser(record):
            keys_to_features = {"id": tf.FixedLenFeature((), dtype=tf.string),
                                "kmer": tf.FixedLenFeature([config["train_eval"]["kmerlength"]], dtype=tf.float32),
                                "codon": tf.FixedLenFeature([config["train_eval"]["codonlength"]], dtype=tf.float32),
                                "minhash": tf.VarLenFeature(dtype=tf.string),
                                "hmmer": tf.VarLenFeature(dtype=tf.string),
                                "label": tf.FixedLenFeature((), dtype=tf.int64)}
            parsed = tf.parse_example(serialized=record, features=keys_to_features)
            return {'id': parsed['id'], 'kmer': parsed['kmer'],
                    'codon': parsed['codon'], 'minhash': parsed['minhash'], 'hmmer': parsed['hmmer']}, parsed['label']
        dataset = dataset.batch(batch)
        dataset = dataset.map(parser)
        if shuffle_buffer_size > 0:
            dataset = dataset.shuffle(shuffle_buffer_size)
        # Batch the dataset
        dataset = dataset.repeat(epochs)
        iterator = dataset.make_one_shot_iterator()
        features, labels = iterator.get_next()
        return features, labels

    # modified input functions for special purposes



    # classift_input_fn = functools.partial(_base_input_fn,
    #     shuffle_buffer_size=0,
    #     batch=config["train_eval"]["eval_batch_size"],
    #     epochs=1,
    #     filenames=classify_files)

    # Define feature_columns

    kmer_feat = tf.feature_column.numeric_column(key='kmer', shape=(config["train_eval"]["kmerlength"]))
    codon_feat = tf.feature_column.numeric_column(key='codon', shape=(config["train_eval"]["codonlength"]))
    # minhashvocab = list(config["split_shred"]["classes"].keys()).append("nohits")
    minhashvocab = ["2","2157","2759","10239", "nohits"]
    minhash_feat = tf.feature_column.categorical_column_with_vocabulary_list(key='minhash', vocabulary_list=minhashvocab)
    hashed_hmm_feat = tf.feature_column.categorical_column_with_hash_bucket(
        key="hmmer", hash_bucket_size=1000)
    embedded_hmm_feat = tf.feature_column.embedding_column(
        categorical_column=hashed_hmm_feat, dimension=6)
    dense_features = [embedded_hmm_feat, codon_feat, kmer_feat]
    all_features = [hashed_hmm_feat, kmer_feat, codon_feat, minhash_feat]

EVAL_INTERVAL = 300
run_config = tf.estimator.RunConfig(save_checkpoints_secs = EVAL_INTERVAL,
                                        keep_checkpoint_max = 3)

    # Model definitions
def create_estimator(modeldir, n_classes):
    dnnlogistic_estimator = tf.estimator.DNNLinearCombinedClassifier(
        model_dir=modeldir,
        n_classes=n_classes,
        weight_column=None,
        linear_feature_columns=[minhash_feat],
        linear_optimizer=tf.train.FtrlOptimizer(
                                         learning_rate=0.1,
                                         l1_regularization_strength=0.001),
        dnn_feature_columns=dense_features,
        dnn_dropout=0.5,
        dnn_activation_fn=tf.nn.relu,
        dnn_hidden_units=[256, 32],
        dnn_optimizer='Adam')
    return dnnlogistic_estimator

def create_DNN_estimator(modeldir, n_classes):
    dnnlogistic_estimator = tf.estimator.DNNClassifier(
        model_dir=modeldir,
        n_classes=n_classes,
        feature_columns=dense_features,
        dropout=0.5,
        activation_fn=tf.nn.relu,
        hidden_units=[256, 32],
        optimizer='Adam',
        config=run_config)
    return dnnlogistic_estimator

def create_log_estimator(modeldir, n_classes):
    logistic_estimator = tf.estimator.LinearClassifier(
    model_dir=modeldir,
    n_classes=n_classes,
    feature_columns=[codon_feat, kmer_feat, embedded_hmm_feat],
    optimizer=tf.train.FtrlOptimizer(
                                     learning_rate=0.1,
                                     l1_regularization_strength=0.001))
    return logistic_estimator


def train_and_eval(train_files, eval_files, modeldir, configpath=vica.CONFIG_PATH):

    labeldict = _create_class_lookup(classdict=config["split_shred"]["classes"],
                                     tf_rec_filenames=train_files + eval_files)
    n_classes = len(config["split_shred"]["classes"])

    train_input_fn = functools.partial(_base_input_fn,
        labeldict=labeldict,
        shuffle_buffer_size=config["train_eval"]["train_shuffle_buffer"],
        batch=config["train_eval"]["train_batch_size"],
        epochs=config["train_eval"]["epochs"],
        filenames=train_files)

    eval_input_fn = functools.partial(_base_input_fn,
        labeldict=labeldict,
        shuffle_buffer_size=0,
        batch= config["train_eval"]["eval_batch_size"],
        epochs=1,
        filenames=eval_files)
    #def my_auc(labels, predictions):
    #    return {'auc': tf.metrics.auc(labels, predictions['probabilities'], curve="PR")}
    my_estimator = create_DNN_estimator(modeldir=modeldir, n_classes=n_classes)
    #my_estimator = tf.estimator.add_metrics(my_estimator, my_auc)
    train_spec = tf.estimator.TrainSpec(input_fn=train_input_fn)
    eval_spec = tf.estimator.EvalSpec(input_fn=eval_input_fn,
                                      start_delay_secs = 60,
                                      throttle_secs = EVAL_INTERVAL)
    t1, t2 = tf.estimator.train_and_evaluate(my_estimator, train_spec, eval_spec)
    return t1, t2
