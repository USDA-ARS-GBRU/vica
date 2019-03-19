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
import ete3
from typing import Tuple
import vica

from sklearn.metrics import confusion_matrix
from sklearn.metrics import accuracy_score
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import f1_score

with open(vica.CONFIG_PATH) as cf:
    config = yaml.safe_load(cf)


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
        labeldict (dict)

    """
    ncbi =  ete3.NCBITaxa()
    class2labels = {}
    for i, val in enumerate(classdict.keys()):
        taxa = ncbi.get_taxid_translator([val])
        class2labels[val] = {"name": taxa[val], "class": i}
    logging.info("Class labels are %s", class2labels)

    idlist = []
    for filename in tf_rec_filenames:
        idlist.append(_ids_from_tfrecords(filename))
    labeldict = {}
    for item in idlist:
        taxid = item.split("|")[1]
        if not taxid in labeldict:
            lineage = ncbi.get_lineage(taxid)
            classtaxid = list(set(class2labels.keys()).intersection(lineage))
            assert len(classtaxid)==1
            labeldict[taxid] = labeldict[classtaxid[0]]["class"]
    return labeldict




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
        dnn_hidden_units=[64, 8],
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
        feature_columns=[codon],
        dropout=0.5,
        activation_fn=tf.nn.relu,
        hidden_units=[64, 8],
        optimizer='Adagrad')
    return dnn_estimator


def base_input_fn(codonlength: int, minhashlength: int, kmerdim: int,
                  labeldict: dict, shuffle: bool, shuffle_buffer_size: int,
                  batch: int, epochs: int, filenames: list):
    """The function for feeding and processing training data

    Tensorflow estimators take a function that processes TFrecord Datasets
    into an example iterator that returns one processed example each time
    it is needed. because an estimator takes a function without arguments
    functools.partial is used to set the parameter values for each application
    (train, evaluate, classify).

    Args:
        codonlength (int): the number of elements in the codon feature set
        minhashlength (int)  the number of elements in the minhash feature set
        kmerdim (int):  the number of elements in the kmer feature set
        labeldict (dict): a list containing the ids of the sequence segments and integer class values
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
    # Batch the dataset
    dataset = dataset.batch(batch)
    # define how to parse the tfrecords
    def parser(record):
        keys_to_features = {"id": tf.FixedLenFeature((), tf.string),
            "kmer": tf.FixedLenFeature([kmerdim], tf.float32),
            "codon": tf.FixedLenFeature([codonlength], tf.float32),
            "minhash": tf.FixedLenFeature([minhashlength], tf.float32),
            "hmmer":tf.FixedLenFeature((),dtype=tf.string)}
        parsed = tf.parse_example(serialized=record, features=keys_to_features)
        recordlabel = labeldict[parsed["id"]]
        return {'id': parsed['id'], 'kmer': parsed['kmer'], 'codon': parsed['codon'], 'minhash': parsed['minhash']}, recordlabel
    dataset = dataset.map(parser)
    if shuffle:
        dataset = dataset.shuffle(shuffle_buffer_size)

    dataset = dataset.repeat(epochs)
    iterator = dataset.make_one_shot_iterator()
    features, labels = iterator.get_next()
    return features
