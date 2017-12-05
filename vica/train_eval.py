#!/usr/bin/env python3
"""train_eval.py: a module to train models and evaluate models from tfrecords
   of features. It uses The Tensorflow 1.3+ datasets api and estimator api"""


import os
import urllib
import tempfile
import time
import datetime
from collections import Counter

import numpy as np
import tensorflow as tf

import vica

with open(configpath) as cf:
    config = yaml.load(cf)

def _featureshape(k):
    """Determine the shape of the features for each feature type including
    for kmers of different lengths."""
    codonlength = 177
    minhashlength= 267
    kmerdim = len(vica.khmer_features.iterate_kmer(k)) - 1
    kmer = tf.feature_column.numeric_column(key='kmer', shape=(kmerdim))
    codon = tf.feature_column.numeric_column(key='codon', shape=(config["train_eval"]["codonlength"]))
    minhash = tf.feature_column.numeric_column(key='minhash', shape=(config["train_eval"]["minhashlength"]))
    return kmerdim, kmer, codon, minhash


def train_input_fn():
    """the function for feeding and processing training data"""
    filenames = tf.placeholder(tf.string, shape=[None])
    dataset = tf.contrib.data.TFRecordDataset(filenames)
    def parser(record):
        keys_to_features = {"id": tf.FixedLenFeature((), tf.string),
            "label": tf.FixedLenFeature((), tf.int64),
            "kmer": tf.FixedLenFeature([kmerdim], tf.float32),
            "codon": tf.FixedLenFeature([config["train_eval"]["codonlength"]], tf.float32),
            "minhash": tf.FixedLenFeature([config["train_eval"]["minhashlength"]], tf.float32)}
        parsed = tf.parse_single_example(record, keys_to_features)
        return {'kmer': parsed['kmer'], 'codon': parsed['codon'], 'minhash': parsed['minhash']}, parsed['label']
    dataset = dataset.map(parser)
    dataset = dataset.shuffle(buffer_size=10000)
    dataset = dataset.batch(32)
    dataset = dataset.repeat(epochs)
    iterator = dataset.make_one_shot_iterator()
    features, labels = iterator.get_next()
    return features, labels

def test_input_fn():
    """    """the function for feeding and processing training data""""""
    filenames = tf.placeholder(tf.string, shape=[None])
    dataset = tf.contrib.data.TFRecordDataset(filenames)
    def parser(record):
        keys_to_features = {"id": tf.FixedLenFeature((), tf.string),
            "label": tf.FixedLenFeature((), tf.int64),
            "kmer": tf.FixedLenFeature([kmerdim], tf.float32),
            "codon": tf.FixedLenFeature([config["train_eval"]["codonlength"]], tf.float32),
            "minhash": tf.FixedLenFeature([config["train_eval"]["minhashlength"]], tf.float32)}
        parsed = tf.parse_single_example(record, keys_to_features)
        return {'kmer': parsed['kmer'], 'codon': parsed['codon'], 'minhash': parsed['minhash']}, parsed['label']
    dataset = dataset.map(parser)
    dataset = dataset.batch(512)
    dataset = dataset.repeat(1)
    iterator = dataset.make_one_shot_iterator()
    features, labels = iterator.get_next()
    return features, labels


# Model definitions
combined_estimator = tf.estimator.DNNLinearCombinedClassifier(
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

dnn_estimator = tf.estimator.DNNClassifier(
    model_dir = modeldir,
    n_classes=n_classes,
    weight_column=None,
    feature_columns=[kmer, codon],
    dropout=0.5,
    activation_fn=tf.nn.relu,
    hidden_units=[256, 32],
    optimizer='Adagrad')

def train(infiles, out,  modeldir, n_classes, configpath):
    """Main training function called by vica_cli trains a Tensorflow model
    returning a modeldir and TFmodel file used by  the tensorflow serving api

    """
    try:
        logging.info("Beginning tensorflow model training. to see results in real-time run 'tensorboard --logdir=path/to/log-directory'")
        with open(configpath, "r") as cf:
            global config
            config = yaml.load(cf)
        kmerdim, kmer, codon, minhash = _featureshape(config["train_eval"]["ksize"])
        filenames = [infiles]
        epochs = config["train_eval"]["epochs"]
        modeldir = modeldir
        global modeldir, epochs, kmerdim, kmer, codon, minhash, n_classes, filenames
        if config["train_eval"]["model"]args.model == "DNN":
            dnn_estimator.train(input_fn={train_input_fn: filenames})
        elif config["train_eval"]["model"] == "DNNLogistic":
            combined_estimator.train(input_fn={train_input_fn: filenames})
    except:
        loggin.exception(" during tensorflow model training the following exception occured:")

def eval(infiles, out,  modeldir, n_classes, configpath):
    """ Main evaluation function called by vica_cli. Load a model from a model directory
    returning a file of predictions. In the fiture it will have other evaluation datd.

    """
    try:
        logging.info("Beginning tensorflow model evaluation. to see results in real-time run 'tensorboard --logdir=path/to/log-directory'")
        with open(configpath, "r") as cf:
            global config
            config = yaml.load(cf)
        kmerdim, kmer, codon, minhash = _featureshape(args.ksize)
        filenames = [infiles]
        epochs = config["train_eval"]["epochs"]
        modeldir = modeldir
        global modeldir, epochs, kmerdim, kmer, codon, minhash, n_classes, filenames
        if config["train_eval"]["model"] == "DNN":
            preds = dnn_estimator.train(input_fn={test_input_fn: filenames})
        elif config["train_eval"]["model"] == "DNNLogistic":
            preds = combined_estimator.train(input_fn={test_input_fn: filenames})
        logging.info("Tensorflow model performance. See also {}.".format(out))
        logging.info(preds)
        if not os.path.exists(out)
            os.mkdir(out)
        predictions = os.path.join(out,"modelpredictions.txt")
        with open(predictions), "w") as outfile:
            csv_writer_instance = csv.writer(outfile, lineterminator='\n')
            for rec in preds:
                plist = rec['probabilities']
                pliststr = [str(x) for x in plist]
                ll = [rec['classes'][0].decode("utf-8"), str(rec['class_ids'][0])]
                ll.extend(pliststr)
                csv_writer_instance.writerow(ll)
    except:
        loggin.exception(" during tensorflow model evaluation the following exception occured:")
