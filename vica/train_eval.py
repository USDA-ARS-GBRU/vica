#!/usr/bin/env python3
"""train_eval.py: a module to train models and evaluate models from tfrecords
   of features. It uses The Tensorflow 1.3+ datasets api and estimator api

"""


import os
import urllib
import tempfile
from collections import Counter
import functools
import logging

import yaml
import numpy as np
import tensorflow as tf
import csv

import vica

with open(vica.CONFIG_PATH) as cf:
    config = yaml.load(cf)

def _featureshape(k=5, codonlength=177, minhashlength=267):
    """Determine the shape of the features for each feature type including
    for kmers of different lengths."""
    kmerdim = len(vica.khmer_features.iterate_kmer(k)) - 1
    kmer = tf.feature_column.numeric_column(key='kmer', shape=(kmerdim))
    codon = tf.feature_column.numeric_column(key='codon', shape=(codonlength))
    minhash = tf.feature_column.numeric_column(key='minhash', shape=(minhashlength))
    return kmerdim, kmer, codon, minhash

global modeldir, epochs, kmerdim, kmer, codon, minhash, n_classes, filenames

def base_input_fn(codonlength, minhashlength, kmerdim, shuffle, shuffle_buffer_size, batch, epochs, filenames):
    """the function for feeding and processing training data"""
    # filenames = tf.placeholder(tf.string, shape=[None])
    dataset = tf.contrib.data.TFRecordDataset(filenames)
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

def test_input_fn():
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
    dataset = dataset.batch(512)
    dataset = dataset.repeat(1)
    iterator = dataset.make_one_shot_iterator()
    features, labels = iterator.get_next()
    return features, labels


# Model definitions
def mk_dnnlogistic_estimator(modeldir, n_classes, minhash, kmer, codon):
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
    """Main training function called by vica_cli trains a Tensorflow model
    returning a modeldir and TFmodel file used by  the tensorflow serving api

    """
    try:
        logging.info("Beginning tensorflow model training. To see results in real-time run 'tensorboard --logdir={}'".format(modeldir))
        with open(configpath, "r") as cf:
            config = yaml.load(cf)
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
        logging.exception("During tensorflow model training the following exception occured:")
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


def eval(infiles, out, modeldir, n_classes, configpath):
    """ Main evaluation function called by vica_cli. Load a model from a model directory
    returning a file of predictions. In the fiture it will have other evaluation datd.

    """
    try:
        logging.info("Beginning tensorflow model evaluation. To see results in real-time run 'tensorboard --logdir={}'".format(modeldir))
        with open(configpath, "r") as cf:
            config = yaml.load(cf)
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
        with open(predictions, "w") as outfile:
            csv_writer_instance = csv.writer(outfile, lineterminator='\n')
            for rec in preds:
                plist = rec['probabilities']
                pliststr = [str(x) for x in plist]
                ll = [rec['classes'][0].decode("utf-8"), str(rec['class_ids'][0])]
                ll.extend(pliststr)
                csv_writer_instance.writerow(ll)
        for key in sorted(results):
            logging.info('{}: {}'.format(key, results[key]))
    except:
        logging.exception("During tensorflow model evaluation the following exception occured:")
