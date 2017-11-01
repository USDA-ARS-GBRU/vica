#!/usr/bin/env python3
'''classifier.py: a module to train models and classify contigs fro mtfrecords
   of features. It uses tensorflow's datasets api and estimator api'''
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os
import urllib

import numpy as np
import tensorflow as tf
import tempfile

tmpdir = tempfile.mkdtemp()
n_classes = 2



# filenames = ["tests/test-data/combined.tfrecord"]
kmer = tf.feature_column.numeric_column(key='kmer', shape=(135,))
codon = tf.feature_column.numeric_column(key='codon', shape=(177,))
minhash = tf.feature_column.numeric_column(key='minhash', shape=(266,))


def dataset_input_fn():
    #filenames = tf.placeholder(tf.string, shape=[None])
    filenames = ["tests/test-data/combined-test.tfrecord"]
    dataset = tf.contrib.data.TFRecordDataset(filenames)
    def parser(record):
        keys_to_features = {"label": tf.FixedLenFeature((), tf.int64),
            "kmer": tf.FixedLenFeature([135], tf.float32),
            "codon": tf.FixedLenFeature([177], tf.float32),
            "minhash": tf.FixedLenFeature([266], tf.float32)}
        parsed = tf.parse_single_example(record, keys_to_features)
        return {'kmer': parsed['kmer'], 'codon': parsed['codon'], 'minhash': parsed['minhash']}, parsed['label']
    dataset = dataset.map(parser)
    dataset = dataset.shuffle(buffer_size=10000)
    dataset = dataset.batch(32)
    dataset = dataset.repeat(1000)
    iterator = dataset.make_one_shot_iterator()
    # `features` is a dictionary in which each value is a batch of values for
    # that feature; `labels` is a batch of labels.
    features, labels = iterator.get_next()
    return features, labels

estimator1 = tf.estimator.DNNClassifier(
    model_dir = "DNNdir",
    n_classes=n_classes,
    # deep settings
    feature_columns=[kmer, codon],
    dropout=0.5,
    activation_fn=tf.nn.relu,
    hidden_units=[80, 20],
    optimizer='Adagrad')

def _print_tfrecords(filename):
    for example in tf.python_io.tf_record_iterator(filename):
        result = tf.train.Example.FromString(example)
        print(result)

estimator = tf.estimator.DNNLinearCombinedClassifier(
    model_dir = "logDNNdir",
    n_classes=n_classes,
    weight_column=None,
    linear_feature_columns=[minhash],
    linear_optimizer='Ftrl',
    dnn_feature_columns=[kmer, codon],
    dnn_dropout=0.5,
    dnn_activation_fn=tf.nn.relu,
    dnn_hidden_units=[80, 20],
    dnn_optimizer='Adagrad')

model1 = estimator.train(input_fn=dataset_input_fn)
