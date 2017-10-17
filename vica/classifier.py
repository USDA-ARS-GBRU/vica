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



filenames = ["tests/test-data/combined.tfrecord"]
densefeatures = tf.feature_column.numeric_column(key='densefeatures', shape=(372,))



def dataset_input_fn():
    #filenames = tf.placeholder(tf.string, shape=[None])
    filenames = ["tests/test-data/combined.tfrecord"]
    dataset = tf.contrib.data.TFRecordDataset(filenames)
    def parser(record):
        keys_to_features = {"label": tf.FixedLenFeature((), tf.int64),
            "densefeatures": tf.FixedLenFeature((),tf.float32),
            "minhashids": tf.VarLenFeature(tf.int64),
            "minhashvalues": tf.VarLenFeature(tf.float32)}
        parsed = tf.parse_single_example(record, keys_to_features)
        return {'densefeatures': parsed["densefeatures"]}, parsed['label']
        #return {'densefeatures': parsed["densefeatures"], 'minhashids': parsed["minhashids"], 'minhashvalues' : parsed["minhashvalues"]}, parsed['label']
    # Use `tf.parse_single_example()` to extract data from a `tf.Example`
    # protocol buffer, and perform any additional per-record preprocessing.
    dataset = dataset.map(parser)
    #dataset = dataset.shuffle(buffer_size=100)
    #dataset = dataset.batch(32)
    #dataset = dataset.repeat(num_epochs)
    iterator = dataset.make_one_shot_iterator()
    # `features` is a dictionary in which each value is a batch of values for
    # that feature; `labels` is a batch of labels.
    features, labels = iterator.get_next()
    return features, labels

estimator = tf.estimator.DNNClassifier(
    model_dir = tmpdir,
    n_classes=n_classes,
    # deep settings
    feature_columns=[densefeatures],
    dropout=0.5,
    activation_fn=tf.nn.relu,
    hidden_units=[80, 20],
    optimizer='Adagrad')

for example in tf.python_io.tf_record_iterator("tests/test-data/combined.tfrecord"):
    result = tf.train.Example.FromString(example)
estimator.train(input_fn=dataset_input_fn)

estimator = tf.estimator.DNNLinearCombinedClassifier(
    model_dir = tmpdir,
    n_classes=n_classes,
    weight_column=None,
    # wide settings
    linear_feature_columns=[],
    linear_optimizer='Ftlr',
    # deep settings
    dnn_feature_columns=[densefeatures],
    dnn_dropout=0.5,
    dnn_activation_fn=tf.nn.relu,
    dnn_hidden_units=[80, 20],
    dnn_optimizer='Adagrad')
