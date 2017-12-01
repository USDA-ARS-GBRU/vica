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

import vica.khmer_features


def _featureshape(k):
    """Determine the shape of the features for each feature type including
    for kmers of different lengths."""
    codonlength = 177
    minhashlength= 267
    kmerdim = len(vica.khmer_features.iterate_kmer(k)) - 1
    kmer = tf.feature_column.numeric_column(key='kmer', shape=(kmerdim))
    codon = tf.feature_column.numeric_column(key='codon', shape=(codonlength))
    minhash = tf.feature_column.numeric_column(key='minhash', shape=(minhashlength))
    return kmerdim, kmer, codon, minhash


def train_eval_input_fn():
    #filenames = tf.placeholder(tf.string, shape=[None])
    filenames = ["vabn.train.k5.tfrecord"]
    dataset = tf.contrib.data.TFRecordDataset(filenames)
    def parser(record):
        keys_to_features = {"id": tf.FixedLenFeature((), tf.string),
            "label": tf.FixedLenFeature((), tf.int64),
            "kmer": tf.FixedLenFeature([kmerdim], tf.float32),
            "codon": tf.FixedLenFeature([177], tf.float32),
            "minhash": tf.FixedLenFeature([267], tf.float32)}
        parsed = tf.parse_single_example(record, keys_to_features)
        return {'kmer': parsed['kmer'], 'codon': parsed['codon'], 'minhash': parsed['minhash']}, parsed['label']
    dataset = dataset.map(parser)
    dataset = dataset.shuffle(buffer_size=10000)
    dataset = dataset.batch(32)
    dataset = dataset.repeat(epochs)
    iterator = dataset.make_one_shot_iterator()
    features, labels = iterator.get_next()
    return features, labels

def old_input_fn():
    filenames = tf.placeholder(tf.string, shape=[None])
    dataset = tf.contrib.data.TFRecordDataset(filenames)
    def parser(record):
        keys_to_features = {"label": tf.FixedLenFeature((), tf.int64),
            "kmer": tf.FixedLenFeature([kmerdim], tf.float32),
            "codon": tf.FixedLenFeature([177], tf.float32),
            "minhash": tf.FixedLenFeature([267], tf.float32)}
        parsed = tf.parse_single_example(record, keys_to_features)
        return {'kmer': parsed['kmer'], 'codon': parsed['codon'], 'minhash': parsed['minhash']}, parsed['label']
    dataset = dataset.map(parser)
    dataset = dataset.shuffle(buffer_size=10000)
    dataset = dataset.batch(32)
    dataset = dataset.repeat(epochs)
    iterator = dataset.make_one_shot_iterator()
    features, labels = iterator.get_next()
    return features, labels

def evaluate_input_fn():
    filenames = ["vabn.test.tfrecord"]
    dataset = tf.contrib.data.TFRecordDataset(filenames)
    def parser(record):
        keys_to_features = {"label": tf.FixedLenFeature((), tf.int64),
            "kmer": tf.FixedLenFeature([kmerdim], tf.float32),
            "codon": tf.FixedLenFeature([177], tf.float32),
            "minhash": tf.FixedLenFeature([267], tf.float32)}
        parsed = tf.parse_single_example(record, keys_to_features)
        return {'kmer': parsed['kmer'], 'codon': parsed['codon'], 'minhash': parsed['minhash']}, parsed['label']
    dataset = dataset.map(parser)
    #dataset = dataset.shuffle(buffer_size=10000)
    dataset = dataset.batch(512)
    dataset = dataset.repeat(1)
    iterator = dataset.make_one_shot_iterator()
    features, labels = iterator.get_next()
    return features, labels

def _print_tfrecords(n, filename):
    for i, example in enumerate(tf.python_io.tf_record_iterator(filename)):
        if i < n:
            result = tf.train.Example.FromString(example)
            return(result)
        else:
            break

def _count_tfrecords(filename):
    for i, example in enumerate(tf.python_io.tf_record_iterator(filename)):
        pass
    return i

def _count_labels(filename):
    """read tfrecords, return a count of"""
    def parser(record):
        keys_to_features = {"label": tf.FixedLenFeature((), tf.int64)}
        parsed = tf.parse_single_example(record, keys_to_features)
        return parsed['label']
    cnt = Counter()
    for i, example in enumerate(tf.python_io.tf_record_iterator(filename)):
        cnt[parser(example)] += 1
    cd = dict(cnt)
    total = sum(cd.values())
    wtdict = {k: v / total for k, v in cd.items()}
    return wtdict

combined_estimator = tf.estimator.DNNLinearCombinedClassifier(
    model_dir = modeldir,
    n_classes=n_classes,
    weight_column=None,
    linear_feature_columns=[minhash],
    linear_optimizer='Ftrl',
    dnn_feature_columns=[kmer, codon],
    dnn_dropout=0.5,
    dnn_activation_fn=tf.nn.relu,
    dnn_hidden_units=[80, 20],
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

def train():
    nclasses = args.nclasses
    kmerdim, kmer, codon, minhash = _featureshape(args.ksize)
    filenames = [args.inputs]
    epochs = args.epochs
    modeldir = args.output
    if args.model == "DNN":
        dnn_estimator.train(input_fn=train_eval_input_fn)
    elif args.model == "DNNLogistic":
        combined_estimator.train(input_fn=train_eval_input_fn)


def main():
    t0 = time.perf_counter()
    args = parser()
    # Set up logging
    logging.basicConfig(filename=args.logfile,
                        level=logging.DEBUG,
                        format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                        datefmt ='%Y-%m-%d %H:%M:%S')
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    # set a format which is simpler for console use
    formatter = logging.Formatter('%(asctime)s: %(levelname)-8s %(message)s',
                                 datefmt='%H:%M:%S')
    # tell the handler to use this format
    console.setFormatter(formatter)
    # add the handler to the root logger
    logging.getLogger('').addHandler(console)
    logging.info('Starting Vica feature selection workflow.')

estimator.train(input_fn=dataset_input_fn)
dnn_estimator.train(input_fn=train_eval_input_fn)
estimator.evaluate(input_fn=evaluate_input_fn)


dnnpreds = dnn_estimator.predict(input_fn=evaluate_input_fn)
with open("dnnpreds.csv", "w") as outfile:
    csv_writer_instance = csv.writer(outfile, lineterminator='\n')
    for rec in dnnpreds:
        ll = [rec['classes'][0].decode("utf-8"), rec['class_ids'][0],str(rec['probabilities'][0]),str(rec['probabilities'][1])]
        csv_writer_instance.writerow(ll)
with open("k5preds.csv", "w") as outfile:
    csv_writer_instance = csv.writer(outfile, lineterminator='\n')
    for rec in k5predictions:
        ll = [rec['classes'][0].decode("utf-8"), rec['class_ids'][0],str(rec['probabilities'][0]),str(rec['probabilities'][1])]
        csv_writer_instance.writerow(ll)
with open("class4preds.csv", "w") as outfile:
    csv_writer_instance = csv.writer(outfile, lineterminator='\n')
    for rec in class4predict:
        plist = rec['probabilities']
        pliststr = [str(x) for x in plist]
        ll = [rec['classes'][0].decode("utf-8"), str(rec['class_ids'][0])]
        ll.extend(pliststr)
        csv_writer_instance.writerow(ll)
