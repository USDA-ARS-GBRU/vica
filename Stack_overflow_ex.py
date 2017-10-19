#!/usr/bin/env python3
'''Example code to reprodce TFRecord import error using hte Datasets api in stack overflow'''
import os
import subprocess
import tensorflow as tf # version 1.3

#CSV

def csv_to_tfrecords(label, csvfile, tfrecordfile):
    '''convert a csv file with 1 text ID column follwed by 270 float columns into a tfrecord file'''
    writer = tf.python_io.TFRecordWriter(tfrecordfile)
    with open(csvfile, 'r') as file:
        for line in file:
            lline = line.split(',')
            ldata = lline[1:]
            densefeatures = [float(i) for i in ldata]
            example = tf.train.Example(features=tf.train.Features(feature={
                "label":
                    tf.train.Feature(int64_list=tf.train.Int64List(
                                     value=[int(label)])),
                "densefeatures":
                    tf.train.Feature(float_list=tf.train.FloatList(
                                     value=densefeatures)),
                }))
            writer.write(example.SerializeToString())
    writer.close()
    print("Successfully converted {} to tfrecord file named {}".format( csvfile, tfrecordfile))

# Create TFRecords
csv_to_tfrecords(label=0,
                 csvfile="tests/test-data/human_herpesvirus_3-kmercodon-sorted.csv",
                 tfrecordfile="tests/test-data/human_herpesvirus_3-test.tfrecord")
csv_to_tfrecords(label=1,
                 csvfile="tests/test-data/med4-kmercodon-sorted.csv",
                 tfrecordfile="tests/test-data/med4-test.tfrecord")

# Concatenate TFRecords with different labels

subprocess.call("cat tests/test-data/human_herpesvirus_3-test.tfrecord tests/test-data/med4-test.tfrecord > tests/test-data/combined-test.tfrecord", shell=True)

def dataset_input_fn():
    '''process Tfrecord file for input into Tensorflow estimator'''
    filenames = ["tests/test-data/combined-test.tfrecord"]
    dataset = tf.contrib.data.TFRecordDataset(filenames)
    def parser(record):
        keys_to_features = {"label": tf.FixedLenFeature((), tf.int64),
            "densefeatures": tf.FixedLenFeature([270],tf.float32)}
        parsed = tf.parse_single_example(record, keys_to_features)
        return {'densefeatures': parsed["densefeatures"]}, parsed['label']
    dataset = dataset.map(parser)
    dataset = dataset.shuffle(buffer_size=10000)
    dataset = dataset.batch(32)
    dataset = dataset.repeat(100)
    iterator = dataset.make_one_shot_iterator()
    # `features` is a dictionary in which each value is a batch of values for
    # that feature; `labels` is a batch of labels.
    features, labels = iterator.get_next()
    return features, labels

# Identify the type of data column
densefeatures = tf.feature_column.numeric_column(key='densefeatures', shape=(270,))

# Create a Tensorflow canned estimator
estimator = tf.estimator.DNNClassifier(
    model_dir="tests/test-data/DNNdir",
    n_classes=2,
    feature_columns=[densefeatures],
    dropout=0.5,
    activation_fn=tf.nn.relu,
    hidden_units=[80, 20],
    optimizer='Adagrad')

# Call the estimator for training
estimator.train(input_fn=dataset_input_fn)
