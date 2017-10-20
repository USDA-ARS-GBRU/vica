#!/usr/bin/env python3
'''tfrecord_maker.py: a module to create tf record files'''
import tensorflow as tf
import subprocess
import itertools
import tempfile
import os
import csv
import argparse

def external_sort(infile, outfile, sep, size=1048576, key=1):
    options = ['sort', '--field-separator=' + sep, '--buffer-size=' + str(size) , '--key=' + str(key), '--output=' + outfile, infile ]
    subprocess.run(options)


def csv_to_tfrecords(kmerfile, codonfile, minhashfile, tfrecordfile, label):
    '''convert csv files of features to ttfrecord files'''
    writer = tf.python_io.TFRecordWriter(tfrecordfile)
    with open(kmerfile, 'r') as kfile, open(codonfile, 'r') as cfile, open(minhashfile, 'r') as mfile:
        dataiter = zip(kfile, cfile, mfile)
        for kline, cline, mline in dataiter:
            kl = kline.split(',')
            cl = cline.split(',')
            ml = mline.split(',')
            assert kl[0] == cl[0] and cl[0] == ml[0], "the ids for the data files do not match. kmer id: {}, codon ID: {}, minhashid: {}".format(kl[0], cl[0], ml[0])
            lab = kl[0]
            kdat = [float(i) for i in kl[1:]]
            cdat = [float(i) for i in cl[1:]]
            mdat = [float(i) for i in ml[1:]]
            example = tf.train.Example(features=tf.train.Features(feature={
                "label":
                    tf.train.Feature(int64_list=tf.train.Int64List(value=[int(label)])),
                "kmer":
                    tf.train.Feature(float_list=tf.train.FloatList(value=kdat)),
                "codon":
                    tf.train.Feature(float_list=tf.train.FloatList(value=cdat)),
                "minhash":
                    tf.train.Feature(float_list=tf.train.FloatList(value=mdat))
                }))
            writer.write(example.SerializeToString())
    writer.close()
    print("Successfully converted data to the tfrecord file: {}".format(tfrecordfile))


def convert_to_tfrecords(dtemp, kmerfile, codonfile, minhashfile,
                         tfrecordfile, label, sort=False):
    if sort:
        ksorted = os.path.join(dtemp, "kmer_sorted.csv")
        csorted = os.path.join(dtemp, "codon_sorted.csv")
        msorted = os.path.join(dtemp, "minhash_sorted.csv")
        external_sort(infile=kmerfile, outfile=ksorted, sep=",")
        external_sort(infile=codonfile, outfile=csorted, sep=",")
        external_sort(infile=minhashfile, outfile=msorted, sep=",")
    else:
        ksorted = kmerfile
        csorted = codonfile
        msorted = minhashfile
    csv_to_tfrecords(kmerfile=ksorted, codonfile=csorted, minhashfile=msorted,
                     tfrecordfile=tfrecordfile, label=label)

def main():

    parser = argparse.ArgumentParser(description='A script to generate a tfrecord file from feature files')
    parser.add_argument('--kmerin', help="A csv of kmer frequencies")
    parser.add_argument('--codonin', help="A csv of codon frequencies")
    parser.add_argument('--minhashin', help="A libsvm file of taxonids and minhash matches")
    parser.add_argument('--outfile', help= "A tfrecord file")
    parser.add_argument('--label', help= "An interger label of the class")
    parser.add_argument('--sort', help= "A flag to sort input files", action="store_true")
    args = parser.parse_args()

    dtemp = tempfile.mkdtemp()
    convert_to_tfrecords(dtemp= dtemp, kmerfile=args.kmerin, codonfile=args.codonin,
                         minhashfile=args.minhashin, tfrecordfile=args.outfile,
                         label=args.label, sort=args.sort)

if __name__ == '__main__':
    main()
