#!/usr/bin/env python3
'''tfrecord_maker.py: a module to create tf record files'''
import tensorflow as tf
import subprocess
import itertools
import tempfile
import os
import csv
import argparse
import numpy as np
import logging

def external_sort(infile, outfile, sep, key=1):
    '''Externally sort and make unique csv files using built-in linux utilities'''
    try:
        sortoptions = ['sort', '-t', sep, '-k', str(key), '-s','-u', '-o', outfile,  infile ]
        p1 = subprocess.run(sortoptions, check=True,)
        return outfile
    except:
        logging.exception("Input files could not be sorterd")

def join(kmerfile, codonfile, minhashfile, dtemp):
    '''Externally join with built-in linux utilities in hte order label, kmers,codons,minhash'''
    kcfile= os.path.join(dtemp, "kcfile.csv")
    mergefile = os.path.join(dtemp, "mergefile.csv")
    try:
        with open(kcfile, 'w') as kcf:
            options = ['join', '-t', ',', '-1', '1', '-2', '1', kmerfile, codonfile]
            subprocess.run(options,  check=True, stdout=kcf)
        with open(mergefile, "w") as  mf:
            options2 = ['join', '-t', ',', '-1', '1', '-2', '1', kcfile, minhashfile]
            subprocess.run(options2,  check=True, stdout=mf)
        os.remove(kcfile)
        return mergefile
    except RuntimeError:
        logging.exception("Could not merge csv files using unix join command")

def count_features(**kwargs):
    '''given a dictionary of fileypes and file locations return a dictionary of feature lengths'''
    featuredict ={}
    for key, val in kwargs.items():
        with open(val, 'r') as f:
            features = len(f.readline().strip().split(",")) - 1
        featuredict[key] = features
    return featuredict


def csv_to_tfrecords(kmerfile, codonfile, minhashfile, mergefile, tfrecordfile, label):
    '''convert csv files of features to ttfrecord files'''
    writer = tf.python_io.TFRecordWriter(tfrecordfile)
    features = count_features(kmers=kmerfile,codons=codonfile,minhash=minhashfile)
    kstart = 1
    kend = features['kmers'] + 1
    cend = kend + features['codons']
    i = 0
    with open(mergefile, 'r') as mergedata:
        for i, lstring in enumerate(mergedata, 1):
            line = lstring.strip().split(",")
            lab = line[0]
            kdat = np.array(line[kstart:kend], dtype='float32')
            cdat = np.array(line[kend:cend], dtype='float32')
            mdat = np.array(line[cend:], dtype='float32')
            example = tf.train.Example(features=tf.train.Features(feature={
                "id":
                    tf.train.Feature(bytes_list=tf.train.BytesList(value=[tf.compat.as_bytes(lab)])),
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
    logging.info("Successfully converted {} records to to TFrecord format".format(i))


def convert_to_tfrecords(dtemp, kmerfile, codonfile, minhashfile,
                         tfrecordfile, label, sort=False):
    if sort:
        ksorted = os.path.join(dtemp, "kmer_sorted.csv")
        csorted = os.path.join(dtemp, "codon_sorted.csv")
        msorted = os.path.join(dtemp, "minhash_sorted.csv")
        mergefile = os.path.join(dtemp, "mergefile.csv")
        external_sort(infile=kmerfile, outfile=ksorted, sep=",")
        external_sort(infile=codonfile, outfile=csorted, sep=",")
        external_sort(infile=minhashfile, outfile=msorted, sep=",")
    else:
        ksorted = kmerfile
        csorted = codonfile
        msorted = minhashfile
    mergefile = join(kmerfile=ksorted, codonfile=csorted, minhashfile=msorted, dtemp=dtemp)
    csv_to_tfrecords(kmerfile=ksorted, codonfile=csorted, minhashfile=msorted,mergefile=mergefile,
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
    # dtemp = '/Users/rivers/Desktop/tfrecordtest'
    convert_to_tfrecords(dtemp= dtemp, kmerfile=args.kmerin, codonfile=args.codonin,
                         minhashfile=args.minhashin, tfrecordfile=args.outfile,
                         label=args.label, sort=args.sort)

if __name__ == '__main__':
    main()
