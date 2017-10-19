#!/usr/bin/env python3
'''tfrecord_maker.py: a module to create tf record files'''
import tensorflow as tf
import subprocess
import itertools
import tempfile
import os
import csv
import itertools
import argparse

def external_sort(infile, outfile, sep, size=1048576, key=1):
    options = ['sort', '--field-separator=' + sep, '--buffer-size=' + str(size) , '--key=' + str(key), '--output=' + outfile, infile ]
    subprocess.run(options)


def python_join(infile1, infile2, outfile, sep):
    with open(outfile, 'w', newline='') as csvfile, open(infile1, 'r') as f1, open(infile2, 'r') as f2:
        csv_writer_instance = csv.writer(csvfile, delimiter=sep, quoting=csv.QUOTE_MINIMAL)
        for line1, line2 in zip(f1, f2):
            l1 = line1.strip().split(',')
            l2 = line2.strip().split(',')
            assert l1[0] == l2[0], "the ids for the codon and kmer data do not match: {} and {}".format(l1, l2)
            ll = l1 + l1[1:]
            csv_writer_instance.writerow(ll)


def combine_dict_csv_to_tfrecords(label, minhashdict, csvfile,  tfrecordfile):
    '''convert a dictionary of minhash data and a sorted csv file into a tfrecord file'''
    writer = tf.python_io.TFRecordWriter(tfrecordfile)
    with open(csvfile, 'r') as file:
        for line in file:
            ll1 = line.split(',')
            ll1lab = ll1[0]
            ll1data = ll1[1:]
            densefeatures = [float(i) for i in ll1data]
            if minhashdict[ll1lab]:
                taxids = sorted(minhashdict[ll1lab])
                values =  [float(minhashdict[ll1lab][k]) for k in taxids]
            else:
                taxids = []
                values = []
            example = tf.train.Example(features=tf.train.Features(feature={
                "label":
                    tf.train.Feature(int64_list=tf.train.Int64List(value=[int(label)])),
                "densefeatures":
                    tf.train.Feature(float_list=tf.train.FloatList(value=densefeatures)),
                "minhashids":
                    tf.train.Feature(int64_list=tf.train.Int64List(value=taxids)),
                "minhashvalues":
                    tf.train.Feature(float_list=tf.train.FloatList(value=values))
                }))
            writer.write(example.SerializeToString())
    writer.close()
    print("Successfully converted minhash data and {} to tfrecord file named {}".format( csvfile,
                                                                  tfrecordfile))


def convert_to_tfrecords(kmerfile, codonfile, minhashfile, tfrecordfile, label):
    dtemp = tempfile.mkdtemp()
    ksorted = os.path.join(dtemp, "kmer_sorted.csv")
    csorted = os.path.join(dtemp, "codon_sorted.csv")
    densefile = os.path.join(dtemp, "dense.csv")
    minhashdict = parse_sketchout(minhashfile)
    external_sort(infile=kmerfile, outfile=ksorted, sep=",")
    external_sort(infile=codonfile, outfile=csorted, sep=",")
    python_join(infile1=ksorted, infile2=csorted, outfile=densefile, sep=",")
    combine_dict_csv_to_tfrecords(label=label, minhashdict=minhashdict,
                                    csvfile=densefile,  tfrecordfile=tfrecordfile)

def main():

    parser = argparse.ArgumentParser(description='A script to generate a tfrecord file from feature files')
    parser.add_argument('--kmerin', help="A csv of kmer frequencies")
    parser.add_argument('--codonin', help="A csv of codon frequencies")
    parser.add_argument('--minhashin', help="A libsvm file of taxonids and minhash matches")
    parser.add_argument('--outfile', help= "A tfrecord file")
    parser.add_argument('--label', help= "An interger label of the class")
    args = parser.parse_args()

    convert_to_tfrecords(kmerfile=args.kmerin, codonfile=args.codonin,
                         minhashfile=args.minhashin, tfrecordfile=args.outfile,
                         label=args.label)

if __name__ == '__main__':
    main()
