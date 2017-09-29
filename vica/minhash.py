#!/usr/bin/env python3
'''minhash.py: a module to run bbtools minhash function on a set of data and
return a sparse vector of classification data'''

import subprocess
import os

import tensorflow as tf

# Constants

jgi_server_url='https://refseq-sketch.jgi-psf.org/sketch'

def send_sketch(file, outfile, tempdir):
    '''Runs bbtools sendsketch.sh on a file of sequences returning a classification for each'''
    options = ["sendsketch.sh",
               "in=" + file,
               "out=" + outfile,
               "address=" + jgi_server_url,
               "mode=sequence",
               "level=3",
               "color=f",
               "overwrite=t",
               "printani=f",
               "completeness=f",
               "score=t",
               "printmatches=t",
               "printgsize=f",
               "printgseqs=f",
               "printtaxname=f",
               "printname0=f",
               "printcontam=f",
               "printunique=t",
               "printnohit=f",
               "printtaxid=t"]
    sendsketchout = subprocess.run(options, stderr=subprocess.PIPE)
    return sendsketchout.stderr.decode('utf-8')
    #return sendsketchout

def parse_sketch(file, outfile, tempdir):
    '''parses bbtools sendsketch output returning python dictionary'''
    try:
        tempdf = {}
        with open(file, 'r') as f:
            for i, line in enumerate(f):
                if line.strip() == '':
                    next
                elif line.startswith("Query:"):
                    ll = line.strip().split("\t")
                    key1 = ll[0].split(":")[1].split()[0]
                    tempdf[key1] = {}
                elif line.startswith("WKID"):
                    next
                else:
                    ll2 = line.strip().split("\t")
                    tempdf[key1][int(ll2[5])] = int(ll2[3])
        return tempdf
    except IOError:
        print("could not parse of the sketch file %s" % (file))

def load_labels(labelfile):
    '''load labels from a two column id, label return a dict of sets'''
    labeldict={}
    with open(labelfile, 'r') as labfile:
        for line in labfile:
            ll = line.strip().split()
            if ll[0] in labeldict:
                continue
            else:
                labeldict[ll[0]] = ll[1]
    return labeldict

def parse_sketch_to_tfrecord(sketchfile, tfrecordfile, labeldict):
    writer = tf.python_io.TFRecordWriter(tfrecordfile)
    try:
        with open(sketchfile, 'r') as f:
            label = None
            taxidlist = []
            matchlist = []
            for line in f:
                if line.strip() == '':
                    # write tf record examples
                    pass
                    # 3 reset lists and label
                    taxidlist = []
                    matchlist = []
                    label = None
                elif line.startswith("Query:"):
                    ll = line.strip().split("\t")
                    key1 = ll[0].split(":")[1].split()[0]
                    label = labeldict[key1]
                elif line.startswith("WKID"):
                    next
                else:
                    ll2 = line.strip().split("\t")
                    taxidlist.append(int(ll2[5]))
                    matchlist.append(int(ll2[3]))
    except IOError:
        print("could not parse of the sketch file %s" % (file))
