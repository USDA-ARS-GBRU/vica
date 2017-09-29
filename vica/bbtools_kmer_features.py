#!/usr/bin/env python3
'''kmer.py: a module to run bbtools kmercountexact function on a set of data and
return a dense vector of classification data'''

import subprocess
import os
import logging
import pandas
import itertools
from Bio import SeqIO
from collections import OrderedDict

min = 4
max = 8

def limit_kmers(min, max, k):
    '''Limit the range of kmer values'''
    assert isinstance(k, int), \
        "kmer value % is not an integer" % (k)
    assert k >= min and k<= max \
        "kmer value % is not between % and %" % (k, min, max)
    continue

def initialize_dict(kmerlist):
    '''create an ordered dict of all possible kmers''''
    kmer_dict={}
    for kmer in kmerlist:
        kmer_dict[kmer]=0
    return OrderedDict(sorted(td.items(), key=lambda t: t[0]))

def iterate_kmer(k=5):
    limit_kmers(min=min, max=max, k=k)
    '''get the list of kmers'''
    bases = ['A','C','T','G']
    kmers = [''.join(p) for p in itertools.product(bases, repeat=k)]
    core_kmer = []
    for kmer in kmers:
        if not str(Seq(kmer).reverse_complement()) in core_kmer:
            core_kmer.append(kmer)
    return core_kmer


def kmer_count_exact(file, filehandle, k=5, kmer_dict):
    '''kmercountexact file of sequences returning a classification for each'''
    limit_kmers(min=min, max=max, k=k)
        options = ["kmercountexact.sh",
               "in=" + file,
               "k=" + k,
               "fastadump=f",
               "out=stdout"]
    kmercountexactout = subprocess.run(options, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    # logging.info(sendsketchout.stderr.decode('utf-8'))
    reader = csv.reader(kmercountexactout.stdout, delimiter="\t")
    for line in reader:
        ll = line.strip().split()
        kmer_dict[ll[0]] = ll[1]
    kmertotal = sum(kmer_dict.values())
    return list(map(lambda x: x/kmertotal, kmer_dict.values()))

def process_fasta(infile, tempdir, out, k):
         seqs = SeqIO.parse(infile)
         for record in seqs:
             # unfinised


    '''process all sequences in a file returning a csv with data'''
    limit_kmers(min,max,k)
    kmer_list = iterate_kmer(k)
    kmer_dict = initialize_dict(kmer_list)
    for

def python_kmer_count(infile, tempdir, out, k):
    #unfinishesd
    pass
