#!/usr/bin/env python


import itertools
import csv
import logging

import khmer
from Bio import SeqIO
from Bio.Seq import Seq

import vica


def iterate_kmer(k):
    """ get the list of tetramers"""
    try:
        bases = ['A','C','T','G']
        kmers = [''.join(p) for p in itertools.product(bases, repeat=k)]
        core_kmer = []
        for kmer in kmers:
            if not str(Seq(kmer).reverse_complement()) in core_kmer:
                core_kmer.append(kmer)
        return core_kmer
    except:
        logging.exception("Could not calculate the list of kmers for k = {}".format(k))

def get_composition(ksize, seq, kmers, norm):
    """ get the composition profile and return a list of kmer counts or normalized kmer counts"""
    try:
        nkmers = 4**ksize
        tablesize = nkmers + 100
        counting_hash = khmer.Countgraph(ksize, tablesize, 1)
        counting_hash.consume(seq)
        composition = [counting_hash.get(kmer) for kmer in kmers]
        if norm == True:
            total = sum(composition)
            nc = []
            for item in composition:
                if item == 0:
                    nc.append(0.0)
                else:
                    nc.append(float(item)/float(total))
                composition = nc
        return composition
    except:
        logging.exception("Could not calculate composition using khmer")


def write_kmers_as_csv(infile, outfile, ksize, kmers):
    try:
        with open(infile, 'r') as f1:
            with open(outfile, 'w') as csvfile:
                mywriter = csv.writer(csvfile, lineterminator='\n')
                header = ["id"]
                header.extend(kmers)
                # mywriter.writerow(header)
                ksize = int(ksize)
                kmers = iterate_kmer(ksize)
                rn = 0
                for record in SeqIO.parse(f1, 'fasta'):
                    rl = [record.id]
                    kmer_frequency = get_composition(ksize,str(record.seq).upper(), kmers, False)
                    kmer_ilr = vica.prodigal.ilr(kmer_frequency)
                    rl.extend(kmer_ilr)
                    mywriter.writerow(rl)
                    rn += 1
                logging.info("Wrote {} kmer records to {}.".format(rn, outfile))
    except:
        logging.exception("Could not write kmer profiles to file")

def run(infile, outfile, ksize):
    kmers = iterate_kmer(ksize)
    logging.info("identifying kmer features with a k of {}".format(ksize))
    write_kmers_as_csv(infile=infile, outfile=outfile, ksize=ksize, kmers=kmers)
