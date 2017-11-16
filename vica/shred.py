#!/usr/bin/env python3

import os
import argparse
import subprocess
import numpy as np
import scipy.stats
import argparse
from pyfaidx import Fasta
import random
import logging

## Functions

def parse_samples(samples):
    '''Takes input value for samples and returns type all, fraction or fixed'''
    try:
        def eq( a, b, eps=0.0001):
            return abs(a - b) <= eps
        if 0 < samples < 1:
            return 'fraction'
        if eq(samples, 1.0):
            return 'all'
        if samples > 1:
            return "set_number"
    except:
        logging.exception("the input parameter 'samples' had an unexpected input")

def writeseq(record, pos, length, handle):
    '''writes a fasta sequence to a file, adding position information to the id'''
    seqlist = []
    label = (">" + record.name + "|pos|" + str(pos) + ".." + str(pos + length) + "\n")
    end = pos + length
    result = record[pos:end].seq
    seqlist = [result[n:n + 60] +'\n' for n in range(0, len(result), 60)]
    handle.write(label)
    handle.writelines(seqlist)

def _calc_seg_len(samplemethod, length, shape, loc, scale):
    '''Returns a length value for a sequence'''
    try:
        if samplemethod == 'fixed':
            seglen = length
        elif samplemethod == 'lognorm':
            seglen = int(scipy.stats.lognorm.rvs(s=shape, loc=loc, scale=scale))
        return seglen
    except:
        logging.exception("Could not calculate sequence length")

def shred_all(inhandle, outhandle, samples, samplemethod, testing, length=5000, shape=1.333, loc=3000, scale=1140):
    if testing:
        seed = 435903852
        logging.debug("Shred.py is using testing mode a with seed of {}".format(seed))
        np.random.seed(seed=seed)
    sample_type = parse_samples(samples)
    if sample_type == 'fraction':
        rw = 0
        for i, record in enumerate(inhandle):
            length_used = 0
            seq_length = len(record)
            while length_used/seq_length < samples:
                seglen = _calc_seg_len(samplemethod, length, shape, loc, scale)
                if seq_length > seglen:
                    pos = np.random.choice(seq_length - seglen)
                    endpos = pos+ seglen
                    if record[pos: seglen + pos].seq.count("N")/seglen < 0.1:
                        writeseq(record, pos, seglen, outhandle)
                        rw += 1
                length_used += seglen
        logging.info("Processed {} input sequences and wrote out {} shredded sequences".format(i, rw))
    elif sample_type =='set_number':
        rw = 0
        for i, record in enumerate(inhandle):
            seq_length = len(record)
            samples_written = 0
            while samples_written < samples:
                seglen = _calc_seg_len(samplemethod, length, shape, loc, scale)
                if seq_length > seglen:
                    pos = np.random.choice(seq_length - seglen)
                    endpos = pos+length
                    if record[pos: seglen + pos].seq.count("N")/seglen < 0.1:
                        writeseq(record, pos, seglen, outhandle)
                        rw += 1
                samples_written += 1

        logging.info("Processed {} input sequences and wrote out {} shredded sequences".format(i, rw))
    elif sample_type =='all':
        rw = 0
        for i, record in enumerate(inhandle):
            seq_length = len(record)
            pos = 0
            seglen = _calc_seg_len(samplemethod, length, shape, loc, scale)
            while pos + seglen < seq_length:
                if record[pos: seglen + pos].seq.count("N")/seglen < 0.1:
                    writeseq(record, pos, seglen, outhandle)
                    rw += 1
                pos += seglen
                seglen = _calc_seg_len(samplemethod, length, shape, loc, scale)
        logging.info("Processed {} input sequences and wrote out {} shredded sequences".format(i, rw))


def parser():
        parser = argparse.ArgumentParser(description='A script to shred a genome file into sizes specified by a fixed length or a lognormal distbution')
        subparsers = parser.add_subparsers(help ="Commands", dest= 'command')
        fixed_parser = subparsers.add_parser('fixed', help="Divide all sequences into segments of fixed length")
        lognorm_parser = subparsers.add_parser('lognorm', help="Divide all sequences into segments with log-normally distibuted length")
        parser.add_argument('--samples', help="The total number of fragments to create per contig, or if between 0 and 1, the proportion of the contigs to sample.",
                                         default = 0.5,
                                         type=float)
        parser.add_argument('--input', help="A fasta file", required=True)
        parser.add_argument('--output', help= "Name of output file written",
                                        type=argparse.FileType('w'), default='-')
        parser.add_argument('--testing', help="Testing mode, use a fixed seed for repodicability",  action="store_true")
        fixed_parser.add_argument('--length', help="The length of the genome subsamples if fixed is selected",
                                               default = 5000,
                                               type=int)
        fixed_parser.set_defaults(whichmethod='fixed')
        lognorm_parser.add_argument('--shape', help="Shape parameter of lognormal distribution, size", default =1.333, type=float)
        lognorm_parser.add_argument('--loc', help= "Offset, or minimum contig length allowed determined from fitting", default = 3000, type=int)
        lognorm_parser.add_argument('--scale', help="Scale parameter of gamma distribution, theta", default = 1140, type=float)
        lognorm_parser.set_defaults(whichmethod='lognorm')
        args = parser.parse_args()
        return args


def main():
    args = parser()
    inhandle = Fasta(args.input, read_ahead=10000)
    if args.whichmethod == 'fixed':
        shred_all(inhandle=inhandle, outhandle=args.output, samples=args.samples, samplemethod=args.whichmethod, testing=args.testing, length=args.length)
    elif args.whichmethod == 'lognorm':
        shred_all(inhandle=inhandle, outhandle=args.output, samples=args.samples, samplemethod=args.whichmethod, testing=args.testing, shape=args.shape, loc=args.loc, scale=args.scale)
if __name__ == '__main__':
    main()
