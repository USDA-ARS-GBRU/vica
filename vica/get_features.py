#!/usr/bin/env python3
'''prodigal.py: a module to call genes with prodigal then count codon usage
and transform into centered log ratio returning values as a CSV'''

import argparse
import tempfile
import os
import vica.shred
import vica.minhash
import vica.khmer_features
import vica.prodigal
import vica.tfrecord_maker
import glob
from pyfaidx import Fasta
import shutil

def main():
    parser = argparse.ArgumentParser(description='A script to take a fasta with contigs and extract fratures')
    parser.add_argument('--input', help="A fasta file", required=True)
    parser.add_argument('--output', help="path to tfrecord file written", required=True)
    parser.add_argument('--nodes', help="path to NCBI taxonomy nodes file", required=True)
    parser.add_argument('--label', help="An integer representing the classification label", type=int, required=True)
    parser.add_argument('--shred', help="A flag to shred data into randomly sampled contigs. Default is false", action="store_true")
    parser.add_argument('--length', help="The length of the genome subsamples if fixed is selected", default = 5000, type=int)
    parser.add_argument('--samples', help="Total number of shreded contigs to create, or if between 0 and 1, the proportion of the genome to sample", default = 0.5, type=float)
    parser.add_argument('--testing', help="Testing mode",  action="store_true")
    parser.add_argument('--minhashlocal', help="Use local minhash database instead of remote server",  action="store_true")
    parser.add_argument('--minhashrefs', help=" a directory containing a local minhash database", default=".")
    parser.add_argument('--minhashtree', help="local minhash tree")
    parser.add_argument('--minhashblacklist', help="local minhash blacklist")
    args = parser.parse_args()

    #create workng directory and file names
    #dtemp = "/Users/rivers/Desktop/tmp"
    dtemp = tempfile.mkdtemp()
    segments = os.path.join(dtemp, "segments.fasta")
    reflist = glob.glob(os.path.join(args.minhashrefs, "*.sketch"))
    if reflist:
        refs = ",".join(reflist)
    minhashout = os.path.join(dtemp,"minhashout.txt")
    kmerout = os.path.join(dtemp,"kmerout.csv")
    codonout = os.path.join(dtemp, "codonout.csv")
    ksize = 4

    # shred gneomes into contigs
    if args.shred:
        inhandle = Fasta(args.input, read_ahead=10000)
        with open( segments, 'w') as outhandle:
            samples_all = vica.shred. shred_all(inhandle=inhandle,
                                                outhandle=outhandle,
                                                samples=args.samples,
                                                samplemethod="fixed",
                                                testing=args.testing,
                                                length=args.length,
                                                shape=1.333,
                                                loc=3000,
                                                scale=1140)
    else:
        segments = args.input

    # Extract minhash features
    if args.minhashlocal:
        vica.minhash.minhashlocal(dtemp=dtemp, infile=segments, outfile=minhashout,
                     nodesfile=args.nodes, localsketchdir=args.minhashrefs,
                     blacklist=args.minhashblacklist, tree=args.minhashtree)
    else:
        vica.minhash.minhashremote(dtemp=dtemp, infile=segments, outfile=minhashout,
                     nodesfile=args.nodes)


    # Extract kmer features
    kmers = vica.khmer_features.iterate_kmer(k=ksize)
    vica.khmer_features.write_kmers_as_csv(infile=segments, outfile=kmerout,
                              ksize=ksize, kmers=kmers)

    # Extract codons
    vica.prodigal.contigs_to_feature_file(infile=segments, outfile=codonout)

    #shutil.copytree(dtemp, "dtempout")

    # Combine data into a Tensorflow TF record file
    vica.tfrecord_maker.convert_to_tfrecords(dtemp=dtemp, kmerfile=kmerout, codonfile=codonout,
             minhashfile=minhashout, tfrecordfile=args.output,
             label=str(args.label), sort=True)


if __name__ == '__main__':
    main()
