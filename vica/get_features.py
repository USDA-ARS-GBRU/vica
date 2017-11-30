#!/usr/bin/env python3
"""prodigal.py: a module to call genes with prodigal then count codon usage
and transform into centered log ratio returning values as a CSV"""

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
import logging
import time
import datetime

def parser():
        parser = argparse.ArgumentParser(description='A script to take a fasta with contigs and extract features to a TFrecord file')
        parser.add_argument('--input', help="A fasta file", required=True)
        parser.add_argument('--output', help="path to tfrecord file written", required=True)
        parser.add_argument('--nodes', help="path to NCBI taxonomy nodes file", required=True)
        parser.add_argument('--label', help="An integer representing the classification label", type=int, required=True)
        parser.add_argument('--shred', help="A flag to shred data into randomly sampled contigs. Default is false", action="store_true")
        parser.add_argument('--length', help="The length of the genome subsamples if fixed is selected", default = 5000, type=int)
        parser.add_argument('--samples', help="Total number of shreded contigs to create, or if between 0 and 1, the proportion of the genome to sample", default = 0.5, type=float)
        parser.add_argument('--testing', help="Testing mode",  action="store_true")
        parser.add_argument('--ksize', help="Kmer size", default='5', choices=['4','5','6','7','8'])
        parser.add_argument('--logfile', help="location of log file", default="vica_feature_selection.log")
        parser.add_argument('--tempdir', help="designate a specific file as a temporary dir", default= None)
        parser.add_argument('--minhashlocal', help="Use local minhash database instead of remote server",  action="store_true")
        parser.add_argument('--minhashrefs', help=" a directory containing a local minhash database", default=".")
        parser.add_argument('--minhashtree', help="local minhash tree")
        parser.add_argument('--minhashblacklist', help="local minhash blacklist")
        args = parser.parse_args()
        return args
def main():
    # Parse command line optiions
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

    # Set up temp directory
    if args.tempdir:
        if not os.path.exists(args.tempdir):
            os.makedirs(args.tempdir)
        dtemp = args.tempdir
    else:
        dtemp = tempfile.mkdtemp()
    logging.info("temporary directory: {}".format(dtemp))
    # Ceate list of minhash files
    reflist = glob.glob(os.path.join(args.minhashrefs, "*.sketch"))
    if reflist:
        refs = ",".join(reflist)
    # Set paths for temporary output files
    segments = os.path.join(dtemp, "segments.fasta")
    minhashout = os.path.join(dtemp,"minhashout.txt")
    kmerout = os.path.join(dtemp,"kmerout.csv")
    codonout = os.path.join(dtemp, "codonout.csv")
    ksize = int(args.ksize)
    logging.info("Set kmer size to {}".format(ksize))
    # Shred gneomes into contigs
    if args.shred:
        logging.info("Shredding sequence data into contigs")
        inhandle = Fasta(args.input, read_ahead=10000)
        with open( segments, 'w') as outhandle:
            samples_all = vica.shred.shred_all(inhandle=inhandle,
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
        logging.info("Processing the input sequences without shredding")

    # Extract minhash features
    s1 = time.perf_counter()
    if args.minhashlocal:
        logging.info("Extacting minhash signatures and identifiing them locally")

        vica.minhash.minhashlocal(dtemp=dtemp, infile=segments, outfile=minhashout,
                     nodesfile=args.nodes, localsketchdir=args.minhashrefs,
                     blacklist=args.minhashblacklist, tree=args.minhashtree)
    else:
        logging.info("Extacting minhash signatures and sending them to a server for identification")
        vica.minhash.minhashremote(dtemp=dtemp, infile=segments, outfile=minhashout,
                     nodesfile=args.nodes)
    s2 = time.perf_counter()
    t1 = s2-s1
    timestring1 =  str(datetime.timedelta(seconds=t1))
    logging.info("Processed Minhash features in: {}".format(timestring1))


    # Extract kmer features
    logging.info("Calculating Kmer features")
    s5 = time.perf_counter()
    kmers = vica.khmer_features.iterate_kmer(k=ksize)
    vica.khmer_features.write_kmers_as_csv(infile=segments, outfile=kmerout,
                              ksize=ksize, kmers=kmers)
    s6 = time.perf_counter()
    t3 = s6 - s5
    timestring3 =  str(datetime.timedelta(seconds=t3))
    logging.info("Processed Kmer features in: {}".format(timestring3))

    # Extract codons
    logging.info("Calculating Codon features")
    s3 = time.perf_counter()
    vica.prodigal.contigs_to_feature_file(infile=segments, outfile=codonout, dtemp=dtemp)
    s4 = time.perf_counter()
    t2 = s4 - s3
    timestring2 =  str(datetime.timedelta(seconds=t2))
    logging.info("Processed Codon features in: {}".format(timestring2))
    #shutil.copytree(dtemp, "dtempout")

    # Combine data into a Tensorflow TF record file
    logging.info("Writing data to the TFrecord file {}".format(args.output))
    s7 = time.perf_counter()
    vica.tfrecord_maker.convert_to_tfrecords(dtemp=dtemp, kmerfile=kmerout, codonfile=codonout,
             minhashfile=minhashout, tfrecordfile=args.output,
             label=str(args.label), sort=True)
    s8 = time.perf_counter()
    t4 = s8-s7
    timestring4 =  str(datetime.timedelta(seconds=t4))
    logging.info("Wrote TFrecord file in: {}".format( timestring4))
    tfinal = time.perf_counter()
    ttot = str(datetime.timedelta(seconds=(tfinal - t0)))
    logging.info("All features processed in: {}".format(ttot) )
    if not args.tempdir:
        shutil.rmtree(dtemp)

if __name__ == '__main__':
    main()
