#!/usr/bin/env python3
"""prodigal.py: a module to call genes with prodigal then count codon usage
and transform into centered log ratio returning values as a CSV"""

import argparse
import tempfile
import os
import glob
import shutil
import logging
import time
import datetime

import pyfaidx
import yaml

import vica



def run(input, output, label, tempdir, minhashlocal=None, ksize=5, configpath=vica.CONFIG_PATH):
    """A command to run all the steps in thefeature selection selection workflow
    1. the selection of minhash features
    2. the selection of codon usage features
    3. the selection of kmer features
    4. writing features to a tfrecord file


    """
    # Begin timing opperation
    t0 = time.perf_counter()
    # read config
    with open(configpath) as cf:
        config = yaml.load(cf)
    # Set up temp directory
    if config["get_features"]["tempdir"]:
        if not os.path.exists():
            os.makedirs(config["get_features"]["tempdir"])
        dtemp = tempdir
    else:
        dtemp = tempfile.mkdtemp()
    logging.info("temporary directory: {}".format(dtemp))


    # Set paths for temporary output files
    segments = os.path.join(dtemp, "segments.fasta")
    minhashout = os.path.join(dtemp,"minhashout.txt")
    kmerout = os.path.join(dtemp,"kmerout.csv")
    codonout = os.path.join(dtemp, "codonout.csv")
    logging.info("Set kmer size to {}".format(ksize))
    # Shred gneomes into contigs

    # Extract minhash features
    s1 = time.perf_counter()
    if minhashlocal:
        logging.info("Extacting minhash signatures and identifing them locally")
        try:
            # Create list of minhash files
            reflist = glob.glob(os.path.join(config["minhash"]["refs"], "*.sketch"))
            if reflist:
                refs = ",".join(reflist)
            vica.minhash.minhashlocal(dtemp=dtemp,
                         infile=input,
                         outfile=minhashout,
                         nodesfile=config["minhash"]["nodesfile"],
                         refs=reflist,
                         blacklist=config["minhash"]["blacklist"],
                         tree=config["minhash"]["tree"],
                         taxfilter= config["get_features"]["taxfilter"],
                         taxfilterlevel=config["get_features"]["taxfilterlevel"],
                         configpath=configpath)
        except:
            logging.exception("vica get_features: during minhash feature selection the following exception occcured:")
            raise SystemExit(1)
    else:
        logging.info("Extacting minhash signatures and sending them to a server for identification")
        vica.minhash.minhashremote(dtemp=dtemp, infile=segments,
            outfile=minhashout,configpath=configpath)
    s2 = time.perf_counter()
    t1 = s2-s1
    timestring1 =  str(datetime.timedelta(seconds=t1))
    logging.info("Processed Minhash features in: {}".format(timestring1))


    # Extract kmer features
    logging.info("Calculating Kmer features")
    s5 = time.perf_counter()
    vica.khmer_features.run(infile=segments, outfile=kmerout, configpath=configpath)

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
    logging.info("Writing data to the TFrecord file {}".format(output))
    s7 = time.perf_counter()
    vica.tfrecord_maker.convert_to_tfrecords(dtemp=dtemp, kmerfile=kmerout, codonfile=codonout,
             minhashfile=minhashout, tfrecordfile=output,
             label=str(label), sort=True)
    s8 = time.perf_counter()
    t4 = s8-s7
    timestring4 =  str(datetime.timedelta(seconds=t4))
    logging.info("Wrote TFrecord file in: {}".format( timestring4))
    tfinal = time.perf_counter()
    ttot = str(datetime.timedelta(seconds=(tfinal - t0)))
    logging.info("All features processed in: {}".format(ttot) )
    if not tempdir:
        shutil.rmtree(dtemp)
