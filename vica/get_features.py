"""A module with functions to run all the steps in the feature
    selection workflow for vica"""


import tempfile
import os
import shutil
import logging
import time
import datetime

import yaml

import vica


def run(infile, output, filtertaxa=False, configpath=vica.CONFIG_PATH):
    """Run all the steps in the feature selection selection workflow.

    This command: 1) selects minhash features, 2) codon usage features, 3)
    kmer features and 4) writes features to a tfrecord file.

    Args:
        infile (str): a fasta file with names in the format
            "tid|<NCBI taxonomy ID>|<optional accession>".
            Example: "tid|1026970|NW_008342263.1
        output (str): a name for the TFrecords file to be generated. It should
            end in ".tfrecords".
        filtertaxa (bool: Should minhash features be filtered for taxa to prevent bleedover in the the test set)
        configpath (str): path to the yaml configuration file.

    Returns:
        None


    """
    # read config
    with open(configpath) as cf:
        config = yaml.safe_load(cf)

    # Set up temp directory
    if config["get_features"]["tempdir"]:
        if not os.path.exists(config["get_features"]["tempdir"]):
            os.makedirs(config["get_features"]["tempdir"])
        dtemp = config["get_features"]["tempdir"]
    else:
        dtemp = tempfile.mkdtemp()
    logging.info("temporary directory: {}".format(dtemp))

    # Set paths for temporary output files
    minhashout = os.path.join(dtemp, "minhashout.txt")
    kmerout = os.path.join(dtemp, "kmerout.csv")
    codonout = os.path.join(dtemp, "codonout.csv")
    transout = os.path.join(dtemp, "translations.fasta")
    hmmerout = os.path.join(dtemp, "hmmerout.json")


    # Extract minhash features
    try:
        logging.info("Extacting minhash signatures and sending them to a server for identification")
        vica.minhash.minhashremote(infile=infile,
                                   outfile=minhashout,
                                   server_url=config["minhash"]["server_url"],
                                   filtertaxa=filtertaxa)
    except:
        logging.exception("vica get_features: during minhash remote feature selection the following exception occurred:")
        raise SystemExit(1)

    logging.info("Processed Minhash features")


    # Extract kmer features
    logging.info("Calculating Kmer features")
    try:
        vica.khmer_features.run(infile=infile, outfile=kmerout, ksize=config["khmer_features"]["ksize"])
    except:
        logging.exception("vica get_features: during kmer feature selection the following exception occurred:")
        raise SystemExit(1)
    logging.info("Processed Kmer features")


    # Extract codon features
    try:
        vica.prodigal.contigs_to_feature_file(infile=infile,
                                              outfile=codonout,
                                              translations=transout,
                                              dtemp=dtemp,
                                              codon_list=config["prodigal"]["codon_list"])
    except:
        logging.exception("vica get_features: during codon feature selection the following exception occurred:")
        raise SystemExit(1)
    logging.info("Processed Codon features")

    # Identify proteins
    logging.info("Finding protein homology with Hmmer")
    try:
        vica.hmmer.get_hmmer_features(dtemp=dtemp,
                                     seqfile=transout,
                                     outfile=hmmerout,
                                     hmmfile=os.path.join(vica.DATA_PATH, config["hmmer"]["hmmer_file"]))
    except:
        logging.exception("vica get_features: during HMM feature selection the following exception occurred:")
        raise SystemExit(1)
    logging.info("Processed protein homology features")


    # Combine data into a Tensorflow TF record file
    logging.info("Writing data to the TFrecord file %s", output)
    s7 = time.perf_counter()
    try:
        vica.tfrecord_maker.convert_to_tfrecords(dtemp=dtemp,
                                                 kmerfile=kmerout,
                                                 codonfile=codonout,
                                                 minhashfile=minhashout,
                                                 hmmerfile=hmmerout,
                                                 tfrecordfile=output,
                                                 sort=True)
    except:
        logging.exception("vica get_features: While creating a TFrecord file the following exception occurred:")
        raise SystemExit(1)
    logging.info("Wrote TFrecord file")

    logging.info("All features processed")
    if not config["get_features"]["tempdir"]:
        shutil.rmtree(dtemp)
