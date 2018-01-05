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


def run(infile, output, label, minhashlocal=None, configpath=vica.CONFIG_PATH):
    """Run all the steps in the feature selection selection workflow.

    This command: 1) selects minhash features, 2) codon usage features, 3)
    kmer features and 4) writes features to a tfrecord file.

    Args:
        infile (str): a fasta file with names in the format
            "tid|<NCBI taxonomy ID>|<optional accession>".
            Example: "tid|1026970|NW_008342263.1
        output (str): a name for the TFrecords file to be generated. It should
            end in ".tfrecords".
        label (int): This is an integer for the taxonomic class used by the
            classifier.  It should be -1 if the true class is unknown or,
            if the class is known it should  begin with 0 and
            increase sequentially for each training class.
        minhashlocal (bool): If True the function uses a local instance of
            the bbtools minhash database and taxonomy files. Locations of
            the reference files should be set in the configuration file.
        configpath (str): path to the yaml configuration file.

    Returns:
        None


    """
    # Begin timing opperation
    t0 = time.perf_counter()
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
    minhashout = os.path.join(dtemp,"minhashout.txt")
    kmerout = os.path.join(dtemp,"kmerout.csv")
    codonout = os.path.join(dtemp, "codonout.csv")

    # Extract minhash features
    s1 = time.perf_counter()
    if minhashlocal:
        logging.info("Extacting minhash signatures and identifying them locally")
        try:
            vica.minhash.minhashlocal(dtemp=dtemp,
                infile=infile,
                outfile=minhashout,
                ref=config["minhash"]["ref"],
                blacklist=config["minhash"]["blacklist"],
                tree=config["minhash"]["tree"],
                taxfilter=config["minhash"]["taxfilter"],
                taxfilterlevel=config["minhash"]["taxfilterlevel"],
                memory=config["minhash"]["memory"],
                nodesfile=config["minhash"]["nodesfile"],
                noncellular=config["minhash"]["noncellular"])
        except:
            logging.exception("vica get_features: during minhash local feature selection the following exception occurred:")
            raise SystemExit(1)
    else:
        try:
            logging.info("Extacting minhash signatures and sending them to a server for identification")
            vica.minhash.minhashremote(dtemp=dtemp,
                infile=infile,
                outfile=minhashout,
                server_url=config["minhash"]["server_url"],
                nodesfile=config["minhash"]["nodesfile"],
                noncellular=config["minhash"]["noncellular"])
        except:
            logging.exception("vica get_features: during minhash remote feature selection the following exception occurred:")
            raise SystemExit(1)

    s2 = time.perf_counter()
    t1 = s2-s1
    timestring1 =  str(datetime.timedelta(seconds=t1))
    logging.info("Processed Minhash features in: {}".format(timestring1))


    # Extract kmer features
    logging.info("Calculating Kmer features")
    s5 = time.perf_counter()
    try:
        vica.khmer_features.run(infile=infile, outfile=kmerout, ksize=config["khmer_features"]["ksize"])
    except:
        logging.exception("vica get_features: during kmer feature selection the following exception occurred:")
        raise SystemExit(1)

    s6 = time.perf_counter()
    t3 = s6 - s5
    timestring3 =  str(datetime.timedelta(seconds=t3))
    logging.info("Processed Kmer features in: {}".format(timestring3))

    # Extract codons
    logging.info("Calculating Codon features")
    s3 = time.perf_counter()
    try:
        vica.prodigal.contigs_to_feature_file(infile=infile, outfile=codonout, dtemp=dtemp, codon_list=config["prodigal"]["codon_list"])
    except:
        logging.exception("vica get_features: during codon feature selection the following exception occurred:")
        raise SystemExit(1)
    s4 = time.perf_counter()
    t2 = s4 - s3
    timestring2 =  str(datetime.timedelta(seconds=t2))
    logging.info("Processed Codon features in: {}".format(timestring2))
    #shutil.copytree(dtemp, "dtempout")

    # Combine data into a Tensorflow TF record file
    logging.info("Writing data to the TFrecord file {}".format(output))
    s7 = time.perf_counter()
    try:
        vica.tfrecord_maker.convert_to_tfrecords(dtemp=dtemp, kmerfile=kmerout, codonfile=codonout,
                 minhashfile=minhashout, tfrecordfile=output,
                 label=str(label), sort=True)
    except:
        logging.exception("vica get_features: While creating a TFrecord file the following exception occurred:")
        raise SystemExit(1)
    s8 = time.perf_counter()
    t4 = s8-s7
    timestring4 =  str(datetime.timedelta(seconds=t4))
    logging.info("Wrote TFrecord file in: {}".format( timestring4))
    tfinal = time.perf_counter()
    ttot = str(datetime.timedelta(seconds=(tfinal - t0)))
    logging.info("All features processed in: {}".format(ttot) )
    if not config["get_features"]["tempdir"]:
        shutil.rmtree(dtemp)
