#!/usr/bin/env python3
"""vica_cli.py a command line interface to Vica: A classifier for indentifing
    highly divergent DNA and RNA viruses
    in metagenomic and metatranscriptomic contigs
    """

import argparse
import sys
import logging
import ast

import yaml

import vica

# Set config path for modules in the event that they are accessed directly and
# not from the vica-cli

def config_logging(logfile, level=logging.DEBUG):
    """Set up logging
    Set up logging for vica_cli.py

    Args:
        logfile (str): the path to the log file
        level (str): the level of logging to record

    Returns:
        None

    """
    logging.basicConfig(filename=logfile,
                        level=level,
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
    logging.info("Starting Vica: software for the identification of highly divergent DNA and RNA viruses")


def parser():
    """Parse the input arguments for the classifier and return an args object"""
    parser = argparse.ArgumentParser(
        description="Vica: A classifier for identifying highly divergent \
        DNA and RNA viruses in metagenomic and metatranscriptomic contigs")
    subparsers = parser.add_subparsers(help ="Commands", dest= 'command')
    classify = subparsers.add_parser(
        'classify',help="Predict viral contigs from fasta or tfrecord files.")
    split = subparsers.add_parser(
        'split',help="Split fasta data into testing and training data \
        subsampled and fragmented to a specified length.")
    get_features = subparsers.add_parser(
        'get_features', help="Extract features from a fasta file of contigs.")
    train = subparsers.add_parser(
        'train', help="Train a classifier model from tfrecord data.")
    evaluate = subparsers.add_parser(
        'evaluate', help="Evaluate the classification performance of a \
        model from test data.")
    # vica classify subparser
    classify.set_defaults(whichmethod='classify')
    classify.add_argument(
         '-i','--infile', help="A fasta (.fasta or .fa) file  or tfrecord file \
        (.tfrecord) of sequences or features to be classified.", required=True)
    classify.add_argument(
        '--out', help="A text file of class predictions for each sequence or \
        example in the tfrecord file.", required=True)
    classify.add_argument(
        '--threshold', help="A probability value to call a contig viral.",
        type=float, default=0.5)
    classify.add_argument(
        '--modeldir', help="A directory containing the Tensorflow model and \
        evaluation data for analysis. If training has been done in this directory \
        previously the results will be added", required=True)
    classify.add_argument(
        '--n_classes', help="the number of classes present in the training \
        data, default 4", default=4)
    classify.add_argument(
        '--logfile',help="A file to record the analysis. If the same log file \
        is given for multiple vica commands all the steps in the workflow will \
        be recorded.", default="vica.log")
    classify.add_argument(
        '--config', help="path to YAML formatted configuration file, default \
        is " + vica.CONFIG_PATH,
        default=vica.CONFIG_PATH)
    # vica split subparser
    split.set_defaults(whichmethod='split')
    split.add_argument(
         '-i','--infile', help="A training fasta processed into BBtools taxonomy format \
        using the scripts in the 'pipelines' directory of vica.", required=True)
    split.add_argument(
        '--out', help="A directory to create containing the testing and \
        training split at the requested taxonomic level for each class.", required=True)
    split.add_argument(
        '--length', help="The length of training fragments.", type=int,
        default=5000)
    split.add_argument(
        '--testfrac', help="The proportion of data to put into the test group.",
        type=float, default=0.1)
    split.add_argument(
        '--split_depth', help="The depth above the leaves at which the testing \
        and training data will be split at. The final distribution of taxonomic \
        ranks sampled will be displayed in the log file.",
        default=7, type=int)
    split.add_argument(
        '--classes', help="The classes to separate data into. This should be a \
        dictionary of NCBI taxonomy identifiers and the number of samples from each class.",
        default='{2: 100000, 2157: 100000, 2759: 100000, 10239: 100000}')
    split.add_argument(
        '--logfile',help="A file to record the analysis. If the same log file \
        is given for multiple Vica commands all the steps in the workflow will \
        be recorded.", default="vica.log")
    split.add_argument(
        '--config', help="path to YAML formatted configuration file, default \
        is " + vica.CONFIG_PATH,
        default=vica.CONFIG_PATH)
    # vica get_features subparser
    get_features.set_defaults(whichmethod='get_features')
    get_features.add_argument(
         '-i','--infile', help="A fasta (.fasta or .fa) file with sequences that have \
        been wrapped to a consistent line length", required=True)
    get_features.add_argument(
        '--out', help="An uncompressed tfrecord file containing the name of \
        the sequence, the class label, and vectors for minhash, codon and \
        5-mer features", required=True)
    get_features.add_argument(
        '--label', help="An integer label for the classification class of \
        training or evaluation data. Needed to for training and test data. \
        for data to be classified use -1.",
        type=int, required=True)
    get_features.add_argument(
        '--minhashlocal', help="A flag to use a local version of a minhash \
        database rather than a remote server. Default is false. relevant minhash \
        paths in the config file should be updated before using --minhashlocal.",
        action="store_true")
    get_features.add_argument(
        '--logfile',help="A file to record the analysis. If the same log file \
        is given for multiple Vica commands all the steps in the workflow will \
        be recorded.", default="vica.log")
    get_features.add_argument(
        '--config', help="path to YAML formatted configuration file, default \
        is " + vica.CONFIG_PATH,
        default=vica.CONFIG_PATH)
    # vica train subparser
    train.set_defaults(whichmethod='train')
    train.add_argument(
         '-i','--infile', help="One or more tfrecord files containing training instances",
         required=True, nargs="+")
    train.add_argument(
        '--out', help="A file containing the Tensorflow model of use in \
        future classification.", required=True)
    train.add_argument(
        '--modeldir', help="A directory containing the Tensorflow model and \
        evaluation data for analysis. If training has been done in this \
        directory previously, the results will be added.", required=True)
    train.add_argument(
        '--n_classes', help="the number of classes present in the training \
        data, default 4", default=4)
    train.add_argument(
        '--logfile',help="A file to record the analysis. If the same log file \
        is given for multiple Vica commands all the steps in the workflow will \
        be recorded.", default="vica.log")
    train.add_argument(
        '--config', help="path to YAML formatted configuration file, default \
        is " + vica.CONFIG_PATH,
        default=vica.CONFIG_PATH)
    # vica evaluate subparser
    evaluate.set_defaults(whichmethod='evaluate')
    evaluate.add_argument(
         '-i','--infile', help="One or more tfrecord files containing testing instances.",
         required=True, nargs="+")
    evaluate.add_argument(
        '--out', help="A directory the predictions and analysis plots  \
        summarizing model performance.", default="modelresults")
    evaluate.add_argument(
        '--modeldir', help="A directory containing the Tensorflow model and \
        evaluation data for analysis. If training has been done in this directory \
        previously the results will be added", required=True)
    evaluate.add_argument(
        '--n_classes', help="the number of classes present in the training \
        data, default 4", default=4)
    evaluate.add_argument(
        '--logfile',help="A file to record the analysis. If the same log file \
        is given for multiple Vica commands all the steps in the workflow will \
        be recorded.", default="vica.log")
    evaluate.add_argument(
        '--config', help="path to YAML formatted configuration file, default \
        is " + vica.CONFIG_PATH,
        default=vica.CONFIG_PATH)
    if len(sys.argv)==1:
        parser.print_help()
        raise SystemExit(1)
    args = parser.parse_args()
    return args

def main():
    """The main entry point for using Vica. it will parse arguments,
    set up logging and run selected subprogram
    """
    try:
        args = parser()
    except:
        print("Could not parse command line arguments.")
        raise SystemExit(1)

    try:
        with open(args.config) as cf:
            config = yaml.safe_load(cf)
    except:
        print("Could not parse the configuration file.")
        raise SystemExit(1)

    try:
        config_logging(args.logfile)
        logging.info("Configuration data loaded from {}:".format(args.config))
        logging.info(config)
    except:
        print("Could not set up logging, exiting.")
        raise SystemExit(1)
    try:
        if args.whichmethod == 'classify':
            vica.train_eval.classify(infile=args.infile,
                out= args.out,
                modeldir= args.modeldir,
                n_classes=args.n_classes,
                configpath= args.config)
        elif args.whichmethod == 'split':
            vica.split_shred.run(infile=args.infile,
                outdir=args.out,
                length=args.length,
                testfrac=args.testfrac,
                split_depth=args.split_depth,
                classes=ast.literal_eval(args.classes))
        elif args.whichmethod == 'get_features':
            vica.get_features.run(infile=args.infile,
                output=args.out,
                label= args.label,
                minhashlocal=args.minhashlocal,
                configpath=args.config)
        elif args.whichmethod == 'train':
            vica.train_eval.train(infiles=args.infile,
                out= args.out,
                modeldir= args.modeldir,
                n_classes= args.n_classes,
                configpath= args.config)
        elif args.whichmethod == 'evaluate':
            vica.train_eval.evaluate(infiles=args.infile,
                out= args.out,
                modeldir= args.modeldir,
                n_classes= args.n_classes,
                configpath= args.config)
    except:
        logging.exception("vica_cli.py: The following exception occurred:")
        raise SystemExit(1)



if __name__ == '__main__':
    main()
