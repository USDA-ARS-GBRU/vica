#!/usr/bin/env python3
"""Vica: A classifier for indentifing highly divergent DNA and RNA viruses
   in metagenomic and metatranscriptomic contigs"""

import argparse
import os
import pkg_resources

def parser():
    """Parse the input arguments for the classifier and return an args object"""
    parser = argparse.ArgumentParser(
        description="Vica: A classifier for indentifing highly divergent \
        DNA and RNA viruses in metagenomic and metatranscriptomic contigs")
    subparsers = parser.add_subparsers(help ="Commands", dest= 'command')
    classify = subparsers.add_parser(
        'predict',help="Predict viral contigs from fasta or tfrecord files.")
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
        '--in', help="A fasta (.fasta or .fa) file  or tfrecord file \
        (.tfrecord) of sequences or features to be classified.", required=True)
    classify.add_argument(
        '--out', help="A text file of class predictions for each sequence or \
        example in the tfrecord file.", required=True)
    classify.add_argument(
        '--threshold', help="A probability value to call a contig viral.",
        type=float, default=0.5)
    classify.add_argument(
        '--config', help="A config file with additional parameters.",
        default=pkg_resources.resource_filename(__name__, "config.yaml"))
    classify.add_argument(
        '--logfile',help="A file to record the analysis. If the same log file \
        is given for multiple vica commands all the setps in the workflow will \
        be recorded.", default="vica.log")
    # vica get_features subparser
    get_features.set_defaults(whichmethod='get_features')
    get_features.add_argument(
        '--in', help="A fasta (.fasta or .fa) file with sequences that have \
        been wrapped to a consistent line length", required=True)
    get_features.add_argument(
        '--out', help="An uncompressed tfrecord file. containing tne name of \
        the sequence, the class label, and vectors for minhash, codon and \
        5-mer features", required=True)
    get_features.add_argument(
        '--label', help="An interger label for the classifcation class of \
        training or evaluation data. Not needed for data to be classified.",
        type=int)
    get_features.add_argument(
        '--minhashlocal', help="A flag to use a local version of a minhash \
        database rather than a remote server. Default is false.",
        action=store_true)
    get_features.add_argument(
        '--config', help="A config file with additional parameters.",
        default=pkg_resources.resource_filename(__name__, "config.yaml"))
    cget_features.add_argument(
        '--logfile',help="A file to record the analysis. If the same log file \
        is given for multiple vica commands all the setps in the workflow will \
        be recorded.", default="vica.log")
    # vica train subparser
    train.set_defaults(whichmethod='train')
    train.add_argument(
        '--in', help="One or more tfrecord files containing training instances",
         required=True, nargs="+")
    train.add_argument(
        '--out', help="A file contianing the Tensorflow model of use in \
        future classification.", required=True)
    train.add_argument(
        '--modeldir', help="A directory contianing the Tensorflow model and \
        evaluation data for analysis. If training has been done in this directory \
        previously the results will be added", required=True)
        # vica evaluate subparser
    evaluate.add_argument(
        '--in', help="One or more tfrecord files containing testing instances",
         required=True, nargs="+")
    evaluate.add_argument(
        '--out', help="A directory the predictions and analysis plots  \
        summarizing model performance.", default="modelresults")
    evaluate.add_argument(
        '--modeldir', help="A directory contianing the Tensorflow model and \
        evaluation data for analysis. If training has been done in this directory \
        previously the results will be added", required=True)
    args = parser.parse_args()
