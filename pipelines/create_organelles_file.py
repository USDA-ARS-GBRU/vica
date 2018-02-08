#!/usr/bin/env python3

import argparse
import subprocess
import tempfile
import logging
import json
import os
import gzip

## Variables
MITO = 'ftp://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/*.1.genomic.fna.gz'
PLASTID = 'ftp://ftp.ncbi.nlm.nih.gov/refseq/release/plastid/*.1.genomic.fna.gz'


def parser():
    parser = argparse.ArgumentParser(description='Create organelle file')
    parser.add_argument('--out', help="location to write file containing organelles")
    args = parser.parse_args()
    return args

def find_organelles():
    """ There are approximately twice as many organelle only genomes as complete
    eukaryotic genomes in refseq. we need to exclude these when building
    the training set so we do not oversample them.
    """
    logging.info("Downloading organelle sequences from Refseq")
    seqidlist = []
    dtemp = tempfile.mkdtemp()
    for item in [MITO, PLASTID]:
        options = ["wget", item, "-P", dtemp,]
        wgetout = subprocess.run(options, stderr=subprocess.PIPE)
    for seqfile in os.listdir(dtemp):
        with gzip.open(os.path.join(dtemp, seqfile), 'rt') as f:
            for line in f:
                if line.startswith(">"):
                    ll = line.split(" ")[0]
                    seqid = ll[1:]
                    seqidlist.append(seqid)
    logging.info("{} organelle sequence accessions saved".format(len(seqidlist)))
    print("{} organelle sequence accessions saved".format(len(seqidlist)))
    #shutil.rmtree(dtemp)
    return seqidlist

def main():
    args = parser()
    seqidlist = find_organelles()
    with open(args.out, 'w') as f:
        json.dump(seqidlist, f, ensure_ascii=False)

if __name__ == '__main__':
    main()
