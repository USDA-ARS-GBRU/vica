import tempfile
import os
import filecmp
import random
import shutil
import yaml


import pandas
import numpy
import pyfaidx
from nose.tools import ok_, eq_
from ete3 import NCBITaxa

import vica.split_shred

with open(vica.CONFIG_PATH) as cf:
    config = yaml.safe_load(cf)

classes = {2: 100,
         2157: 100,
         2759:100,
         10239: 100}

def test_split_init():
    data = vica.split_shred.Split(
        fasta_file='tests/test-data/test_contigs.fasta',
            split_depth=5,
        classes=classes,
        testfrac=0.2)

def test_find_organelles():
    data = vica.split_shred.Split(
        fasta_file='tests/test-data/test_contigs.fasta',
        split_depth=5,
        classes=classes,
        testfrac=0.2)
    seqidlist = data._find_organelles()
    ok_(len(seqidlist) > 5000)

def test_test_or_train():
    data = vica.split_shred.Split(
        fasta_file='tests/test-data/test_contigs.fasta',
        split_depth=5,
        classes=classes,
        testfrac=0.2)
    ncbi = NCBITaxa()
    subtree = ncbi.get_topology([10239])
    test_subtrees, train_subtrees = data._test_or_train(subtree)
    ok_(len(test_subtrees)>0)
    ok_(len(train_subtrees)>0)

# def _process_samples():
