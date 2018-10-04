import tempfile
import os
import filecmp
import random
import shutil
import yaml
import collections
import glob

from Bio import SeqIO


import pandas
import numpy
import pyfaidx
from nose.tools import ok_, eq_
from ete3 import NCBITaxa

import vica.split_shred

with open(vica.CONFIG_PATH) as cf:
    config = yaml.safe_load(cf)

classes = {2: 1000,
         2157: 1000,
         2759:1000,
         10239: 1000}

def test_split_init():
    data = vica.split_shred.Split(
        fasta_file='tests/test-data/test_contigs.fasta',
            split_depth='order',
        classes=classes,
        testfrac=0.2)

def test_find_organelles():
    data = vica.split_shred.Split(
        fasta_file='tests/test-data/test_contigs.fasta',
        split_depth='order',
        classes=classes,
        testfrac=0.2)
    seqidlist = data._find_organelles()
    ok_(len(seqidlist) > 5000)


def test_test_or_train():
    data = vica.split_shred.Split(
        fasta_file='tests/test-data/test_contigs.fasta',
        split_depth='genus',
        classes=classes,
        testfrac=0.2)
    subtree = data.tax_instance.get_topology([10239])
    test_subtrees, train_subtrees = data._test_or_train(subtree)
    ok_(len(test_subtrees)>0)
    ok_(len(train_subtrees)>0)

def test_assign_samples_attribute():
    data = vica.split_shred.Split(
        fasta_file='tests/test-data/test_contigs.fasta',
        split_depth='order',
        classes=classes,
        testfrac=0.2)
    subtree = data.tax_instance.get_topology([10239])
    data._assign_samples_attribute(n=1000, depth='order', nodelist=[subtree])
    eq_(subtree.samples,1000)

def test_add_samples_feature_to_test_train_nodes():
    data = vica.split_shred.Split(
        fasta_file='tests/test-data/test_contigs.fasta',
        split_depth='order',
        classes=classes,
        testfrac=0.2)
    subtree = data.tax_instance.get_topology([2])
    test = [subtree&str(57723), subtree&str(1783272)]
    train = [subtree&str(1783270), subtree&str(1224)]
    n = 1000
    data._add_samples_feature_to_test_train_nodes(n, test, train)
    print(test[0].samples)
    print(train[0].samples)
    print(test[1].samples)
    print(train[1].samples)
    eq_(test[0].samples, 192)
    eq_(test[1].samples, 8)
    eq_(train[0].samples, 34)
    eq_(train[1].samples, 766)

def test_add_samples_to_children():
    data = vica.split_shred.Split(
        fasta_file='tests/test-data/test_contigs.fasta',
        split_depth='order',
        classes=classes,
        testfrac=0.2)
    subtree = data.tax_instance.get_topology([97050])
    subtree.add_features(samples = 10000)
    data._add_samples_feature_to_children(subtree)
    children = subtree.get_children()
    for child in children:
        print(child.samples)
        ok_(child.samples > 11)
        ok_(child.samples < 30)

def test_propagate_samples_feature_from_nodes_to_leaves():
    data = vica.split_shred.Split(
        fasta_file='tests/test-data/test_contigs.fasta',
        split_depth='genus',
        classes=classes,
        testfrac=0.2)
    subtree = data.tax_instance.get_topology([31989])
    subtree.add_features(samples = 10000000)
    data._propagate_samples_feature_from_nodes_to_leaves(subtree)
    leaves = subtree.get_leaves()
    for leaf in leaves:
        ok_(leaf.samples > 0)

def test_calculate_tax_composition():
    data = vica.split_shred.Split(
        fasta_file='tests/test-data/test_contigs.fasta',
        split_depth='order',
        classes=classes,
        testfrac=0.2)
    subtree = data.tax_instance.get_topology([31989])
    children = subtree.get_children()
    compdict = data._calculate_tax_composition(children)
    # expected = collections.Counter({'genus': 128, 'species': 39, 'no rank': 5})
    # ok_(compdict, expected)

def test_split_test_train_nodes():
    data = vica.split_shred.Split(
        fasta_file='tests/test-data/test_contigs.fasta',
        split_depth='order', # or less?
        classes=classes,
        testfrac=0.5)
    data.split_test_train_nodes()
    for key in data.test_subtrees:
        for node in data.test_subtrees[key]:
            ok_(node.samples >= 0)
    # for node in data.train_leaves:
    #     ok_(node.samples >= 0)

def test_writeseq():
    dtemp = tempfile.mkdtemp()
    data = vica.split_shred.Split(
        fasta_file='tests/test-data/test_contigs.fasta',
        split_depth='order',
        classes=classes,
        testfrac=0.5)
    with open(os.path.join(dtemp,"test.fasta"), 'w') as outfile:
        for record in data.pyfaidx_obj:
            data._writeseq(record, 1, 60, outfile)
        SeqIO.parse(os.path.join(dtemp,"test.fasta"), 'fasta')
    shutil.rmtree(dtemp)

# def test_select_random_segment():
#     data = vica.split_shred.Split(
#         fasta_file='tests/test-data/test_contigs.fasta',
#         split_depth=2,
#         classes=classes,
#         testfrac=0.5)
#     self.

def test_select_random_segment_and_write():
    dtemp = tempfile.mkdtemp()
    data = vica.split_shred.Split(
        fasta_file='tests/test-data/test_contigs.fasta',
        split_depth='order',
        classes=classes,
        testfrac=0.5)
    data.split_test_train_nodes()
    data._select_fragments_and_write(dtemp, seq_length=5000, test=False)
    data._select_fragments_and_write(dtemp, seq_length=5000, test=True)
    for filename in glob.glob(os.path.join(dtemp, "*/*.fasta")):
        ok_(os.path.getsize(filename) > 0)
    shutil.rmtree(dtemp)

def test_write_sequence_data():
    dtemp = tempfile.mkdtemp()
    data = vica.split_shred.Split(
        fasta_file='tests/test-data/test_contigs.fasta',
        split_depth='order',
        classes=classes,
        testfrac=0.5)
    data.split_test_train_nodes()
    data.write_sequence_data(dtemp, overwrite=True, seq_length=5000)
    for filename in glob.glob(os.path.join(dtemp, "*/*.fasta")):
        ok_(os.path.getsize(filename) > 0)
    shutil.rmtree(dtemp)
