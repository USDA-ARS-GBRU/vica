import tempfile
import os
import filecmp
import csv
import shutil

import numpy as np
import nose
import yaml
from Bio import SeqIO


import vica

with open(vica.CONFIG_PATH) as cf:
    config = yaml.safe_load(cf)

codon_list = config["prodigal"]["codon_list"]


def test_call_genes():
    td = tempfile.mkdtemp()
    outfile1 = os.path.join(td,"genes.fasta")
    trans = os.path.join(td, "translations.faa")
    vica.prodigal._call_genes(infile="tests/test-data/2testseqs.fasta",
                         outfile=outfile1, translations=trans)
    nose.tools.ok_(filecmp.cmp("tests/test-data/calledgenes.fasta", outfile1, shallow=False))
    shutil.rmtree(td)

def test_gene_to_codon():
    genestring = 'ATCGATCGCATCGACTAGCCACGTACGACGACTGACGTCG'
    genedict = ['ATC', 'GAT', 'CGC', 'ATC', 'GAC', 'TAG', 'CCA', 'CGT', 'ACG',
                'ACG', 'ACT', 'GAC', 'GTC']
    result = vica.prodigal._gene_to_codon(genestring)
    nose.tools.eq_(result, genedict)

def test_codon_to_dict():
    genestring = 'ATCGATCGCATCGACTAGCCACGTACGACGACTGACGTCG'
    d0 = {'CGT': 1, 'GAT': 1, 'ATC': 2, 'ACT': 1, 'ACG': 2, 'GAC': 2, 'GTC': 1, 'CGC': 1, 'CCA': 1, 'TAG': 1}
    d1 = {'ACG': 1, 'CAC': 1, 'GTA': 1, 'GCA': 1, 'AGC': 1, 'TCG': 3, 'CTG': 1, 'ACT': 1, 'ATC': 1, 'CGA': 2}
    d2 = {'TCG': 1, 'CGA': 2, 'CGT': 1, 'GCC': 1, 'CTA': 1, 'GAC': 2, 'TGA': 1, 'CAT': 1, 'ACG': 1, 'TAC': 1}
    r0 = vica.prodigal._codon_to_dict(genestring, 0)
    r1 = vica.prodigal._codon_to_dict(genestring, 1)
    r2 = vica.prodigal._codon_to_dict(genestring, 2)
    nose.tools.eq_(r0, d0)
    nose.tools.eq_(r1, d1)
    nose.tools.eq_(r2, d2, msg="{} did not match expected value of {}".format(r2, d2))


def _parse_prodigal_id_from_biopython():
    stripped_id = '>NC_000'
    result = vica.prodigal._parse_prodigal_id_from_biopython('>NC_000_001')
    nose.tools.eq_(result, stripped_id)

def almost_equal(value_1, value_2, accuracy = 10**-8):
    '''a Function to compare tuples of float values'''
    return abs(value_1 - value_2) < accuracy

def test_count_dict_to_ilr_array():
    d0 = {'CGT': 1, 'GAT': 1, 'ATC': 2, 'ACT': 1, 'ACG': 2, 'GAC': 2, 'GTC': 1,
          'CGC': 1, 'CCA': 1, 'TAG': 1}
    expected = np.array([-8.82897217e-16, -2.70205994e-17,  4.36278823e-16, -2.60685790e-15,
        1.36645246e-15,  7.79053701e-16, -4.46755806e-15, -2.72346331e-15,
        1.42899667e-15, -6.07576770e-16,  3.48032867e-15, -4.90782530e-15,
       -8.82516983e-15, -5.41953679e+00,  3.62108033e-01,  3.40141269e-01,
        3.20688264e-01,  3.03340730e-01, -5.17993711e+00,  5.47455887e-01,
        5.21978708e-01,  4.98767926e-01, -5.01410629e+00,  6.87051686e-01,
       -4.84071799e+00,  8.46906011e-01,  8.16099053e-01,  7.87455022e-01,
        7.60753826e-01,  7.35804259e-01,  7.12439410e-01,  6.90512925e-01,
       -4.85674553e+00,  8.13093153e-01,  7.90184510e-01, -5.44861146e+00,
        9.16125926e-01,  8.92326380e-01,  8.69732171e-01,  8.48253988e-01,
        8.27811128e-01,  8.08330491e-01,  7.89745696e-01,  7.71996335e-01,
       -5.47898717e+00,  8.74342355e-01,  8.55933085e-01,  8.38283085e-01,
        8.21346327e-01, -4.74940324e+00,  8.98378702e-01, -4.67531415e+00,
        9.69651246e-01, -5.29348079e+00,  1.04827636e+00,  1.02972133e+00,
        1.01181176e+00,  9.94514562e-01,  9.77798841e-01])

    result = vica.prodigal.count_dict_to_ilr_array(d0, codon_list)
    nose.tools.ok_(all(almost_equal(*values) for values in zip(expected, result)))

def test_dsum():
    d0 = {'CGT': 1, 'GAT': 1, 'ATC': 2, 'ACT': 1, 'ACG': 2, 'GAC': 2, 'GTC': 1,
          'CGC': 1, 'CCA': 1, 'TAG': 1}
    d1 = {'ACG': 1, 'CAC': 1, 'GTA': 1, 'GCA': 1, 'AGC': 1, 'TCG': 3, 'CTG': 1,
          'ACT': 1, 'ATC': 1, 'CGA': 2}
    expected = {'CGT': 1, 'GTA': 1, 'CCA': 1, 'ACT': 2, 'AGC': 1, 'TAG': 1,
                'CGA': 2, 'TCG': 3, 'GCA': 1, 'CGC': 1, 'ATC': 3, 'GAT': 1,
                'GTC': 1, 'CAC': 1, 'ACG': 3, 'GAC': 2, 'CTG': 1}
    result = vica.prodigal.dsum(d0, d1)
    nose.tools.eq_(result, expected)

def test_count_codon_in_gene():
    testfasta = 'tests/test-data/2testseqs.fasta'
    expected = {0: {'GAC': 53, 'TAC': 11, 'AGA': 18, 'CAG': 56, 'CGT': 29,
                    'ACC': 49, 'ATG': 46, 'GCC': 93, 'AAC': 33, 'AAA': 34,
                    'ATT': 28, 'GCA': 46, 'TAT': 22, 'GAG': 41, 'TCC': 33,
                    'GCG': 95, 'GAA': 49, 'GTG': 51, 'CAC': 17, 'TCG': 79,
                    'CGG': 51, 'ACG': 38, 'TTC': 42, 'TGG': 27, 'TGT': 14,
                    'ACT': 14, 'TGA': 18, 'CTG': 55, 'GCT': 30, 'AAG': 34,
                    'GTT': 31, 'CCA': 34, 'CAT': 29, 'CAA': 27, 'CTA': 8,
                    'TTT': 23, 'GGA': 36, 'GTA': 10, 'GGG': 53, 'TCA': 25,
                    'GAT': 48, 'TGC': 36, 'CTC': 22, 'ACA': 18, 'CCT': 16,
                    'CTT': 25, 'CCC': 41, 'CCG': 75, 'GGC': 63, 'TAA': 5,
                    'GGT': 43, 'CGA': 38, 'TCT': 25, 'AGG': 36, 'CGC': 78,
                    'ATC': 85, 'AAT': 19, 'GTC': 44, 'TTG': 44, 'ATA': 9,
                    'AGT': 13, 'AGC': 51, 'TAG': 13, 'TTA': 4},
                1: {'CTT': 21, 'GAC': 31, 'TAC': 8, 'AGA': 45, 'CAG': 35,
                    'CGT': 35, 'ACC': 33, 'ATG': 44, 'GCC': 54, 'AAC': 38,
                    'CGA': 63, 'ATT': 17, 'GCA': 51, 'TAT': 8, 'GAG': 16,
                    'TCC': 35, 'GCG': 94, 'GAA': 36, 'TGA': 49, 'ACT': 17,
                    'CAC': 22, 'TCG': 67, 'CGG': 93, 'ACA': 23, 'TTC': 22,
                    'TGG': 64, 'TCA': 65, 'AGC': 36, 'GTG': 35, 'CTG': 43,
                    'CCC': 34, 'AAG': 29, 'GTT': 29, 'CCA': 57, 'CAT': 29,
                    'AAT': 19, 'CTA': 8, 'TTT': 27, 'GGA': 28, 'GTA': 12,
                    'GGG': 40, 'GCT': 29, 'GAT': 27, 'TGC': 50, 'CCT': 57,
                    'CTC': 13, 'CCG': 68, 'GGC': 59, 'TAA': 7, 'GGT': 40,
                    'AAA': 29, 'TCT': 26, 'AGG': 56, 'CGC': 96, 'ATC': 47,
                    'ACG': 41, 'GTC': 23, 'TGT': 33, 'TTG': 53, 'ATA': 10,
                    'AGT': 7, 'TAG': 8, 'CAA': 37, 'TTA': 5},
                2: {'GAC': 42, 'AGA': 22, 'CAG': 67, 'CGC': 68, 'CAA': 41,
                    'ATG': 16, 'GCC': 80, 'TGT': 34, 'CGA': 85, 'ATT': 15,
                    'GCA': 39, 'TAT': 13, 'GAG': 20, 'TCC': 23, 'GCG': 80,
                    'GAA': 38, 'TGA': 44, 'CAC': 39, 'TCG': 33, 'CGG': 80,
                    'ACA': 25, 'TTC': 45, 'TGG': 35, 'TCA': 28, 'ACT': 16,
                    'GTG': 17, 'CTG': 41, 'CCC': 35, 'TAA': 5, 'GTT': 38,
                    'CCA': 37, 'CAT': 49, 'ACC': 28, 'CTA': 14, 'TTT': 22,
                    'GGA': 40, 'GTA': 14, 'AGG': 25, 'GGG': 55, 'AGT': 6,
                    'GAT': 85, 'TGC': 62, 'CTC': 36, 'CCT': 31, 'CTT': 38,
                    'AAC': 31, 'CCG': 53, 'GGC': 99, 'GGT': 58, 'AAA': 36,
                    'TCT': 21, 'TAC': 7, 'CGT': 37, 'ATC': 35, 'AAT': 21,
                    'TTG': 21, 'GTC': 46, 'GCT': 42, 'ATA': 17, 'AGC': 35,
                    'TAG': 10, 'AAG': 21, 'TTA': 6, 'ACG': 30}}
    records = SeqIO.parse(testfasta, 'fasta')
    record = next(records)
    result = vica.prodigal.count_codon_in_gene(record)
    nose.tools.eq_(result, expected)

def test_count_codons():
    dt  = tempfile.mkdtemp()
    print(dt)
    outfile = os.path.join(dt, "codons.csv")
    expected =  'tests/test-data/codons.csv'
    genefile = 'tests/test-data/calledgenes.fasta'
    seqs = SeqIO.parse(genefile, 'fasta')
    with open(outfile, 'w') as csvfile:
        csv_writer_instance = csv.writer(csvfile, lineterminator='\n')
        vica.prodigal.count_codons(seqio_iterator= seqs, csv_writer_instance=csv_writer_instance, codon_list=codon_list)
    # nose.tools.ok_(filecmp.cmp(expected, outfile, shallow=False))
    shutil.rmtree(dt)


def test_contigs_to_feature_file():
    infile = "tests/test-data/2testseqs.fasta"
    dt = tempfile.mkdtemp()
    outfile = os.path.join(dt, "codons.csv")
    trans = os.path.join(dt, "translations.faa")
    expected =  'tests/test-data/codons.csv'
    vica.prodigal.contigs_to_feature_file(infile=infile, outfile=outfile, translations=trans, dtemp=dt, codon_list=codon_list)
    nose.tools.ok_(filecmp.cmp(expected, outfile, shallow=False))
    shutil.rmtree(dt)
