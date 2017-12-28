import vica.tfrecord_maker
import tempfile
import os
import filecmp
import yaml
from nose.tools import ok_, eq_
import difflib
import sys

def test_external_sort():
    td = tempfile.mkdtemp()
    outfile1 = os.path.join(td,"sortout.txt")
    vica.tfrecord_maker.external_sort(infile="tests/test-data/sort_test_in.csv",
        outfile=outfile1,
        sep = ",",
        key = 1)
    ok_(filecmp.cmp("tests/test-data/sort_test_out.csv", outfile1, shallow=False))

def test_join():
    td = tempfile.mkdtemp()
    outfile1 = os.path.join(td, "mergefile.csv")
    vica.tfrecord_maker.join(kmerfile="tests/test-data/join_1.csv",
                             codonfile="tests/test-data/join_2.csv",
                             minhashfile="tests/test-data/join_3.csv",
                             dtemp=td)
    ok_(filecmp.cmp("tests/test-data/join_out.csv", outfile1, shallow=False))

def test_count_features():
    outdict = vica.tfrecord_maker.count_features(j1="tests/test-data/join_1.csv",
                                       j2="tests/test-data/join_2.csv",
                                       j3="tests/test-data/join_3.csv")
    eq_(outdict, {"j1":2, "j2":3, "j3":3})

# def test_csv_to_tfrecords():
#     vica.tfrecord_maker._csv_to_tfrecords(kmerfile=,
#         codonfile=,
#         minhashfile=,
#         mergefile=
#         tfrecordfile=)
