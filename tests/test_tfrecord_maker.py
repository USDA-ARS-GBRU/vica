import tempfile
import os
import filecmp

from nose.tools import ok_, eq_

import vica.tfrecord_maker


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

def test_csv_to_tfrecords():
    td = tempfile.mkdtemp()
    print(td)
    outfile1 = os.path.join(td, "2testseqs.tfrecord")
    vica.tfrecord_maker._csv_to_tfrecords(kmerfile="tests/test-data/5mers.csv",
        codonfile="tests/test-data/codons.csv",
        minhashfile="tests/test-data/minhashcsvtest.csv",
        mergefile="tests/test-data/mergefile.csv",
        tfrecordfile=outfile1,
        label=0)
    # ordering of fields in TFrecords is stochastic so just compare size
    refsize = os.stat("tests/test-data/2testseqs.tfrecord").st_size
    eq_(refsize, os.stat(outfile1).st_size)

def test_convert_to_tfrecords():
    td = tempfile.mkdtemp()
    print(td)
    outfile1 = os.path.join(td, "2testseqs.tfrecord")
    vica.tfrecord_maker.convert_to_tfrecords(dtemp= td,
        kmerfile="tests/test-data/5mers.csv",
        codonfile="tests/test-data/codons.csv",
        minhashfile="tests/test-data/minhashcsvtest.csv",
        tfrecordfile=outfile1,
        label=0,
        sort=True)
    # ordering of fields in TFrecords is stochastic so just compare size
    refsize = os.stat("tests/test-data/2testseqs.tfrecord").st_size
    eq_(refsize, os.stat(outfile1).st_size)
