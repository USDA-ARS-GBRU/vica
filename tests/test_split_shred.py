import tempfile
import os
import filecmp
import random
import shutil

import pandas
import numpy
import pyfaidx
from nose.tools import ok_, eq_
from ete3 import NCBITaxa

import vica.split_shred


classes = {2: "Bacteria",
         2157: "Archaea",
         2759: "Eukaryota",
         10239: "Viruses"}

def test_read_data_fasta():
    pyfaidx_obj = vica.split_shred._read_data("tests/test-data/2testseqs.fasta")
    ok_(isinstance(pyfaidx_obj, pyfaidx.Fasta))

def test_read_data_bfgz():
        td = tempfile.mkdtemp()
        # copy file so that pyfaidx must generate a new index
        shutil.copy("tests/test-data/2testseqs.fasta.gz", td)
        infile = os.path.join(td, "2testseqs.fasta.gz")
        pyfaidx_obj = vica.split_shred._read_data(infile)
        ok_(isinstance(pyfaidx_obj, pyfaidx.Fasta))
        shutil.rmtree(td)

def test_shuffle_keys():
    ddict = {'a':5, 'b':6, 'c':7, 'd':8, 'f':6, 'h':10, 'j':11}
    shuffled_list = ['d', 'h', 'f', 'b', 'a', 'c', 'j']
    # setting random seed is not stable across python versions
    # for that reason we only compare the content of the lists not the order
    random.seed(a=123456)
    random_list = vica.split_shred._shuffle_keys(ddict)
    eq_(sorted(shuffled_list), sorted(random_list))

def test_profile_sequences():
    seqobj = vica.split_shred._read_data("tests/test-data/bbtools-taxheader.fa")
    ncbiobj = NCBITaxa()
    splitlevel = 'genus'
    # setting random seed is not stable across python versions
    # for that reason we re-sort the daa frames
    random.seed(a=123456)
    df = vica.split_shred._profile_sequences(seqobj=seqobj,
        ncbiobj=ncbiobj,
        splitlevel = splitlevel,
        classes= classes)
    df = df.sort_index()
    df2 = pandas.read_csv('tests/test-data/test_profile_sequences_data.csv', index_col=0).sort_index()
    ok_(df.equals(df2))

def test_split_levels():
    testfrac =0.5
    df = pandas.read_csv('tests/test-data/test_profile_sequences_data.csv', index_col=0)
    numpy.random.seed(seed=123456)
    cd = vica.split_shred._split_levels(testfrac=testfrac, df=df, classes=classes)
    print(cd)
    expected = {2: {'test': {1378, 212743, 114248, 1765682, 92793, 946234,
                190972, 160798}, 'train': {1017280, 33986, 28067, 1484898,
                467084, 54066, 1938003, 702}, 'total': 16}, 2759: {'test':
                {1955842, 5931, 12967},
                'train': {1633384, 4783, 148959}, 'total': 6}}

    eq_(cd, expected)

def test_writeseq():
    recordid = 'tid|2340|NZ_MPQT01000037.1'
    data = pyfaidx.Fasta('tests/test-data/bbtools-taxheader.fa')
    td = tempfile.mkdtemp()
    outfile = os.path.join(td,'test_writeseq.fa')
    with open(outfile, 'w') as f:
        vica.split_shred._writeseq(data[recordid], pos=0, length=4000, handle=f)
    ok_(filecmp.cmp(outfile, 'tests/test-data/test_writeseq.fa'))

def test_select_random_segment():
    seqobj = pyfaidx.Fasta('tests/test-data/bbtools-taxheader_badexamples.fa')
    name = 'good_seq'
    length = 500
    numpy.random.seed(seed=123456)
    seq = vica.split_shred._select_random_segment(seqobj=seqobj,
        name=name,
        length=length,
        tries=10,
        ns= 0.1)
    eq_(seq, 3121)

def test_select_random_segment_ns():
    seqobj = pyfaidx.Fasta('tests/test-data/bbtools-taxheader_badexamples.fa')
    name = 'too_many_ns'
    length = 500
    numpy.random.seed(seed=123456)
    seq = vica.split_shred._select_random_segment(seqobj=seqobj,
        name=name,
        length=length,
        tries=10,
        ns= 0.1)
    eq_(seq, None)

def test_select_random_segment_short():
    seqobj = pyfaidx.Fasta('tests/test-data/bbtools-taxheader_badexamples.fa')
    name = 'too_short'
    length = 500
    numpy.random.seed(seed=123456)
    seq = vica.split_shred._select_random_segment(seqobj=seqobj,
        name=name,
        length=length,
        tries=10,
        ns= 0.1)
    eq_(seq, None)

def _read_taxid_from_fasta():
    td = tempfile.mkdtemp()
    # copy file so that pyfaidx must generate a new index
    outdir=os.path.join(td, 'ex_test_train_split_dir')
    outfile= os.path.join(outdir,'test','test_taxids.txt')
    resultsfile = ("tests/test-data/ex_test_train_split_dir/test/test_taxids.txt")
    shutil.copytree("tests/test-data/ex_test_train_split_dir", outdir)
    recs = vica.split_shred._read_taxid_from_fasta(outdir)
    ok_(filecmp.cmp(outfile, resultsfile))
    eq_(recs, 62)

# def _process_samples():
