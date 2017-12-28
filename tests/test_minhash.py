import vica.minhash
import tempfile
import os
import filecmp
import yaml
from nose.tools import ok_, eq_

from ete3 import NCBITaxa

with open(vica.CONFIG_PATH) as cf:
    config = yaml.safe_load(cf)

def test_send_sketch():
    td = tempfile.mkdtemp()
    outfile1 = os.path.join(td,"sendsketchout.txt")
    vica.minhash._send_sketch(infile="tests/test-data/2testseqs.fasta",
                             outfile=outfile1,
                             server_url=config["minhash"]["server_url"])
    ok_(filecmp.cmp("tests/test-data/testsketch1.txt", outfile1, shallow=False))


def test_compare_sketch():
    td = tempfile.mkdtemp()
    print(td)
    #test comparesketch without any filtering
    outfile2 = os.path.join(td,"comparesketchout.txt")
    vica.minhash._compare_sketch(infile="tests/test-data/2testseqs.fasta",
                              outfile=outfile2,
                              ref="tests/test-data/2testseqs_ref.sketch",
                              blacklist= "tests/test-data/blacklist_refseq_species_300.sketch",
                              tree="tests/test-data/tree.taxtree.gz",
                              taxfilter=None,
                              taxfilterlevel=None,
                              memory="-Xmx1g")

def test_compare_sketch_with_taxfilter():
    td = tempfile.mkdtemp()
    print(td)
    #test comparesketch without any filtering
    outfile3 = os.path.join(td,"comparesketchout3.txt")
    vica.minhash._compare_sketch(infile="tests/test-data/2testseqs.fasta",
                              outfile=outfile3,
                              ref="tests/test-data/2testseqs_ref.sketch",
                              blacklist= "tests/test-data/blacklist_refseq_species_300.sketch",
                              tree="tests/test-data/tree.taxtree.gz",
                              taxfilter="246200",
                              taxfilterlevel="species",
                              memory="-Xmx1g")
    ok_(filecmp.cmp("tests/test-data/testsketch3.txt", outfile3, shallow=False))


def test_parse_sendsketch():
    dout = vica.minhash._parse_sendsketch("tests/test-data/testsketch1.txt")
    expected = {'NC_003911.12': {246200: 181.8}, 'NC_005072.1': {190047: 259.8}}
    eq_(dout, expected)

def test_parse_comparesketch():
    dout = vica.minhash._parse_comparesketch("tests/test-data/testsketch2.txt")
    expected = {'NC_003911.12': {246200: 181.8}, 'NC_005072.1': {59919: 259.9}}
    eq_(dout, expected)

def test_find_key():
    testdict = {'NC_003911.12': {246200: 181.8}, 'NC_005072.1': {190047: 259.8}}
    val = {246200: 181.8}
    output = vica.minhash._find_key(testdict,val)
    eq_('NC_003911.12', output)

def test_pick_higher_level():
    ncbi = NCBITaxa()
    higher_level = vica.minhash._pick_higher_level(9606, ncbi)
    ok_(7711, higher_level)

def test_raise_taxdict_level():
    ncbi = NCBITaxa()
    testdict = {246200: 181.8, 190047: 259.8}
    higher_level = vica.minhash._raise_taxdict_level(testdict, ncbi)
    ok_({1224: 181.8, 1117: 259.8}, higher_level)
