import vica.minhash
import tempfile
import os
import filecmp
import yaml
from nose.tools import ok_, eq_
import difflib
import sys

from ete3 import NCBITaxa

with open(vica.CONFIG_PATH) as cf:
    config = yaml.safe_load(cf)

vallist = [0, 976, 1090, 1117, 1224, 1239, 1297, 2833, 2836, 2870, 3035, 3041, 4761, 4890, 5204, 5302, 5747, 5794, 6029, 6040, 6073, 6157, 6217, 6231, 6340, 6447, 6656, 6657, 6843, 7568, 7586, 7711, 7712, 7735, 10190, 10197, 10205, 10213, 10214, 10219, 10226, 10229, 10232, 12333, 12429, 12877, 12884, 27563, 28889, 28890, 29000, 29258, 31291, 32066, 33209, 33310, 33313, 33467, 35237, 35268, 35325, 35493, 39759, 40117, 42241, 43120, 51516, 51967, 57723, 61985, 62680, 65842, 66780, 67810, 67812,   67814, 67817, 67818, 67819, 67820, 68297, 68525, 69815, 74015, 74152, 74201, 89593, 91989, 95818, 104731, 110814, 115360, 115365, 115369, 120557, 131221, 134625, 142182, 142187, 147099, 147537, 147538, 157124, 192989, 200783, 200795, 200918, 200930, 200938, 200940, 201174, 203682, 203691, 204428, 214504, 221216, 256845, 265317, 310840, 363464, 419944, 422282, 431837, 431838, 439488, 447827, 451344, 451459, 451507, 451827, 451828, 451866, 452284, 456828, 508458, 544448, 552364, 569577, 640293, 651137, 686617, 743724, 743725, 877183, 928852, 1031332, 1052197, 1052815, 1104542, 1134404, 1137986, 1154675, 1154676, 1215728, 1264859, 1293497, 1312402, 1379697, 1383058, 1448051, 1448933, 1448937, 1462422, 1462430, 1492816, 1618330, 1618338, 1618339, 1618340, 1619053, 1655434, 1696033, 1703755, 1704031, 1706441, 1714266, 1729712, 1752708, 1752716, 1752717, 1752718, 1752719, 1752720, 1752721, 1752722, 1752723, 1752724, 1752725, 1752726, 1752727, 1752728, 1752729, 1752730, 1752731, 1752732, 1752733, 1752734, 1752735, 1752736, 1752737, 1752738, 1752739, 1752740, 1752741, 1798710, 1801616, 1801631, 1802339, 1817796, 1817797, 1817798, 1817799, 1817800, 1817801, 1817802, 1817803, 1817804, 1817805, 1817806, 1817807, 1817808, 1817809, 1817810, 1817811, 1817812, 1817898, 1817899, 1817900, 1817901, 1817902, 1817903, 1817904, 1817905, 1817906, 1817907, 1817908, 1817909, 1817910, 1817911, 1817912, 1817913, 1817914, 1817915, 1817916, 1817917, 1817918, 1817919, 1817920, 1817921, 1817922, 1817923, 1817924, 1817925, 1817926, 1819803, 1853220, 1855361, 1910928, 1913637, 1913638, 1915410, 1922347, 1930617, 1936271, 1936272, 1936987, 2013583]


def test_send_sketch():
    td = tempfile.mkdtemp()
    outfile1 = os.path.join(td,"sendsketchout.txt")
    vica.minhash._send_sketch(infile="tests/test-data/2testseqs.fasta",
                             outfile=outfile1,
                             server_url=config["minhash"]["server_url"])
    ok_(filecmp.cmp("tests/test-data/testsketch1.txt", outfile1, shallow=False))


def test_compare_sketch():
    td = tempfile.mkdtemp()
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
    with open(outfile3, "r") as a:
        with open("tests/test-data/testsketch3.txt", 'r') as b:
            al = a.readlines()
            bl = b.readlines()
            cd = cd = difflib.context_diff(al, bl)
    sys.stdout.writelines(cd)
    ok_(filecmp.cmp("tests/test-data/testsketch3.txt", outfile3, shallow=True))


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
    eq_(7711, higher_level)

def test_raise_taxdict_level():
    ncbi = NCBITaxa()
    testdict = {246200: 181.8, 190047: 259.8}
    higher_level = vica.minhash._raise_taxdict_level(testdict, ncbi)
    eq_({1224: 181.8, 1117: 259.8}, higher_level)

def test_get_feature_list():
    testlist = vica.minhash._get_feature_list("vica/data/phylum.txt",
        vica.minhash.config["minhash"]["noncellular"])
    eq_(testlist, vallist)

def test_dict_to_csv():
    testdict = {'NC_003911.12': {246200: 181.8}, 'NC_005072.1': {190047: 259.8}}
    td = tempfile.mkdtemp()
    outfile = os.path.join(td,"minhashcsvtest.csv")
    vica.minhash._dict_to_csv(testdict, vallist, outfile)
    ok_(filecmp.cmp("tests/test-data/minhashcsvtest.csv", outfile, shallow=False))

# def test_minhashlocal():
#         pass
#         outfile3 = os.path.join(td,"comparesketchout3.txt")
#         vica.minhash._compare_sketch(infile="tests/test-data/2testseqs.fasta",
#                                   outfile=outfile3,
#                                   ref="tests/test-data/2testseqs_ref.sketch",
#                                   blacklist= "tests/test-data/blacklist_refseq_species_300.sketch",
#                                   tree="tests/test-data/tree.taxtree.gz",
#                                   taxfilter="246200",
#                                   taxfilterlevel="species",
#                                   memory="-Xmx1g")
