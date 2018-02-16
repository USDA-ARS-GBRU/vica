import vica
import tempfile
import os
import filecmp
import shutil
from nose.tools import ok_ , eq_

kmer_list = ['AAAA', 'AAAC', 'AAAT', 'AAAG', 'AACA', 'AACC', 'AACT', 'AACG',
 'AATA', 'AATC', 'AATT', 'AATG', 'AAGA', 'AAGC', 'AAGT', 'AAGG', 'ACAA',
 'ACAC', 'ACAT', 'ACAG', 'ACCA', 'ACCC', 'ACCT', 'ACCG', 'ACTA', 'ACTC',
 'ACTG', 'ACGA', 'ACGC', 'ACGT', 'ACGG', 'ATAA', 'ATAC', 'ATAT', 'ATAG',
 'ATCA', 'ATCC', 'ATCT', 'ATCG', 'ATTA', 'ATTC', 'ATTG', 'ATGA', 'ATGC',
 'ATGG', 'AGAA', 'AGAC', 'AGAG', 'AGCA', 'AGCC', 'AGCT', 'AGCG', 'AGTA',
 'AGTC', 'AGTG', 'AGGA', 'AGGC', 'AGGG', 'CAAA', 'CAAC', 'CAAG', 'CACA',
 'CACC', 'CACG', 'CATA', 'CATC', 'CATG', 'CAGA', 'CAGC', 'CAGG', 'CCAA',
 'CCAC', 'CCAG', 'CCCA', 'CCCC', 'CCCG', 'CCTA', 'CCTC', 'CCGA', 'CCGC',
 'CCGG', 'CTAA', 'CTAC', 'CTAG', 'CTCA', 'CTCC', 'CTCG', 'CTTA', 'CTTC',
 'CTGA', 'CTGC', 'CGAA', 'CGAC', 'CGCA', 'CGCC', 'CGCG', 'CGTA', 'CGTC',
 'CGGA', 'CGGC', 'TAAA', 'TAAC', 'TACA', 'TACC', 'TATA', 'TATC', 'TAGA',
 'TAGC', 'TCAA', 'TCAC', 'TCCA', 'TCCC', 'TCTC', 'TCGA', 'TCGC', 'TTAA',
 'TTAC', 'TTCA', 'TTCC', 'TTTC', 'TTGC', 'TGAC', 'TGCA', 'TGCC', 'TGTC',
 'TGGC', 'GAAC', 'GACC', 'GATC', 'GAGC', 'GCAC', 'GCCC', 'GCGC', 'GTAC',
 'GTCC', 'GGCC']

def almost_equal(value_1, value_2, accuracy = 10**-8):
    '''a Function to compare tuples of float values'''
    return abs(value_1 - value_2) < accuracy

def test_iterate_kmer():
    expected = kmer_list
    results = vica.khmer_features.iterate_kmer(4)
    eq_(results, expected, msg="The fuction itrate Khmer did not return the \
        expected codon table for k=4")

def test_get_composition():
    seq='GTGAAACATTCGGATTTCGATATTGTCGTGATCGGGGCCGGACATGCCGGCGCCGAGGCTGCACATGCTGCGGCACGCATGGGAATGCGTACTGCCTTGGTTTCCCTGTCCGAACGCGACATTGGCGTGATGTCCTGTAA'
    expected1 = [0, 2, 1, 0, 1, 1, 0, 1, 1, 1, 0, 3, 0, 0, 0, 1, 1, 0, 5, 2, 1,
                 0, 0, 0, 0, 0, 1, 1, 4, 0, 0, 0, 0, 1, 0, 2, 1, 0, 2, 0, 2, 2,
                 0, 4, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 2, 1, 0, 0, 1, 1, 0,
                 3, 0, 1, 3, 0, 2, 2, 2, 0, 0, 1, 1, 1, 0, 1, 4, 1, 2, 0, 0, 0,
                 0, 0, 1, 0, 0, 0, 3, 3, 2, 3, 3, 1, 1, 0, 3, 5, 0, 0, 1, 0, 0,
                 1, 0, 0, 0, 3, 0, 2, 0, 1, 1, 0, 1, 1, 2, 3, 0, 0, 1, 3, 5, 1,
                 1, 0, 1, 0, 2, 1, 1, 1, 3, 1]

    r1 = vica.khmer_features.get_composition(ksize=4, seq=seq, kmers=kmer_list, norm=False)
    eq_(r1,expected1, msg="The un-normalized composition did not match it's expected result")

    expected2 = [0, 0.014598540145985401, 0.0072992700729927005, 0, 0.0072992700729927005, 0.0072992700729927005, 0, 0.0072992700729927005, 0.0072992700729927005, 0.0072992700729927005, 0, 0.021897810218978103, 0, 0, 0, 0.0072992700729927005, 0.0072992700729927005, 0, 0.0364963503649635, 0.014598540145985401, 0.0072992700729927005, 0, 0, 0, 0, 0, 0.0072992700729927005, 0.0072992700729927005, 0.029197080291970802, 0, 0, 0, 0, 0.0072992700729927005, 0, 0.014598540145985401, 0.0072992700729927005, 0, 0.014598540145985401, 0, 0.014598540145985401, 0.014598540145985401, 0, 0.029197080291970802, 0.0072992700729927005, 0, 0, 0, 0.0072992700729927005, 0.0072992700729927005, 0, 0, 0.0072992700729927005, 0, 0, 0.0072992700729927005, 0.014598540145985401, 0.0072992700729927005, 0, 0, 0.0072992700729927005, 0.0072992700729927005, 0, 0.021897810218978103, 0, 0.0072992700729927005, 0.021897810218978103, 0, 0.014598540145985401, 0.014598540145985401, 0.014598540145985401, 0, 0, 0.0072992700729927005, 0.0072992700729927005, 0.0072992700729927005, 0, 0.0072992700729927005, 0.029197080291970802, 0.0072992700729927005, 0.014598540145985401, 0, 0, 0, 0, 0, 0.0072992700729927005, 0, 0, 0, 0.021897810218978103, 0.021897810218978103, 0.014598540145985401, 0.021897810218978103, 0.021897810218978103, 0.0072992700729927005, 0.0072992700729927005, 0, 0.021897810218978103, 0.0364963503649635, 0, 0, 0.0072992700729927005, 0, 0, 0.0072992700729927005, 0, 0, 0, 0.021897810218978103, 0, 0.014598540145985401, 0, 0.0072992700729927005, 0.0072992700729927005, 0, 0.0072992700729927005, 0.0072992700729927005, 0.014598540145985401, 0.021897810218978103, 0, 0, 0.0072992700729927005, 0.021897810218978103, 0.0364963503649635, 0.0072992700729927005, 0.0072992700729927005, 0, 0.0072992700729927005, 0, 0.014598540145985401, 0.0072992700729927005, 0.0072992700729927005, 0.0072992700729927005, 0.021897810218978103, 0.0072992700729927005]

    r2 = vica.khmer_features.get_composition(ksize=4, seq=seq, kmers=kmer_list, norm=True)
    ok_(all(almost_equal(*values) for values in zip(r2, expected2)), "The normalized composition did not match it's expected result")


def test_write_kmers_as_csv():
    infile = 'tests/test-data/2testseqs.fasta'
    expected = 'tests/test-data/4mers.csv'
    dtemp = tempfile.mkdtemp()
    print(dtemp)
    outfile = os.path.join(dtemp,"outfile.csv")
    vica.khmer_features._write_kmers_as_csv(infile=infile, outfile=outfile, ksize=4, kmers=kmer_list)
    ok_(filecmp.cmp(outfile, expected))
    shutil.rmtree(dtemp)
