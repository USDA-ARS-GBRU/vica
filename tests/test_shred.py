import vica.shred
import tempfile
import os
import filecmp
import nose
from pyfaidx import Fasta


def test_parse_samples():
    samples = 0.5
    pout = vica.shred.parse_samples(samples)
    assert filecmp.cmp(pout == 'fraction')
    samples = 0.9999999998
    pout = vica.shred.parse_samples(samples)
    assert filecmp.cmp(pout == 'all')
    samples = 1
    pout = vica.shred.parse_samples(samples)
    assert filecmp.cmp(pout == 'all')
    samples = 10
    pout = vica.shred.parse_samples(samples)
    assert filecmp.cmp(pout == 'set_number')
    samples = 0
    nose.tools.asert_raises(ValueError,vica.shred.parse_samples, samples)
    samples = -1
    nose.tools.asert_raises(ValueError,vica.shred.parse_samples, samples)

def test_write_seq():
    seqhandle = Fasta('tests/test-data/2testseqs.fasta')
    tf = tempfile.NamedTemporaryFile(mode='w')
    vica.shred.writeseq(record=seqhandle,pos=1, length=1000,handle=tf)
    assert filecmp.cmp(tf, 'tests/test-data/shredwrite.Fasta', shallow=False)

def test_shred_all():
    pass
