import vica.minhash
import tempfile
import os
import filecmp
import yaml
from nose.tools import ok_

with open(vica.CONFIG_PATH) as cf:
    config = yaml.safe_load(cf)

def test_send_sketch():
    td = tempfile.mkdtemp()
    outfile1 = os.path.join(td,"sendsketchout.txt")
    vica.minhash._send_sketch(infile="tests/test-data/2testseqs.fasta",
                             outfile=outfile1,
                             server_url=config["minhash"]["server_url"])
    ok_(filecmp.cmp("tests/test-data/testsketch1.txt", outfile1, shallow=False))
