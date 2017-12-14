import vica.minhash
import tempfile
import os
import filecmp
import yaml

with open(vica.CONFIG_PATH) as cf:
    config = yaml.safe_load(cf)

def test_send_sketch():
    td = tempfile.mkdtemp()
    outfile1 = os.path.join(td,"sendsketchout.txt")
    vica.minhash._send_sketch(infile="tests/test-data/2testseqs.fasta",
                             outfile=outfile1,
                             server_url=config["minhash"]["server_url"])
    assert filecmp.cmp("tests/test-data/testsketch1.txt", outfile1, shallow=False)
