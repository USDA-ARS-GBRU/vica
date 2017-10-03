import vica.minhash
import tempfile
import os
import filecmp
import nose

def test_send_sketch():
    td = tempfile.mkdtemp()
    outfile1 = os.path.join(td,"sendsketchout.txt")
    vica.minhash.send_sketch(infile="tests/test-data/2testseqs.fasta",
                             outfile=outfile1)
    assert filecmp.cmp("tests/test-data/testsketch1.txt", outfile1, shallow=False)
