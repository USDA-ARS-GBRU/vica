
import tempfile
import os
import filecmp
import yaml
import shutil

from nose.tools import eq_ # ok_,

import vica.minhash

with open(vica.CONFIG_PATH) as cf:
    config = yaml.safe_load(cf)


def test_send_sketch():
    expected = '{\n   "Name": "NC_003911.12 Ruegeria pomeroyi DSS-3, first 7000nt",\n   "DB": "RefSeq",\n   "SketchLen": 191,\n   "Seqs": 1,\n   "Bases": 7000,\n   "gSize": 6326,\n   "GC": "0.607",\n   "Ruegeria pomeroyi DSS-3": {\n      "taxName": "Ruegeria pomeroyi DSS-3",\n      "WKID": "100.0000",\n      "KID": "0.1490",\n      "ANI": "100.000",\n      "Complt": "0.149",\n      "Contam": "0.000",\n      "Score": 197.13017,\n      "Matches": 33,\n      "Unique": 27,\n      "TaxID": 246200,\n      "gSize": 4544804,\n      "gSeqs": 2,\n      "GC": "0.641"\n   }\n}{\n   "Name": "NC_005072.1 Prochlorococcus marinus MED4 first 7000nt",\n   "DB": "RefSeq",\n   "SketchLen": 191,\n   "Seqs": 1,\n   "Bases": 7000,\n   "gSize": 6969,\n   "GC": "0.293",\n   "Prochlorococcus marinus subsp. pastoris str. CCMP1986": {\n      "taxName": "Prochlorococcus marinus subsp. pastoris str. CCMP1986",\n      "WKID": "100.0000",\n      "KID": "0.4121",\n      "ANI": "100.000",\n      "Complt": "0.412",\n      "Contam": "0.000",\n      "Score": 333.99326,\n      "Matches": 65,\n      "Unique": 0,\n      "TaxID": 59919,\n      "gSize": 1652527,\n      "gSeqs": 1,\n      "GC": "0.308"\n   }\n}\n'
    dataout = vica.minhash._send_sketch(infile="tests/test-data/2testseqs.fasta",
                                        server_url=config["minhash"]["server_url"])
    eq_(dataout, expected)


def test_parse_sendsketch():
    dataout = vica.minhash._send_sketch(infile="tests/test-data/2testseqs.fasta",
                                        server_url=config["minhash"]["server_url"])
    dout = vica.minhash._parse_sendsketch(dataraw = dataout)
    expected = {'NC_003911.12 Ruegeria pomeroyi DSS-3, first 7000nt': 2, 'NC_005072.1 Prochlorococcus marinus MED4 first 7000nt': 2}
    eq_(dout, expected)


def test_minhashremote():
    td = tempfile.mkdtemp()
    outfile = os.path.join(td, "minhashremote.json")
    vica.minhash.minhashremote(infile="tests/test-data/2testseqs.fasta",
                               outfile=outfile,
                               server_url=vica.minhash.config["minhash"]["server_url"],
                               filtertaxa=False)
    #ok_(filecmp.cmp("tests/test-data/minhashremote.json", outfile, shallow=True))
    shutil.rmtree(td)
