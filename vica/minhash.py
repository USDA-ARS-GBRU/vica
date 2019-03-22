"""A module to run bbtools minhash functions sequences data and
return a file of tab delimited classification data

"""

import subprocess
import os
import csv

import logging
from ete3 import NCBITaxa
import yaml
import json

import vica

with open(vica.CONFIG_PATH) as cf:
    config = yaml.safe_load(cf)

def _send_sketch(infile, server_url):
    """Runs bbtools sendsketch.sh on a file of sequences returning a
        classification for each

    Calculates minhash sketches if (k=32, 24) for each sequence in a fasta
    sends the sketches to a server and returns the classification data to
    a file.

    Args:
        infile (str): A multi-sequence fasta file for which to calulate
            minhash taxonomy.
        server_url (str): a URL for the minhash server

    Returns:
        (str): raw string of json sent by _send_sketch

    """

    options = ["sendsketch.sh",
               "in=" + infile,
               "address=" + server_url,
               "mode=sequence",
               "level=3",
               "k=32,24",
               "records=1",
               "format=json",
               "score=t",
               "fixjunk"]
    sendsketchout = subprocess.run(options, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    dataraw = sendsketchout.stdout.decode('utf-8')
    logging.info(sendsketchout.stderr.decode('utf-8'))
    return dataraw

def _taxid_2_taxclass(taxid, classdict, taxinstance):
    lineage = taxinstance.get_lineage(taxid)
    classtaxid = list(set(classdict.keys()).intersection(lineage))
    assert len(classtaxid) == 1
    return classtaxid[0]

def _parse_sendsketch(dataraw, cutoff=100):
    """Parse json string into a dict with

    """

    json_str = "".join(dataraw.splitlines())
    dec = json.JSONDecoder()
    pos = 0
    datadict = {}
    ncbi = NCBITaxa()
    while not pos == len(str(json_str)):
        j, json_len = dec.raw_decode(str(json_str)[pos:])
        pos += json_len
        for key, val in j.items():
            if key == "Name":
                name = val
            elif key not in ["DB", "SketchLen", "Seqs", "Bases", "gSize", "File"]:
                score = val["Score"]
                if score > cutoff:
                    taxid = val["TaxID"]
                    classid = _taxid_2_taxclass(taxid=taxid,
                                                classdict=config["split_shred"]["classes"],
                                                taxinstance=ncbi)
                    datadict[name] = classid
                    continue
    return datadict

def minhashremote(infile, outfile, server_url):
    """Runs bbtools sendeketch.sh on a file of sequences returning a
        classification for each

    Calculates minhash sketches if (k=24, 31) for each sequence in a fasta
    and uses a remote server then returns the classification data to
    a file.

    Args:
        infile (str): A multi-sequence fasta file for which to calulate
            minhash taxonomy.
        outfile (str): a path to write the file retuned from the minhash
            server
        server_url (str): a URL for the minhash server

    Returns:

    """
    logging.info("Using BBtools Sendsketch.sh to send minhash sketches to the server {}".format(server_url))
    dataraw = _send_sketch(infile=infile, server_url=server_url)
    logging.info("Parsing results file from BBtools Sendsketch.sh")
    datadict = _parse_sendsketch(dataraw)
    print(datadict)
    with open(outfile, 'w') as ofile:
        json.dump(datadict, ofile)
