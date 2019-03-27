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
               "records=10",
               "format=json",
               "score=t",
               "fixjunk"]
    sendsketchout = subprocess.run(options, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    dataraw = sendsketchout.stdout.decode('utf-8')
    logging.info(sendsketchout.stderr.decode('utf-8'))
    return dataraw

def _taxid_2_taxclass(taxid, classdict, taxinstance):
    try:
        lineage = taxinstance.get_lineage(taxid)
        classtaxid = list(set(classdict.keys()).intersection(lineage))
        if not len(classtaxid) == 1:
            logging.info("Could not assign taxid %s to a higher taxonomic level", taxid)
            return None
        return classtaxid[0]
    except:
        logging.info("could not assign %s to a class", taxid)
        return None
#>>> Stopped here
def _same_clade_as_query(hit, query, taxinstance, taxtree, level):
    try:
        lca = taxtree.get_common_ancestor[str(hit), str(query)]

        rank = taxinstance.get_rank([int(lca.name)]).decode("utf-8")
        if rank in level:
            return True
        else:
            return False
    except:
        logging.info("there was an issue parsing the lowest common ancestor for query %s", query)
        return False


def _parse_sendsketch(dataraw, cutoff=100, filtertaxa=False):
    """Parse json string into a dict with

    """

    json_str = "".join(dataraw.splitlines())
    dec = json.JSONDecoder()
    pos = 0
    datadict = {}
    ncbi = NCBITaxa()
    taxtree = ncbi.get_topology([1])
    while not pos == len(str(json_str)):
        j, json_len = dec.raw_decode(str(json_str)[pos:])
        pos += json_len
        for key, val in j.items():
            if key == "Name":
                name = val
                query_taxid = name.split("|")[1]
            elif key not in ["DB", "SketchLen", "Seqs", "Bases", "gSize", "File"]:
                score = val["Score"]
                taxid = val["TaxID"]
                if filtertaxa:
                    filter = _same_clade_as_query(taxid, query_taxid,
                                                  ncbi, taxtree,
                                                  config["minhash"]["taxfilterlevel"])
                else:
                    filter = False
                if score > cutoff and not filter:
                    classid = _taxid_2_taxclass(taxid=taxid,
                                                classdict=config["split_shred"]["classes"],
                                                taxinstance=ncbi)
                    datadict[name] = classid
                    continue
    return datadict

def minhashremote(infile, outfile, server_url, filtertaxa=False ):
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
        filtertaxa (bool): should minhash hits in the same group as the query be removed?

    Returns:

    """
    logging.info("Using BBtools Sendsketch.sh to send minhash sketches to the server {}".format(server_url))
    dataraw = _send_sketch(infile=infile, server_url=server_url)
    logging.info("Parsing results file from BBtools Sendsketch.sh")
    datadict = _parse_sendsketch(dataraw, taxfilte=taxfilter)
    with open(outfile, 'w') as ofile:
        json.dump(datadict, ofile)
