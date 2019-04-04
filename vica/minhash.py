"""A module to run bbtools minhash functions sequences data and
return a file of tab delimited classification data

"""

import subprocess

import logging
from ete3 import NCBITaxa
import yaml
import json

import vica

with open(vica.CONFIG_PATH) as cf:
    config = yaml.safe_load(cf)

def _send_sketch(infile, server_url, level=3):
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
               "level=" + str(level),
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
    """Takes taxid and returns classifier class of that id

        Args:
            taxid (int): the NCBI TaxID
            classdict (dict): dictionary where the keys are the classes used by the classifier
            taxinstance (obj): an ete3.NCBITaxa instance

        Returns:
            int: the class of the taxid
    """
    #try:
    lineage = taxinstance.get_lineage(taxid)
    classtaxid = list(set(classdict.keys()).intersection(lineage))
    if not len(classtaxid) == 1:
        logging.info("Could not assign taxid %s to a higher taxonomic level", taxid)
        return None
    return classtaxid[0]
    #except Exception:
    #    logging.info("could not assign %s to a class, classtaxid is %s" % (taxid, classtaxid))
    #    return None



def _same_clade_as_query(hit, query):
    """Check if the hit taxid is in the same clade as the query taxid

        Args:
            hit (int): the taxid of the minhash hit
            query (int): the taxid of the query sequence
        Returns:
            bool: True if taxid is in hte sam clase as the query, false otherwise

    """
    if str(hit) == str(query):
        return True
    return False


def _parse_sendsketch(dataraw: str, cutoff: float=100., filtertaxa: bool=False) -> dict:
    """Parse json string from minhashremote into a dict of sequences and
        the highest scoring class they hit to.

        Args:
            dataraw (str): the string of multiple json objects returned by _send_sketch
            cutoff (float): the minimum score to record a hits
            filtertaxa (bool): should hits be filtered to remove hits in the \
                        same genera (by default) as the query. This prevents \
                        information bleedover in the test set
        Returns:
            dict: a Dictionary containing the {query: classid, ...}

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
            try:
                if key == "Name":
                    name = val
                    query_taxid = name.split("|")[1]
                elif key not in ["DB", "SketchLen", "Seqs", "Bases", "gSize", "file"]:
                    score = val["Score"]
                    taxid = val["TaxID"]
                    if filtertaxa:
                        filter_them = _same_clade_as_query(taxid, query_taxid)
                    else:
                        filter_them = False
                    if score > cutoff and filter_them==False:
                        classid = _taxid_2_taxclass(taxid=taxid,
                                                    classdict=config["split_shred"]["classes"],
                                                    taxinstance=ncbi)
                        datadict[name] = classid
                        continue
            except:
                logging.info("error parsing minhash sample, continuing")
    return datadict

def minhashremote(infile, outfile, server_url, filtertaxa=False):
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
        none

    """
    logging.info("Using BBtools Sendsketch.sh to send minhash sketches to the server %s", server_url)
    dataraw = _send_sketch(infile=infile, server_url=server_url)
    logging.info("Parsing results file from BBtools Sendsketch.sh")
    datadict = _parse_sendsketch(dataraw, filtertaxa=filtertaxa)
    with open(outfile, 'w') as ofile:
        json.dump(datadict, ofile)
