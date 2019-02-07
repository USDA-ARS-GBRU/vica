"""A module to run bbtools minhash functions sequences data and
return a file of tab delimited classification data

"""

import subprocess
import os
import csv

import logging
from ete3 import NCBITaxa
import yaml

import vica

with open(vica.CONFIG_PATH) as cf:
    config = yaml.safe_load(cf)

def _send_sketch(infile, outfile, server_url):
    """Runs bbtools sendsketch.sh on a file of sequences returning a
        classification for each

    Calculates minhash sketches if (k=32, 24) for each sequence in a fasta
    sends the sketches to a server and returns the classification data to
    a file.

    Args:
        infile (str): A multi-sequence fasta file for which to calulate
            minhash taxonomy.
        outfile (str): a path to write the file retuned from the minhash
            server
        server_url (str): a URL for the minhash server

    Returns:
        (str) The standard output from BBtools comparesketch.sh

    """

    options = ["sendsketch.sh",
               "in=" + infile,
               "out=" + outfile,
               "address=" + server_url,
               "mode=sequence",
               "level=3",
               "color=f",
               "overwrite=t",
               "k=32,24",
               "printani=f",
               "completeness=f",
               "score=t",
               "printmatches=t",
               "printgsize=f",
               "printgseqs=f",
               "printtaxname=f",
               "printname0=f",
               "printcontam=f",
               "printunique=t",
               "printnohit=f",
               "printtaxid=t"]
    sendsketchout = subprocess.run(options, stderr=subprocess.PIPE)
    return sendsketchout.stderr.decode('utf-8')
    #return sendsketchout



def _parse_sendsketch(file):
    """Parses bbtools sendsketch output returning python dictionary.

    Args:
        file (str): a text file created by BBtools sendsketch.sh

    Returns:
        (dict): A dictionary with the ID as a key and as a value, a dict
        with taxid: score for each hit identified

    See Also:
        vica.minhash._parse_comparesketch

    """
    try:
        tempdf = {}
        with open(file, 'r') as f:
            for line in f:
                if line.strip() == '':
                    next
                elif line.startswith("Query:"):
                    ll = line.strip().split("\t")
                    key1 = ll[0].split(":")[1].split()[0]
                    tempdf[key1] = {}
                elif line.startswith("WKID"):
                    next
                elif line.startswith("No hits."):
                    tempdf[key1]['0'] = 0
                else:
                    ll2 = line.strip().split("\t")
                    tempdf[key1][int(ll2[5])] = float(ll2[2])
        return tempdf
    except RuntimeError:
        logging.error("could not parse sketch file {}".format(file))



def _get_feature_list(nodesfile, noncellular):
    """Takes a NCBI taxonomy nodes.dmp file and a dict with high level
       noncellular categories and returns a list of taxids at the selected level.

    """
    with open(nodesfile, 'r') as nodes:
        phylist =[0]
        for line in nodes:
            linelist = line.strip().split()
            if linelist[4] in ["superphylum", "phylum", "subphylum"]:
                phylist.append(int(linelist[0]))
            feat_list = phylist + list(noncellular)
    return sorted(feat_list)


def _raise_taxdict_level(taxdict, taxlist, taxinstance):
    """Takes a dict in the form {taxid1: score1, taxid2: score2, ...} and
       returns dict with scores summed at the phylum level for cellular
       organisms or top ICVT level for viruses.

    """
    # pick highest scoring taxid
    print(taxdict)
    hi_score, hi_score_taxa = max(zip(taxdict.values(), taxdict.keys()))
    print(hi_score)
    print(hi_score_taxa)
    lineage = taxinstance.get_lineage(hi_score_taxa)
    inters = list(set(lineage).intersection(taxlist))
    print(inters)
    phyid = inters[0]
    newdict = {phyid: hi_score}
    return newdict




def _dict_to_csv(sketchdict, taxlist, outfile):
    """Takes a dictionary of data parsed from a minhash sketch file and returns a
    CSV file with the scores and each top level taxa category.

    """
    newdict = {}
    ncbi = NCBITaxa()
    for key, item in sketchdict.items():
        newdict[key] = _raise_taxdict_level(taxdict=item, taxlist=taxlist, taxinstance=ncbi)
    with open(outfile, 'w') as csvhandle:
        csv_writer_instance = csv.writer(csvhandle, lineterminator='\n')
        for key, item in newdict.items():
            line = []
            line.append(key)
            for i in taxlist:
                if i in item:
                    line.append(item[i])
                else:
                    line.append(0.)
            csv_writer_instance.writerow(line)


def minhashremote(dtemp, infile, outfile, server_url, nodesfile, noncellular):
    """Runs bbtools sendeketch.sh on a file of sequences returning a
        classification for each

    Calculates minhash sketches if (k=24, 31) for each sequence in a fasta
    and uses a remote server then returns the classification data to
    a file.

    Args:
        dtemp (str): a temporary to write intermediate files
        infile (str): A multi-sequence fasta file for which to calulate
            minhash taxonomy.
        outfile (str): a path to write the file retuned from the minhash
            server
        server_url (str): a URL for the minhash server
        nodesfile (str): a file in NCBI 'taxdump' nodes format containing
            the phyla super phyla and subphyla that should be used as
            classification categories for cellular organisms. A filtered
            version of the nodes files is in the package's data directory.
        noncellular (dict): a dictionary of taxid: names pairs containing
            the high level classifications for viruses.

    Returns:
        (str): The standard output from BBtools sendsketch.sh

    Notes:
        This function is prefered over `vica.minhash.minhashlocal` unless
        you need to sue taxfilter files to exclude training data when
        evaluating the performance of a custom trained classifier.

    See Also:
        vica.minhash.minhashlocal

    """
    sketchfile = os.path.join(dtemp,"sketchout.txt")
    logging.info("Using BBtools Sendsketch.sh to send minhash sketches to the server {}".format(server_url))
    sresult = _send_sketch(infile=infile, outfile=sketchfile, server_url=server_url)
    logging.info(sresult)
    logging.info("Parsing results file from BBtools Sendsketch.sh")
    sketchdict = _parse_sendsketch(sketchfile)
    taxlist = _get_feature_list(nodesfile=os.path.join(vica.DATA_PATH, nodesfile), noncellular=noncellular)
    _dict_to_csv(sketchdict, taxlist=taxlist, outfile=outfile)
