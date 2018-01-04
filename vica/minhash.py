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

    Calculates minhash sketches if (k=24, 31) for each sequence in a fasta
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



def _compare_sketch(infile, outfile, ref, blacklist, tree, taxfilter, taxfilterlevel, memory):
    """Runs bbtools compareketch.sh on a file of sequences returning a classification for each

    Calculates minhash sketches if (k=24, 31) for each sequence in a fasta
    and uses a local sketch database, then returns the classification data to
    a file.

    Args:
        infile (str): A multi-sequence fasta file for which to calculate
            minhash taxonomy.
        outfile (str): a path to write the file retuned from the minhash
            server
        ref (str):  A directory containing the .sketch reference files
        blacklist (str): A file containing blacklisted sketches which are to
            common to be informative
        tree (str): a path to a BBtools formatted taxonomy archive.
        taxfilter (str): a path to a file containing NCBI taxonomy IDs
            (one ID per line) to exclude from use in classification.
        taxfileter level (str): A taxonomic rank to use for taxonomic
            filtering. For example if a the taxid 562 (E. coli) is in the
            taxfiter file and the taxfilterlevel= 'genus' then all
            references in the genus Escherichia will be excluded.

            ["species", "genus", "family", "order", "class", "phylum", "kingdom"]

        memory (str): a Java style specification of the memory available
            for the process. for Example 16GB: "-Xmx16g". If not specified
            BBtools will autodetect.

    Returns:
        (str) The standard output from BBtools comparesketch.sh

    """
    options = ["comparesketch.sh",
               "in=" + infile,
               "out=" + outfile,
               "ref=" + ref,
               "blacklist=" + blacklist,
               "tree=" + tree,
               "mode=sequence",
               "k=31,24",
               "level=3",
               "color=f",
               "overwrite=t",
               "printani=f",
               "completeness=f",
               "score=t",
               "printmatches=t",
               "printgsize=f",
               "printgseqs=f",
               "printtaxname=f",
               "printname0=t",
               "printcontam=f",
               "printunique=t",
               "printnohit=f",
               "printtaxid=t",
               ]
    if memory:
        options.append(memory)
    if taxfilter:
        options.append("taxfilter=" + taxfilter)
    if taxfilterlevel:
        options.append("taxfilterlevel=" + taxfilterlevel)
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

def _parse_comparesketch(file):
    """Parses bbtools comparesketch output returning python dictionary.

    Args:
        file (str): a text file created by BBtools comparesketch.sh

    Returns:
        (dict): A dictionary with the ID as a key and as a value, a dict
        with taxid: score for each hit identified

    See Also:
        vica.minhash._parse_sendsketch
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


def _find_key(input_dict, value):
    """Finds key in a dictionary given a value

    Args:
        input_dict (dict): A dictionary
        value (str, int): a value to search for

    Returns:
        (str, int):  the key for the given entry

    """
    return next((k for k, v in input_dict.items() if v == value), None)

def _pick_higher_level(taxid, taxinstance):
    """ Returns the high level taxonomy class that a taxid belongs to.

    Takes a taxid and an ete3 taxonomy instance and returns a higher
    level taxid.

    Args:
        taxid (str): An NCBI taxonomy id
        taxinstance (obj): an ETE3  NCBI taxonomy object

    Returns:
        (int): The NCBI taxonomy id of the higher taxonomic rank

    """
    try:
        lineage = taxinstance.get_lineage(taxid)
        rank = taxinstance.get_rank(lineage)
        cellularlist = [2, 2157, 2759]
        noncellularlist = [12884, 10239]
        if set(lineage).intersection(cellularlist): # is it archaea, bacteria or euk?
            if "superphylum" in rank.values():
                hightax = _find_key(rank, "superphylum")
            elif "phylum" in rank.values():
                hightax = _find_key(rank, "phylum")
            elif "subphylum" in rank.values():
                hightax = _find_key(rank, "subphylum")
            else:
                hightax = 1
        elif set(lineage).intersection(noncellularlist): #is it virus or viroid?
            for key in rank:
                if key in noncellular:
                    hightax = key
        else:
            hightax = 1
        return hightax
    except:
        if taxid == 0:
            return 0


def _raise_taxdict_level(taxdict, taxinstance):
    """Takes a dict in the form {taxid1: score1, taxid2: score2, ...} and
       returns dict with scores summed at the phylum level for cellular
       oganisms or top ICVT level for viruses.

    """
    newdict ={}
    for key, item in taxdict.items():
        phyid = _pick_higher_level(taxid=key, taxinstance=taxinstance)
        if phyid in newdict:
            newdict[phyid] = newdict[phyid] + item
        else:
            newdict[phyid] = item
    return newdict



def _get_feature_list(nodesfile, noncellular):
    """Takes a NCBI taxonomy nodes.dmp file and a dict with high level
       noncellular categories and returns a list of taxids at the selected level.

    """
    with open(nodesfile, 'r') as nodes:
        phylist =[0]
        for line in nodes:
            ll = line.strip().split()
            if ll[4] in ["superphylum", "phylum", "subphylum"]:
                phylist.append(int(ll[0]))
        for key in noncellular:
                phylist.append(int(key))
    return sorted(phylist)


def _dict_to_csv(sketchdict, taxlist, outfile):
    """Takes a dictionary of data parsed from a minhash sketch file and returns a
    CSV file with the scores and each top level taxa category.

    """
    newdict = {}
    ncbi = NCBITaxa()
    for key, item in sketchdict.items():
        newdict[key]= _raise_taxdict_level(taxdict=item, taxinstance=ncbi)
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


def minhashlocal(dtemp, infile, outfile, ref, blacklist, tree, taxfilter, taxfilterlevel, memory, nodesfile, noncellular):
    """Runs bbtools compareketch.sh on a file of sequences returning a
        classification for each

    Calculates minhash sketches if (k=24, 31) for each sequence in a fasta
    and uses a local sketch database, then returns the classification data to
    a file.

    Args:
        dtemp (str): a temporary to write intermediate files
        infile (str): A multi-sequence fasta file for which to calculate
            minhash taxonomy.
        outfile (str): a path to write the file retuned from the minhash
            server
        ref (str):  A directory containing the .sketch reference files
        blacklist (str): A file containing blacklisted sketches which are to
            common to be informative
        tree (str): a path to a BBtools formatted taxonomy archive.
        taxfilter (str): a path to a file containing NCBI taxonomy IDs
            (one ID per line) to exclude from use in classification.
        taxfileter level (str): A taxonomic rank to use for taxonomic
            filtering. For example if a the taxid 562 (E. coli) is in the
            taxfiter file and the taxfilterlevel= 'genus' then all
            references in the genus Escherichia will be excluded.

            ["species", "genus", "family", "order", "class", "phylum", "kingdom"]
        memory (str): a Java style specification of the memory available
            for the process. for Example 16GB: "-Xmx16g". If None, BBtools
            will autodetect.
        nodesfile (str): a file in NCBI 'taxdump' nodes format containing
            the phyla super phyla and subphyla that should be used as
            classification categories for cellular organisms. A filtered
            version of the nodes files is in the package's data directory.
            This is the name of the file in the data directory not the path
            to the file.
        noncellular (dict): a dictionary of taxid: names pairs containing
            the high level classifications for viruses.

    Returns:
        (str): The standard output from BBtools comparesketch.sh

    Notes:
        Before running comparesketch.sh data must be prepared for bbtools
        comparesketch.sh. Pipelines for downloading and preparing the data
        are available in the pipeline directory. Refseq and NCBI taxonomy
        data must be downloaded and processed before running.

         the function `vica.minhash.minhashremote` creates the same
        classification files slightly faster without requiring data
        proprocessing.  The main reason to run run minhash locally is if
        you need to pass taxfilter files to it. This is primarily used
        for excluding training data when evaluating the performance of
        a custom trained classifier.

    See Also:
        vica.minhash.minhashremote

    """
    sketchfile = os.path.join(dtemp,"sketchout.txt")
    logging.info("Using BBtools Comparesketch.sh to identify matching taxa in the local database {}".format(str(ref)))
    cresult = _compare_sketch(infile=infile,
        outfile=sketchfile,
        ref= ref,
        blacklist=blacklist,
        tree=tree,
        taxfilter=taxfilter,
        taxfilterlevel=taxfilterlevel,
        memory=memory)
    logging.info(cresult)
    logging.info("Parsing results file from BBtools Comparesketch.sh")
    sketchdict = _parse_comparesketch(sketchfile)
    taxlist = _get_feature_list(nodesfile=os.path.join(vica.DATA_PATH, nodesfile), noncellular=noncellular)
    _dict_to_csv(sketchdict, taxlist=taxlist, outfile=outfile)

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
