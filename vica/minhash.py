#!/usr/bin/env python3
"""minhash.py: a module to run bbtools minhash function on a set of data and
return a file of tab delimited classification data"""

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
    """Runs bbtools sendsketch.sh on a file of sequences returning a classification for each"""
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
    """Runs bbtools sendsketch.sh on a file of sequences returning a classification for each"""
    options = ["comparesketch.sh",
               "in=" + infile,
               "out=" + outfile,
               "ref=" + ref,
               "blacklist=" + blacklist,
               "tree=" + tree,
               "taxfilter=" + taxfilter,
               "taxfilterlevel=" + taxfilterlevel,
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
               memory]
    sendsketchout = subprocess.run(options, stderr=subprocess.PIPE)
    return sendsketchout.stderr.decode('utf-8')
    #return sendsketchout

def _parse_sendsketch(file):
    """parses bbtools sendsketch output returning python dictionary"""
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
    """parses bbtools comparesketch output returning python dictionary"""
    try:
        tempdf = {}
        with open(file, 'r') as f:
            for line in f:
                if line.strip() == '':
                    next
                elif line.startswith("Query:"):
                    ll = line.strip().split("\t")
                    key1 = ll[6].split(":")[1].strip()
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
    return next((k for k, v in input_dict.items() if v == value), None)

def _pick_higher_level(taxid, taxinstance):
    """take a taxid and an ete3 taxonomy instance and returns a higher level taxid"""
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
    """takes a dict in the form {taxid1: score1, taxid2: score2, ...} and
       returns dict with scores summed at the phylum level for cellular oganisms or top ICVT level for viruses"""
    newdict ={}
    for key, item in taxdict.items():
        phyid = _pick_higher_level(taxid=key, taxinstance=taxinstance)
        if phyid in newdict:
            newdict[phyid] = newdict[phyid] + item
        else:
            newdict[phyid] = item
    return newdict



def _get_feature_list(nodesfile, noncellular):
    """takes a NCBI taxonomy nodes.dmp file and a dict with high level
       noncellular categories and returns a list of taxids at the selected level"""
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
    sketchfile = os.path.join(dtemp,"sketchout.txt")
    logging.info("Using BBtools Sendsketch.sh to send minhash sketches to the server {}".format(server_url))
    sresult = _send_sketch(infile=infile, outfile=sketchfile, server_url=server_url)
    logging.info(sresult)
    logging.info("Parsing results file from BBtools Sendsketch.sh")
    sketchdict = _parse_sendsketch(sketchfile)
    taxlist = _get_feature_list(nodesfile=os.path.join(vica.DATA_PATH, nodesfile), noncellular=noncellular)
    _dict_to_csv(sketchdict, taxlist=taxlist, outfile=outfile)
