#!/usr/bin/env python3
'''minhash.py: a module to run bbtools minhash function on a set of data and
return a file of tab delimited classification data'''

import subprocess
import argparse
import os
from ete3 import NCBITaxa
import csv
import tempfile
import shutil
import logging

# Constants

jgi_server_url='https://refseq-sketch.jgi-psf.org/sketch'

def send_sketch(infile, outfile):
    '''Runs bbtools sendsketch.sh on a file of sequences returning a classification for each'''
    options = ["sendsketch.sh",
               "in=" + infile,
               "out=" + outfile,
               "address=" + jgi_server_url,
               "mode=sequence",
               "level=9",
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



def compare_sketch(infile, outfile, ref, blacklist, tree):
    '''Runs bbtools sendsketch.sh on a file of sequences returning a classification for each'''
    options = ["comparesketch.sh",
               "in=" + infile,
               "out=" + outfile,
               "ref=" + ref,
               "blacklist=" + blacklist,
               "tree=" + tree,
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
               "printname0=f",
               "printcontam=f",
               "printunique=t",
               "printnohit=f",
               "printtaxid=t"]
    sendsketchout = subprocess.run(options, stderr=subprocess.PIPE)
    return sendsketchout.stderr.decode('utf-8')
    #return sendsketchout

def parse_sketchout(file):
    '''parses bbtools sendsketch output returning python dictionary'''
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

superkingdom = {2:"Bacteria", 2157: "Archaea", 2759: "Eukaryota", 12884: "Viroids", 10239: "Viruses"}

noncellular = {39759: "Deltavirus",
               35237: "dsDNA viruses, no RNA stage",
               35325: "dsRNA viruses",
               1910928: "Genomoviridae",
               35268: "Retro-transcribing viruses",
               12877: "Satellites",
               29258: "ssDNA viruses",
               439488: "ssRNA viruses",
               1714266: "Virus families not assigned to an order",
               686617: "unassigned viruses",
               451344: "unclassified archaeal viruses",
               12333: "unclassified bacterial viruses",
               1922347: "unclassified RNA viruses",
               552364: "unclassified virophages",
               12429: "unclassified viruses",
               12884: "Viroids"}

def _find_key(input_dict, value):
    return next((k for k, v in input_dict.items() if v == value), None)

def pick_higher_level(taxid, taxinstance):
    try:
        '''take a taxid and an ete3 taxonomy instance and returns a higher level taxid'''
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
            for key, val in rank.items():
                if key in noncellular:
                    hightax = key
        else:
            hightax = 1
        return hightax
    except:
        if taxid == 0:
            return 0


def raise_taxdict_level(taxdict, taxinstance):
    '''takes a dict in the form {taxid1: score1, taxid2: score2, ...} and
       returns dict with scores summed at the phylum level for cellular oganisms or top ICVT level for viruses'''
    newdict ={}
    for key, item in taxdict.items():
        phyid = pick_higher_level(taxid=key, taxinstance=taxinstance)
        if phyid in newdict:
            newdict[phyid] = newdict[phyid] + item
        else:
            newdict[phyid] = item
    return newdict



def get_feature_list(nodesdmpfile, noncellular):
    '''takes a NCBI taxonomy nodes.dmp file and a dict with high level
       noncellular categories and returns a list of taxids at the selected level'''
    with open(nodesdmpfile, 'r') as nodes:
        phylist =[0]
        for line in nodes:
            ll = line.strip().split()
            if ll[4] in ["superphylum", "phylum", "subphylum"]:
                phylist.append(int(ll[0]))
        for key in noncellular:
                phylist.append(int(key))
    return sorted(phylist)


def dict_to_csv(sketchdict, taxlist, outfile):
    newdict = {}
    ncbi = NCBITaxa()
    for key, item in sketchdict.items():
        newdict[key]= raise_taxdict_level(taxdict=item, taxinstance=ncbi)
    with open(outfile, 'w') as csvhandle:
        csv_writer_instance = csv.writer(csvhandle, linrminator='\n')
        for key, item in newdict.items():
            line = []
            line.append(key)
            for i in taxlist:
                if i in item:
                    line.append(item[i])
                else:
                    line.append(0.)
            csv_writer_instance.writerow(line)

def minhashlocal(dtemp, infile, outfile, nodesfile, localsketchdir, blacklist, tree):
    sketchfile = os.path.join(dtemp,"sketchout.txt")
    refs = os.path.join(localsketchdir, "*.sketch")
    cresult = compare_sketch(infile=infile, outfile=sketchfile, ref=refs, blacklist=blacklist, tree=tree)
    logging.debug(cresult)
    logging.info("Parsing results file from BBtools Comparesketch.sh")
    sketchdict = parse_sketchout(sketchfile)
    phylumlist = get_phylum_list(nodes)
    dict_to_csv(sketchdict, phylumlist, outfile=outfile)

def minhashremote(dtemp, infile, outfile, nodesfile):
    sketchfile = os.path.join(dtemp,"sketchout.txt")
    logging.info("Using BBtools Sendsketch.sh to send minhash sketches to the server {}".format(jgi_server_url))
    sresult = send_sketch(infile=infile, outfile=sketchfile)
    logging.debug(sresult)
    logging.info("Parsing results file from BBtools Sendsketch.sh")
    sketchdict = parse_sketchout(sketchfile)
    taxlist = get_feature_list(nodesdmpfile=nodesfile, noncellular=noncellular)
    dict_to_csv(sketchdict, taxlist=taxlist, outfile=outfile)

def main():

    parser = argparse.ArgumentParser(description='A script to generate minhash compositions usung BBtools Sendsketch')
    parser.add_argument('--input', help="A multi-sequence fasta file")
    parser.add_argument('--output', help= "A tabular sketch file")
    parser.add_argument('--nodes', help= "An ncbi nodes taxonomy file")
    parser.add_argument('--local', action='store_true')
    parser.add_argument('--localsketchdir', help= "Directory containing sketches, if local flag is used", default=None)
    parser.add_argument('--blacklist', help= "Blacklisted kmers, if local flag is used", default=None)
    parser.add_argument('--tree', help= "Taxonomy tree, if local flag is used", default=None)

    args = parser.parse_args()
    dtemp = tempfile.mkdtemp()
    if args.local:
        minhashlocal(dtemp=dtemp, infile=args.input, outfile=args.output,
                     nodesfile=args.nodes, localsketchdir=args.localsketchdir,
                     blacklist=args.blacklist, tree=args.tree)
    else:
        minhashremote(dtemp=dtemp, infile=args.input, outfile=args.output,
                     nodesfile=args.nodes)

    shutil.rmtree(dtemp)


if __name__ == '__main__':
    main()
