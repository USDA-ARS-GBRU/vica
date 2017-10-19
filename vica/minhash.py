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
                else:
                    ll2 = line.strip().split("\t")
                    tempdf[key1][int(ll2[5])] = float(ll2[2])
        return tempdf
    except IOError:
        print("could not parse of the  line {} of sketch file {}".format(i, file))


def taxid_to_phylum(taxdict, taxinstance):
    '''takes a dict in the form {taxid1: score1, taxid2: score2, ...} and
       returns dict with scores summed at the phylum level'''
    newdict ={}
    for key, item in taxdict.items():
        lineage = taxinstance.get_lineage(key)
        rank = taxinstance.get_rank(lineage)
        for key2, item2 in rank.items():
            if item2 in ["superphylum", "phylum", "subphylum"]:
                phyid = key2
                if phyid in newdict:
                    newdict[phyid] = newdict[phyid] + item
                else:
                    newdict[phyid] = item
                break
    return newdict



def get_phylum_list(nodesdmpfile):
    '''takes a NCBI taxonomy nodes.dmp file and returns a list of taxids at the superphylum, phylum or subphylum level'''
    with open(nodesdmpfile, 'r') as nodes:
        phylist =[]
        for line in nodes:
            ll = line.strip().split()
            if ll[4] in ["superphylum", "phylum", "subphylum"]:
                phylist.append(int(ll[0]))
    return sorted(phylist)


def dict_to_csv(sketchdict, phylumlist, csvfile):
    newdict = {}
    ncbi = NCBITaxa()
    for key, item in sketchdict.items():
        newdict[key]= taxid_to_phylum(item, ncbi)
    with open(csvfile, 'w') as csvhandle:
        csv_writer_instance = csv.writer(csvhandle)
        for key, item in newdict.items():
            line = []
            line.append(key)
            for i in phylumlist:
                if i in item:
                    line.append(item[i])
                else:
                    line.append(0.)
            csv_writer_instance.writerow(line)


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
    outfile = os.path.join(dtemp,"sketchout.txt")
    if args.local:
        refs = os.path.join(os.path.abspath(args.localsketchdir), "*.sketch")
        compare_sketch(infile=args.input, outfile=outfile, ref=refs, blacklist=args.blacklist, tree=args.tree)
    else:
        send_sketch(infile=args.input, outfile=outfile)

    sketchdict = parse_sketchout(outfile)
    phylumlist = get_phylum_list(args.nodes)
    dict_to_csv(sketchdict, phylumlist, args.output)

    shutil.rmtree(dtemp)


if __name__ == '__main__':
    main()
