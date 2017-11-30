#!/usr/bin/env python3
"""organelle_splitter.py: a script to randomly sleect reads from a file """

import argparse
from pyfaidx import Fasta
import random




#load fasta records and create index
def _read_data(file):
    """read a fasta and optionally an equivelantly named faidx and return a pyfaidx handle"""
    genes = Fasta(file, read_long_names=True, read_ahead=10000)
    return genes

def _shuffle_keys(ddict):
    """take the pyfaidx file and return a suffled list of the keys"""
    keylist = []
    for key in ddict.keys():
        keylist.append(key)
    kls = random.shuffle(keylist)
    return keylist

def _is_organelle(string):
    """parses the name and returns a tuple of lists conatining the taxid and rank from lowest to highest order"""
    ll = string.split(" ")
    name = 'll[0]'
    taxstringlist = ll[1].split(",")
    td = {}
    for item in taxstringlist:
        pair = item.split('=')
        td[pair[0]] = pair[1]
    if "organelle" in td:
        return True
    else:
        return False

def _write_record(record,handle):
    seqlist = []
    label = (">" + record.long_name + "\n")
    seqlist=(["".join(x)+"\n" for x in zip(*[iter(str(record[:].seq))] * 60)])
    handle.write(label)
    handle.writelines(seqlist)


# iterate through records
def split_records(seqrecord, keylist, organellefile, nuclearfile):
    with open(organellefile, 'w') as ofile, open(nuclearfile, 'w') as nfile:
        for key in keylist:
            record = seqrecord[key]
            org = _is_organelle(record.long_name)
            if org:
                _write_record(record, ofile)
            else:
                _write_record(record, nfile)

def split_all(genomesfile, organellefile, nuclearfile):
    seqrecord = _read_data(genomesfile)
    genekeys = _shuffle_keys(seqrecord)
    split_records(seqrecord, genekeys, organellefile, nuclearfile)

def main():

    parser = argparse.ArgumentParser(description='A script to generate a split a reftree fasta file into test and train as specific taxonomic levels')
    parser.add_argument('--infile', help="A fasta file to split")
    parser.add_argument('--organelleout', help="the test file in fasta format")
    parser.add_argument('--nuclearout', help="the train file in fasta format")
    args = parser.parse_args()

    split_all(genomesfile=args.infile, organellefile=args.organelleout, nuclearfile=args.nuclearout)


if __name__ == '__main__':
    main()
