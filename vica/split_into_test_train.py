#!/usr/bin/env python3
'''test_train_splitter.py'''

import argparse
from pyfaidx import Fasta
import random




#load fasta records and create index
def _read_data(file):
    '''read a fasta or bgzf fasta and optionally an equivelantly named faidx and return a pyfaidx handle'''
    genes = Fasta(file)#, read_long_names=True)
    return genes


def _shuffle_keys(ddict):
    '''take the pyfaidx file and return a suffled list of the keys'''
    keylist = []
    for key in ddict.keys():
        keylist.append(key)
    kls = random.shuffle(keylist)
    return keylist



def _parse_fancyname(string):
    '''parses the name and returns a tuple of lists conatining the taxid and rank from lowest to highest order'''
    ll = string.split(" ")
    name = 'll[0]'
    taxstring = ll[1].split(",")[0].split("=")[1]
    taxlist = taxstring.split("/")
    taxidlist = []
    taxranklist = []
    for item in taxlist:
        ll3 = item.split(":")
        taxid = int(ll3[0])
        if not ll3[1] == '':
            rank = int(ll3[1])
        else:
            rank = None
        taxidlist.append(taxid)
        taxranklist.append(rank)
    taxidlist.reverse()
    taxranklist.reverse()
    return zip(taxidlist, taxranklist)


def pick_taxid_at_level(level, string):
    '''picks a taxid at the selected level or then next highest ranked level'''
    taxzip = _parse_fancyname(string)
    for taxid, rank in taxzip:
        if rank:
            if rank > level:
                continue
            elif rank <= level:
                return taxid
        else:
            continue
    else:
        return None


def assign_to_test(train_seq_tot, test_seq_tot, testpct):
    '''decide if a new taxid whould go to test or train'''
    tot = float(test_seq_tot) + float(train_seq_tot)
    if tot > 0:
        currentpct = float(test_seq_tot)/(float(test_seq_tot) + float(train_seq_tot))
    else: currentpct = 0
    if currentpct <= testpct:
        v = True
    else:
        v = False
    return v

def _write_record(record,handle):
    seqlist = []
    label = (">" + record.long_name + "\n")
    seqlist=(["".join(x)+"\n" for x in zip(*[iter(str(record[:].seq))] * 60)])
    handle.write(label)
    handle.writelines(seqlist)


# iterate through records
def split_records(genes, genekeys, testout, trainout, testpct, taxlevel):
    train_set = set()
    test_set = set()
    test_record_count = 0
    train_record_count = 0
    train_seq_tot = 0
    test_seq_tot = 0
    with open(testout, "w") as test, open(trainout, "w") as train:
        for key in genekeys:
            record = genes[key]
            sequence = genes[key][:]
            taxid = pick_taxid_at_level(taxlevel, record.long_name)
            if taxid:
                if taxid in train_set:
                    _write_record(record, train)
                    length = sequence.end
                    train_seq_tot += length
                    train_record_count += 1
                elif taxid in test_set:
                    _write_record(record, test)
                    length = sequence.end
                    test_seq_tot += length
                    test_record_count += 1
                else:
                    if assign_to_test(train_seq_tot, test_seq_tot, testpct):
                        _write_record(record, test)
                        length = sequence.end
                        test_seq_tot += length
                        test_record_count += 1
                    else:
                        _write_record(record, train)
                        length = sequence.end
                        train_seq_tot += length
                        train_record_count += 1
            else:
                continue
    print("Wrote, {} sequences and {} nucleotides the test file {}.".format(test_record_count, test_seq_tot, testout))
    print("Wrote, {} sequences and {} nucleotides the train file {}.".format(train_record_count, train_seq_tot, trainout))
    print("The requested percent of test data was {:.2%}. {:.2%} percent of sequences and {:.2%} percent of nucleotides are in the test dataset".format(testpct, float(test_record_count)/float(test_record_count +train_record_count),
          float(test_seq_tot)/float(test_seq_tot + train_seq_tot)))


def split_all(genesfile, testout, trainout, testpct, taxlevel):
    genes = _read_data(genesfile)
    genekeys = _shuffle_keys(genes)
    split_records(genes, genekeys, testout, trainout, testpct, taxlevel)

def main():

    parser = argparse.ArgumentParser(description='A script to generate a split a reftree fasta file into test and train as specific taxonomic levels')
    parser.add_argument('--infile', help="A fasta file to split")
    parser.add_argument('--testout', help="the test file in fasta format")
    parser.add_argument('--trainout', help="the train file in fasta format")
    parser.add_argument('--testpct', help="the percent of data to keep in the test set")
    parser.add_argument('--taxlevel', help= "the taxonomic level to select on 1-24 [ 24= species, 16= family]")
    args = parser.parse_args()

    split_all(genesfile=args.infile, testout=args.testout, trainout=args.trainout, testpct=float(args.testpct), taxlevel=int(args.taxlevel))

if __name__ == '__main__':
    main()
