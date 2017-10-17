#!/usr/bin/env python3
'''prodigal.py: a module to call genes with prodigal then count codon usage
and transform into centered log ratio returning values as a CSV'''

import subprocess
import os
import numpy as np
from collections import defaultdict
import math
from Bio import SeqIO
import tempfile
import csv
import shutil
import scipy.linalg
import argparse
import scipy.stats
codon_list = ["TTT", "TCT", "TAT", "TGT",
              "TTC", "TCC", "TAC", "TGC",
              "TTA", "TCA", "TAA", "TGA",
              "TTG", "TCG", "TAG", "TGG",
              "CTT", "CCT", "CAT", "CGT",
              "CTC", "CCC", "CAC", "CGC",
              "CTA", "CCA", "CAA", "CGA",
              "CTG", "CCG", "CAG", "CGG",
              "ATT", "ACT", "AAT", "AGT",
              "ATC", "ACC", "AAC", "AGC",
              "ATA", "ACA", "AAA", "AGA",
              "ATG", "ACG", "AAG", "AGG",
              "GTT", "GCT", "GAT", "GGT",
              "GTC", "GCC", "GAC", "GGC",
              "GTG", "GCG", "GAG", "GGG"]



def clr(composition):
    '''calcualtes a centered log ratio transformation from a list of values'''
    a = np.array(composition)
    am =np.ma.masked_equal(a, 0)
    gm = scipy.stats.mstats.gmean(am)
    clrm = am/gm
    clrarray = np.ma.getdata(clrm)
    return list(clrarray)


def ilr(composition):
    '''claculates isometric log ratio transformation from list of values'''
    clrlen= len(composition)
    clrarray = clr(composition)
    hmat = scipy.linalg.helmert(clrlen)
    ilrmat = np.inner(clrarray, hmat)
    return list(ilrmat)


def call_genes(infile, outfile):
    '''Runs prodigal calling genes'''
    options = ["prodigal",
               "-i", infile,
               "-p", "meta",
               "-d", outfile]
    callgenesout = subprocess.run(options, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    return callgenesout.stderr.decode('utf-8')


def _gene_to_codon(genestring):
    '''Converts a DNA sequence string to a list of codons'''
    try:
        if len(genestring)>=3:
            f1 = [genestring[i:i+3] for i in range(0, len(genestring), 3)]
            if not len(f1[-1]) == 3:
                f1 = f1[:-1]
            return f1
    except ValueError:
        print("Warning: could not convert gene sequence to a list for codon counting")
        return []

def _codon_to_dict(genestring, offset):
    '''counts codons in a gene string, with a reading frame offest returning
       codon counts as a dict.'''
    assert offset in [0,1,2], "Offset must be 0, 1, or 2"
    framen = _gene_to_codon(genestring[offset:])
    cdict = {}
    for codon in framen:
        if not codon in cdict:
            cdict[codon] = 1
        else:
            cdict[codon] += 1
    return cdict

def get_id_list(fasta):
    '''extract the ids from a fasta file'''
    idlist = []
    with open(fasta, 'r') as f:
        for line in f:
            if line.startswith(">"):
                idlist.append(line.strip().split()[0])
    return idlist

def _parse_prodigal_id_from_biopython(id):
    '''strips off prodigal gene annotations and returns the id as it was in the contig file'''
    return '_'.join(str(id).split('_')[:-1])

def count_dict_to_clr_array(count_dict, codon_list):
    '''Takes a dictionary of counts where the key is the upper case codon,
       orders them by codon, and performs a clr transformation returning a list'''
    output_list = []
    for i in codon_list:
        if i in count_dict:
            output_list.append(count_dict[i])
        else:
            output_list.append(0)
    return clr(output_list)

def count_dict_to_ilr_array(count_dict, codon_list):
    '''Takes a dictionary of counts where the key is the upper case codon,
       orders them by codon, and performs a ilr transformation returning a list'''
    output_list = []
    for i in codon_list:
        if i in count_dict:
            output_list.append(count_dict[i])
        else:
            output_list.append(0)
    return ilr(output_list)

def dsum(*dicts):
    '''add up values in two dicts returning their sum'''
    ret = defaultdict(int)
    for d in dicts:
        for k, v in d.items():
            ret[k] += v
    return dict(ret)

def count_codon_in_gene(record, cdict={}):
    '''takes a biopython sequence record and optionally a defaultdict and
       returns a defaultdict with the counts for the three codon frames adding
       them to the existing default dict if one was supplied.'''
    seq = str(record.seq)
    d1 = {}
    d2 = {}
    for i in range(3):
        d1[i] = _codon_to_dict(genestring=seq, offset=i)
    for i in range(3):
        if i in cdict:
            d2[i] = dsum(cdict[i], d1[i])
        else:
            d2[i] = d1[i]
    return d2


def count_codons(seqio_iterator, csv_writer_instance):
    '''Count codons in from sequences in a BioIO seq iterator, and write to a csv handle'''

    def record_line(id, codon_dict, csv_writer_instance):
        '''combine id and codon data from the three frames, writing to csv handle'''
        l0 = count_dict_to_ilr_array(codon_dict[0], codon_list)
        l1 = count_dict_to_ilr_array(codon_dict[1], codon_list)
        l2 = count_dict_to_ilr_array(codon_dict[2], codon_list)
        id_and_data = [id]
        id_and_data.extend(list(np.concatenate((l0, l1, l2))))
        csv_writer_instance.writerow(id_and_data)

    print("running count_codons")
    last_base_id = None
    codon_dict = {}
    for record in seqio_iterator:
        base_id = _parse_prodigal_id_from_biopython(record.id)
        if base_id == last_base_id:
                codon_dict = count_codon_in_gene(record=record, cdict=codon_dict)
            # print("in loop 1")
        elif base_id is not last_base_id:
            if codon_dict != {}:
                record_line(id=last_base_id, codon_dict=codon_dict, csv_writer_instance=csv_writer_instance)
            codon_dict =count_codon_in_gene(record=record, cdict={})
            last_base_id = base_id
    if codon_dict != {}:
        record_line(id=base_id, codon_dict=codon_dict, csv_writer_instance=csv_writer_instance)


def contigs_to_feature_file(infile, outfile):
    '''for each contig in a file, count codons and write to csv'''
    dtemp = tempfile.mkdtemp()
    genefile= os.path.join(dtemp, "genes.fasta")
    call_genes(infile, genefile)
    seqs = SeqIO.parse(genefile, 'fasta')
    with open(outfile, 'w') as csvfile:
        csv_writer_instance = csv.writer(csvfile)
        count_codons(seqio_iterator= seqs, csv_writer_instance=csv_writer_instance)
    shutil.rmtree(dtemp)


def main():

    parser = argparse.ArgumentParser(description='A script to generate codon use frequency from Prodigal')
    parser.add_argument('--input', help="A multi-sequence fasta file")
    parser.add_argument('--output', help= "An output file of the clr transformed codon usage for frames 1, 2, and 3, in csv format")
    args = parser.parse_args()
    contigs_to_feature_file(infile=args.input, outfile=args.output)

if __name__ == '__main__':
    main()
