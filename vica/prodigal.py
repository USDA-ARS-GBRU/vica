"""prodigal.py:  a module with functions to call genes with Prodigal,
   count codon usage, calculate centered log ratio and isometric log ration
   transformations and return values as a CSV."""


import subprocess
import os
import logging
import csv
import yaml

import numpy as np
import scipy.linalg
import scipy.stats
from Bio import SeqIO
from collections import defaultdict


def clr(composition):
    """Calculates a centered log-ratio transformation from a list of values.

    Args:
        composition (list): a list of integers of floats containing the
            compositional data

    Returns:
        a list with the centered log-ratio transformed values

    References:
        Aitchison, J. (John), 2003. The statistical analysis of
        compositional data. Blackburn Press.

    """
    with np.errstate(divide='ignore', invalid='ignore'):
        a = np.array(composition)
        am =np.ma.masked_equal(a, 0)
        gm = scipy.stats.mstats.gmean(am)
        clrm = am/gm
        clrarray = np.ma.getdata(clrm)
    return list(clrarray)


def ilr(composition):
    """Calculates a isometric log-ratio transformation from a list of values.

    Args:
        composition (list): a list of integers of floats containing the
            compositional data

    Returns:
        a list with the isometric log-ratio transformed values. The
        length is len(composition - 1).

    References:
        Aitchison, J. (John), 2003. The statistical analysis of
        compositional data. Blackburn Press.

    """
    with np.errstate(divide='ignore', invalid='ignore'):
        clrlen= len(composition)
        clrarray = clr(composition)
        hmat = scipy.linalg.helmert(clrlen)
        ilrmat = np.inner(clrarray, hmat)
    return list(ilrmat)


def _call_genes(infile, outfile):
    """Runs Prodigal, calling genes.

    The function runs Prodigal, prokaryotic gene calling software, in
       metagenomic mode, saving its output.

    Args:
         infile (str): a multi-sequence fasta to call genes on.
         outfile (str): a Fasta file containing the called genestring

    Returns:
         (str): the Standard output of Prodigal

    References:
        Hyatt, D., Chen, G.-L., Locascio, P.F., Land, M.L., Larimer, F.W.,
        Hauser, L.J., 2010. Prodigal: prokaryotic gene recognition and
        translation initiation site identification.
        BMC Bioinformatics 11, 119. doi:10.1186/1471-2105-11-119

        https://github.com/hyattpd/Prodigal

    """
    options = ["prodigal",
               "-i", infile,
               "-p", "meta",
               "-d", outfile]
    callgenesout = subprocess.run(options, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    return callgenesout.stderr.decode('utf-8')


def _gene_to_codon(genestring):
    """Converts a DNA sequence string to a list of codons

    Args:
        genestring (str): A DNA sequence

    Returns:
        (list): A list containing the codons of the sequence

    """
    try:
        if len(genestring)>=3:
            f1 = [genestring[i:i+3] for i in range(0, len(genestring), 3)]
            if not len(f1[-1]) == 3:
                f1 = f1[:-1]
            return f1
    except:
        logger.exception("Warning: could not convert gene sequence to a list for codon counting")
        return []

def _codon_to_dict(genestring, offset):
    """counts codons in a gene string, with a reading frame offest returning
       codon counts as a dict.

    Args:
        genestring (str): A DNA sequence
        offset (int):  the starting point of the sequence, used to shift
            the reading frame

    Returns:
        (list): A list containing the codons of the sequence in the
            selected reading frame

    """
    assert offset in [0,1,2], "Offset must be 0, 1, or 2"
    framen = _gene_to_codon(genestring[offset:])
    cdict = {}
    for codon in framen:
        if not codon in cdict:
            cdict[codon] = 1
        else:
            cdict[codon] += 1
    return cdict



def _parse_prodigal_id_from_biopython(idval):
    """strips off prodigal gene annotations and returns the ID as it was in
        the contig file

    Args:
        idval (str): the ID value returned by Prodigal

    Returns: The ID value as fed to Prodigal

    """
    return '_'.join(str(idval).split('_')[:-1])

def count_dict_to_clr_array(count_dict, codon_list):
    """ Converts a count dictionary to a CLR list

    Takes a dictionary of counts where the key is the upper case codon,
    orders them by codon, and performs a centered log-ratio transformation
    returning a list.

    Args:
        count_dict (dict): a dictionary where codon is the key and the
            value is the count
        codon_list (list): A lexicographically sorted list of codons

    Returns:
        (list):  A vector of centered, log-ratio transformed values in
            ordered by the lexicographically sorted codons they correspond to.

    """
    output_list = []
    for i in codon_list:
        if i in count_dict:
            output_list.append(count_dict[i])
        else:
            output_list.append(0)
    return clr(output_list)

def count_dict_to_ilr_array(count_dict, codon_list):
    """ Converts a count dictionary to a ILR list

    Takes a dictionary of counts where the key is the upper case codon,
    orders them by codon, and performs a isometric log-ratio transformation
    returning a list.

    Args:
        count_dict (dict): a dictionary where codon is the key and the
            value is the count
        codon_list (list): A lexicographically sorted list of codons

    Returns:
        (list):  A vector of isometric log-ratio transformed values in
            ordered by the lexicographically sorted codons they correspond to.
            The length is len(codon_list - 1).

    """
    output_list = []
    for i in codon_list:
        if i in count_dict:
            output_list.append(count_dict[i])
        else:
            output_list.append(0)
    return ilr(output_list)

def dsum(*dicts):
    """Add up values in multiple dicts returning their sum.

    Args:
        *dicts (*awks): Dictionaries to summed

    Returns:
        (dict): a Dict with the summed values

    """
    ret = defaultdict(int)
    for d in dicts:
        for k, v in d.items():
            ret[k] += v
    return dict(ret)

def count_codon_in_gene(record, cdict={}):
    """Counts codons for all three frames in a gene.

   Takes a biopython sequence record and optionally a dict and
   returns a dict with the counts for the three codon frames adding
   them to the existing cdict if one was supplied.

   Args:
       record (obj): A Biopython sequence record object
       cdict (dict): A dictionary containing count data to be added to

   Returns:
       (dict): Counts of codons for the record, added to the optionally
           supplied count dictionary

    """
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


def count_codons(seqio_iterator, csv_writer_instance, codon_list):
    """Count codons from sequences in a BioIO seq iterator and
           write to a csv handle.

    Args:
        seqio_iterator (obj): A Biopython SeqIO iterator object
        csv_writer_instance (obj): A csv module file handle
        codon_list (list): a lexicographically sorted list of codons

    Returns:
        (int): the number of records writen

    """

    def record_line(idval, codon_dict, csv_writer_instance):
        """Combine ID and codon data from the three frames, writing to csv handle

        Args:

        """
        l0 = count_dict_to_ilr_array(codon_dict[0], codon_list)
        l1 = count_dict_to_ilr_array(codon_dict[1], codon_list)
        l2 = count_dict_to_ilr_array(codon_dict[2], codon_list)
        id_and_data = [idval]
        id_and_data.extend(list(np.concatenate((l0, l1, l2))))
        csv_writer_instance.writerow(id_and_data)

    last_base_id = None
    codon_dict = {}
    lc = 0
    for record in seqio_iterator:
        base_id = _parse_prodigal_id_from_biopython(record.id)
        if base_id == last_base_id:
                codon_dict = count_codon_in_gene(record=record, cdict=codon_dict)
        elif base_id is not last_base_id:
            if codon_dict != {}:
                record_line(idval=last_base_id, codon_dict=codon_dict, csv_writer_instance=csv_writer_instance)
                lc +=1
            codon_dict =count_codon_in_gene(record=record, cdict={})
            last_base_id = base_id
    if codon_dict != {}:
        record_line(idval=base_id, codon_dict=codon_dict, csv_writer_instance=csv_writer_instance)
        lc += 1
    return lc


def contigs_to_feature_file(infile, outfile, dtemp, codon_list):
    """for each contig in a file, count codons and write to csv"""
    genefile= os.path.join(dtemp, "genes.fasta")
    cgout = _call_genes(infile, genefile)
    logging.debug("From prodigal gene caller:")
    logging.debug(cgout)
    seqs = SeqIO.parse(genefile, 'fasta')
    with open(outfile, 'w') as csvfile:
        csv_writer_instance = csv.writer(csvfile, lineterminator='\n')
        lc = count_codons(seqio_iterator= seqs, csv_writer_instance=csv_writer_instance, codon_list= codon_list)
        logging.info("Wrote {} examples to the temporary file".format(lc))
