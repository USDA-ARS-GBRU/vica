"""A module with functions to calculate kmer frequency for short kmers and
    transform those into isometric log ratio compositions"""

import itertools
import csv
import logging

import khmer
from Bio import SeqIO
from Bio.Seq import Seq

import vica

__all__ = ["iterate_kmer", "get_composition", "run"]

def iterate_kmer(k):
    """Create a list of Kmers.

    Creates a list of kmers that have only the first kmer in the reverse
    complement pair. For example, in a list of 3-mers 'ATG' would be present
    but its reverse complement 'CAT' would not.

    Args:
        k (int): the kmer size between 4-8

    Returns:
        list: A list with lexographically sorted kmers

    """
    try:
        bases = ['A','C','T','G']
        kmers = [''.join(p) for p in itertools.product(bases, repeat=k)]
        core_kmer = []
        for kmer in kmers:
            if not str(Seq(kmer).reverse_complement()) in core_kmer:
                core_kmer.append(kmer)
        return core_kmer
    except:
        logging.exception("Could not calculate the list of kmers for k = {}".format(k))

def get_composition(ksize, seq, kmers, norm):
    """Calculate the kmer compsosition from a sequence.

    Count the kmers in a sequence and return a list of kmer counts or
    normalized kmer counts

    Args:
        ksize (int): Kmer size between 4-8
        seq (str): a string representing a DNA sequence
        kmers (list): a list of lexographically sorted kmers
        norm (bool): Should values be normalized by dividing by total

    Returns:
        list: A list of float values with the counts or proportions for
            each kmer in the list 'kmers'.

        """
    try:
        nkmers = 4**ksize
        tablesize = nkmers + 100
        counting_hash = khmer.Countgraph(ksize, tablesize, 1)
        counting_hash.consume(seq)
        composition = [counting_hash.get(kmer) for kmer in kmers]
        if norm == True:
            total = sum(composition)
            nc = []
            for item in composition:
                if item == 0:
                    nc.append(0.0)
                else:
                    nc.append(float(item)/float(total))
                composition = nc
        return [float(x) for x in composition]
    except:
        logging.exception("Could not calculate composition using khmer")


def _write_kmers_as_csv(infile, outfile, ksize, kmers):
    """Calculate ilr transformed kmer compositions for sequences in a fasta
       file.

   Takes a multi-sequence fasta file and a list of kmers and calculates the
   isometric log-ratio transformed kmer composition for each sequence,
   writing a CSV file with the data.

   Args:
       infile (str): a Fasta file
       outfile (str): a path to a CSV output file
       ksize (int): the kmer size, 4-8
       kmers (list): A list of the kmers to count

   Returns:
       None

   References:
       Aitchison, J. (John), 2003. The statistical analysis of
       compositional data. Blackburn Press.

    """

    try:
        with open(infile, 'r') as f1:
            with open(outfile, 'w') as csvfile:
                mywriter = csv.writer(csvfile, lineterminator='\n')
                header = ["id"]
                header.extend(kmers)
                # mywriter.writerow(header)
                ksize = int(ksize)
                kmers = iterate_kmer(ksize)
                rn = 0
                for record in SeqIO.parse(f1, 'fasta'):
                    rl = [record.id]
                    kmer_frequency = get_composition(ksize,str(record.seq).upper(), kmers, False)
                    kmer_ilr = vica.prodigal.ilr(kmer_frequency)
                    rl.extend(kmer_ilr)
                    mywriter.writerow(rl)
                    rn += 1
                logging.info("Wrote {} kmer records to {}.".format(rn, outfile))
    except:
        logging.exception("Could not write kmer profiles to file")

def run(infile, outfile, ksize):
    """Calculate the ilr transformed kmer composition for each sequence
        in a fasta file.

   Takes a multi-sequence fasta file and a kmer size and calculates the
   isometric log-ratio (ilr) transformed kmer composition for each sequence,
   writing a CSV file with the data.

   Args:
       infile (str): a Fasta file
       outfile (str): a path to a CSV output file
       ksize (int): the kmer size, 4-8

   Returns:
       None

   References:
       Aitchison, J. (John), 2003. The statistical analysis of
       compositional data. Blackburn Press.

    """

    kmers = iterate_kmer(ksize)
    logging.info("identifying kmer features with a k of {}".format(ksize))
    _write_kmers_as_csv(infile=infile, outfile=outfile, ksize=ksize, kmers=kmers)
