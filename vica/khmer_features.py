"""A module with functions to calculate kmer frequency for short kmers and
    transform those into isometric log ratio compositions"""

import itertools
import csv
import logging

import khmer
import skbio.stats.composition
from Bio import SeqIO
from Bio.Seq import Seq
import scipy


import vica

__all__ = ["iterate_kmer", "get_composition", "run"]


def iterate_kmer(k, rc=True):
    """Create a list of Kmers.

    Creates a list of kmers. If rc=True, only the first kmer in the reverse
    complement pair is kept. For example, in a list of 3-mers 'ATG' would be
    present but its reverse complement 'CAT' would not. 'CAT' would be added to
    the total for 'ATG'

    Args:
        k (int): the kmer size, 3,5,7
        rc (bool): Should reverse complement primers be used? Default True
    Returns:
        list: A list with lexographically sorted kmers

    """
    try:
        bases = ['A','C','G','T']
        kmers = [''.join(p) for p in itertools.product(bases, repeat=k)]
        if rc:
            core_kmer = []
            for kmer in kmers:
                if not str(Seq(kmer).reverse_complement()) in core_kmer:
                    core_kmer.append(kmer)
            return core_kmer
        return kmers
    except:
        logging.exception("Could not calculate the list of kmers for k = {}".format(k))

def get_composition(ksize, seq, kmers):
    """Calculate the kmer composition from a sequence.

    Count the kmers in a sequence and return a list of kmer counts

    Args:
        ksize (int): Kmer size [3,5,7]
        seq (str): a string representing a DNA sequence
        kmers (list): a list of lexographically sorted kmers

    Returns:
        list: A list of float values with the counts or proportions for
            each kmer in the list 'kmers'.

        """
    try:
        nkmers = 4**ksize
        tablesize = nkmers + 100
        counting_hash = khmer.Countgraph(ksize, tablesize, 4)
        counting_hash.set_use_bigcount(True)
        counting_hash.consume(seq)
        composition = [counting_hash.get(kmer) for kmer in kmers]
        return [float(x) for x in composition]
    except:
        logging.exception("Could not calculate composition using khmer")


def _write_kmers_as_csv(infile, outfile, ksize, kmers, ):
    """Calculate isometric log ratio transformed transformed kmer compositions for sequences in a fasta
       file.

   Takes a multi-sequence fasta file and a list of kmers and calculates the
   isometric log-ratio transformed kmer composition for each sequence,
   writing a CSV file with the data. Uses the Gram-Schmidt orthonormal basis.

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
    # The length id non-redundant kmer vectors for each k
    try:
        with open(infile, 'r') as f1:
            with open(outfile, 'w', buffering=16777216) as csvfile:
                mywriter = csv.writer(csvfile, lineterminator='\n')
                header = ["id"]
                header.extend(kmers)
                recnum = 0
                for record in SeqIO.parse(f1, 'fasta'):
                    rl = [record.id]
                    kmer_frequency = get_composition(ksize=int(ksize),
                                                     seq=str(record.seq).upper(),
                                                     kmers=kmers)
                    kmer_z = skbio.stats.composition.multiplicative_replacement(kmer_frequency)
                    kmer_ilr = skbio.stats.composition.ilr(kmer_z)
                    rl.extend(kmer_ilr)
                    mywriter.writerow(rl)
                    recnum += 1
                logging.info("Wrote {} kmer records to {}.".format(recnum, outfile))
    except:
        logging.exception("Could not write kmer profiles to file")

def run(infile, outfile, ksize, rc=True):
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

       J. J. Egozcue., “Isometric Logratio Transformations for
       Compositional Data Analysis” Mathematical Geology, 35.3 (2003)

    """

    kmers = iterate_kmer(k=ksize, rc=rc)
    logging.info("identifying kmer features with a k of {}".format(ksize))
    _write_kmers_as_csv(infile=infile, outfile=outfile, ksize=ksize, kmers=kmers)
