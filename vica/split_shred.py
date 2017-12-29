"""A module to process genomic training and test data formatted
    by BBtools.

    The module contains functions to split data into testing and training
    datasets at a desired taxonomic split level (family, order, etc.).
    It then selects a user-specified number of random contigs a desired length
    from the data.  It also returns a file of test taxids that are used to
    build a minhash database that excludes the taxa in the test dataset.
"""

import os
import logging
import random

import pandas
import numpy
import pyfaidx
import ete3
import yaml

import vica

config = yaml.safe_load(vica.CONFIG_PATH)

def _read_data(file):
    """read a fasta or bgzf fasta and optionally an equivelantly named faidx
        and return a pyfaidx handle.

    Args:
        file (str):A fasta file or blocked gzip format fasta file with names
            in the format "tid|<NCBI taxonomy ID>|<optional accession>"".
            Example: "tid|1026970|NW_008342263.1"

    Returns:
        (obj): A pyfaidx Fasta file object

    """
    seqobj = pyfaidx.Fasta(file, read_ahead=10000)#, read_long_names=True)
    return seqobj


def _shuffle_keys(ddict):
    """take the pyfaidx file and return a suffled list of the keys.

    """
    keylist = []
    for key in ddict.keys():
        keylist.append(key)
    random.shuffle(keylist)
    return keylist


def _profile_sequences(seqobj, ncbiobj, splitlevel, classes):
    """Collect data on all sequences in the reference db.
    Valid split levels are: 'species', 'genus', 'family', 'order','class',
        'phylum', 'kingdom'.

    """
    datadict = {}
    keylist = _shuffle_keys(seqobj)
    for key in keylist:
        try:
            rec = key.strip().split("|")
            tid = rec[1]
            length = len(seqobj[key])
            revlinlist = ncbiobj.get_lineage(tid)[::-1]
            rankdict = ncbiobj.get_rank(revlinlist)
            sltaxon = None
            cltaxon = None
            for item in revlinlist:
                if rankdict[item] == splitlevel:
                    sltaxon = item
                if item in classes:
                    cltaxon = item
            if sltaxon and cltaxon:
                datadict[key] = [tid, sltaxon, cltaxon, length]
        except Exception:
            logging.exception("An error occured while profiling the sequence {} in the reference database. Coninuing with the next sequence.".format(str(key)))
    df = pandas.DataFrame.from_dict(datadict, orient="index", dtype='int64')
    df.columns = ["taxid", "taxlevelid", "classid", "length"]
    df.astype(dtype={"taxid":"int64","taxlevelid": "int64","classid":"int64","length":"int64"})
    return df

def _split_levels(testfrac, df, classes):
    """Split the taxa at the selected level into test and train"""
    cd = {}
    for taxid in classes:
         dff = df[df.classid==taxid]
         if len(dff)>0:
             classids =set(dff['taxlevelid'])
             clength = len(classids)
             test = round(clength*testfrac)
             testids = set(numpy.random.choice(a=list(classids), size = test,
                 replace=False))
             trainids = set(classids) - testids
             cd[taxid]={'test':testids,'train': trainids, 'total':clength}
         else:
             logging.warning("the class {} was not present in the sequence file being rocessed".format(classes[taxid]))
    return cd

def _writeseq(record, pos, length, handle):
    """writes a fasta sequence to a file, adding position information
        to the id"""
    seqlist = []
    label = (">" + record.name + "|pos|" + str(pos) + ".." +
        str(pos + length) + "\n")
    end = pos + length
    result = record[pos:end].seq
    seqlist = [result[n:n + 60] +'\n' for n in range(0, len(result), 60)]
    handle.write(label)
    handle.writelines(seqlist)

def _select_random_segment(seqobj, name, length, tries=10, ns= 0.1):
    """Select a random fragment from sequence checking for short length and N's"""
    seq_length = len(seqobj[name])
    # verify the contig is long enough
    if seq_length > length:
        i = 0
        while i < tries:
            # select a random starting point
            pos = numpy.random.choice(seq_length - length)
            # calculate an end position
            endpos = pos + length
            # verify the reads is less than 10% N's
            if seqobj[name][pos: endpos].seq.count("N")/length < ns:
                return pos
            i += 1
    return None

def _read_taxid_from_fasta(outdir):
    """read taxids from completed fastas and write a file with the taxids"""
    species_set = set()
    indir = os.path.join(outdir,"test")
    for testfile in os.listdir(indir):
        with open(os.path.join(indir,testfile), "r") as testf:
            for line in testf:
                if line.startswith(">"):
                    speciesid = line.strip().split("|")[1]
                    species_set.add(speciesid)
    species_list = list(species_set)
    outfile = os.path.join(indir,"test_taxids.txt")
    with open(outfile, "w") as outf:
        outf.write("\n".join(species_list))
    return len(species_list)


def _process_examples(exampletype, n_per_class, cd, outdir, length, df, seqobj):
    """Sample genomes from the test/train split in each class, writing the
        desired number of fragments out to a directory.

    """
    tot_recs = 0
    for classid in cd:
        outfile = os.path.join(outdir, exampletype, str(classid) + ".fasta")
        with open(outfile, 'w') as outhandle:
            recs_written = 0
            #calculate the samples required to get an even
            # sampling from each taxa level
            rec_per_level = round(n_per_class/cd[classid]['total'])
            # For each selcted level identify all the unique species
            for ltid in cd[classid][exampletype]:
                ttdf = df[df.taxlevelid==ltid]
                species = list(set(ttdf["taxid"]))
                # select N  species at random (N = rec_per_level)
                ttvect = numpy.random.choice(a=species,size=rec_per_level)
                # for each species selected
                for item in ttvect:
                    #get records for the genome fragments from the species
                    speciesdf = ttdf[ttdf.taxid==item]
                    # Select genome fragment based on size
                    contigid = numpy.random.choice(a=speciesdf.index.values,
                        p=speciesdf["length"]/sum(speciesdf["length"]),size=1)[0]
                    pos = _select_random_segment(seqobj,contigid,length)
                    if pos:
                        _writeseq(seqobj[contigid], pos, length, outhandle)
                        recs_written += 1
                        tot_recs += 1
                        if tot_recs % 50000 == 0:
                            logging.info("{} total records processed".format(tot_recs))
        logging.info("Wrote {} fragmented sequences to the {} directory in the class {}".format(recs_written, exampletype, classid))
    return tot_recs




# >>> Left off here<<<
def _select_contigs(n_per_class, cd, outdir, length, df, seqobj):
    """select contigs for testing and training for each classifier class with
       even sampling in each taxon at the selected taxonomic level.
       Writes to the output directory.

    """
    # create output directory
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    # make output directories
    testdir = os.path.join(outdir,"test")
    if not os.path.exists(testdir):
        os.mkdir(testdir)
    traindir = os.path.join(outdir,"train")
    if not os.path.exists(traindir):
        os.mkdir(traindir)
    # for each test and train:
    testcount = _process_examples(exampletype='test', n_per_class=n_per_class,
        cd=cd, outdir=outdir, length=length, df=df, seqobj=seqobj)
    traincount = _process_examples(exampletype='train', n_per_class=n_per_class,
        cd=cd, outdir=outdir, length=length, df=df, seqobj=seqobj)
    logging.info("Wrote a total of {} testing and {} training fragmented sequences.".format(testcount, traincount))



def run(fastafile, outdir, length=5000, n_per_class=100000,
              testfrac =0.1, splitlevel="family",
              classes={2: "Bacteria",
                       2157: "Archaea",
                       2759: "Eukaryota",
                       10239: "Viruses"},
              configpath=vica.CONFIG_PATH):
    """Select randomly sampled sequence fragments for each class put them into test and train data sets.

    This selects a target number of sequences for each model training class.
    The length of the sequence fragments and the taxonomic level at which to
    split the test and train data can be selected.

    Args:
        fastafile (str): A fasta file or blocked gzip format fasta file with
            names in the format "tid|<NCBI taxonomy ID>|<optional accession>".
            Example: "tid|1026970|NW_008342263.1"
        outdir (str): A directory to write the output data
        length (int): the length of the sequence fragments sampled
        n_per_class (int): the number ot test and train sequences to sample.
            Sampling is probabilistic so this is a target not a gaurantee.
        testfrac (float): The fraction of the data to be in the test set
        splitlevel (str): the taxonomic level to split at:

            [“species”, “genus”, “family”, “order”, “class”, “phylum”, “kingdom”]

        classes (dict):  A dictionary of the classes to train the model on in
            the format  {taxid1: taxname1, taxid2: name2}
        configpath (str): A yaml configuration file

    Returns:
        None

    """
    try:
        global config
        config = yaml.safe_load(configpath)
        # Read data as pyfaidx object
        seqobj = _read_data(fastafile)
        ncbi = ete3.NCBITaxa()
        df = _profile_sequences(seqobj, ncbi, splitlevel, classes)
        cd = _split_levels(testfrac=testfrac, df=df, classes=classes)
        _select_contigs(n_per_class=n_per_class, cd=cd, outdir=outdir,length=length,df=df, seqobj=seqobj)
        testtaxa = _read_taxid_from_fasta(outdir=outdir)
        logging.info("Wrote {} NCBI taxomomy ids to the file 'test_taxids.txt'. This file isused to exclude test taxa from minhash during training".format(testtaxa))
    except:
        logging.exception("vica.split_train logged the following exception:")
