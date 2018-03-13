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
import collections
import shutil
import tempfile
import gzip
import subprocess
import io
import json

import pandas
import numpy
import pyfaidx
import ete3
import yaml

import vica

with open(vica.CONFIG_PATH) as cf:
        config = yaml.safe_load(cf)

class Split:
    """A class to create taxonomically distributed training data.

    The split class is used to sample refseq data formatted by BBtools. The
    methods hold summary dat necessary for sampling. for faster loading a
    precomputed summary file can be passed at initialization.

    Attributes:
        pyfaidx_obj (obj): An instance of the Pyfaidx class from the reference
            fasta file
        tax_instance (obj): An instance of the ETE NCBITaxa class
        profile (dict): a dict with the pyfaidx key as its key and the taxid
            and length as values
        test_subtree (dict): a dict of class with value list of subtree ids
        train_subtree (dict): a list of class with value list of subtree ids
        depth (int): the depth of the NCBI taxonomy tree at which to sample
        composition (dict): a summary of the number of samples taken at each
            taxonomic level for a selected sampling split_depth
        classes (list): a list of the lasses to be sampled
        n_samples (list): a list of the samples to be taken for each class.

    """

    def __init__(self, fasta_file, split_depth, classes, testfrac):
        """Initialization of the Split class

        Args:
            fasta_file (str): the path to the reference fasta file to use
            split_depth (int): the tree depth at which to split data between
                test and train.
            classes (dict): a dict of NCBI Taxonomy ids as keys and the number of
                samples to sample for each class
            testfrac (float): The proportion of the data to use for testing the model

        Returns:
            None
        """
        logging.info("Initializing vica.split_shred.Split object")
        logging.info("loading pyfaidx index, or creating index if not present")
        self.pyfaidx_obj =pyfaidx.Fasta(fasta_file, read_ahead=100)
        logging.info("Loading ete3 NCBI taxonomy data object")
        self.tax_instance = ete3.NCBITaxa()
        self.pruned_tree = None
        logging.info("Profiling sequences taxonomically")
        self.profile = self.set_profile(fasta_file)
        self.test_subtrees = None
        self.train_subtrees = None
        self.depth = split_depth
        self.composition = {}
        self.classes = classes
        self.testfrac = testfrac


    def _find_organelles(self):
        """ There are approximately twice as many organelle only genomes as complete
        eukaryotic genomes in refseq. we need to exclude these when building
        the training set so we do not oversample them.
        """
        organellepath = os.path.join(vica.DATA_PATH, config['split_shred']['organellefile'])
        with open(organellepath, 'r') as data_file:
            seqidlist = json.load(data_file)
        logging.info("{} organelle sequence accessions loaded".format(len(seqidlist)))
        return set(seqidlist)



    def set_profile(self, fasta_file):
        """load the profile from a file or create it from the fasta file

        Args:
            fasta_file
            use_profile_file (str): the string to a three column tab delimited file with
                key, taxonomy id and length for each contig in the reference fasta
            save_profile_file (str): location to save the  contig in the r
                reference data for faster profiling in the future
        Returns:
            profile (dict): Format: {tax_id: {seq_id: seq_length}, ...}
        """
        pdict = {}
        index_file = fasta_file + ".fai"
        # if profile not given, create profile
        logging.info("Profiling reference fasta file")
        organelles = self._find_organelles()
        errorcount = 0
        with open(index_file, 'r') as f:
            for n, line in enumerate(f):
                try:
                    if n % 1000000 == 0:
                        logging.info("processed {:,d} records".format(n))
                    rec = line.strip().split("\t")
                    namelist = rec[0].split("|")
                    seq_id = namelist[2]
                    tax_id = namelist[1]
                    length = int(rec[1])
                    if seq_id not in organelles:
                        if tax_id in pdict:
                            pdict[tax_id][seq_id] = length
                        else:
                            pdict[tax_id]={}
                            pdict[tax_id][seq_id] = length
                except:
                    errorcount += 1
                    continue
        if errorcount > 0:
            logging.info("{:,d} sequence entries out of {:,d} had errors".format(errorcount, n))
        return pdict


    def _test_or_train(self, tree):
        """Split internal nodes at the selected path into test and train given a subtree.

        Args:
            tree (obj): An ETE3 TreeNode object
            testfrac(float): the fraction of the data to use for testing

        Returns:
            (tuple): A list of test internal nodes and list of training internal nodes

        """
        top_nodes = []
        for node in tree.traverse():
            cn, depth = node.get_farthest_leaf()
            if depth == self.depth:
                top_nodes.append(node)
        testn = round(len(top_nodes) * self.testfrac)
        test_subtrees = list(numpy.random.choice(a=top_nodes, size = testn,
            replace=False))
        train_subtrees = list(set(top_nodes) - set(test_subtrees))
        return test_subtrees, train_subtrees


    def _assign_samples_attribute(self, n, nodelist):
        """for each node in nodelist, set split n samples evenly among them,
            randomly assigning modulo. Add a 'samples' attribute to each
            recording their own n.

        Args:
            n (int): the number of samples to take
            nodelist (list): a list of nodes to set 'samples' attribute on

        """
        whole_samples = n // len(nodelist)
        modulo = n % len(nodelist)
        if modulo != 0:
            nodes_with_extra_sample  = numpy.random.choice(a=nodelist, size = modulo,
                replace=False)
            for node in  nodes_with_extra_sample:
                node.add_features(samples=whole_samples + 1)
            normal_nodes = list(set(nodelist)- set(nodes_with_extra_sample))
            for node in normal_nodes:
                node.add_features(samples=whole_samples)
        else:
            for node in nodelist:
                node.add_features(samples=whole_samples)

    def _add_samples_feature_to_test_train_nodes(self, n, test_subtrees, train_subtrees):
        """split n samples up evenly among the test and train subtrees

        Args:
            n (int): the number of samples to take
            test_subtrees (list): a list of test nodes
            train_subtrees (list): a list of train nodes

        """
        testn = round(n * self.testfrac)
        trainn = n - testn
        self._assign_samples_attribute(testn, test_subtrees)
        self._assign_samples_attribute(trainn, train_subtrees)


    def _add_samples_feature_to_children(self, node):
        """for node with the 'samples' attribute add a 'samples' attribute to each
            child the number of samples evenly between them and randomly assigning
            modulo.

        Args:
            node (obj): a node containing the samples attribute and children

        """
        children = node.get_children()
        if len(children) > 0:
            self._assign_samples_attribute(node.samples, children)


    def _propagate_samples_feature_from_nodes_to_leaves(self, node):
        """Given a node at a split N samples across all subnodes down to the leaves

        """
        for node in node.traverse():
                self._add_samples_feature_to_children(node)



    def _calculate_tax_composition(self, taxalist):
        """Take a list of taxa nodes and return a tally of the
            taxonomic rank of the nodes.

            Args:
                taxalist(list}: A list of ETE3 Tree node objects
            Returns:
                (dict): a Counter dictionary of taxonomic ranks and the counts of those ranks
        """
        rankdict = collections.Counter()
        for node in taxalist:
            try:
                rank = self.tax_instance.get_rank([int(node.name)])
                currank = rank[int(node.name)]
                if currank == 'no_rank':
                    lineage_list = self.tax_instance.get_lineage(int(node.name))[::-1]
                    full_rank_dict = self.tax_instance.get_rank(lineage_list)
                    for i in lineage_list:
                        uprank = full_rank_dict[i]
                        if uprank != 'no_rank':
                            rankdict['below_' + uprank] += 1
                            break
                else:
                    rankdict[currank] += 1
            except ValueError:
                logging.exception("there is no node name for this node {}, skipping".format(node))
                continue
        return rankdict


    def split_test_train_nodes(self):
        """Split nodes into a test set and train set at the selected node
            height. Assign the number of samples to take to each node. Record
            the taxonomic level distribution for the selected height.

        """
        # Reset values
        self.test_subtrees = {}
        self.train_subtrees = {}
        self.composition = collections.Counter()
        # Prune the NCBI taxonomy tree to the leaves in the training dataset
        leaves = []
        for tax_id in self.profile:
            leaves.append(tax_id)
        self.pruned_tree = self.tax_instance.get_topology(list(set(leaves)),intermediate_nodes=True) # create ncbi topology tree
        # For each classification class process the data
        for key in self.classes:
            subtree = self.pruned_tree&str(key)
            # Split subtrees into test and train sets
            self.test_subtrees[key], self.train_subtrees[key] = self._test_or_train(subtree)
            n = self.classes[key]
            # split the samples among the subtrees
            self._add_samples_feature_to_test_train_nodes(n, self.test_subtrees[key], self.train_subtrees[key])
            # record the taxonomic levels the sampling occurred at
            comp_counter = self._calculate_tax_composition(self.test_subtrees[key] + self.train_subtrees[key])
            self.composition = self.composition + comp_counter
            # propagate the samples down to the leaves
            for node in self.test_subtrees[key] + self.train_subtrees[key]:
                self._propagate_samples_feature_from_nodes_to_leaves(node)




    def _writeseq(self, record, pos, length, handle):
        """writes a fasta sequence to a file, adding position information
            to the id.

        """
        seqlist = []
        label = (">" + record.name + "|pos|" + str(pos) + ".." +
            str(pos + length) + "\n")
        end = pos + length
        result = record[pos:end].seq
        seqlist = [result[n:n + 60] +'\n' for n in range(0, len(result), 60)]
        handle.write(label)
        handle.writelines(seqlist)

    def _select_random_segment(self, seqobj, name, length, tries=10, ns= 0.1):
        """Select a random fragment from sequence checking for short length and N's.

        """
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

    def _select_fragments_and_write(self, basedir, seq_length, test=True):
        """Select the random fragments and fasta files to a directory.
        """
        if test:
            subtrees =self.test_subtrees
            subdir = os.path.join(basedir, "test")
        else:
            subtrees = self.train_subtrees
            subdir = os.path.join(basedir, "train")
        if not os.path.exists(subdir):
            os.makedirs(subdir)
        for key in self.classes:
            subtree = subtrees[key]
            writefile = os.path.join(subdir, str(key) + ".fasta")
            with open(writefile, 'w') as outfile:
                for item in subtree:
                    for leaf in item.iter_leaves():
                        # make a list of seq_ids for the taxon
                        seq_ids=list(self.profile[leaf.name].keys())
                        full_seq_ids = ["tid|" + leaf.name + "|" + i for i in seq_ids]
                        # make an array of lengths for the contigs
                        length_array=numpy.array(list(self.profile[leaf.name].values()))
                        # select the contigs to sample based on their length.
                        sampling_list = numpy.random.choice(a=full_seq_ids,
                                            size=leaf.samples,
                                            replace=True,
                                            p=length_array/sum(length_array))
                        for i in sampling_list:
                            pos = self._select_random_segment(seqobj=self.pyfaidx_obj,
                                    name=i,
                                    length=seq_length,
                                    tries=10,
                                    ns= 0.1)
                            if pos:
                                self._writeseq(self.pyfaidx_obj[i],
                                    pos=pos,
                                    length=seq_length,
                                    handle=outfile)

    def write_sequence_data(self, directory, overwrite=False, seq_length=5000):
        """Write the training and test data to a directory

        Args:
            directory (str): The directory to write the training and test data
            overwrite (bool): Should the directory be overwritten if it exists
            seq_length (int): the length to sample the training space

        Returns:
            None, a directory is written

        """
        # if not self.test_subtrees and self.train_subtrees:
        #     raise Exception("Please run the method 'split_test_train_nodes() before writing sequence data")
        if os.path.exists(directory) and overwrite == False:
            raise Exception("File or directory {} exists and overwrite is set to False".format(directory))
        if os.path.exists(directory) and overwrite == True:
            if os.path.isfile(directory):
                os.remove(directory)
            if os.path.isdir(directory):
                shutil.rmtree(directory)
            os.makedirs(os.path.join(directory, "train"))
            os.makedirs(os.path.join(directory, "test"))
        else:
            os.makedirs(os.path.join(directory, "train"))
            os.makedirs(os.path.join(directory, "test"))
        self._select_fragments_and_write(basedir=directory, seq_length=seq_length, test=True )
        self._select_fragments_and_write(basedir=directory, seq_length=seq_length, test=False)



def run(infile, outdir, length, testfrac,
    split_depth, classes):
    """Run sample selection workflow

    Args:
        infile (str):A fasta containing the reference database (usually RefSseq)
        outdir (str): The directory to write the training data
        length (int): the length of the training fragments
        testfrac (float): the fraction of the data to put into the test set.
        split_depth (int): the depth above the leaf to split the trst and train data at.
        classes (dict):a dictionary of classes and the number of samples to collect from each class
    Returns:
        None, a directory is written

    """
    logging.info("Creating Split object.")
    data = vica.split_shred.Split(
        fasta_file=infile,
        split_depth=split_depth,
        classes=classes,
        testfrac=testfrac)
    logging.info("Dividing testing and training nodes.")
    data.split_test_train_nodes()
    logging.info("Writing data to the output directory.")
    data.write_sequence_data(outdir, overwrite=True, seq_length=length)
    logging.info("The distribution of taxonomic levels for split depth {} is {}.".format(data.depth,data. composition))
