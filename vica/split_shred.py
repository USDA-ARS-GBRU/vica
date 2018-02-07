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
import urllib
import gzip
import pandas
import numpy
import pyfaidx
import ete3
import yaml

import vica



class Split:
    """A class to create taxonomically distributed training data.

    The split class is used to sample refseq data formatted by BBtools. The
    methods hold summary dat necessary for sampling. for faster loading a
    precomputed summary file can be passed at initialization.

    Attributes:
        pyfaidx_obj (obj): An instance of the Pyfaidx class from the reference
            fasta file
        tax_instance (obj): An instance of the ETE NCBITaxa class
        rand_key_list (list): the keys in the Pyfaidx reference, randomized
        profile (dict): a dict with the yfaidx key as its key and the taxid
            and length as values
        test_subtree (dict): a dict of class with value list of subtree ids
        train_subtree (dict): a list of class with value list of subtree ids
        test_leaves (dict): a dict of class with value list of leaf ids
        train_leaves (dict): a list of class with value list of leaf ids
        depth (int): the depth of the NCBI taxonomy tree at which to sample
        composition (dict): a summary of the number of samples taken at each
            taxonomic level for a selected sampling split_depth
        classes (list): a list of the lasses to be sampled
        n_samples (list): a list of the samples to be taken for each class.

    """

    def __init__(self, fasta_file, profile_file, split_depth, classes, testfrac):
        """Initialization of the Split class

        Args:
            fasta_file (str): the path to the reference fasta file to use
            profile (str): the path to the summary data for the reference
                file (optional)
            split_depth (int): the tree depth at which to split data between
                test and train.
            classes (dict): a dict of NCBI Taxonomy ids as keys and the number of
                samples to sample for each class
            testfrac (float): The proportion of the data to use for testing the model

        Returns:
            None
        """

        self.pyfaidx_obj =pyfaidx.Fasta(fasta_file, read_ahead=100)
        self.tax_instance = ete3.NCBITaxa()
        self.pruned_tree = None
        self.rand_key_list = []
        for  key in self.pyfaidx_obj:
            self.rand_key_list.append(key)
            random.shuffle(self.rand_key_list)
        if profile_file:
            self.profile = set_profile(profile_file)
        else:
            self.profile = set_profile()
        self.test_subtrees = None
        self.train_subtrees = None
        self.test_leaves = None
        self.train_leaves = None
        self.depth = split_depth
        self.composition = {}
        self.classes = classes
        self.testfrac = None


    def _find_organelles():
        """ There are approximately twice as many organelle only genomes as complete
        eukaryotic genomes in refseq. we need to exclude these when building
        the training set so we do not oversample them.
        """
        seqidlist = []
        mito = 'ftp://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/*.1.genomic.fna.gz'
        plast = 'ftp://ftp.ncbi.nlm.nih.gov/refseq/release/plastid/*.1.genomic.fna.gz'
        dtemp = tempfile.mkdtemp()
        r = urllib.request.urlopen()
        for item in [mito, plastid]:
            options = ["wget", item, "-P", dtemp,]
        wgetout = subprocess.run(options, stderr=subprocess.PIPE)
        for seqfile in os.listdir(dtemp):
            with gzip.open(seqfile, 'rb') as f:
                for line in f:
                    if line.startswith(">"):
                        ll = line.split(" ")[0]
                        seqid = ll[1:]
                        seqidlist.append(seqid)
        shutil.rmtree(dtemp)
        return set(seqidlist)

        #return sendsketchout




    def set_profile(self, profile_file=None):
        """load the profile from a file or create it from the fasta file

        Args:
            profile_file (str): the string to a three column tab delimited file with
                key, taxonomy id and length for each contig in the reference fasta

        Returns:
            profile (dict): Format: {tax_id: {seq_id: seq_length}, ...}
        """
        pdict = {}
        try:
            with open(profile_file, 'r') as f:
                for line in f:
                    ll = f.strip().split('\t')
                    seq_id = ll[0]
                    tax_id = ll[1]
                    length = ll[2]
                    pdict[tax_id][seq_id] = length
            self.profile = pdict
        except:
            seqidset = _find_organelles()
            for seq_id in self.pyfaidx_obj:
                try:
                    rec = seq_id.strip().split("|")
                    tax_id = rec[1]
                    accession = rec[3]
                    length = len(self.pyfaidx_obj[seq_id])
                    with open(profile_file, 'w') as f:
                        if accession not in seqidset:
                            pdict[tax_id][seq_id] = length
                            f.writeline([seq_id, tax_id, length ])
                except Exception:
                    logging.exception("An error occurred while profiling the sequence {} in the reference database. Continuing with the next sequence.".format(str(key)))
                self.profile = pdict

    def _test_or_train(tree):
        """Split internal nodes at the selected path into test and train given a subtree.

        Args:
            tree (obj): An ETE3 TreeNode object
            testfrac(float): the fraction of the data to use for testing

        Returns:
            (tuple): A list of test internal nodes and list of training internal nodes

        """
        top_nodes = []
        top_test_nodes =  []
        top_train_nodes = []
        for node in tree.traverse():
            cn, depth = node.get_farthest_leaf()
            if depth == self.depth:
                top_nodes.append(node)
        testn = round(len(top_nodes) * self.testfrac)
        test_subtrees = numpy.random.choice(a=top_nodes, size = testn,
            replace=False)
        train_subtrees = list(set(top_nodes) - set(top_test_nodes))
        return test_subtrees, train_subtrees


    def _assign_samples_attribute(n, nodelist):
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
                node.add_feature(samples=whole_samples + 1)
            normal_nodes = list(set(children)- set(nodes_with_extra_sample))
            for node in normal_nodes:
                node.add_feature(samples=whole_samples)
        else:
            for node in nodelist:
                node.add_feature(samples=whole_samples)

    def _add_samples_feature_to_test_train_nodes(n, test_subtrees, train_subtrees):
        """split n samples up evenly among the test and train subtrees

        Args:
            n (int): the number of samples to take
            test_subtrees (list): a list of test nodes
            train_subtrees (list): a list of train nodes

        """
        testn = round(n * self.testfrac)
        trainn = n - testn
        _assign_samples_attribute(testn, test_subtrees)
        _assign_samples_attribute(trainn, train_subtrees)


    def _add_samples_feature_to_children(node):
        """for node with the 'samples' attribute add a 'samples' attribute to each
            child the number of samples evenly between them and randomly assigning
            modulo.

        Args:
            node (obj): a node containing the samples attribute and children

        """
        children = node.get_children()
        _assign_samples_attribute(node.samples, children)


    def _propagate_samples_feature_from_nodes_to_leaves(node):
        """Given a node at a split N samples across all subnodes down to the leaves

        """
        for node in tree.traverse():
            if not node.is_leaf():
                _add_samples_feature_to_children(node)

    def _calculate_tax_composition(taxalist):
        """Take a list of taxa nodes and return a tally of the
            taxonomic rank of the nodes.

            Args:
                taxalist(list}: A list of ETE3 Tree node objects
            Returns:
                (dict): a Counter dictionary of taxonomic ranks and the counts of those ranks
        """
        rankdict = collections.Counter()
        for node in taxalist:
            rank = self.tax_instance.get_rank([a.name])
            currank = rank[int(a.name)]
            if currank == 'no_rank':
                lineage_list = self.tax_instance.get_lineage(int(a.name))[::-1]
                full_rank_dict = self.tax_instance.get_rank(lineage_list)
                for i in lineage_list:
                    uprank = full_rank_dict[i]
                    if uprank != 'no_rank':
                        rankdict['below_' + uprank] += 1
                        break
            else:
                rankdict[currank] += 1
        return rankdict


    def split_test_train_nodes():
        """Split nodes into a test set and train set at the selected node
            height. Assign the number of samples to take to each node. Record
            the taxonomic level distribution for the selected height.

        """
        # Reset values
        self.test_subtrees = {}
        self.train_subtrees = {}
        self.test_leaves = {}
        self.train_leaves = {}
        self.composition = Counter()
        # Prune the NCBI taxonomy tree to the leaves in the training dataset
        leaves = []
        for tax_id in self.profile:
            leaves.append(tax_id)
        self.pruned_tree = self.tax_instance.get_topology(Set(leaves)) # create ncbi topology tree
        # For each classification class process the data
        for key in self.classes:
            subtree = self.pruned_tree&str(key)
            # Split subtrees into test and train sets
            self.test_subtrees[key], self.train_subtrees[key] = _test_or_train(subtree)
            n = self.classes[key]
            # split the samples among the subtrees
            _add_samples_feature_to_test_train_nodes(n, self.test_subtrees[key], self.train_subtrees[key])
            # record the taxonomic levels the sampling occurred at
            self.composition + _calculate_tax_composition(self.test_subtrees[key] + self.train_subtrees[key])
            # propagate the samples down to the leaves
            for node in self.test_subtrees[key] + self.train_subtrees[key]:
                _propagate_samples_feature_from_nodes_to_leaves(node)



    def _writeseq(record, pos, length, handle):
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

    def _select_random_segment(seqobj, name, length, tries=10, ns= 0.1):
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

    def _select_fragments_and_write(basedir, seq_length, test=True):
        """Select the random fragments and fasta files to a directory.
        """
        if test:
            leaf_list =self.test_leaves
            subdir = os.path.join(basedir,"test")
        else:
            leaf_list = self.train_leaves
            subdir = os.path.join(basedir,"train")
        for key in self.classes:
            r_leaf_list= random.shuffle(leaf_list[key])
            with open(os.path.join(subdir, str(key)), 'w') as outfile:
                for node in r_leaf_list:
                    # make a list of seq_ids for the taxon
                    seq_ids=list(self.profile[node.name].keys())
                    # make an array of lengths for the contigs
                    length_array=array(list(self.profile[node.name].values()))
                    # select the contigs to sample based on their length.
                    sampling_list = numpy.random.choice(a=seq_ids,
                                        size=node.samples,
                                        replace=True,
                                        p=length_array/sum(length_array)  )
                    for i in sampling_list:
                        pos = _select_random_segment(seqobj=self.pyfaidx_obj,
                                name=i,
                                length=seq_length,
                                tries=10,
                                ns= 0.1)
                        _writeseq(self.pyfaidx_obj,
                            pos=pos,
                            length=seq_length,
                            outfile=outfile)

    def write_sequence_data(directory, overwrite=False, seq_length=5000):
        """Write the training and test data to a directory

        Args:
            directory (str): The directory to write the training and test data
            overwrite (bool): Sould the directory be overwritten if it exists
            seq_length (int): the length to sample the training space

        Returns:
            None, a directory is written

        """
        if self.test_leaves or self.train_leaves == None:
            raise Exception("Please run the method 'split_test_train_nodes() before writing sequence data")
        if os.path.exists(directory) and overwrite == False:
            raise Exception("file or directory {} exists and overwrite is set to False".format(directory))
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
        _select_fragments_and_write(basedir=directory, seq_length=seq_length, test=True )
        _select_fragments_and_write(basedir=directory, seq_length=seq_length, test=False)
