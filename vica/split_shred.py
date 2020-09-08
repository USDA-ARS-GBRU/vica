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
    object holds summary data necessary for sampling and methods for sampling.

    Attributes:
        pyfaidx_obj (obj): An instance of the Pyfaidx class from the reference
            fasta file
        tax_instance (obj): An instance of the ETE NCBITaxa class
        pruned_tree(obj): An ete3.Tree object containing just the nodes with genomic data
        profile (dict): a dict with the pyfaidx key as its key and the taxid
            and length as values
        test_subtrees (dict): a dict of classes with a list of subtrees
        train_subtrees (dict): a list of classes with a list of subtrees
        depth (int): the taxonomic level of the NCBI taxonomy tree at which to sample e.g. 'order'
        composition (dict): a summary of the number of subtree nodes taken at each
            taxonomic level. On level is targeted but thih it the actual distribution.
        classes (dict): a dict of NCBI Taxonomy ids as keys and the number of
                samples to sample for each class
        testfrac (int): the fraction of samples to add to the test set
        ranks (dict): a dictionary mapping NCBI ranks to their sequencial order
        iranks (dict): the inverse of ranks


    """

    def __init__(self, fasta_file, split_depth, classes, testfrac):
        """Initialization of the Split class

        Args:
            fasta_file (str): the path to the reference fasta file to use
            split_depth (int): the taxonomic level at which to split data between
                test and train.
            classes (dict): a dict of NCBI Taxonomy ids as keys and the number of
                samples to sample for each class
            testfrac (float): The proportion of the data to use for testing the model

        Returns:
            None
        """
        logging.info("Initializing vica.split_shred.Split object")
        logging.info("loading pyfaidx index, or creating index if not present")
        self.pyfaidx_obj = pyfaidx.Fasta(fasta_file, read_ahead=100)
        logging.info("Loading ete3 NCBI taxonomy data object")
        self.tax_instance = ncbi = NCBITaxa(dbfile=config["minhash"]["dbfile"])
        self.pruned_tree = None
        logging.info("Profiling sequences taxonomically")
        self.profile = self.set_profile(fasta_file)
        self.test_subtrees = None
        self.train_subtrees = None
        self.depth = split_depth
        self.composition = {}
        self.classes = classes
        self.testfrac = testfrac
        self.ranks= {'no rank':None, 'superkingdom':0, 'kingdom':1, 'subkingdom':2,
                     'superphylum':3,'phylum':4, 'subphylum':5, 'superclass':6, 'class':7,
                     'subclass':8, 'infraclass':9,'superorder':10, 'order': 11,
                     'suborder':12, 'infraorder':13,'parvorder':14,'superfamily':15,
                     'family':16, 'subfamily':17, 'tribe':18,'subtribe':19,'genus':20,
                     'subgenus':21,'species group':22,'species subgroup':23,'species':24,
                     'subspecies':25,'varietas':26,'forma':27}
        self.iranks = {v: k for k, v in self.ranks.items()}



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
            for n, line in enumerate(f, 1):
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
        """Split internal nodes at the selected split_depth into test and train given a subtree.

        Args:
            tree (obj): An ETE3 TreeNode object
            testfrac(float): the fraction of the data to use for testing

        Returns:
            (tuple): A list of test internal nodes and list of training internal nodes

        """
        target_level = self.ranks[self.depth]
        top_nodes = []
        # iterate across all leaves
        for node in tree.iter_leaves():
            # get ordered lineage list
            ranklist = node.lineage
            # get dict mapping taxids to rank names
            rankdict = self.tax_instance.get_rank(node.lineage)
            # write id to top_nodes if below or at the target_taxa
            for id in node.lineage:
                named_rank = rankdict[id]
                if named_rank not in self.ranks:
                    continue
                num_rank = self.ranks[named_rank]
                if named_rank == 'no rank':
                    continue
                if num_rank < target_level:
                    continue
                elif num_rank == target_level:
                    top_nodes.append(tree&str(id))
                    break
                elif num_rank > target_level:
                    top_nodes.append(tree&str(id))
                    break
        unique_top_nodes = list(set(top_nodes))
        test_subtrees, train_subtrees = self._list_to_test_or_train(unique_top_nodes)
        return test_subtrees, train_subtrees

    def _list_to_test_or_train(self, top_nodes):
        """Split a list of internal nodes into test and train.

        Args:
            top_nodes(list):

        Returns:
            (tuple): A list of test internal nodes and list of training internal nodes

        """
        testn = round(len(top_nodes) * self.testfrac)
        test_subtrees = list(numpy.random.choice(a=top_nodes, size=testn,
            replace=False))
        train_subtrees = list(set(top_nodes) - set(test_subtrees))
        # check for  lists with no training or test data and move one node in if needed
        if len(train_subtrees) < 1 and len(test_subtrees) > 1:
            train_subtrees.append(test_subtrees[0])
            test_subtrees = test.subtrees[1:]
        elif len(test_subtrees) < 1 and len(train_subtrees) > 1:
            test_subtrees.append(train_subtrees[0])
            train_subtrees = train.subtrees[1:]
        return test_subtrees, train_subtrees

    def _weight(self, rank, target_rank):
        """Create a weight for the node

        The weight diminishes exponentially as the distance from the target level
        to the node level increases.

        Args:
            rank(str): the rank of the node, e.g. order
            target_rank(str): the targeted rank for dividing the data, e.g. genus
        """
        r0 = self.ranks[target_rank]
        r = self.ranks[rank]
        if not r:
            r = r0 + 2 # for node with no rank assume we went down 2 levels
        dif = r - r0
        level_penalty = 0.5
        weight = 1/(2**(level_penalty * dif))
        return weight

    def _assign_samples_attribute(self, n, depth, nodelist):
        """for each node in nodelist, set split n samples among them,
            with less samples given to nodes below the desired taxonomic level.
            Add a 'samples' attribute to each recording their own n.

        Args:
            n (int): the number of samples to take
            depth (str): the taxonomic level of the target or parent node e.g. 'order'
            nodelist (list): a list of nodes to set 'samples' attribute on

        """

        try:
            wvect = []
            for node in nodelist:
                wvect.append(self._weight(node.rank, depth))
            warray = numpy.array(wvect)
            x = n/sum(wvect)
            samplevect = list(x * warray)
            # round to whole numbers and give any node at least 1 sample
            samplelist = []
            for item in samplevect:
                if item < 1.0:
                    samplelist.append(1)
                else:
                    samplelist.append(int(round(item)))

            for i, node in enumerate(nodelist):
                node.add_features(samples=samplelist[i])

        except ZeroDivisionError:
            logging.warning("node list was empty")

    def _add_samples_feature_to_test_train_nodes(self, n, test_subtrees, train_subtrees):
        """split n samples up among the test and train subtrees

        Args:
            n (int): the number of samples to take
            test_subtrees (list): a list of test nodes
            train_subtrees (list): a list of train nodes

        """

        testn = round(n * self.testfrac)
        trainn = n - testn
        if testn < len(test_subtrees):
            logging.warning("The number of test samples requested, {} is less than the number of taxa at the requested taxonomic level, {}".format(testn, len(test_subtrees) ))
        self._assign_samples_attribute(testn, self.depth, test_subtrees)
        if trainn < len(train_subtrees):
            logging.warning("The number of training samples requested, {} is less than the number of taxa at the requested taxonomic level, {}".format(trainn, len(train_subtrees) ))
        self._assign_samples_attribute(trainn, self.depth, train_subtrees)

    def _add_samples_feature_to_children(self, node):
        """for node with the 'samples' attribute add a 'samples' attribute to each
            child the number of samples evenly between them and randomly assigning
            modulo.

        Args:
            node (obj): a node containing the samples attribute and children

        """
        children = node.get_children()
        noderank = node.rank
        # if the node has no rank estimate its position from the number of steps to a
        # parent node with a rank
        if noderank == 'no rank':
            no_rank_count = 0
            upstream_rank = None
            for id in node.lineage[::-1][1:]:
                ancestor = self.pruned_tree&str(id)
                if ancestor.rank == 'no rank':
                    no_rank_count += 1
                else:
                    upstream_rank = ancestor.rank
                    break
            nodeval = 1 * no_rank_count + (self.ranks[upstream_rank])
            if nodeval > 27:
                nodeval == 27
            noderank = self.iranks[nodeval]

        if len(children) > 0:
            self._assign_samples_attribute(node.samples, noderank, children)


    def _propagate_samples_feature_from_nodes_to_leaves(self, node):
        """Given a node at a split N samples across all subnodes down to the leaves

        """
        for node in node.traverse():
            try:
                self._add_samples_feature_to_children(node)
            except:
                logging.warning("Could not add samples feature to node {}. Continuing.".format(node))
                continue



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

                if node.rank == 'no_rank':
                    lineage_list = node.lineage[::-1]
                    full_rank_dict = self.tax_instance.get_rank(lineage_list)
                    for i in lineage_list:
                        uprank = full_rank_dict[i]
                        if uprank != 'no_rank':
                            rankdict['below_' + uprank] += 1
                            break
                else:
                    rankdict[node.rank] += 1
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
        try:
            self.pruned_tree = self.tax_instance.get_topology(list(set(leaves)),intermediate_nodes=True) # create ncbi topology tree
        except KeyError:
            logging.error("A taxonomy id in the sequence database could not be \
             found in the taxonomy database used by ETE3. Please update the \
             taxonomy database by removing the directory ~/.etetoolkit") \
        # For each classification class process the data
        for key in self.classes:
            logging.info("splitting testing and training nodes for class %s", key)
            subtree = self.pruned_tree&str(key)
            # Split subtrees into test and train sets
            self.test_subtrees[key], self.train_subtrees[key] = self._test_or_train(subtree)
            n = self.classes[key]
            # split the samples among the subtrees
            logging.info("Assigning initial sample values to top nodes in the class %s", key)
            self._add_samples_feature_to_test_train_nodes(n, self.test_subtrees[key], self.train_subtrees[key])
            # record the taxonomic levels the sampling occurred at
            comp_counter = self._calculate_tax_composition(self.test_subtrees[key] + self.train_subtrees[key])
            logging.info("the subnode level distribution for %s is %s" % (key, str(comp_counter)))
            self.composition = self.composition + comp_counter
            # propagate the samples down to the leaves
            logging.info(" propagating sample values down the tree for the class {}".format(key))
            for node in self.test_subtrees[key] + self.train_subtrees[key]:
                logging.info("propagating for node {}".format(node.sci_name))
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
            with open(writefile, 'w', 16777216) as outfile:
                for item in subtree:
                    for leaf in item.iter_leaves():
                        # make a list of seq_ids for the taxon
                        try: # in case leaf.name is not in self.profile
                            seq_ids=list(self.profile[leaf.name].keys())

                            full_seq_ids = ["tid|" + leaf.name + "|" + i for i in seq_ids]
                            # make an array of lengths for the contigs
                            length_array=numpy.array(list(self.profile[leaf.name].values()))
                            # select the contigs to sample based on their length.
                            logging.debug("writing {} segments for {}".format(leaf.samples, leaf.sci_name))
                            sampling_list = numpy.random.choice(a=full_seq_ids,
                                                size=int(leaf.samples),
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
                        except KeyError:
                            logging.info('leaf.name '+leaf.name+' is not in self.profile')

    def write_sequence_data(self, directory, overwrite=False, seq_length=5000, shuffle=True):
        """Write the training and test data to a directory optionally shuffle it

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
        if shuffle:
            for ctype in ["train", "test"]:
                for file in os.listdir(os.path.join(directory, ctype)):
                    tempfastafile = os.path.join(directory, ctype, "temp.fa")
                    fullpath = os.path.join(directory, ctype, file)
                    shutil.move(fullpath, tempfastafile)
                    runshuffle(tempfastafile, fullpath)
                    os.remove(tempfastafile)


def runshuffle(infile, outfile):
    """Runs bbtools shuffle2.sh to randomize the training segments

    Args:
         infile (str): a multi-sequence fasta
         outfile (str): a shuffled multi-sequence fasta

    Returns:
         (str): shuffle output

    """
    options = ["shuffle2.sh",
               "in=" + infile,
               "out=" + outfile]
    callgenesout = subprocess.run(options, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    return callgenesout.stderr.decode('utf-8')

def run(infile, outdir, length, testfrac,
    split_depth, classes):
    """Run sample selection workflow

    Args:
        infile (str):A fasta containing the reference database (usually RefSseq)
        outdir (str): The directory to write the training data
        length (int): the length of the training fragments
        testfrac (float): the fraction of the data to put into the test set.
        split_depth (int): The desired taxonomic level to split the test and train data
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
    data.write_sequence_data(os.path.abspath(outdir), overwrite=True, seq_length=length, shuffle=True)
    logging.info("The distribution of taxonomic levels for split depth {} is {}.".format(data.depth, data.composition))
