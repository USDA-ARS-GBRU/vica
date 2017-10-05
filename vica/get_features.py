#!/usr/bin/env python3
'''prodigal.py: a module to call genes with prodigal then count codon usage
and transform into centered log ratio returning values as a CSV'''

import subprocess
import os
import vica.shred
import vica.minhash
import vica.khmer_features
import vica.prodigal
import vica.tfrecord_maker
import glob

def main():
	parser = argparse.ArgumentParser(description='A script to take a fasta with contigs and extract fratures')
	parser.add_argument('--input', help="A fasta file")
	parser.add_argument('--output', help="path to tfrecord file written")
	parser.add_argument('--label', help="An integer representing the classification label", type='int')
	parser.add_argument('--shread', help="A flag to shred data into rangomly sampled contigs. Default is false", action="store_true")
	parser.add_argument('--length', help="The length of the genome subsamples if fixed is selected", default = 5000, type='int')
	parser.add_argument('--samples', help="Total number of shreded contigs to create, or if between 0 and 1, the proportion of the genome to sample", default = 0.5, type='float')
	parser.add_argument('--testing', help="Testing mode",  action="store_true")
	parser.add_argument('--minhashlocal', help="Use local minhash database instead of remote server",  action="store_true")
	parser.add_argument('--minhashrefs', help="local minhash database")
	parser.add_argument('--minhashtree', help="local minhash tree")
	parser.add_argument('--minhashblacklist', help="local minhash blacklist")
	args = parser.parse_args()

	#create workng directory and file names
	dtemp = tempfile.mkdtemp()
	segments = os.path.join(dtemp, "segments.fasta")
	refs = ",".join(glob.glob(os.path.join(args.minhashrefs, "*.sketch")))
	minhashout = os.path.join(dtemp,"minhashout.txt")
	kmerout = os.path.join(dtemp,"kmerout.csv")
	codonout = os.path.join(dtemp, "codonout.csv")

	# shred gneomes into contigs
	if args.shred:
		samples_frags = vica.shred.shred(fasta=args.input, shred="fixed",
						 samples=args.samples, length=args.length,
						 test=args.test)
	if samples_frags:
		SeqIO.write(samples_frags, segments, "fasta")
	else:
		segments = args.input

	# Extract minhash features
	if args.minhashlocal:
		vica.minhash.compare_sketch(infile=segments, outfile=minhashout,
						ref=refs, blacklist=args.blacklist,
						tree=args.tree)
	else:
		vica.minhash.send_sketch(infile=segments, outfile=minhashout)

	# Extract kmer features
	kmers = vica.khmer_features.iterate_kmer(ksize)
	vica.khmer_featureswrite_kmers_as_csv(infile=segments, outfile=kmerout,
							  ksize=4, kmers=kmers)
	# Extract codons
	vica.prodigal.contigs_to_feature_file(infile=segments, outfile=codonout)

	# Combine data into a Tensorflow TF record file
	convert_to_tfrecords(kmerfile=kmerout, codonfile=codonout,
			 minhashfile=minhashout, tfrecordfile=args.output,
			 label=args.label)

if __name__ == '__main__':
	main()
