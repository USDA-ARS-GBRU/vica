#!/bin/bash
# Setup_refseq_test_train.sh -- a script to download and format all data for
# training and evaluating the performance of a classifier
module load bbtools


# run a modified version of the fetch taxonomy script from bbtools
./fetchTaxonomy.sh
#fetch and process RefSeq data  using bbtools
./fetchRefSeq.sh


#Make a blacklist of kmers occuring in at least 300 different species.
time sketchblacklist.sh -Xmx63g in=$datadir/sorted.fa.bgzip.gz prepasses=1 tree=auto taxa taxlevel=species ow out=blacklist_refseq_species_300.sketch mincount=300 k=19,24

#Generate 31 sketch files, with one sketch per species.
time bbsketch.sh -Xmx63g in=$datadir/sorted.fa.bgzip.gz out=taxa#.sketch mode=taxa tree=auto accession=null gi=null files=31 ow unpigz minsize=400 prefilter autosize blacklist=blacklist_refseq_species_300.sketch k=19,24

#A query such as contigs.fa can now be compared to the new reference sketches like this:
# comparesketch.sh in=contigs.fa k=31,24 tree=auto taxa*.sketch blacklist=blacklist_refseq_species_300.sketch

#On NERSC systems, you can then set the default path to nt by pointing /global/projectb/sandbox/gaag/bbtools/refseq/current at the path to the new sketches.
#Then you can use the default set of nt sketches like this:
#comparesketch.sh in=contigs.fa refseq tree=auto
#That command automatically adds the default path to the sketches, the blacklist, and the correct values for K.
