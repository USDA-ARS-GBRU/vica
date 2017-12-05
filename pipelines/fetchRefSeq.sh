#!/bin/bash
#$ -swd
#$ -M adam.rivers@ars.usda.gov
#$ -m abe
## Specify a run time of 11 hours (or 39600 seconds)
#$ -l h_rt=12:00:00
#$ -l ram.c=120G
#Fetches and sketches RefSeq.
#Be sure taxonomy is updated first!
#To use this script outside of NERSC,
#add "taxpath=/path/to/taxonomy_directory/" to each BBTools command

#Ensure necessary executables are in your path
module load bbtools
module load pigz
module load samtools

# Fetch RefSeq
time wget -nv ftp://ftp.ncbi.nlm.nih.gov/refseq/release/complete/*genomic.fna.gz

#Concatenate into a single file
time cat *genomic.fna.gz > all.fa.gz

#Optionally, delete the old files
rm *genomic.fna.gz

#Rename by taxID by looking up gi numbers or accessions
time gi2taxid.sh -Xmx63g in=all.fa.gz out=renamed.fa.gz tree=auto table=auto accession=auto zl=6

# Remove all.fa.gz
rm all.fa.gz

#Sort by taxonomy and compress with bgzip for random access by pyfaidx
#This makes sketching by taxa use much less memory because sketches can be written to disk as soon as they are finished.
time sortbyname.sh -Xmx63g in=renamed.fa.gz out=stdout.fa taxa tree=auto gi=ignore fastawrap=255 minlen=60 | pbgzip -c -t 16 -6 > sorted.fa.bfgz


#remove all.fa.gz and reanmed.fa.gz
rm renamed.fa.gz
rm all.fa.gz
