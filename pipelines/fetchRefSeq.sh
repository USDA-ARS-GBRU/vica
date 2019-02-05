#!/bin/bash
#SBATCH --job-name=refseqformat
#SBATCH --output=refseqformat.out
#SBATCH --error=refseqformat.err
#SBATCH --time=48:00:00
#SBATCH -p <queue-name>
#SBATCH -N 1
#SBATCH -n 40

#Written by Brian Bushnell
#Last updated July 19, 2018

#Modified  by Adam Rivers  on February 4, 2019

#Fetches and renames RefSeq. Compresses frreads using
#Be sure the taxonomy server is updated first, or run with local taxonomy data!
#To use this script outside of NERSC when using local taxonomy data,
#add "taxpath=/path/to/taxonomy_directory/" to each BBTools command


#Ensure necessary executables are in your path
module load bbtools
module load pigz
module load pbgzip
module laod samtools


# Set this based on where the taxonomy ws downloaded and processed
TAXTREE=/path/to/taxonomy/

#Fetch RefSeq and stream to gi2taxid for reformatting
wget -q -O - ftp://ftp.ncbi.nlm.nih.gov/refseq/release/complete/*.genomic.fna.gz \
| gi2taxid.sh -Xmx16g in=stdin.fa.gz out=renamed.fa.gz unpigz \
taxpath=$TAXTREE table=null accession=null server ow

#Sort by taxonomy and stream to pbgzip for compression which can be indexed by samtools
#This makes sketching by taxa use much less memory because sketches can be written to disk as soon as they are finished.
sortbyname.sh -Xmx1275g in=renamed.fa.gz memmult=0.33 out=stdout.fa \
 taxpath=$TAXTREE gi=ignore fastawrap=1023 minlen=60 readbufferlen=2 \
readbuffers=1 | pbgzip -c -n 32 -6 > sorted.fa.gz

if [ $? -eq 0 ]
then
  rm renamed.fa.gz
fi

#index with samtools
samtools faidx sorted.fa.gz
