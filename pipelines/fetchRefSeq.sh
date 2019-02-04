#!/bin/bash

#Written by Brian Bushnell
#Last updated July 19, 2018

#Modified  by Adam Rivers  on January 30, 2019

#Fetches and renames RefSeq.
#Be sure the taxonomy server is updated first, or run with local taxonomy data!
#To use this script outside of NERSC when using local taxonomy data,
#add "taxpath=/path/to/taxonomy_directory/" to each BBTools command


#Ensure necessary executables are in your path
#module load bbtools
module load pigz
module laod bbtools

# Set this based on where the taxonomy ws downloaded and processed
TAXTREE=taxonomy/tree.taxtree.gz

#Fetch RefSeq
#time wget -nv ftp://ftp.ncbi.nlm.nih.gov/refseq/release/complete/*genomic.fna.gz
wget -q -O - ftp://ftp.ncbi.nlm.nih.gov/refseq/release/complete/*genomic.fna.gz \
| gi2taxid.sh -Xmx2g in=stdin.fa.gz out=renamed.fa.gz pigz=32 unpigz zl=9 \
tree=$TAXTREE table=null accession=null server ow

#Sort by taxonomy.
#This makes sketching by taxa use much less memory because sketches can be written to disk as soon as they are finished.
sortbyname.sh -Xmx500g in=renamed.fa.gz memmult=0.33 out=sorted.fa.gz zl=8 pigz=32 taxa tree=$TAXTREE gi=ignore fastawrap=1023 minlen=60 readbufferlen=2 readbuffers=1

if [ $? -eq 0 ]
then
  rm renamed.fa.gz
fi
