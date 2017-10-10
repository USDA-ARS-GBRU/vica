#!/usr/bin/env python3
'''minhash.py: a module to run bbtools minhash function on a set of data and
return a file of tab delimited classification data'''

import subprocess
import os

# Constants

jgi_server_url='https://refseq-sketch.jgi-psf.org/sketch'

def send_sketch(infile, outfile):
    '''Runs bbtools sendsketch.sh on a file of sequences returning a classification for each'''
    options = ["sendsketch.sh",
               "in=" + infile,
               "out=" + outfile,
               "address=" + jgi_server_url,
               "mode=sequence",
               "level=3",
               "color=f",
               "overwrite=t",
               "printani=f",
               "completeness=f",
               "score=t",
               "printmatches=t",
               "printgsize=f",
               "printgseqs=f",
               "printtaxname=f",
               "printname0=f",
               "printcontam=f",
               "printunique=t",
               "printnohit=f",
               "printtaxid=t"]
    sendsketchout = subprocess.run(options, stderr=subprocess.PIPE)
    return sendsketchout.stderr.decode('utf-8')
    #return sendsketchout



def compare_sketch(infile, outfile, ref, blacklist, tree):
    '''Runs bbtools sendsketch.sh on a file of sequences returning a classification for each'''
    options = ["comparesketch.sh",
               "in=" + infile,
               "out=" + outfile,
               "ref=" + ref,
               "blacklist=" + blacklist,
               "tree=" + tree,
               "k=31,24",
               "level=3",
               "color=f",
               "overwrite=t",
               "printani=f",
               "completeness=f",
               "score=t",
               "printmatches=t",
               "printgsize=f",
               "printgseqs=f",
               "printtaxname=f",
               "printname0=f",
               "printcontam=f",
               "printunique=t",
               "printnohit=f",
               "printtaxid=t"]
    sendsketchout = subprocess.run(options, stderr=subprocess.PIPE)
    return sendsketchout.stderr.decode('utf-8')
    #return sendsketchout

def main():

    parser = argparse.ArgumentParser(description='A script to generate minhash compositions usung BBtools Sendsketch')
    parser.add_argument('--input', help="A multi-sequence fasta file")
    parser.add_argument('--output', help= "A tabular sketch file")
    parser.add_argument('--local', action='store_true')
    parser.add_argument('--localsketchdir', help= "Directory containing sketches, if local flag is used", default=None)
    parser.add_argument('--blacklist', help= "Blacklisted kmers, if local flag is used", default=None)
    parser.add_argument('--tree', help= "Taxonomy tree, if local flag is used", default=None)

    args = parser.parse_args()

    if args.local:
        refs = os.path.join(os.path.abspath(args.localsketchdir), "*.sketch")
        compare_sketch(infile=args.input, outfile=args.output, ref=refs, blacklist=args.blacklist, tree=args.tree)
    else:
        send_sketch(infile=args.outfile, outfile=args.output)

if __name__ == '__main__':
    main()
