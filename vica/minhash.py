#!/usr/bin/env python3
'''minhash.py: a module to run bbtools minhash function on a set of data and
return a file of tab delimited classification data'''

import subprocess
import os

# Constants

jgi_server_url='https://refseq-sketch.jgi-psf.org/sketch'

def send_sketch(file, outfile):
    '''Runs bbtools sendsketch.sh on a file of sequences returning a classification for each'''
    options = ["sendsketch.sh",
               "in=" + file,
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


def main():

    parser = argparse.ArgumentParser(description='A script to generate minhash compositions usung BBtools Sendsketch')
    parser.add_argument('--input', help="A multi-sequence fasta file")
    parser.add_argument('--output', help= "A tabular sketch file")
    parser.add_argument(])
    args = parser.parse_args()
    send_sketch(args.outfile, args.output)

if __name__ == '__main__':
    main()
