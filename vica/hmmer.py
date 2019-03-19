"""A module to run hmmer's hmmsearch against a custom database parse the reports


"""

import subprocess
import os

import logging
import yaml
import json

import vica

with open(vica.CONFIG_PATH) as cf:
    CONFIG = yaml.safe_load(cf)

def _run_hmmsearch(hmmfile, seqfile, tblout):
    """Run hmmsearch and report the results

    Args:
        hmmfile (str): a binary hmmfile
        seqfile (str): a path to write the file retuned from the minhash
            server
        reportfile (str): a URL for the minhash server

    Returns:
        (str) The standard output from BBtools comparesketch.sh

    """

    options = ["hmmsearch",
               "--tblout",
               tblout,
               "-T", "30",
               "--noali",
               "--cpu", "2",
               hmmfile,
               seqfile]
    sendsketchout = subprocess.run(options, stdout= subprocess.PIPE, stderr=subprocess.PIPE)
    return sendsketchout.stderr.decode('utf-8')
    #return sendsketchout


def _parse_tblout(tblout):
    sampledict = {}
    with open(tblout, "r") as ifile:
        for line in ifile:
            if not line.startswith('#'):
                linelist = line.strip().split()
                hmm = linelist[2]
                source = linelist[0]
                source_nt = "".join(source.split("_")[:-1])
                if source_nt in sampledict:
                    sampledict[source_nt].append(hmm)
                else:
                    sampledict[source_nt] = [hmm]
    return sampledict


def _get_tokens(hmmfile):
    """reads hmmfile and returns a list of all the possible hmms

    """
    options = ["hmmstat",
               hmmfile]
    hmmstatout = subprocess.run(options, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    namelist = []
    for rec in  hmmstatout.stdout.decode('utf-8').splitlines():
        try:
            if rec == '': # remove blank lines
                continue
            if not rec.startswith('#'):
                linelist = rec.split()
                namelist.append(linelist[1])
        except Exception:
            raise

    return namelist






def _write_file(sampledict, report_file):
    with open(report_file, "w") as ofile:
        json.dump(sampledict, ofile)

def get_hmmer_features(dtemp, hmmfile, seqfile, outfile):
    """ Run and parde hmm results

    """
    tblout = os.path.join(dtemp, "tblout.txt")
    _run_hmmsearch(hmmfile, seqfile, tblout)
    sampledict = _parse_tblout(tblout)
    linelist = _get_tokens(hmmfile)
    datadict = {"tokens":linelist, "data":sampledict}
    _write_file(datadict, outfile)
