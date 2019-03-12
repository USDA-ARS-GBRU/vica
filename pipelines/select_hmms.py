#!/usr/bin/env python
"""A script to take the input of HMMer DOM report files and return a list of the
minimum number of hmms requires to achieve X% coverage. This is used for
selecting the hmms to use as features it reduces the search time required to
classify and train models"
"""

from collections import Counter
import argparse


def myparser():
    """Parse the input arguments for the classifier and return an args object"""
    parser = argparse.ArgumentParser(
        description="takes a list of HMMER DOM report files and returns smallest \
        subset the subset of hmms which are needed to achieve the selected read coverage")
    parser.add_argument("reports", nargs="+", help="report files to process")
    parser.add_argument("--out", required=True, help="the file to write")
    parser.add_argument("--cutoff", help="the cutoff value to use for coverage \
                        default 0.99", default=0.99, type=float)
    args = parser.parse_args()
    return args

def read_dom_files(argv):
    """Reads one or more dom files and parses the results
    returning a count dictionary and a dictionary of hmm models and the matching sequences

    """
    cnt = Counter()
    modelset = {}
    sourcelist = []
    for arg in argv:
        with open(arg, 'r') as fhandle:
            psource = None
            phmm = None
            for line in fhandle:
                if not line.startswith("#"):
                    ll = line.strip().split()
                    hmm = ll[3]
                    source = ll[0]
                    source_nt = "".join(source.split("_")[:-1])
                    cnt[hmm] += 1
                    sourcelist.append(source_nt)
                    if source == psource and hmm == phmm:
                        continue
                    if hmm in modelset:
                        modelset[hmm].append(source_nt)
                    else:
                        modelset[hmm] = [source_nt]
                    psource = source
                    phmm = hmm
    return modelset, len((set(sourcelist)))


# modellist = list(modelset.items())

def _check_list(testset, modellist):
    """ Given a set of reads matching, runs through an ordered test set
        identifying the set that creates the largest union with the testset

    Args:
        testset (tuple): a set object containing test sequences matched
        modellist (list): a list containing tuples in the form
        (hmm,[sequences,...]) sorted from largest to smallest number of
        matching an hmm

    Returns:
        n: the model test index of the item in modellist that most expands the set

    """
    bestlen = 0
    bestindex = None
    for i, item in enumerate(modellist):
        if len(item[1]) < bestlen:
            return bestindex
        tset = testset.union(set(item[1]))
        ulen = len(tset) - len(testset)
        # print("i: {}, ulen: {}, bestlen: {}".format(i, ulen, bestlen))
        if ulen > bestlen:
            bestlen = ulen
            bestindex = i
    return bestindex


def estimate_fraction_covered2(modelset, tot, cutoff=.99):
    """given a modellist return a list containing the hmms ranked with
        cumulative coverage
    Args:

    """
    #Sort the moeltest
    modellist = list(modelset.items())
    modellist.sort(key=lambda x: -len(x[1]))
    datalist = []
    testset = set(modellist.pop(0)[1])
    cov = 0
    i = 1
    while cov < cutoff:
        num = _check_list(testset, modellist)
        # if not n:
        #     return datalist
        next_element = modellist.pop(num)
        testset = testset.union(set(next_element[1]))
        cov = len(testset)/tot
        datalist.append((i, next_element[0], cov))
        i += 1

    return datalist

def main():
    """Run the cli

    """
    args = myparser()
    modelset, num = read_dom_files(args.reports)
    datalist = estimate_fraction_covered2(modelset=modelset, tot=num, cutoff=args.cutoff)
    with open(args.out, "w") as fhandle:
        fhandle.writelines("%s\n" % item[1] for item in datalist)

if __name__ == '__main__':
    main()
