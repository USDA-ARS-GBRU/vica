Vica: Software to identify highly divergent DNA and RNA viruses and phages in microbiomes
=========================================================================================
.. image:: https://travis-ci.org/USDA-ARS-GBRU/vica.svg?branch=master
    :target: https://travis-ci.org/USDA-ARS-GBRU/vica

.. image:: https://codecov.io/gh/USDA-ARS-GBRU/vica/branch/master/graph/badge.svg
    :target: https://codecov.io/gh/USDA-ARS-GBRU/vica

.. image:: https://readthedocs.org/projects/vica/badge/?version=latest
    :target: http://vica.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

.. image:: https://api.codacy.com/project/badge/Grade/f39e8359ea334739842bba35e596cfdc
    :target: https://www.codacy.com/app/arivers/vica?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=USDA-ARS-GBRU/vica&amp;utm_campaign=Badge_Grade


Authors
-------
* Adam R. Rivers, US Department of Agriculture, Agricultural Research Service
* Qingpeng Zhang, US Department of Energy, Joint Genome Institute
* Susannah G. Tringe, US Department of Energy, Joint Genome Institute

Introduction
------------

Vica is designed to identify highly divergent viruses and phage representing new
families or orders in assembled metagenomic and metatranscriptomic data. Vica
does this by combining information from across the spectrum of composition
to homology. The current version of Vica uses three feature sets (5-mers,
codon usage in all three frames, and minhash sketches from long kmers (k=24,31).
The classifier uses a jointly trained deep neural network and logistic model
implemented in Tensorflow. The software is designed to identify  both DNA
and RNA viruses and phage in metagenomes and metatranscriptomes.

Models
------

The current leases does not include trained models but we will be adding them
in the future to allow for the rapid identification of viruses without model training.

Usage
-----

This package can classify assembled data and train new classification models.
Most users will only use the classification functionality in Vica. We will provide
trained models for classifying contigs in future releases. classification can be
easily invoked with the command::

   vica classify -infile contigs.fasta -out classifications.txt -modeldir modeldir

The package also has a suite of tools to prepare data, train and evaluate new
classification models. Many of the workflows for doing this can be evoked with
the same sub-command interface::

   vica split
   vica get_features
   vica train
   vica evaluate

For details see the Tutorial.

Requirements
------------

The package relies on a number of python dependencies that are resolved when
the package is installed with PIP.

The non-python dependencies are:

- Bbtools > v37.75- https://jgi.doe.gov/data-and-tools/bbtools/
- Prodigal > v2.6.3 - https://github.com/hyattpd/Prodigal
- GNU Coreutils - http://www.gnu.org/software/coreutils/coreutils.html

Documentation
-------------
Documentation for the package is at http://vica.readthedocs.io/en/latest/

Package availability
--------------------
- PyPi: https://pypi.python.org/pypi/vica
- Github: https://github.com/USDA-ARS-GBRU/vica


Copyright information
---------------------

ViCA Copyright (c) 2018, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required
approvals from the U.S. Dept. of Energy).  All rights reserved.

If you have questions about your rights to use or distribute this software,
please contact Berkeley Lab's Innovation & Partnerships Office at  IPO@lbl.gov
referring to " Viral Classification Algorithm Using Supervised Learning (LBNL
Ref 2017-125)."

NOTICE.  This software was developed under funding from the U.S. Department of
Energy.  As such, the U.S. Government has been granted for itself and others
acting on its behalf a paid-up, nonexclusive, irrevocable, worldwide license in
the Software to reproduce, prepare derivative works, and perform publicly and
display publicly.  The U.S. Government is granted for itself and others acting
on its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the
Software to reproduce, prepare derivative works, distribute copies to the
public, perform publicly and display publicly, and to permit others to do so.
