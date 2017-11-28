Vica: A classifier for the identification of RNA and RNA viruses in complex metagenomic samples
===============================================================================================

Authors
-------
* Adam R. Rivers, US Department of Agriculture, Agricultural Research Service
* Qingpeng Zhang, US Department of Energy, Joint Genome Institute
* Susannah Tringe, US Department of Energy, Joint Genome Institute

Introduction
------------

Vica is designed to identify highly divergent viruses and phage representing new
families or orders in assembled metagenomic and metatranscriptomic data. Vica
does this by combining information from across the spectrum of composition
to homology. The current version of Vica uses three feature sets (5-mers,
codon usage, and minhash sketches from long kmers (k=24,31). A Deep neural
network was trained on this data to identify virus fragment in both DNA
and RNA datasets.

Usage
-----

This package can classify assembled data and train new classification models.
Most users will only use the classification functionality in Vica. We provide
trained models for classifying contigs. This can be
easily invoked with the command::

   vica classify [options]

The package also has a suite of tools to train evaluate and use new
classification models. many of the workflows for doing this can be evoked with
the same sub command interface. For details see the Tutorial.

Requirements
------------

The package relies on a number of python dependencies that are resolved when
the package is installed with PIP.

The non-python dependencies are:

- Bbtools - https://jgi.doe.gov/data-and-tools/bbtools/
- Prodigal - https://github.com/hyattpd/Prodigal
- GNU Coreutils (standard on Mac and Linux systems) - http://www.gnu.org/software/coreutils/coreutils.html
