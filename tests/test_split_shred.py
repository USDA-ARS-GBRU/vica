import tempfile
import os
import filecmp
import random
import shutil

import pandas
import numpy
import pyfaidx
from nose.tools import ok_, eq_
from ete3 import NCBITaxa

import vica.split_shred


classes = {2: 1000,
         2157: 1000,
         2759:1000,
         10239: 1000}



# def _process_samples():
