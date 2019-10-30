
"""This script takes a list of tfecord files are creates a plot of the
taxonomic distribution of the sampled segments

"""

# Before running this script: Set up Xvfb graphics software so that QIIME2 can generate figures.
# Xvfb creates a virtual X session. This needs to be run in the background.
# If you encounter an error saying the session already exists, then use a
# different session number, e.g. "Xvfb :2"
# > Xvfb :1 -screen 0 800x600x16 &

# Add a Display variable to the local environment variables
# > export DISPLAY=:1.0


import tensorflow as tf
tf.enable_eager_execution()
import glob
from collections import Counter

import yaml
from ete3 import Tree, NCBITaxa, NodeStyle, TreeStyle
import vica

TFREC_DIR="*.tfrecord"

with open(vica.CONFIG_PATH) as cf:
    config = yaml.safe_load(cf)

def parser(record):
    keys_to_features = {"id": tf.FixedLenFeature((), dtype=tf.string),
                        "kmer": tf.FixedLenFeature([config["train_eval"]["kmerlength"]], dtype=tf.float32),
                        "codon": tf.FixedLenFeature([config["train_eval"]["codonlength"]], dtype=tf.float32),
                        "minhash": tf.VarLenFeature(dtype=tf.string),
                        "hmmer": tf.VarLenFeature(dtype=tf.string),
                        "label": tf.FixedLenFeature((), dtype=tf.int64)}
    return tf.parse_single_example(record, keys_to_features)


datacnt = Counter()
fnames2 = [f for f in glob.glob(TFREC_DIR)]
raw_dataset = tf.data.TFRecordDataset(fnames2)
parsed_dataset = raw_dataset.map(parser)

for parsed_record in parsed_dataset: # .take(1000):
    testdict = parsed_record
    taxid = (testdict['id']).numpy().decode('utf8').split("|")[1]
    if taxid in datacnt:
        datacnt[taxid] += 1
    else:
        datacnt[taxid] = 1




ncbi = NCBITaxa(dbfile=config["minhash"]["dbfile"])
tree = ncbi.get_topology(list(datacnt))

for n in tree.traverse():
    nstyle = NodeStyle()
    if n.name in datacnt:
        counts = datacnt[n.name]
        size = 1 + 2*((counts/3.14159)**(0.5))
        nstyle["size"] = size
        n.set_style(nstyle)

circular_style = TreeStyle()
circular_style.mode = "c" # draw tree in circular mode
circular_style.scale = 20

tree.render("traintree.pdf", w=183,  units="mm", tree_style=circular_style)

#tree.show(tree_style=ts)
