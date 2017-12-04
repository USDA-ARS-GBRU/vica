
import pkg_resources

# Set data path
DATA_PATH = pkg_resources.resource_filename('vica', 'data/')
CONFIG_PATH = pkg_resources.resource_filename('vica', 'data/config_default.yml')

from .classify import *
from .get_features import *
from .khmer_features import *
from .minhash import *
from .prodigal import *
from .shred import *
from .split_shred import *
from .tfrecord_maker import *
from .vica_cli import *
# from .train_eval import *
