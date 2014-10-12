# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function
try:
    # Python 3
    from configparser import ConfigParser
except ImportError:
    # Python 2
    from ConfigParser import ConfigParser
import os
from collections import OrderedDict
# the followings are necessary for parsing genecoder.ini
from genecoder.lab.codec import Coder_Cyclic, Coder_Convolution, coder_Iwadare, coder_SelfOrthogonal
from genecoder.lab.poly_gf4 import GF4
try:
    Coder_Cyclic
    Coder_Convolution
    coder_Iwadare
    coder_SelfOrthogonal
    GF4
except:
    pass

# === config ===
INI_FILE = os.path.join(os.path.dirname(__file__), 'genecoder.ini')

# === global variables ===
# CODERS = { coder_name: coder_object }
CODERS = OrderedDict


def load_coders():
    '''load coders from ini file.

    [coder]
    coder_name = (describe, n, k, coder object)

    '''
    if not os.path.exists(INI_FILE):
        raise Exception('{0} is not found.'.format(INI_FILE))
    config = ConfigParser()
    config.read(INI_FILE)

    section_name = 'coder'
    if section_name not in config.sections():
        raise Exception(('resource error. there is not the section "coder".'
                         'sections:{0}').format(config.sections()))
    a_coder_name_list = config.options(section_name)

    global CODERS
    CODERS = OrderedDict()
    for a_coder_name in a_coder_name_list:
        CODERS[a_coder_name] = eval(config.get(section_name, a_coder_name))

# === main ===
load_coders()
