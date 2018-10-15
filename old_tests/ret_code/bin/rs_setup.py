#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os

#-- Retrieval System requires the subdirectory 'input'
#   to be present in its main directory
RS_INPUT_DIR = 'input'

try:
    import signaturesimulator as ss
    ss_dir_path = os.path.dirname(os.path.realpath(ss.__file__))
    srf2a_file = os.path.join(ss_dir_path,'data','srf','s2a.srf')
except ImportError as exc:
    msg = "------------------------------"
    msg += '\n'
    msg += " ! ! ! Signature Simulator installation not found on system ! ! !"
    msg += '\n'
    msg += " ! ! ! Retrieval System cannot be setup properly            ! ! !"
    msg += '\n'
    msg += "------------------------------"
    print(msg)
    sys.exit(1)

#-- S2 Surface response files must be available to
#   the Retrieval System in a predefined location
for sfile in ['s2a.srf']:
    src = os.path.join(ss_dir_path, 'data', 'srf', sfile)
    dst = os.path.join(RS_INPUT_DIR, sfile)
    #-- os.symlink does not work like 'ln -s'
    if os.path.exists(dst):
        os.remove(dst)
    os.symlink(src, dst)

msg =  "-"*30
msg += '\n'
msg += " Retrieval System Setup properly finished, Retrieval System ready for use!"
msg += '\n'
msg += "-"*30
print(msg)
