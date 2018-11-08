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

#-- required input files for the Retrieval System:
#   - S2 Surface response file(s)
#   - solar spectrum response file
for sfile in ['s2a.srf', 'ASTMG173.csv']:
    src = os.path.join(ss_dir_path, 'data', 'srf', sfile)
    dst = os.path.join(RS_INPUT_DIR, sfile)
    if not os.path.exists(src):
        msg = "------------------------------"
        msg += '\n'
        msg += " Retrieval System Setup failed: file ***{}*** provided by signature simulator was not found.".format(src)
        msg += '\n'
        msg += "------------------------------"
        print(msg)
        sys.exit(0)
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
