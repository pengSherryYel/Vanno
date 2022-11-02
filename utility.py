#!/usr/bin/env python
# author: sherry peng
# mail: xue.peng@helmholtz-muenchen.de
# date: 2021.12.6

import os
from subprocess import Popen


def mkdirs(dirname):
    '''
    makedirs
    '''
    if not os.path.exists(dirname):
        os.makedirs(dirname)


def checkEnv(sft):
    '''
    check software in the PATH
    '''
    cmd = "which %s" % sft
    status = Popen(cmd, shell=True)
    status.wait()
    if status:
        print("%s exist" % sft)
    else:
        print("Please add %s in your PATH" % sft)

def load_phrogs_des(phrogs):
    d = {}
    with open(phrogs) as f:
        for line in f:
           t = line.strip("\n").split("\t")
           pid,color,des,cati = t
           new_id = "phrog_%s"%pid
           d[new_id] = "Des:%s;Catigory:%s"%(des,cati)
    return d
