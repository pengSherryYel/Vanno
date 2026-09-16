#!/usr/bin/env python
# author: sherry peng
# mail: xue.peng@helmholtz-muenchen.de
# date: 2021.12.6

import os
from subprocess import Popen
import pandas as pd


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


def read_kegg_anno(ko_list_file):
    ## load kegg annotation
    f = pd.read_csv(ko_list_file, compression='gzip', header=0, sep='\t', quotechar='"')
    keggD = dict(zip(f['knum'], f['definition']))
    return keggD


def read_vog_anno(vog_annotation_file):
    ## load vog annotation
    f = pd.read_csv(vog_annotation_file, compression='gzip', header=0, sep='\t', quotechar='"')
    vogD_fc = dict(zip(f['#GroupName'], f['FunctionalCategory']))
    vogD_des = dict(zip(f['#GroupName'], f['ConsensusFunctionalDescription']))
    return vogD_fc, vogD_des