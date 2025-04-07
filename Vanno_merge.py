#!/usr/bin/python
# coding: utf-8

import argparse
import pandas as pd

anno = argparse.ArgumentParser(description="merge multi vanno summary")
def splitarg(string):
    return [i for i in string.split(",")]
anno.add_argument('-i', type=splitarg, required=True, help="path to for summary(sep must be comma). eg file1,file2")
anno.add_argument('-o', type=str, default="Vanno_merged.tsv", help="path to for output")

args = anno.parse_args()

fileList = args.i

def mergeVanno(fileList,optfile):
    filet = []
    for infile in fileList:
        print(infile)
        df = pd.read_csv(infile,sep="\t",header=0,index_col=0)
        #print(df)
        filet.append(df)
    finalDf = pd.concat(filet).fillna("-")
    print(finalDf)
    finalDf.to_csv(optfile, sep="\t")

#mergeVanno(["../../Vanno/prophage/Vanno_summary.tsv","../../Vanno/virus/Vanno_summary.tsv"])
mergeVanno(args.i, args.o)
