#!/usr/bin/env python3

from Bio import SeqIO
from subprocess import Popen, PIPE, STDOUT
import argparse
import sys
import logging
import os
import re
from threading import Thread
from collections import defaultdict
from utility import mkdirs, checkEnv, read_kegg_anno, read_vog_anno
from fmtparser import parse_hhsuit_hhr, parse_hmmsearch_opt, parse_m8_fmt_opt
import pandas as pd

################
## default #####
sys.path.append('/home/viro/xue.peng/script/module_annotation/Vanno/')

script_dir = str(os.path.dirname(os.path.abspath(__file__)))
wkdir = str(os.getcwd())

anno = argparse.ArgumentParser(description="function annotation for virus. current support: kegg, vog, pfam, phrog, swissprot. Usage: python annotation.py -i test.faa -k -v -p -r")
anno.add_argument('--version', action='version', version='Vanno v1.1')

##require
anno.add_argument('-i', type=str, required=True, help='input faa file')

##optional
anno.add_argument('-o', type=str, default="./Vanno_opt",
                  help="path to deposit output folder and temporary files, will create if doesn't exist [default= working directory]")
anno.add_argument('-t', type=int, default='1',
                  help='number of threads, each software occupies 1 CPU [default=1, max of 1 CPU per scaffold]')
#anno.add_argument('-virome', action='store_true',
#                  help='use this setting if dataset is known to be comprised mainly of viruses. More sensitive to viruses, less sensitive to false identifications [default=off]')
#anno.add_argument('-no_plot', action='store_true',
#                  help='suppress the generation of summary plots [default=off]')
anno.add_argument('-d', type=str, default=str(script_dir) + '/databases/',
                  help='path to original "databases" directory that contains .HMM files (if moved from default location)')

anno.add_argument('-k', '--kegg',action='store_true',dest="kegg",default=False, help="run kegg")
anno.add_argument('-kc', '--keggC',type=float, default="1e-5", dest="kc", help="kegg creteria. discard the not meet this creteria")
anno.add_argument('-kf', '--keggF',action='store_true',dest="kf",default=False, help="force rerun kegg")

anno.add_argument('-v', '--vog',action='store_true',dest="vog",default=False, help="run vog")
anno.add_argument('-vc', '--vogC',type=float, default="1e-5", dest="vc", help="vogdb creteria. discard the not meet this creteria")
anno.add_argument('-vf', '--vogF',action='store_true',dest="vf",default=False, help="force rerun vog")

anno.add_argument('-p', '--pfam',action='store_true',dest="pfam",default=False, help="run pfam")
anno.add_argument('-pc', '--pfamC',type=float, default="1e-5", dest="pc", help="pfam creteria. discard the not meet this creteria")
anno.add_argument('-pf', '--pfamF',action='store_true',dest="pf",default=False, help="force rerun pfam")

anno.add_argument('-r', '--phrog',action='store_true',dest="phrog",default=False, help="run phrog")
anno.add_argument('-rs', '--phrog_mode',choices=['hmmsearch','mmseqs','hhblits'],dest="phrog_mode",default="hhblits", help="run phrog use mmseqs|hmmsearch|hhblits")
anno.add_argument('-rc', '--phrogC',type=float, default="1e-5", dest="rc", help="phrog creteria. discard the not meet this creteria")
anno.add_argument('-rp', '--phrogProb',type=float, default="95", dest="rp", help="phrog creteria when using hhblits. discard the not meet this creteria")

anno.add_argument('-rf', '--phrogF',action='store_true',dest="rf",default=False, help="force rerun phrog")

anno.add_argument('-u', '--uniprot',action='store_true',dest="uniprot",default=False, help="run uniprot(default swiss-prot)")
anno.add_argument('-ud', '--uniprotDB',choices=['sprot', 'trembl', 'all'],default="sprot", help="run uniprot using sprot(swiss-prot); trembl  or all (sprot+trembl)")
anno.add_argument('-uc', '--uniprotC',type=float, default="1e-5", dest="uc", help="uniprot creteria. discard the not meet this creteria")
anno.add_argument('-uf', '--uniprotF',action='store_true',dest="uf",default=False, help="force rerun uniprot")

anno.add_argument('-b', '--pdb',action='store_true',dest="pdb",default=False, help="run pdb")
anno.add_argument('-bc', '--pdbC',type=float, default="1e-5", dest="bc", help="pdb creteria. discard the not meet this creteria")
anno.add_argument('-bf', '--pdbF',action='store_true',dest="bf",default=False, help="force rerun pdb")

#anno.add_argument('-l',type=str, nargs=1, default='1000',
#                  help='length in basepairs to limit input sequences [default=1000, can increase but not decrease]')
#anno.add_argument('-m', type=str, nargs=1, default=str(vibrant_path) + '/files/',
#                  help='path to original "files" directory that contains .tsv and model files (if moved from default location)')


args = anno.parse_args()
thread=args.t
input_faa=args.i
outputD = args.o
phrog_mode = args.phrog_mode
print("Results will be store at %s"%outputD)

######################
## perpare the dir
######################
if not os.path.exists(str(outputD)):
    Popen('mkdir -p ' + str(outputD) + ' 2>/dev/null', shell=True)
    print('mkdir -p', str(outputD))
logging.basicConfig(filename=os.path.join(str(outputD)+'anno.log'), level=logging.INFO, format='%(message)s')

## check database file
databases=args.d
kegg_db=os.path.join(databases,"kegg/KEGG_profiles_prokaryotes.HMM")
kegg_anno_file=os.path.join(databases,"kegg/ko_list.gz")

pfam_db=os.path.join(databases,"pfam/Pfam-A.hmm")

vog_db=os.path.join(databases,"vog/VOGDB_all.HMM")
vog_anno_file=os.path.join(databases,"vog/vog.annotations.tsv.gz")

if phrog_mode == "hmmsearch":
    phrog_db=os.path.join(databases,"phrog/all_phrogs.hmm")
    phrog_db_anno=os.path.join(databases,"phrog/phrog_annot.tsv")
elif phrog_mode == "mmseqs":
    phrog_db=os.path.join(databases,"phrog/phrogs_profile_db")
    phrog_db_anno=os.path.join(databases,"phrog/phrog_annot.tsv")
elif phrog_mode == "hhblits":
    phrog_db=os.path.join(databases,"phrog/phrogs_hhm.ffdata")
    phrog_db_anno=os.path.join(databases,"phrog/phrog_annot.tsv")


## unipriot and pdb unfinished, so comment out
# if args.uniprotDB == "sprot":
#     uniprot_db=os.path.join(databases,"uniprot_sprot.fasta")
# elif args.uniprotDB == "trembl":
#     uniprot_db=os.path.join(databases,"uniprot_trembl.fasta")
# elif args.uniprotDB == "all":
#     uniprot_db=os.path.join(databases,"uniprot_trembl_sprot.merge.fasta")
# else:
#     print("Wrong parameter")

# pdb_db=os.path.join(databases,"pdb_seqres.txt")
# pdb_db_anno=os.path.join(databases,"pdb_seqres.header.anno.txt")


def checkdb(db_file):
    if not os.path.exists(db_file):
        print("check %s! this file is not exist!"%db_file)
        exit()
    else:
        logging.info("using db_file: %s"%db_file)

if args.kegg:
    checkdb(kegg_db)

if args.pfam:
    checkdb(pfam_db)

if args.vog:
    checkdb(vog_db)

if args.phrog:
    checkdb(phrog_db)

if args.pdb:
    checkdb(pdb_db)


## run hmmer to annotation
thread = args.t

## step1 : split input faa file
print("Using %s threads" % thread)


#################################
#### Search Function ############
#################################

def runHmmsearch(inputfile, prefix, wd, hmmModel, otherPara="-T 40 --cpu 1"):
    '''
    Aim: run hmmer search for a pfam hmm database

    Usage: runHmmsearch(inputfile,prefix,wd,hmmModel,otherPara="--cpu 1")
        inputfile: a protein set from a metabin or a genome
        prefix: sample ientifier, will be used as a prefix to the output.
        wd: work path where put the result
        hmmModel: pfam *.hmm file
        otherPara: parameter for run the hmmer search.
            default: "--cpu 1"

    Return: output file path (*.tblout)
    '''

    checkEnv("hmmsearch")
    mkdirs(wd)
    cmd = "hmmsearch --noali {4} -o {2}/{1}.hmmsearch.out --tblout {2}/{1}.hmmsearch.tblout {3} {0}".format(
        inputfile, prefix, wd, hmmModel, otherPara)
    print("RUN command: %s\n" % cmd)
    obj = Popen(cmd, shell=True, stdout=PIPE, stderr=STDOUT)
    [logging.info(line.rstrip()) for line in obj.stdout]
    obj.wait()
    print("hmmsearch done!")
    return "%s/%s.hmmsearch.tblout" % (wd, prefix)


def runPhmmer(inputfile, prefix, wd, dbseq, otherPara="-T 40 --cpu 1"):
    '''
    Aim: run phmmer search for a sequence database

    Usage: runPhmmer(inputfile,prefix,wd,dbseq,otherPara="--cpu 1")
        inputfile: a protein set from a metabin or a genome
        prefix: sample ientifier, will be used as a prefix to the output.
        wd: work path where put the result
        dbseq: sequence database(not gzip) file
        otherPara: parameter for run the hmmer search.
            default: "--cpu 1"

    Return: output file path (*.tblout)
    '''

    checkEnv("phmmer")
    mkdirs(wd)
    cmd = "phmmer --noali {4} -o {2}/{1}.phmmer.out --tblout {2}/{1}.phmmer.tblout {0} {3}".format(
        inputfile, prefix, wd, dbseq, otherPara)
    print("RUN command: %s\n" % cmd)
    obj = Popen(cmd, shell=True, stdout=PIPE, stderr=STDOUT)
    [logging.info(line.rstrip()) for line in obj.stdout]
    obj.wait()
    print("phmmer done!")
    return "%s/%s.phmmer.tblout" % (wd, prefix)


def runHHblits(inputfile, prefix, wd, dbseq, otherPara="-n 1 -cpu 1 -e 0.001 -E 1e-5-p 90"):
    checkEnv("hhblits")
    mkdirs(wd)
    dbseq = os.path.dirname(os.path.realpath(dbseq))
    cmd  = "ffindex_from_fasta -s {2}/{1}.queryindex.multifasta.ff{{data,index}} {0} ".format(inputfile, prefix, wd, dbseq, otherPara)
    cmd2 = r"hhblits_omp -i {2}/{1}.queryindex.multifasta -d {3}/phrogs -o {2}/{1}.res.hhr -blasttab {2}/{1}.res.m8 {4} &&\
     tr -cd '\11\12\15\40-\176' < {2}/{1}.res.m8.ffdata > {2}/{1}.res.m8.plain_out.tsv".format(inputfile, prefix, wd, dbseq, otherPara)
    
    print("RUN command: %s\n" % cmd)
    obj = Popen(cmd, shell=True, stdout=PIPE, stderr=STDOUT)
    [logging.info(line.rstrip()) for line in obj.stdout]
    obj.wait()

    print("RUN command: %s\n" % cmd2)
    obj = Popen(cmd2, shell=True, stdout=PIPE, stderr=STDOUT)
    [logging.info(line.rstrip()) for line in obj.stdout]
    obj.wait()

    print("HHblits done!")
    return "%s/%s.res.hhr.ffdata.filtered_besthit.tsv" % (wd, prefix)


def runMMseqs(inputfile, prefix, wd, dbseq, otherPara="-s 6.5"):
    checkEnv("mmseqs")
    mkdirs(wd)
    dbseq = os.path.realpath(dbseq)
    cmd = "mmseqs createdb {0} {2}/{1}.query_seq ".format(inputfile, prefix, wd)
    cmd2 = "mmseqs search {2}/{1}.query_seq {3} {2}/{1}.results_mmseqs  ./tmp {4}".format(inputfile, prefix, wd, dbseq, otherPara)
    cmd3 = "mmseqs createtsv {3} {2}/{1}.query_seq {2}/{1}.results_mmseqs {2}/{1}.results.tsv".format(inputfile, prefix, wd, dbseq, otherPara)

    print("RUN command: %s\n" % cmd)
    obj = Popen(cmd, shell=True, stdout=PIPE, stderr=STDOUT)
    [logging.info(line.rstrip()) for line in obj.stdout]
    obj.wait()

    print("RUN command: %s\n" % cmd2)
    obj = Popen(cmd2, shell=True, stdout=PIPE, stderr=STDOUT)
    [logging.info(line.rstrip()) for line in obj.stdout]
    obj.wait()

    print("RUN command: %s\n" % cmd3)
    obj = Popen(cmd3, shell=True, stdout=PIPE, stderr=STDOUT)
    [logging.info(line.rstrip()) for line in obj.stdout]
    obj.wait()
    print("mmseqs done!")
    return "%s/%s.results.tsv" % (wd, prefix)



#################################
#### Merge Function ##############
#################################
def oneStepRun(inputfile, prefix, wd, db, outD, otherPara, creteria, hhblits_prob_cutoff=90, force=False, program="hmmsearch"):
    print("###### %s begin ######"%prefix)
    hmm_outPath = ""
    if program == "hmmsearch":
        hmm_outPath = "%s/%s.hmmsearch.tblout" % (wd, prefix)
        hmmModel = db
        if not os.path.exists(hmm_outPath) or force:
            hmm_outPath = runHmmsearch(inputfile, prefix, wd, hmmModel, otherPara=otherPara)
        else:
            print("Skip %s the running part, cause output file found!"%prefix)
        annoD = parse_hmmsearch_opt(hmm_outPath, creteria=creteria, reverse=False)

    # elif program == "phmmer":
    #     hmm_outPath = "%s/%s.phmmer.tblout" % (wd, prefix)
    #     dbseq = db
    #     if not os.path.exists(hmm_outPath) or force:
    #         hmm_outPath = runPhmmer(inputfile, prefix, wd, dbseq, otherPara=otherPara)
    #     else:
    #         print("Skip %s the running part, cause output file found!"%prefix)
    #     annoD = parse_hmmsearch_opt(hmm_outPath, creteria=creteria, reverse=True)

    # elif program == "mmseqs":
    #     mmseq_outPath = "%s/%s.mmseqs.results.tsv" % (wd, prefix)
    #     dbseq = db
    #     if not os.path.exists(mmseq_outPath) or force:
    #         mmseq_outPath = runMMseqs(inputfile, "%s_mmseqs"%prefix, "%s_mmseqs"%wd, dbseq, otherPara=otherPara)
    #     else:
    #         print("Skip %s the running part, cause output file found!"%prefix)
    #     annoD = parse_hmmsearch_opt(mmseq_outPath, creteria=1e-5, reverse=True)

    elif program == "hhblits":
        hhblits_outPath = "%s/%s.res.hhr.ffdata" % (wd, prefix)
        dbseq = db
        if not os.path.exists(hhblits_outPath) or force:
            hhblits_outPath = runHHblits(inputfile, prefix, wd, dbseq, otherPara=otherPara)
        else:
            print("Skip %s the running part, cause output file found!"%prefix)
        annoD = parse_hhsuit_hhr(hhblits_outPath, evalue_cutoff=creteria, prob_cutoff=hhblits_prob_cutoff)


    else:
        print("Wrong program parameter! please use hmmsearch|phmmer|mmseqs|hhblits")
        exit()

    outD[prefix] = annoD
    print("###### %s end ######"%prefix)


def split_dict_for_pandas(indict):
    outd = defaultdict(dict)
    for db,annos in indict.items():
        if db == "pfam":
            for query, values in annos.items():
                for target, accession in values.items():
                    ## because the pfam fmt is different. where descibe is in query_name col; and target id in accession
                    acc = "%s_acc"%db
                    des = "%s_des"%db
                    #print(db,query,name)
                    outd[des][query]= target
                    outd[acc][query]= accession
        else:
            ## for vog, kegg, phrogs using hmmer will add descrption below
            for query, values in annos.items():
                for target, info in values.items():
                    acc = "%s_acc"%db
                    outd[acc][query]= target

    return outd


###########################
#### main Programe ########
###########################

outD = {}
###########################  Run/Parse KEGG hmmsearch #########################
if args.kegg:
    keggOptD=os.path.join(outputD,"kegg")
    argsL = [input_faa, "kegg", keggOptD, kegg_db, outD]
    kwargsD = {"otherPara":"-T 40 --cpu %s"%(thread),
                "creteria":args.kc,
                "force":args.kf,
                "program":"hmmsearch"}

    keggt = Thread(target=oneStepRun,args=argsL, kwargs=kwargsD)
    keggt.start()
    #runHmmsearch(input_faa, "kegg", keggOptD, kegg_db, otherPara="-T 40 --cpu %s"%(thread))

###########################  Run/Parse VOG hmmsearch ##########################
if args.vog:
    vogOptD=os.path.join(outputD,"vog")
    argsL = [input_faa, "vog", vogOptD, vog_db, outD]
    kwargsD = {"otherPara":"-T 40 --cpu %s"%(thread),
                "creteria":args.vc,
                "force":args.vf,
                "program":"hmmsearch"}

    vogt = Thread(target=oneStepRun,args=argsL, kwargs=kwargsD)
    vogt.start()
    #runHmmsearch(input_faa, "vog", vogOptD, vog_db, otherPara="-T 40 --cpu %s"%(thread))


###########################  Run/Parse pfam hmmsearch ##########################
if args.pfam:
    pfamOptD=os.path.join(outputD,"pfam")
    argsL = [input_faa, "pfam", pfamOptD, pfam_db, outD]
    kwargsD = {"otherPara":"-T 40 --cpu %s"%(thread),
                "creteria":args.pc,
                "force":args.pf,
                "program":"hmmsearch"}

    pfamt = Thread(target=oneStepRun,args=argsL, kwargs=kwargsD)
    pfamt.start()
    #runHmmsearch(input_faa, "pfam", pfamOptD, pfam_db, otherPara="-T 40 --cpu %s"%(thread))

###########################  Run/Parse PHROG hmmsearch|mmseqs|hhblits ##########################
if args.phrog:

    if args.phrog_mode == "hmmsearch":
        phrogOptD=os.path.join(outputD,"phrog")
        argsL = [input_faa, "phrog", phrogOptD, phrog_db, outD]
        kwargsD = {"otherPara":"-T 40 --cpu %s"%(thread),
                    "creteria":args.rc,
                    "force":args.rf,
                    "program":args.phrog_mode}

    elif args.phrog_mode == "mmseqs":
        phrogOptD=os.path.join(outputD,"phrog_mmseqs")
        argsL = [input_faa, "phrog_mmseqs", phrogOptD, phrog_db, outD]
        kwargsD = {"otherPara":"-s 7 --threads %s"%(thread),
                    "creteria":args.rc,
                    "force":args.rf,
                    "program":args.phrog_mode}

    elif args.phrog_mode == "hhblits":
        phrogOptD=os.path.join(outputD,"phrog_hhblits")
        argsL = [input_faa, "phrog_hhblits", phrogOptD, phrog_db, outD]
        kwargsD = {"otherPara":"-n 1 -e 0.001 -E 1e-5 -p 90 -cpu %s "%(thread),
                    "creteria":args.rc,
                    "hhblits_prob_cutoff":args.rp,  ## only for hhblits
                    "force":args.rf,
                    "program":args.phrog_mode}
    else:
        print("phrogs search mode should be hmmsearch|hhblits|mmseqs")
    phrogt = Thread(target=oneStepRun,args=argsL, kwargs=kwargsD)
    phrogt.start()

###########################  Run/Parse Uniprot(Swiss-prot) phmmer ##########################
if args.uniprot:
    uniprotOptD=os.path.join(outputD,"uniprot")
    argsL = [input_faa, "uniprot", uniprotOptD, uniprot_db, outD]
    kwargsD = {"otherPara":"-T 40 --cpu %s"%(thread),
                "creteria":args.uc,
                "force":args.uf,
                "program":"phmmer"}
    uniprott = Thread(target=oneStepRun,args=argsL, kwargs=kwargsD)
    uniprott.start()

###########################  Run/Parse Uniprot(Swiss-prot) phmmer ##########################
if args.pdb:
    pdbOptD=os.path.join(outputD,"pdb")
    argsL = [input_faa, "pdb", pdbOptD, pdb_db, outD]
    kwargsD = {"otherPara":"-T 40 --cpu %s"%(thread),
               "creteria":args.bc,
               "force":args.bf,
               "program":"phmmer"}
    pdbt = Thread(target=oneStepRun,args=argsL, kwargs=kwargsD)
    pdbt.start()

###########################
#### fetch result #########
###########################
if args.kegg:
    keggt.join()

if args.vog:
    vogt.join()

if args.pfam:
    pfamt.join()

if args.phrog:
    phrogt.join()

if args.uniprot:
    uniprott.join()

if args.pdb:
    pdbt.join()


##################################
## fmt function annotation
##################################

summaryFile = os.path.join(outputD,"Vanno_summary.tsv")
fmt_outD = split_dict_for_pandas(outD)
res_df = pd.DataFrame.from_dict(fmt_outD)
res_df = res_df.fillna("NA")

##################################
## add annotaion for the results
##################################
## add kegg annotation
if args.kegg:
    keggD = read_kegg_anno(kegg_anno_file)
    res_df["kegg_des"] = res_df["kegg_acc"].map(keggD)  


## add vog annotation
if args.vog:
    vogD_fc, vogD_des = read_vog_anno(vog_anno_file)
    res_df["vog_FunctionalCategory"] = res_df["vog_acc"].map(vogD_fc)
    res_df["vog_FunctionalDescription"] = res_df["vog_acc"].map(vogD_des)  


## add phrog annotation
if args.phrog:
    phrog_db_anno_df = pd.read_csv(phrog_db_anno,sep="\t",names=["phrog_ori","color","phrog_annot","phrog_category"])
    phrog_db_anno_df["phrog_acc"] = ["phrog_%s"%i for i in phrog_db_anno_df.phrog_ori]
    phrog_db_anno_df_sub = phrog_db_anno_df.loc[:,["phrog_annot","phrog_category","phrog_acc"]]
    if args.phrog_mode == "hmmsearch":
        res_df = res_df.reset_index().merge(phrog_db_anno_df_sub,left_on="phrog_acc",right_on="phrog_acc",how="left")
    elif args.phrog_mode  == "hhblits":
        res_df = res_df.reset_index().merge(phrog_db_anno_df_sub,left_on="phrog_hhblits_acc",right_on="phrog_acc",how="left")

## add pdb annotation
if args.pdb:
    pdb_db_anno_df = pd.read_csv(pdb_db_anno,sep="\t",names=["pdb_id","pdb_annot"])
    res_df = res_df.merge(pdb_db_anno_df,left_on="pdb_des",right_on="pdb_id",how="left")


## select coloumn to save
print(res_df.head(3))

header=[
    "index",
    "phrog_acc","phrog_annot","phrog_category",
    "pfam_acc", "pfam_des",
    "kegg_acc", "kegg_des",
    "vog_acc", "vog_FunctionalCategory","vog_FunctionalDescription",
    ""
]
sorted_header = [i for i in header if i in res_df.columns]
res_df.loc[:,sorted_header].replace("NA","").to_csv(summaryFile,index=True,sep="\t")


if __name__ == "__main__":
    print("Vanno done!")
