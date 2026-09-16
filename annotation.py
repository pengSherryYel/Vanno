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
                  help='number of threads, each occupies 1 CPU [default=1, max of 1 CPU per scaffold]')
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
anno.add_argument('-rs', '--phrog_mode',choices=['hmmsearch','mmseqs'],dest="phrog_mode",default="mmseqs", help="run phrog use mmseqs|hmmsearch")
anno.add_argument('-rc', '--phrogC',type=float, default="1e-5", dest="rc", help="phrog creteria. discard the not meet this creteria")
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
print(thread)



##########################
#### Function ############
##########################

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


def runMMseqs(inputfile, prefix, wd, dbseq, otherPara="-s 7"):
    checkEnv("mmseqs")
    mkdirs(wd)
    dbseq = os.path.realpath(dbseq)
    cmd = "mmseqs createdb {0} {2}/{1}.target_seq ".format(inputfile, prefix, wd)
    cmd2 = "mmseqs search {3} {2}/{1}.target_seq {2}/{1}.results_mmseqs ./tmp {4}".format(inputfile, prefix, wd, dbseq, otherPara)
    cmd3 = "mmseqs createtsv {3} {2}/{1}.target_seq {2}/{1}.results_mmseqs {2}/{1}.results.tsv".format(inputfile, prefix, wd, dbseq, otherPara)

    print("RUN command: %s\n" % cmd)
    obj = Popen(cmd, shell=True, stdout=PIPE, stderr=STDOUT)
    [logging.info(line.rstrip()) for line in obj.stdout]
    obj.wait()

    print("RUN command: %s\n" % cmd2)
    # obj = Popen(cmd2, shell=True)
    obj = Popen(cmd2, shell=True, stdout=PIPE, stderr=STDOUT)
    [logging.info(line.rstrip()) for line in obj.stdout]
    obj.wait()

    print("RUN command: %s\n" % cmd3)
    obj = Popen(cmd3, shell=True, stdout=PIPE, stderr=STDOUT)
    [logging.info(line.rstrip()) for line in obj.stdout]
    obj.wait()
    print("mmseqs done!")
    return "%s/%s.results.tsv" % (wd, prefix)


def load_hmmsearch_opt(hmmsearch_opt, creteria=1e-5, reverse=False):
    '''
    Aim: parse the hmmersearch output. this file contain multiple columns. is the output from -tblout parameter
    reverse: sometimes the query id and ref id are in different place. if false first query, second ref, verse visa.
    Return: dict.  d[refname][annoacc] = description
    '''
    print("loading hmmsearch output")
    annoD = defaultdict(dict)
    annoMinD = defaultdict(dict)
    tmp = {}
    with open(hmmsearch_opt) as f:
        for line in f:
            if not line.startswith("#"):
                t = re.split("\s+", line.strip("\n"))
                if not reverse:
                    target_name, target_accession, query_name, accession, Evalue, score, bias, bst_Evalue, bst_score, bst_bias,\
                        exp, reg, clu, ov, env, dom, rep, inc, *description_of_target = t
                elif reverse:
                    ## because use phmmer, the target and query are change
                    query_name, accession, target_name, target_accession, Evalue, score, bias, bst_Evalue, bst_score, bst_bias,\
                        exp, reg, clu, ov, env, dom, rep, inc, *description_of_target = t
                accession = accession.split(".")[0]
                # print(target_name,Evalue,bst_Evalue)
                if float(Evalue) <= float(creteria) and float(bst_Evalue) <= float(creteria):
                    annoD[target_name][accession] = query_name
                    ## add minEvalue select this only fetch a target a ref with smallest value
                    if target_name not in tmp:
                        tmp[target_name] = [accession,query_name,Evalue]
                    else:
                        if float(Evalue) <= float(tmp[target_name][2]):
                            tmp[target_name] = [accession,query_name,Evalue]
        ##format tmp
        for key,values in tmp.items():
            accession,query_name,Evalue = values
            annoMinD[key][accession] = query_name

        #return annoD
        return annoMinD


def oneStepRun(inputfile, prefix, wd, db, outD, otherPara,creteria, force=False, program="hmmsearch"):
    print("###### %s begin ######"%prefix)
    hmm_outPath = ""
    if program == "hmmsearch":
        hmm_outPath = "%s/%s.hmmsearch.tblout" % (wd, prefix)
        hmmModel = db
        if not os.path.exists(hmm_outPath) or force:
            hmm_outPath = runHmmsearch(inputfile, prefix, wd, hmmModel, otherPara=otherPara)
        else:
            print("Skip %s the running part, cause output file found!"%prefix)
        annoD = load_hmmsearch_opt(hmm_outPath, creteria=creteria, reverse=False)

    elif program == "phmmer":
        hmm_outPath = "%s/%s.phmmer.tblout" % (wd, prefix)
        dbseq = db
        if not os.path.exists(hmm_outPath) or force:
            hmm_outPath = runPhmmer(inputfile, prefix, wd, dbseq, otherPara=otherPara)
        else:
            print("Skip %s the running part, cause output file found!"%prefix)
        annoD = load_hmmsearch_opt(hmm_outPath, creteria=creteria, reverse=True)

    elif program == "mmseqs":
        mmseq_outPath = "%s/%s.mmseqs.results.tsv" % (wd, prefix)
        dbseq = db
        if not os.path.exists(mmseq_outPath) or force:
            mmseq_outPath = runMMseqs(inputfile, "%s_mmseqs"%prefix, "%s_mmseqs"%wd, dbseq, otherPara=otherPara)
        else:
            print("Skip %s the running part, cause output file found!"%prefix)
        annoD = load_hmmsearch_opt(mmseq_outPath, creteria=1e-5, reverse=True)

    else:
        print("Wrong program parameter! please use hmmsearch|phmmer|mmseqs")
        exit()

    outD[prefix] = annoD
    print("###### %s end ######"%prefix)

def split_dict_for_pandas(indict):
    outd = defaultdict(dict)
    for db,annos in indict.items():
        for query, values in annos.items():
            for accession,name in values.items():
                #print(db,query,name)

                des = "%s_des"%db
                acc = "%s_acc"%db
                outd[des][query]= name
                outd[acc][query]= accession
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

###########################  Run/Parse PHROG hmmsearch ##########################
if args.phrog:
    
    phrogOptD=os.path.join(outputD,"phrog")
    argsL = [input_faa, "phrog", phrogOptD, phrog_db, outD]
    if args.phrog_mode == "hmmsearch":
        kwargsD = {"otherPara":"-T 40 --cpu %s"%(thread),
                    "creteria":args.rc,
                    "force":args.rf,
                    "program":args.phrog_mode}
    elif args.phrog_mode == "mmseqs":
        kwargsD = {"otherPara":"-s 7 --threads %s"%(thread),
                    "creteria":args.rc,
                    "force":args.rf,
                    "program":args.phrog_mode}
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

summaryFile = os.path.join(outputD,"Vanno_summary.tsv")
fmt_outD = split_dict_for_pandas(outD)
res_df = pd.DataFrame.from_dict(fmt_outD)
res_df = res_df.fillna("NA")
print(res_df)

##################################
## add annotaion for the results
##################################
## add kegg annotation
if args.kegg:
    keggD = read_kegg_anno(kegg_anno_file)
    res_df["kegg_defination"] = res_df["kegg_des"].map(keggD)  


## add vog annotation
if args.vog:
    vogD_fc, vogD_des = read_vog_anno(vog_anno_file)
    res_df["vog_FunctionalCategory"] = res_df["vog_des"].map(vogD_fc)
    res_df["vog_FunctionalDescription"] = res_df["vog_des"].map(vogD_des)  



## add phrog annotation
if args.phrog:
    phrog_db_anno_df = pd.read_csv(phrog_db_anno,sep="\t",names=["phrog_ori","color","phrog_annot","phrog_category"])
    phrog_db_anno_df["phrogID"] = ["phrog_%s"%i for i in phrog_db_anno_df.phrog_ori]
    phrog_db_anno_df_sub = phrog_db_anno_df.loc[:,["phrog_annot","phrog_category","phrogID"]]
    res_df = res_df.reset_index().merge(phrog_db_anno_df_sub,left_on="phrog_des",right_on="phrogID",how="left")

## add pdb annotation
if args.pdb:
    pdb_db_anno_df = pd.read_csv(pdb_db_anno,sep="\t",names=["pdb_id","pdb_annot"])
    res_df = res_df.merge(pdb_db_anno_df,left_on="pdb_des",right_on="pdb_id",how="left")

## select coloumn to save
res_df.to_csv(summaryFile,index=True,sep="\t")


if __name__ == "__main__":
    print("Vanno done!")
