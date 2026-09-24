
import pandas as pd
import os
from collections import defaultdict
import re

## IMPORTANT: all output parser should return dict. dict format: {query_id:{target_id:import_info}}. 
## import_info is the one I think I need to keep. it can be empty "".

######################################
### format hhsuit hhr into plain text
######################################
def parse_hhsuit_hhr(hhr_file, evalue_cutoff=1e-5, prob_cutoff=90):
    '''
    Only select the best, because this is for the function annotation
    '''

    annoD = defaultdict(dict)
    header = [ "Query", "Target", "Probab", "E-value", "Score", "Aligned_cols", "Identities", "Similarity", "Sum_probs", "Template_Neff"]
    total_record_num = 0
    with open(hhr_file, "r" ) as f:
        resD = {}

        ## parse the hhr file
        for line in f:
            line = line.strip("\x00")
            if line.startswith("Query"):
                _, query_sname = line.split("         ")
                query_sid, query_anno = query_sname.split(" ",1)
                
            if line.startswith("No "):
                idx=line.replace(" ","_") 
                
                target_db_line = f.readline().strip(">")
                target_sid, _ = target_db_line.split(" ",1)
                
                result_line = f.readline().strip().split("  ")
                tmpd = {i.split("=")[0]: i.split("=")[1] for i in result_line}
                tmpd["Query"]=query_sid
                tmpd["Target"]=target_sid

                # print(tmpd)
                resD[total_record_num] = tmpd
                total_record_num +=1
        
        ## fmt hhr
        hhr_df = pd.DataFrame.from_dict(resD, orient='index')
        final_hhr_df = hhr_df.loc[:, header]
        print("Using creteria %s and Prob %s filtering the hhblits results"%(evalue_cutoff, prob_cutoff))

        ## filter hhr by e-value and prob
        final_hhr_df_filtered = final_hhr_df[(final_hhr_df['E-value'].astype(float) <= evalue_cutoff) & (final_hhr_df['Probab'].astype(float) >= prob_cutoff)]
        final_hhr_df_filtered.to_csv(hhr_file + ".filtered.tsv", sep="\t", index=False)
        # print(final_hhr_df_filtered.head(10))

        ## Because the some protein will have more than one protein meet creteria
        ## then select the best hit for annotation 
        besthitD = {}
        for i in final_hhr_df_filtered.index:
            query = final_hhr_df_filtered.loc[i,"Query"]
            Target = final_hhr_df_filtered.loc[i,"Target"]
            Probab = final_hhr_df_filtered.loc[i,"Probab"]
            if query not in besthitD:
                besthitD[query] = [i, Target, Probab]
            else:
                if Probab > besthitD[query][-1]:
                    besthitD[query] = [i, Target, Probab]

        besthit_index = [value[0] for key,value in besthitD.items()]        
        final_hhr_df_filtered_besthit = final_hhr_df_filtered.loc[besthit_index,:]
        # print(final_hhr_df_filtered_besthit.head(10))
        final_hhr_df_filtered_besthit.to_csv(hhr_file + ".filtered_besthit.tsv", sep="\t", index=False)

        ##format besthitD to suit for the rule. {target_id: query_id: important_info}
        for query,values in besthitD.items():
            bh_index, target, Probab = values
            annoD[query][target] = Probab

    return annoD



######################################
### format hmmsearch_opt into dict
######################################
def parse_hmmsearch_opt(hmmsearch_opt, creteria=1e-5, reverse=False):
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
                    ## add minEvalue select 
                    ## Because one query gene can match to several hit in the target database, only fetch a DB target with smallest value
                    if target_name not in tmp:
                        tmp[target_name] = [accession,query_name,Evalue]
                    else:
                        if float(Evalue) <= float(tmp[target_name][2]):
                            tmp[target_name] = [accession,query_name,Evalue]

        ##format tmp
        ## because pfam will have annotaion at query_name column
        ## So this basic fmt is queryid:{targetid:import info}
        for target,values in tmp.items():
            accession,query_name,Evalue = values
            annoMinD[target][query_name] = accession

        return annoMinD


######################################
### format blast m8 into dict
######################################
def parse_m8_fmt_opt(m8_input, criteria=1e-5):
    """
    Parses blastp / mmseqs2 blast-style m8 (-outfmt 6) layout structures.
    Filters by E-value and maps query sequences to their top hit reference IDs.
    """
    d = {}
    if os.path.exists(m8_input):
        with open(m8_input) as f:
            for line in f:
                parts = line.strip("\n").split("\t")
                if len(parts) < 11:
                    continue
                query, ref, _, _, _, _, _, _, _, _, evalue = parts[:11]
                if float(evalue) <= criteria:
                    if query not in d or float(evalue) < float(d[query][-1]):
                        d[query] = [ref, evalue]
    return d



######################################
### format blast m8 into dict
######################################
def parse_pfam_dat(pfam_dfile):
    store_d = {}
    with gzip.open(pfam_dfile,"rt") as f:
        for l in f:
            # parse file
            if l.strip() == "# STOCKHOLM 1.0":
                acc_d = {}

            if l.strip().startswith("#=GF"):
                res = l.strip().strip("#=GF ").split("   ",1)
                # print(res)
                if len(res) == 2:
                    category, des = res
                else:
                    # print(l)
                    category, des = ["NA", "NA"]
                    category = res[0]
                acc_d[category] = des

            if l.strip() == "//":

                store_d[acc_d["AC"]] = acc_d

        ## manual add two can not read：PF01846.26, PF00167.25
        store_d["PF01846.26"]["ID"] = "FF"
        store_d["PF00167.25"]["ID"] = "FGF"

    pfam_dat_df = pd.DataFrame.from_dict(store_d, orient='index').reset_index().drop(["index"],axis=1)
    # print(pfam_dat_df.head(10))
    return pfam_dat_df.to_dict(orient='index')