# Vanno
annotation phage genome using hmmer3

## Install
1. Download database (KEGG, VOG, Pfam, uniprot, phrog)

```
sh envpre_db_download.sh 
```


## Usage
```
python annotation.py -i test.faa -k -v -p -r -u
```
WARINING:  Run with swissprot db (`-u`) will take long time, suggest to run other db first. 

### Usage - Paramater
```
optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
  -i I                  input faa file
  -o O                  path to deposit output folder and temporary files, will create if doesn't exist [default= working
                        directory]
  -t T                  number of threads, each occupies 1 CPU [default=1, max of 1 CPU per scaffold]
  -d D                  path to original "databases" directory that contains .HMM files (if moved from default location)
  -k, --kegg            run kegg
  -kc KC, --keggC KC    kegg creteria. discard the not meet this creteria
  -kf, --keggF          force rerun kegg
  -v, --vog             run vog
  -vc VC, --vogC VC     vogdb creteria. discard the not meet this creteria
  -vf, --vogF           force rerun vog
  -p, --pfam            run pfam
  -pc PC, --pfamC PC    pfam creteria. discard the not meet this creteria
  -pf, --pfamF          force rerun pfam
  -r, --phrog           run phrog
  -rc RC, --phrogC RC   phrog creteria. discard the not meet this creteria
  -rf, --phrogF         force rerun phrog
  -u, --uniprot         run uniprot(default swiss-prot)
  -ud {sprot,trembl,all}, --uniprotDB {sprot,trembl,all}
                        run uniprot using sprot(swiss-prot); trembl or all (sprot+trembl)
  -uc UC, --uniprotC UC
                        uniprot creteria. discard the not meet this creteria
  -uf, --uniprotF       force rerun uniprot
```

## Input & Output

### Input
amino acid sequences

## Output
Main output: Vanno_summary.tsv. It merge all the annotation results. Please check the example file example/Vanno_opt/Vanno_summary.tsv

These folder will store the results of each db.
* kegg
* pfam
* phrog
* uniprot
* vog


## Example
For more detail, please check the example in example folder.
