#!/usr/bin/bash

## value
db_dir=${1:-"./databases"}

## db version
kegg_version=${2:-'2022-02-01'}
vog_version=${3:-'latest'}
pfam_version=${4:-'Pfam35.0'}
uniprot_version=${5:-"sprot"} ## sprot(swissprot), trembl, all(both)
phrog_version=${6:-"v3"}
################################

mkdir -p $db_dir
fullpath=`realpath $0`
script_dir=`dirname $fullpath `
## based on vibarnt software to pick the vog related with virus(detail below)
vog_only_virus_vibrant="$script_dir/VIBRANT_vog_profiles.txt"
echo $vog_only_virus_vibrant


################################

## kegg db hmm
function download_kegg_hmm(){
    version=${1:-'2022-02-01'}
    db_dir=${2:-'.'}

    mkdir -p $db_dir/kegg/$version  && cd $db_dir/kegg/$version

    wget ftp://ftp.genome.jp/pub/db/kofam/archives/$version/ko_list.gz
    wget ftp://ftp.genome.jp/pub/db/kofam/archives/$version/profiles.tar.gz
    tar -zxvf profiles.tar.gz
    cat profiles/prokaryote.hal |xargs -i cat profiles/{} >KEGG_profiles_prokaryotes.HMM
    ln -s `pwd`/KEGG_profiles_prokaryotes.HMM ..
    cd -  && echo "kegg done with version: $version"

}


## vog db
function download_vog_hmm(){
     version=${1:-'latest'}
     db_dir=${2:-'.'}

     mkdir -p $db_dir/vog/$version && cd $db_dir/vog/$version

     wget https://fileshare.csb.univie.ac.at/vog/$version/vog.hmm.tar.gz &
     wget https://fileshare.csb.univie.ac.at/vog/$version/vog.annotations.tsv.gz
     wait
     ## only pick the virus related vog, this were from vibrant.
     ## retaining profiles that had at least one significant hit to any of the 15,238 NCBI-acquired viruses using BLASTp.
     mkdir vog_hmm && tar -zxvf vog.annotations.tsv.gz vog_hmm
     cat $vog_only_virus_vibrant|sed 's/\r//;s/$/.hmm/'|xargs -i cat vog_hmm/{} >VOGDB_phage.HMM
     ln -s `pwd`/VOGDB_phage.HMM ..
     cd - && echo "vog done with version: $version"
}

## pfam db
function download_pfam_hmm(){
     version=${1:-'Pfam35.0'}
     db_dir=${2:-'.'}

     mkdir -p $db_dir/pfam/$version && cd $db_dir/pfam/$version

     wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/$version/Pfam-A.hmm.gz
     wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/$version/Pfam-A.hmm.dat.gz
     gunzip Pfam-A.hmm.gz
     ln -s `pwd`/Pfam-A.hmm ..
     cd - && echo "pfam done with version: $version"
}

## PHROG db
function download_phrog_hmm(){
    version=${1:-'v3'}
    db_dir=${2:-'.'}

    mkdir -p $db_dir/PHROG/$version && cd $db_dir/PHROG/$version

    ## data download from http://millardlab.org/2021/11/21/phage-annotation-with-phrogs/
    wget http://s3.climb.ac.uk/ADM_share/all_phrogs.hmm.gz
    ## offical website HMM is provied hhm file, which is used for hhsearch
    #wget https://phrogs.lmge.uca.fr/downloads_from_website/HMM_phrog.tar.gz
    wget https://phrogs.lmge.uca.fr/downloads_from_website/phrog_annot_$version.tsv
    gunzip all_phrogs.hmm.gz
    ln -s `pwd`/all_phrogs.hmm ..
    cd - && echo "PHROG done"
}

function download_uniprot_seq(){
    version=${1:-""} ## sprot, trembl, all
    db_dir=${2:-'.'}

    mkdir -p $db_dir/uniprot/$version && cd $db_dir/uniprot/$version
    if [ $version == "all" || $version == "sprot" ];then
        wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
        wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz
        gunzip uniprot_sprot.fasta.gz uniprot_sprot.dat.gz
    elif [ $version == "all" || $version == "trembl" ];then
        wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz
        wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.dat.gz
        gunzip uniprot_trembl.fasta.gz uniprot_trembl.dat.gz
    elif [ $version == "all" ];then
        cat uniprot_sprot.fasta uniprot_trembl.fasta >uniprot_trembl_sprot.merge.fasta
    fi
}

function download_pdb_seq(){
    version=${1:-""}
    db_dir=${2:-'.'}

    mkdir -p $db_dir/pdb && cd $db_dir/pdb
    wget https://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt.gz
    gunzip pdb_seqres.txt.gz
    grep \> pdb_seqres.txt |cut -d " " -f 1,4-200|sed 's/>//g' |sed 's/  /\t/' >pdb_seqres.header.anno.txt
}

##########################
## main ##
##########################
download_kegg_hmm $kegg_version $db_dir
download_vog_hmm $vog_version $db_dir
download_pfam_hmm $pfam_versiom $db_dir
download_uniprot_seq $uniprot_version $db_dir
download_phrog_hmm $phrog_version $db_dir
