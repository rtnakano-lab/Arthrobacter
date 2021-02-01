#!/bin/bash

# scripts to extract rhizobial genes from Arthrobacter pan-genme from Shimasaki et al., 2021
#
# originally by Tomohisa Shimasaki
# shimasaki.tomohisa.45c@st.kyoto-u.ac.jp

db="/Users/shimasakitomohisa/Desktop/github/2021_mBio/Rhizobia/db"
list="/Users/shimasakitomohisa/Desktop/github/2021_mBio/Rhizobia/list"
query="/Users/shimasakitomohisa/Desktop/github/2021_mBio/Rhizobia/query"
data="/Users/shimasakitomohisa/Desktop/github/2021_mBio/Rhizobia/data"
output="/Users/shimasakitomohisa/Desktop/github/2021_mBio/Rhizobia/output"
OG="/Users/shimasakitomohisa/Desktop/GenomeAnalysis/OrthoFinder/Results_May30/Orthogroup_Sequences"

## extract target DNA seq 3841
list_3841=$(cat ${list}/3842_list.txt)

for i in $list_3841
do
seqkit grep -nrp $i ${db}/Rlv3841.fasta  >> ${query}/seq_list_3831.fasta
done

## translate to  amino acid seq
transeq -table 1 ${query}/seq_list_3831.fasta ${query}/protein_list_3841.fasta

## extract target DNA seq 1021
grep -f ${list}/1021_gene_list.txt ${list}/1021_protein_code.txt > ${list}/protein_list_1021.txt

list_1021=$(cut -f 1 ${list}/protein_list_1021.txt)

for i in $list_1021
do
seqkit grep -nrp $i ${db}/1021_proteome.fasta >> ${query}/protein_list_1021.fasta
done

##Run blast 3841
blastp -num_threads 8 -outfmt "6 qseqid sseqid qcovs" -evalue 1e-5  -max_target_seqs 1 -db ${db}/Arthrobacter_all -query ${query}/protein_list_3841.fasta -out ${output}/Rlv3841.txt

##extract OG ID 3841
Art_3841_list=$(cut -f 2 ${output}/Rlv3841.txt)

for i in $Art_3841_list
do
grep -x \>$i -rl ${OG} >> ${output}/OG_list_3841.txt
done

##Run blast 3841
blastp -num_threads 8 -outfmt "6 qseqid sseqid qcovs" -evalue 1e-5  -max_target_seqs 1 -db ${db}/Arthrobacter_all -query ${query}/protein_list_3841.fasta -out ${output}/Rlv3841.txt

##extract OG ID 3841
Art_3841_list=$(awk  '$3 >=50' ${output}/Rlv3841.txt | cut -f 2)

for i in $Art_3841_list
do
grep -x \>$i -rl ${OG} >> ${output}/OG_list_3841.txt
done

cut -f 9 -d "/"  ${output}/OG_list_3841.txt | sed 's/.fa//g' | sort | uniq > ${list}/OG_lD_3841.txt

##Run blast 1021
blastp -num_threads 8 -outfmt "6 qseqid sseqid qcovs" -evalue 1e-5  -max_target_seqs 1 -db ${db}/Arthrobacter_all -query ${query}/protein_list_1021.fasta -out ${output}/1021.txt

##extract OG ID 3841
Art_1021_list=$(awk  '$3 >=50' ${output}/1021.txt | cut -f 2)

for i in $Art_1021_list
do
grep -x \>$i -rl ${OG} >> ${output}/OG_list_1021.txt
done

cut -f 9 -d "/"  ${output}/OG_list_1021.txt | sed 's/.fa//g' | sort | uniq > ${list}/OG_lD_1021.txt