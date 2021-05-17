#!/bin/bash

# scripts to extract rhizobial genes from Arthrobacter pan-genme from Shimasaki et al., 2021
#
# originally by Tomohisa Shimasaki
# shimasaki.tomohisa.45c@st.kyoto-u.ac.jp

#directories
db="/Volumes/Tomohisa/GenomeAnalysis/BLAST/Result/Rhizobia/db"
list="/Volumes/Tomohisa/GenomeAnalysis/BLAST/Result/Rhizobia/list"
query="/Volumes/Tomohisa/GenomeAnalysis/BLAST/Result/Rhizobia/query"
data="/Volumes/Tomohisa/GenomeAnalysis/BLAST/Result/Rhizobia/data"
output="/Volumes/Tomohisa/GenomeAnalysis/BLAST/Result/Rhizobia/output"
OG="/Volumes/Tomohisa/GenomeAnalysis/OrthoFinder/Results_May30/Orthogroup_Sequences"
table="/Volumes/Tomohisa/GenomeAnalysis/BLAST/Result/Rhizobia/Table"

## extract target DNA seq 3841
list_3841=$(cat ${list}/3842_list.txt)

for i in $list_3841
do
seqkit grep -nrp $i ${db}/Rlv3841.fasta  >> ${query}/seq_list_3831.fasta
done

## translate to  amino acid seq
transeq -table 1 ${query}/seq_list_3831.fasta ${query}/protein_list_3841.fasta

##Run blast 3841
blastp -num_threads 8 -outfmt "6 qseqid sseqid qcovs" -evalue 1e-5  -max_target_seqs 1 -db ${db}/Arthrobacter_all -query ${query}/protein_list_3841.fasta -out ${output}/Rlv3841_e5.txt

##blastp -num_threads 8 -outfmt "6 qseqid sseqid qcovs" -evalue 1e-15  -max_target_seqs 1 -db ${db}/Arthrobacter_all -query ${query}/protein_list_3841.fasta -out ${output}/Rlv3841_e15.txt

##extract OG ID 3841
Art_3841_list=$(awk  '$3 >=50' ${output}/Rlv3841_e5.txt | cut -f 2)

for i in $Art_3841_list
do
grep -x \>$i -rl ${OG} >> ${output}/OG_list_3841_e5.txt
done

cut -f 8 -d "/"  ${output}/OG_list_3841_e5.txt | sed 's/.fa//g' | sort | uniq > ${list}/OG_lD_3841_e5.txt

## extract target DNA seq 1021
grep -f ${list}/1021_gene_list.txt ${list}/1021_protein_code.txt > ${list}/protein_list_1021.txt

list_1021=$(cut -f 1 ${list}/protein_list_1021.txt)

for i in $list_1021
do
seqkit grep -nrp $i ${db}/1021_proteome.fasta >> ${query}/protein_list_1021.fasta
done

##Run blast 1021
blastp -num_threads 8 -outfmt "6 qseqid sseqid qcovs" -evalue 1e-5  -max_target_seqs 1 -db ${db}/Arthrobacter_all -query ${query}/protein_list_1021.fasta -out ${output}/1021_e5.txt

##blastp -num_threads 8 -outfmt "6 qseqid sseqid qcovs" -evalue 1e-15  -max_target_seqs 1 -db ${db}/Arthrobacter_all -query ${query}/protein_list_1021.fasta -out ${output}/1021_e15.txt

##extract OG ID 1021
Art_1021_list=$(awk  '$3 >=50' ${output}/1021_e5.txt | cut -f 2)

for i in $Art_1021_list
do
grep -x \>$i -rl ${OG} >> ${output}/OG_list_1021_e5.txt
done

cut -f 8 -d "/"  ${output}/OG_list_1021_e5.txt | sed 's/.fa//g' | sort | uniq > ${list}/OG_lD_1021_e5.txt

##extract rhizobial genes itentidied in Arthrobacter
##1021
genes_1021=$(cut -f 1 ${output}/1021_e5.txt | cut -f 2 -d "_")

for i in $genes_1021
do
grep  $i ${list}/1021_protein_code.txt  >> ${output}/identified_genes_1021_e5.txt
done

##3841
genes_3841=$(cut -f 1 ${output}/Rlv3841_e5.txt | cut -f 1 -d "_")

for i in $genes_3841
do
grep  $i ${list}/3841_protein_code.txt  >> ${output}/identified_genes_3841_e5.txt
done


##make table data 

##1021
Pro_ID_1021=$(awk  '$3 >=50' ${output}/1021_e5.txt | cut -f 1)
Gene_ID_1021=$(awk  '$3 >=50' ${output}/1021_e5.txt | cut -f 2)


for i in $Pro_ID_1021
do 

a=$(grep $i ${db}/1021_proteome.fasta)
echo -e ${i}"\t"${a} >> ${table}/Prot_Table_1021_e5.txt

done

for j in $Gene_ID_1021
do
　
b=$(grep -x \>$j -rl ${OG}  |cut -f 8 -d "/" | sed 's/.fa//g')

echo -e ${j}"\t"${b} >> ${table}/OG_Table_1021_e5.txt

done

R
#directories
list="/Volumes/Tomohisa/GenomeAnalysis/BLAST/Result/Rhizobia/list/"
output="/Volumes/Tomohisa/GenomeAnalysis/BLAST/Result/Rhizobia/output/"
table="/Volumes/Tomohisa/GenomeAnalysis/BLAST/Result/Rhizobia/Table/"

OG <-  read.table(file=paste(table,"OG_Table_1021_e5.txt", sep=""), stringsAsFactors=F, sep="\t", check.names=F)
Pro <-  read.table(file=paste(table,"Prot_Table_1021_e5.txt", sep=""), stringsAsFactors=F, sep="\t", check.names=F, quote = "")
Base <- read.table(file=paste(output,"1021_e5.txt", sep=""), stringsAsFactors=F, sep="\t", check.names=F)

idx <- match(Base$V1, Pro$V1)
Base$Pro <- Pro$V2[idx]

idx <- match(Base$V2, OG$V1)
Base$OG <- OG$V2[idx]

write.table(Base, file=paste(table, "1021_Table.txt", sep=""), quote=F, sep="\t", row.names=T, col.names=T)

q()
