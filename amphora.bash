#!/bin/bash
# bash script for runnung AMPHORA

# Ryohei Thomas Nakano, PhD; nakano@mpipz.mpg.de
# Tomohisa Shimasaki; shimasaki.tomohisa.45c@st.kyoto-u.ac.jp

# exits whenever a function returns 1
set -e
# exits if unset variables are used
set -o nounset

# log function
log() {
        echo $(date -u)": "$1 >> $logfile
}

# paths
data="/biodata/dep_psl/grp_psl/ThomasN/arthrobacter"
wd="/netscratch/dep_psl/grp_psl/ThomasN/amphora/200528arthrobacter2"
logfile="${wd}/log.txt"

script_dir="/netscratch/dep_psl/grp_psl/ThomasN/tools/Phyla_AMPHORA/Scripts"
pearl_path="/usr/bin"
fasttree_path="/netscratch/dep_psl/grp_psl/ThomasN/tools/bin"

# initialization
rm -r -f ${wd}
mkdir -p ${wd}/temp
mkdir -p ${wd}/concatenate

# concatenate
log "Concatenating fasta files"
cat ${data}/ORFs_AA/*.faa > ${wd}/temp/pep_all.faa

# run marker scanner
log "Running MarkerScanner"
cd ${wd}
${pearl_path}/perl ${script_dir}/MarkerScanner.pl -Phylum 7 ${wd}/temp/pep_all.faa

# run marker aligner
log "Running MarkerAlignTrim"
${pearl_path}/perl ${script_dir}/MarkerAlignTrim.pl -OutputFormat fasta &>> ./output.txt 


# list of genomes and genes
log "Concatenating alignments"
genome_list=$(cut -f 1 ${data}/mapping.txt | grep -v "Sample\ name")
amphora_list=$(cd ${wd} && ls *.aln | sed s/\.aln//g)


for genome_id in ${genome_list}; do

	echo ">${genome_id}" >> ${wd}/concatenate/amphora.msa

	for gene in ${amphora_list}; do 

		n=$(grep -c ">" ${wd}/${gene}.aln)
		n_u=$(grep ">" ${wd}/${gene}.aln | sed s/_gene.*//g | sort -u | wc -l)
		n_tot=$(echo ${genome_list} | awk '{print NF}')

		# use genes found in all genomes as a single copy
		if [[ ${n} -eq ${n_tot}  && ${n_u} -eq ${n_tot} ]]; then
			
			awk "/${genome_id}/ {flag=1;next} />/{flag=0} flag {print}" ${wd}/${gene}.aln >> ${wd}/concatenate/amphora.msa

		fi

	done
done

# generate tree
log "Infering an approximately-maximum-likelihood phylogenetic tree"
${fasttree_path}/FastTree ${wd}/concatenate/amphora.msa > ${wd}/concatenate/amphora_tree.nwk

# clean up
rm -r -f ${wd}/temp

log "Done"
