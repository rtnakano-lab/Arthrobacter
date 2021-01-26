#!/bin/bash
#qiime2 script for Shimasaki et al., 2021
#orininally from Tomohisa Shimasaki
#shimasaki.tomohisa.45c@st.kyoto-u.ac.jp

##directory
dir="/Users/shimasakitomohisa/Desktop/DataAnalysis/Qiime2/2021_mBio"

## activate Qiime2
conda activate qiime2-2020.2

## load fastq data
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path ${dir}/list/manifest.tsv \
  --output-path ${dir}/output/paired-end-demux.qza \
  --input-format PairedEndFastqManifestPhred33V2
  
 ## convert qza to qzv 
qiime demux summarize \
  --i-data ${dir}/output/paired-end-demux.qza  \
  --o-visualization ${dir}/output/paired-end-demux.qzv
  
  ## sequence quality control and feature table construction by DADA2
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs ${dir}/output/paired-end-demux.qza \
  --p-trim-left-f 20 \
  --p-trim-left-r 20 \
  --p-trunc-len-f 220 \
  --p-trunc-len-r 220 \
  --o-table  ${dir}/output/table.qza \
  --o-representative-sequences ${dir}/output/rep-seqs.qza \
  --o-denoising-stats ${dir}/output/denoising-stats.qza \
  --p-n-threads 0 \
  --verbose

## FeatureTable and FeatureData summarises
qiime feature-table summarize \
  --i-table ${dir}/output/table.qza \
  --o-visualization ${dir}/output//table.qzv \
  --m-sample-metadata-file ${dir}/list/metadata.tsv
  
qiime feature-table tabulate-seqs \
  --i-data ${dir}/output/rep-seqs.qza \
  --o-visualization ${dir}/output/rep-seqs.qzv

## taxonomy assingnment
qiime feature-classifier classify-sklearn \
  --i-classifier ${dir}/list/gg-13-8-99-515-806-nb-classifier.qza\
  --i-reads ${dir}/output/rep-seqs.qza \
  --o-classification ${dir}/output/taxonomy.qza

## convert qza to qzv 
qiime metadata tabulate \
  --m-input-file ${dir}/output/taxonomy.qza \
  --o-visualization ${dir}/output/taxonomy.qzv

## remove sequences assiged to mitochondria and chloroplast
qiime taxa filter-table \
  --i-table ${dir}/output/table.qza \
  --i-taxonomy ${dir}/output/taxonomy.qza \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-table ${dir}/output/table-no-mitochondria-no-chloroplast.qza

qiime feature-table summarize \
  --i-table ${dir}/output/table-no-mitochondria-no-chloroplast.qza \
  --o-visualization ${dir}/output/table-no-mitochondria-no-chloroplast.qzv \
  --m-sample-metadata-file ${dir}/list/metadata.tsv

##Visualize taxonomy bar plot
qiime taxa barplot \
  --i-table ${dir}/output/table-no-mitochondria-no-chloroplast.qza \
  --i-taxonomy ${dir}/output/taxonomy.qza \
  --m-metadata-file ${dir}/list/metadata.tsv \
  --o-visualization ${dir}/output/taxa-bar-plots-no-mitochondria-no-chloroplast.qzv

##extract ASV table
qiime tools export \
  --input-path ${dir}/output/table-no-mitochondria-no-chloroplast.qza \
  --output-path ${dir}/output/feature-table
biom convert \
  --to-tsv -i ${dir}/output/feature-table/feature-table.biom \
  -o ${dir}/data/ASV-table.tsv

rm ${dir}/output/feature-table/feature-table.biom
rmdir ${dir}/output/feature-table

## tree genetarion for phylogenetic diversity analysis
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences ${dir}/output/rep-seqs.qza \
  --o-alignment ${dir}/output/aligned-rep-seqs.qza \
  --o-masked-alignment ${dir}/output/masked-aligned-rep-seqs.qza \
  --o-tree ${dir}/output/unrooted-tree.qza \
  --o-rooted-tree ${dir}/output/rooted-tree.qza

## rarefaction analysis
qiime diversity alpha-rarefaction \
  --i-table ${dir}/output/table-no-mitochondria-no-chloroplast.qza \
  --i-phylogeny ${dir}/output/rooted-tree.qza \
  --p-max-depth 15000 \
  --p-steps 11 \
  --m-metadata-file ${dir}/list/metadata.tsv \
  --o-visualization ${dir}/output/alpha-rarefaction.qzv
  
  ## Alpha and beta diversity analysis
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny ${dir}/output/rooted-tree.qza \
  --i-table ${dir}/output/table-no-mitochondria-no-chloroplast.qza \
  --p-sampling-depth 6000 \
  --m-metadata-file ${dir}/list/metadata.tsv \
  --output-dir ${dir}/output/core-metrics-results-6000 

## Export UniFrac distance
qiime tools export \
    --input-path ${dir}/output/core-metrics-results-6000/weighted_unifrac_distance_matrix.qza \
    --output-path ${dir}/data/weighted_unifrac_distance_matrix

qiime tools export \
    --input-path ${dir}/output/core-metrics-results-6000/unweighted_unifrac_distance_matrix.qza \
    --output-path ${dir}/data/unweighted_unifrac_distance_matrix


