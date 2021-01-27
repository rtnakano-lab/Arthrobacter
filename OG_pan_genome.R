# R script for pan-genome analysis using orthologous groups from shimasaki et al., 2021
#orininally from Tomohisa Shimasaki
#shimasaki.tomohisa.45c@st.kyoto-u.ac.jp

# cleanup
rm(list = ls())
options(warn=-1)

#load library
library(ggplot2,      quietly=T, warn.conflicts=F)
library(dplyr,        quietly=T, warn.conflicts=F)
library(reshape2,     quietly=T, warn.conflicts=F)

# directories
data              <- "/Users/shimasakitomohisa/Desktop/DataAnalysis/comparative_genome_analysis/data/"
list              <- "/Users/shimasakitomohisa/Desktop/DataAnalysis/comparative_genome_analysis/list/"
Fig               <- "/Users/shimasakitomohisa/Desktop/DataAnalysis/comparative_genome_analysis/fig/"
results           <- "/Users/shimasakitomohisa/Desktop/DataAnalysis/comparative_genome_analysis/results/"

# load data
OG_tab    <- read.table(file=paste(data,"Orthogroups.GeneCount.txt", sep=""), header=T, stringsAsFactors=F, row.names=1, sep="\t", check.names=F)
design     <- read.table(file=paste(list,"strain_list.txt", sep=""), sep="\t", header=T, check.names=F, stringsAsFactors=F) 

# make presence/absence table
func.tab  <- (OG_tab>0)*1

# correlation analysis
cor <- 1-cor(func.tab)

# PCoA 
pcoa     <- cmdscale(cor, k=2, eig=T)　
points   <- as.data.frame(pcoa$points, check.names=F, stringsAsFactors=F)
eig      <- as.data.frame(pcoa$eig, check.names=F, stringsAsFactors=F)

p1 <- (100 * eig[1,1]/sum(eig))
p2 <- (100 * eig[2,1]/sum(eig))

##drawing figures
idx           <- match(rownames(points), design$Strain)
table        <- data.frame(points, origin=design$Origin[idx], lineage=design$sublineages[idx], check.names=F, stringsAsFactors=F)

origin       <- data.frame(name   = c("Plant", "Soil", "Other"),
                           shapes = c(18, 16, 8),
                           col    = c("#07346c","#008000", "#fdc23e"),
                           check.names=F, stringsAsFactors=F)
lineage      <- data.frame(name   = c("A", "B", "C"),
                           col = c("#ef7318","#8fc93e", "#10677f"),
                           check.names=F, stringsAsFactors=F)

table$origin    <- factor(table$origin,   levels=origin$name)
table$lineage   <- factor(table$lineage,  levels=lineage$name)

## PCoA
plot1 <- ggplot(table, aes(x=V1, y=V2, color=lineage, shape=origin)) +
  geom_point(alpha=.9, size=4) +
  theme_classic(base_size = 20) +
  scale_shape_manual(values = origin$shapes) +
  scale_colour_manual(values = lineage$col) +
  labs(x = "PC1 (23.0%)",
       y = "PC2 (15.7%)") +
  theme(axis.text=element_text(size = 15, colour = "black"),
        legend.text=element_text(size=15),
        legend.title=element_blank())

ggsave(plot1, file=paste(Fig, "PCoA_OG.pdf", sep=""), width=6, height=5, bg="transparent")

## genome size comparison
genome_size  <- data.frame(size=colSums(func.tab))
idx          <- match(rownames(genome_size), design$Strain)
genome_size$origin     <- design$Origin[idx]                 
genome_size$sublineage <- design$sublineages[idx]     

##sublineage
size_lineage <- ggplot(genome_size, aes(x=sublineage, y=size, color=sublineage)) +
  geom_boxplot(size=0.8, outlier.shape = NA, width=0.8) +
  scale_x_discrete(limit=c("A","B","C")) +
  geom_jitter(size=3, alpha = 0.8, position=position_jitter(width=0.15)) +
  theme_classic(base_size = 20) +
  scale_y_continuous(limits = c(1500, 5000)) +
  scale_color_manual(values=lineage$col) +
  theme(axis.ticks=element_line(colour = "black"), 
        axis.text.x =element_text(size = 20, colour = "black"),
        axis.text.y =element_text(size = 15, colour = "black"),
        legend.position = "none", 
        axis.title.x = element_blank()) +
  ylab("Number of Orthogroups")
size_lineage
ggsave(size_lineage, file=paste(Fig, "OG_size_lineages.pdf", sep=""), width=6, height=5, bg="transparent")

##origin
size_origin <- ggplot(genome_size, aes(x=origin, y=size, color=origin)) +
  geom_boxplot(size=0.8, outlier.shape = NA, width=0.8) +
  scale_x_discrete(limit=c("Plant","Soil","Other")) +
  geom_jitter(size=3, alpha = 0.8, position=position_jitter(width=0.15)) +
  theme_classic(base_size = 20) +
  scale_y_continuous(limits = c(1500, 5000)) +
  scale_color_manual(values=origin$col) +
  theme(axis.ticks=element_line(colour = "black"), 
        axis.text.x =element_text(size = 20, colour = "black"),
        axis.text.y =element_text(size = 15, colour = "black"),
        legend.position = "none", 
        axis.title.x = element_blank()) +
  ylab("Number of Orthogroups")

ggsave(size_origin, file=paste(Fig, "OG_size_origin.pdf", sep=""), width=6, height=5, bg="transparent")

## similarity analysis
out           <- matrix(0, nrow = 99 , ncol = 99)
colnames(out) <- colnames(func.tab)
rownames(out) <- colnames(func.tab)
func.tab      <- as.data.frame(func.tab, check.names=F, stringsAsFactors=F)

list <- colnames(func.tab)

for (i in list){
  
  idx_A  <-  func.tab[i]
  idx_A  <-  idx_A %>%
    filter(idx_A[1]>0)
  list_A  <- rownames(idx_A)
  
  for (j in list){
    idx_B  <-  func.tab[j]
    idx_B  <-  idx_B %>%
      filter(idx_B[1]>0)
    list_B  <- rownames(idx_B)
    
    similarity <- length(unique(intersect(list_A, list_B)))/length(unique(c(list_A, list_B)))*100
    
    out[i,j] = similarity
    
    print(c(i,j,similarity))
  }}

out[upper.tri(out)] <- NA
diag(out) <- NA
    
sim_table     <- melt(out)
sim_table     <- sim_table[!is.na(sim_table$value),]

sim_strain  　<- sim_table[sim_table$value>99,]
list_strain   <- unique(sim_strain$Var1)      

## remove high similarity genomes
idx <- match(rownames(genome_size), list_strain)
r_genome_size    <- genome_size %>% filter(!rownames(genome_size) %in% list_strain)

## wilcox
ori_wilcox <-pairwise.wilcox.test(r_genome_size$size, r_genome_size$origin, p.adj="bonferroni", exact=F)
sink(paste(results, "wilcox_KEGG_origin.txt", sep=""))
print(ori_wilcox)
sink()

lineage_wilcox <-pairwise.wilcox.test(r_genome_size$size, r_genome_size$sublineage, p.adj="bonferroni", exact=F)
sink(paste(results, "wilcox_KEGG_lineage.txt", sep=""))
print(lineage_wilcox)
sink()
