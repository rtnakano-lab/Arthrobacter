# R script for bacterial community analysis from shimasaki et al., 2021
#orininally from Tomohisa Shimasaki
#shimasaki.tomohisa.45c@st.kyoto-u.ac.jp
        
# clean up
rm(list=ls())
options(warn=-1)

# load library
library(vegan,     quietly=T, warn.conflicts=F)
library(stringr,   quietly=T, warn.conflicts=F)
library(ggplot2,   quietly=T, warn.conflicts=F)
library(reshape2,  quietly=T, warn.conflicts=F)

# directories
data <- "/Users/shimasakitomohisa/Desktop/DataAnalysis/Qiime2/2021_mBio/data/"
list <- "/Users/shimasakitomohisa/Desktop/DataAnalysis/Qiime2/2021_mBio/list/"
Fig  <- "/Users/shimasakitomohisa/Desktop/DataAnalysis/Qiime2/2021_mBio/Fig/"

# load data
w_UniFrac         <- as.matrix(read.table(file=paste(data,"weighted_distance-matrix.tsv", sep=""), sep="\t", header=T, check.names=F, stringsAsFactors=F, row.names=1))
u_UniFrac         <- as.matrix(read.table(file=paste(data,"unweighted_distance-matrix.tsv", sep=""), sep="\t", header=T, check.names=F, stringsAsFactors=F, row.names=1))
design            <- read.table(file=paste(list,"metadata.tsv", sep=""),  sep="\t", header=F, check.names=F, stringsAsFactors=F)

# refine metadata
design$rep       <- factor(str_sub(design$V1, -1))
design$group     <- str_replace(design$V2, ".*_", "")
design$treatment <- str_replace(design$V2, "_.*", "")

# comparison of weighted unifrac distance
Uni_table            <- melt(w_UniFrac)
idx_a                <- match(Uni_table$Var1, design$V1)
idx_b                <- match(Uni_table$Var2, design$V1)
Uni_table$sample     <- design$V2[idx_a]
Uni_table$group1      <- design$group[idx_a]
Uni_table$group2      <- design$group[idx_b]
Uni_table$treatment1  <- design$treatment[idx_a]
Uni_table$treatment2  <- design$treatment[idx_b]


# santhopine
uni_sto             <- Uni_table[Uni_table$treatment1=="Sto" & Uni_table$treatment2=="Nt" & Uni_table$group1 %in% c("Mock", "High")  & Uni_table$group2 %in% c("Bulk", "Root"),]      
uni_sto$group1      <- str_replace(uni_sto$group1, "High", "1000")
uni_sto$group1      <- str_replace(uni_sto$group1, "Mock", "0")
uni_sto$group2      <- str_replace(uni_sto$group2, "Root", "ES")
uni_sto$ID          <- paste(uni_sto$group1, uni_sto$group2, sep="vs" )


Fig_S1A <- ggplot(uni_sto , aes(x=ID, y=value)) +
  geom_boxplot(outlier.shape=NA, position=position_dodge()) +
  geom_jitter(size=2, position=position_jitter(width=0.15), alpha=.8) +
  scale_x_discrete(limits=c("0vsBulk", "1000vsBulk","0vsES", "1000vsES")) +
  scale_y_continuous(limits = c(0, 0.4)) +
  theme_classic(base_size = 20) +
  theme(axis.ticks=element_line(colour = "black"), 
        axis.text.x =element_text(size= 20, colour = "black", angle = 45, hjust = 1),
        axis.text.y =element_text(size= 25, colour = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size= 25, colour = "black"),
        legend.position = "none") +
  ylab("Weighted UniFrac distance")

ggsave(Fig_S1A, file=paste(Fig, "weigted_uni_sto_Fig_S3A.pdf", sep=""), width=6, height=5, bg="transparent")

# nicotine
uni_nic             <- Uni_table[Uni_table$treatment1=="Nic" & Uni_table$treatment2=="Nt" & Uni_table$group1 %in% c("Mock", "High")  & Uni_table$group2 %in% c("Bulk", "Root"),]      
uni_nic$group1      <- str_replace(uni_nic$group1, "High", "1000")
uni_nic$group1      <- str_replace(uni_nic$group1, "Mock", "0")
uni_nic$group2      <- str_replace(uni_nic$group2, "Root", "ES")
uni_nic$ID          <- paste(uni_nic$group1, uni_nic$group2, sep="vs" )


Fig_S1B <- ggplot(uni_nic , aes(x=ID, y=value)) +
  geom_boxplot(outlier.shape=NA, position=position_dodge()) +
  geom_jitter(size=2, position=position_jitter(width=0.15), alpha=.8) +
  scale_x_discrete(limits=c("0vsBulk", "1000vsBulk","0vsES", "1000vsES")) +
  scale_y_continuous(limits = c(0, 0.4)) +
  theme_classic(base_size = 20) +
  theme(axis.ticks=element_line(colour = "black"), 
        axis.text.x =element_text(size= 20 , colour = "black", angle = 45, hjust = 1),
        axis.text.y =element_text(size= 25, colour = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size= 25, colour = "black"),
        legend.position = "none") +
  ylab("Weighted UniFrac distance")

ggsave(Fig_S1B, file=paste(Fig, "weigted_uni_nic_Fig_S3B.pdf", sep=""), width=6, height=5, bg="transparent")

# dual metabolites
uni_dual             <- Uni_table[Uni_table$treatment1=="dual" & Uni_table$treatment2=="Nt" & Uni_table$group1 %in% c("Mock", "Low")  & Uni_table$group2 %in% c("Bulk", "Root"),]      
uni_dual$group1      <- str_replace(uni_dual$group1, "Low", "Dual500")
uni_dual$group1      <- str_replace(uni_dual$group1, "Mock", "0")
uni_dual$group2      <- str_replace(uni_dual$group2, "Root", "ES")
uni_dual$ID          <- paste(uni_dual$group1, uni_dual$group2, sep="vs" )


Fig_S5C <- ggplot(uni_dual , aes(x=ID, y=value)) +
  geom_boxplot(outlier.shape=NA, position=position_dodge()) +
  geom_jitter(size=2, position=position_jitter(width=0.15), alpha=.8) +
  scale_x_discrete(limits=c("0vsBulk", "Dual500vsBulk","0vsES", "Dual500vsES")) +
  scale_y_continuous(limits = c(0, 0.4)) +
  theme_classic(base_size = 20) +
  theme(axis.ticks=element_line(colour = "black"), 
        axis.text.x =element_text(size= 20 , colour = "black", angle = 45, hjust = 1),
        axis.text.y =element_text(size= 25, colour = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size= 25, colour = "black"),
        legend.position = "none") +
  ylab("Weighted UniFrac distance")

ggsave(Fig_S5C, file=paste(Fig, "weigted_uni_dual_Fig_S7.pdf", sep=""), width=6, height=5, bg="transparent")


# ANOVA and post hoc TukeyHSD
#santhopine
aov_sto <- aov(value ~ ID, uni_sto)
sink(paste(data, "ANOVA_sto.txt", sep=""))
print(anova(aov_sto))
sink()

Tukey_sto <- TukeyHSD(aov_sto)
sink(paste(data, "TukeyHSD_sto.txt", sep=""))
print(Tukey_sto)
sink()

# nicotine
aov_nic <- aov(value ~ ID, uni_nic)
sink(paste(data, "ANOVA_nic.txt", sep=""))
print(anova(aov_nic))
sink()

Tukey_nic <- TukeyHSD(aov_nic)
sink(paste(data, "TukeyHSD_nic.txt", sep=""))
print(Tukey_nic)
sink()

# dual
aov_dual <- aov(value ~ ID, uni_dual)
sink(paste(data, "ANOVA_dual.txt", sep=""))
print(anova(aov_dual))
sink()

Tukey_dual <- TukeyHSD(aov_dual)
sink(paste(data, "TukeyHSD_dual.txt", sep=""))
print(Tukey_dual)
sink()


# permanova test for bacterial community composition
# santhopine
idx             <- match(design$V1[design$treatment=="Sto"], colnames(w_UniFrac))
sto_uni_table   <- w_UniFrac[idx,idx]
idx         <- design$treatment=="Sto"
meta_sto        <- design[idx,c("V1","group")]
identical(colnames(sto_uni_table), meta_sto$V1)

per_sto <-adonis(sto_uni_table~group, meta_sto)
sink(paste(data, "Adonis_sto.txt", sep=""))
print(per_sto)
sink()

# nicotine
idx             <- match(design$V1[design$treatment=="Nic"], colnames(w_UniFrac))
nic_uni_table   <- w_UniFrac[idx,idx]
idx         <- design$treatment=="Nic"
meta_nic        <- design[idx,c("V1","group")]
identical(colnames(nic_uni_table), meta_nic$V1)

per_nic <-adonis(nic_uni_table~group, meta_nic)
sink(paste(data, "Adonis_nic.txt", sep=""))
print(per_nic)
sink()

# Tobacco
idx             <- match(design$V1[design$treatment=="Nt"], colnames(w_UniFrac))
nt_uni_table   <- w_UniFrac[idx,idx]
idx         <- design$treatment=="Nt"
meta_nt        <- design[idx,c("V1","group")]
identical(colnames(nt_uni_table), meta_nt$V1)

per_nt <-adonis(nt_uni_table~group, meta_nt)
sink(paste(data, "Adonis_nt.txt", sep=""))
print(per_nt)
sink()
