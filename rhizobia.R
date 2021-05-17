# R script for pan-genome analysis using orthologous groups from shimasaki et al., 2021
# orininally from Tomohisa Shimasaki
# shimasaki.tomohisa.3e@kyoto-u.ac.jp


# cleanup
rm(list = ls())
options(warn=1)

#load library
library(ggplot2,      quietly=T, warn.conflicts=F)
library(dplyr,        quietly=T, warn.conflicts=F)
library(reshape2,     quietly=T, warn.conflicts=F)

# directories
list       <- "/Volumes/Tomohisa/GenomeAnalysis/BLAST/Result/Rhizobia/list/"
data       <- "/Volumes/Tomohisa/GenomeAnalysis/BLAST/Result/Rhizobia/data/"
fig        <- "/Volumes/Tomohisa/GenomeAnalysis/BLAST/Result/Rhizobia/fig/"
results    <- "/Volumes/Tomohisa/GenomeAnalysis/BLAST/Result/Rhizobia/results/"

# load data
OG_tab        <-  read.table(file=paste(data,"Orthogroups.GeneCount.txt", sep=""), header=T, stringsAsFactors=F, sep="\t", check.names=F, row.names=1)
OG_1021       <-  read.table(file=paste(list,"OG_lD_1021_e5.txt", sep=""), sep="\t", check.names=F,stringsAsFactors=F)   
OG_3841       <-  read.table(file=paste(list,"OG_lD_3841_e5.txt", sep=""), sep="\t", check.names=F,stringsAsFactors=F)  
design        <-  read.table(file=paste(list,"strain_list.txt", sep=""), header=T,sep="\t", check.names=F,stringsAsFactors=F)

##make table of rhizobial genes 
OG_tab[OG_tab >1] <-1
rhizo_tab     <- as.matrix(OG_tab[unique(c(OG_1021$V1,OG_3841$V1)),])
rhizo_tab     <- as.data.frame(t(rhizo_tab))

write.table(rhizo_tab, file=paste(results, "rhizio_OG_table_e5.txt", sep=""), quote=F, sep=",", row.names=T, col.names=T)


ratio   <- apply(rhizo_tab, 1, function(x){
       ratio <- sum(x)/ncol(rhizo_tab)
       return(ratio)
})

idx <- match(design$Strain, rownames(rhizo_tab))
design$ratio <- ratio[idx]


origin       <- data.frame(name   = c("Plant", "Soil", "Other"),
                           shapes = c(18, 16, 8),
                           col    = c("#8fc93e","#ef7318","#10677f"),
                           check.names=F, stringsAsFactors=F)
lineage      <- data.frame(name   = c("A", "B", "C"),
                           col = c("#e23635","#187b3b","#33418d"),
                           check.names=F, stringsAsFactors=F)

design$Origin  <- factor(design$Origin, levels=origin$name)

## drawing fig
Fig5A <- ggplot(design, aes(x=sublineages, y=ratio, color=sublineages)) +
  geom_boxplot(size=1.2, outlier.shape = NA, width=0.8) +
  scale_x_discrete(limit=c("A","B","C")) +
  scale_y_continuous(limits = c(0.1, 1)) +
  geom_jitter(size=1.5, alpha = 0.8, position=position_jitter(width=0.15))+
  scale_color_manual(values=lineage$col) +
  theme_classic(base_size = 20) +
  labs( x= "",y = "Presence ratio") +
  theme(axis.text.x=element_text(size = 20, colour = "black"),
        axis.text.y=element_text(size = 20, colour = "black"),
        axis.title.y=element_text(size = 22, colour = "black"),
        legend.position = "top",
        legend.title=element_blank())
Fig5A
ggsave(Fig5A, file=paste(fig, "box_plot_lineage_e5.pdf", sep=""), width=5, height=5, bg="transparent")


Fig5B <- ggplot(design, aes(x=Origin, y=ratio, color=Origin)) +
  geom_boxplot(size=1.2, outlier.shape = NA, width=0.8) +
  geom_jitter(size=1.5, alpha = 0.8, position=position_jitter(width=0.15))+
  scale_color_manual(values=origin$col) +
  scale_y_continuous(limits = c(0.1, 1)) +
  theme_classic(base_size = 20) +
  labs( x= "",y = "Presence ratio") +
  theme(axis.text.x=element_text(size = 20, colour = "black"),
        axis.text.y=element_text(size = 20, colour = "black"),
        axis.title.y=element_text(size = 22, colour = "black"),
        legend.position = "top",
        legend.title=element_blank())
Fig5B
ggsave(Fig5B, file=paste(fig, "box_plot_origin_e5.pdf", sep=""), width=5, height=5, bg="transparent")


##remove high similarity genomes
out           <- matrix(0, nrow = 99 , ncol = 99)
colnames(out) <- colnames(OG_tab)
rownames(out) <- colnames(OG_tab)

list <- colnames(OG_tab)

for (i in list){
  
  idx_A  <-  OG_tab[i]
  idx_A  <-  idx_A %>%
    filter(idx_A[1]>0)
  list_A  <- rownames(idx_A)
  
  for (j in list){
    idx_B  <-  OG_tab[j]
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
idx <- match(design$Strain, list_strain)
r_design    <- design %>% filter(!design$Strain %in% list_strain)

## wilcox
ori_wilcox <-pairwise.wilcox.test(design$ratio, design$Origin, p.adj="bonferroni", exact=F)
sink(paste(results, "wilcox_origin_e30.txt", sep=""))
print(ori_wilcox)
sink()

lineage_wilcox <- pairwise.wilcox.test(design$ratio, design$sublineages, p.adj="bonferroni", exact=F)
sink(paste(results, "wilcox_lineage_e30.txt", sep=""))
print(lineage_wilcox)
sink()
