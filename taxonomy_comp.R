# R script for drwaing taxonomic composition figs from shimasaki et al., 2021
#orininally from Tomohisa Shimasaki
#shimasaki.tomohisa.45c@st.kyoto-u.ac.jp


# clean up
rm(list=ls())
options(warn=-1)

#load library
library(stringr,   quietly=T, warn.conflicts=F)
library(ggplot2,   quietly=T, warn.conflicts=F)
library(dplyr,     quietly=T, warn.conflicts=F)
library(reshape2,  quietly=T, warn.conflicts=F)
library(scales,      quietly=T, warn.conflicts=F)
library(colourpicker,  quietly=T, warn.conflicts=F)

# directories
data <- "/Users/shimasakitomohisa/Desktop/DataAnalysis/Qiime2/2021_mBio/data/"
list <- "/Users/shimasakitomohisa/Desktop/DataAnalysis/Qiime2/2021_mBio/list/"
Fig  <- "/Users/shimasakitomohisa/Desktop/DataAnalysis/Qiime2/2021_mBio/Fig/"

# load data
order       <- t(as.matrix(read.table(file=paste(data,"level-4.csv", sep=""), sep=",", header=T, check.names=F, row.names=1, stringsAsFactors=F)))
genus       <- t(as.matrix(read.table(file=paste(data,"level-6.csv", sep=""), sep=",", header=T, check.names=F, row.names=1, stringsAsFactors=F)))
ASV         <- as.matrix(read.table(file=paste(data,"ASV-table.tsv", sep=""), sep="\t", header=T, check.names=F, stringsAsFactors=F, row.names=1))
design      <- read.table(file=paste(list,"metadata.tsv", sep=""),  sep="\t", header=F, check.names=F, stringsAsFactors=F)

# refine metadata
design$rep       <- factor(str_sub(design$V1, -1))
design$group     <- str_replace(design$V2, ".*_", "")
design$treatment <- str_replace(design$V2, "_.*", "")

## order level taxonomic composition
## normalization  
ra_order <- apply(order, 2, function(x) x/sum(x))
ra_order100 <-  apply(ra_order, c(1,2), function(x){
  a <- x*100
  return(a)
})
ra_order100 <- data.frame(ra_order100, stringsAsFactors=F)

##taxon average top 20
idx <- design$V1[design$treatment!="KUAS"]
order_selected <- ra_order100[,colnames(ra_order100)%in%idx]
idx <-apply(order_selected, 1, sum)
taxon_av <- data.frame(Av=idx, stringsAsFactors=F)
taxon_av_sort <- as.matrix(taxon_av[order(taxon_av$Av,decreasing = TRUE),,drop = FALSE])

idx_a <- rownames(taxon_av_sort)[1:20]
idx_b <- rownames(taxon_av_sort)[21:356]
taxon_20_table <- data.frame(t(rbind(order_selected[idx_a,], o__Other=apply(order_selected[idx_b,], 2, sum))))

##taxon design
idx <- str_split_fixed(colnames(taxon_20_table), pattern = "__", n=5)
taxon_design <- data.frame(Taxon=colnames(taxon_20_table), phyla=idx[,3], order=idx[,5])
a <- taxon_design[order(taxon_design$phyla),,drop = FALSE]
taxon_design_modi <- read.table(file=paste(list,"taxon_design_modi.txt", sep=""), sep="\t", header=T, check.names=F, stringsAsFactors=F)

##colours
Other <- colorRampPalette(c("#004C00", "#EEFFEE"))
Proteo <- colorRampPalette(c("#004695", "#DDEEFF"))
Actino <- colorRampPalette(c("#C70032", "#FADBDA"))
Bacte  <- colorRampPalette(c("#DB913D", "#FFFFDD"))
Firm   <- c("#8B008B")
taxon_design_modi$colour <- c(Other(6), Proteo(8), Actino(3), Bacte(3), Firm)

##order metadata
idx         <- match(design$V1[design$treatment!="KUAS"], rownames(taxon_20_table))
order_table <- cbind(design[design$treatment!="KUAS",], taxon_20_table[idx,])
order_table <- melt(table) 

idx         <- match(order_table$variable,taxon_design_modi$Taxon)
order_table$order <- taxon_design_modi$order[idx]

##Tobacco
table_tobacco <- order_table[order_table$treatment=="Nt",]
table_tobacco$order    <- factor(table_tobacco$order,    levels=taxon_design_modi$order)
table_tobacco$group    <- factor(table_tobacco$group,    levels=c("Bulk","RS","RP", "Root"))

Fig_1A <- ggplot(table_tobacco, aes(x= rep, y = value, fill = order)) +
  geom_bar(stat = "identity", position = "fill", colour = "black") +
  facet_grid(.~group) +
  theme_minimal(base_size = 20) +
  scale_fill_manual(values=taxon_design_modi$colour) +
  scale_y_continuous(labels = percent, expand = c(0,0)) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 15),
        strip.text  = element_text(size=20),
        legend.position = 'none') +
  labs(x= "", y = "Relative abundance")

ggsave(Fig_1A, file=paste(Fig, "taxa_bar_plot_tobacco_Fig_1A.pdf", sep=""), width=6, height=5, bg="transparent")

##Tobacco
table_tobacco <- order_table[order_table$treatment=="Nt",]
table_tobacco$order    <- factor(table_tobacco$order,    levels=taxon_design_modi$order)
table_tobacco$group    <- factor(table_tobacco$group,    levels=c("Bulk","RS","RP", "Root"))

Fig_1B <- ggplot(table_tobacco, aes(x= rep, y = value, fill = order)) +
  geom_bar(stat = "identity", position = "fill", colour = "black") +
  facet_grid(.~group) +
  theme_minimal(base_size = 20) +
  scale_fill_manual(values=taxon_design_modi$colour) +
  scale_y_continuous(labels = percent, expand = c(0,0)) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 15),
        strip.text  = element_text(size=20),
        legend.position = 'none') +
  labs(x= "", y = "Relative abundance")

ggsave(Fig_1B, file=paste(Fig, "taxa_bar_plot_tobacco_Fig_1B.pdf", sep=""), width=6, height=5, bg="transparent")

##santhopine
table_sto <- order_table[order_table$treatment=="Sto",]
table_sto$order    <- factor(table_sto$order,    levels=taxon_design_modi$order)
table_sto$group 　　<- str_replace(table_sto$group , "Mock", "0")
table_sto$group 　　<- str_replace(table_sto$group , "Low", "50")
table_sto$group 　　<- str_replace(table_sto$group , "Middle", "250")
table_sto$group 　　<- str_replace(table_sto$group , "High", "1000")
table_sto$order    <- factor(table_sto$order,    levels=taxon_design_modi$order)
table_sto$group    <- factor(table_sto$group,    levels=c("0","50","250","1000"))

Fig_1C <- ggplot(table_sto, aes(x= rep, y = value, fill = order)) +
  geom_bar(stat = "identity", position = "fill", colour = "black") +
  facet_grid(.~group) +
  theme_minimal(base_size = 20) +
  scale_fill_manual(values=taxon_design_modi$colour) +
  scale_y_continuous(labels = percent, expand = c(0,0)) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 15),
        strip.text  = element_text(size=20),
        legend.position = 'none') +
  labs(x= "", y = "Relative abundance")
                  
ggsave(Fig_1C, file=paste(Fig, "taxa_bar_plot_sto_Fig_1C.pdf", sep=""), width=6, height=5, bg="transparent")

##Nicotine 
table_nic <- order_table[order_table$treatment=="Nic",]
table_nic$order <- factor(table_nic$order, levels=taxon_design_modi$order)
table_nic$group 　　<- str_replace(table_nic$group , "Mock", "0")
table_nic$group 　　<- str_replace(table_nic$group , "Low", "50")
table_nic$group 　　<- str_replace(table_nic$group , "Middle", "250")
table_nic$group 　　<- str_replace(table_nic$group , "High", "1000")
table_nic$order    <- factor(table_nic$order,    levels=taxon_design_modi$order)
table_nic$group    <- factor(table_nic$group,    levels=c("0","50","250","1000"))

Fig_1D <- ggplot(table_nic, aes(x=rep, y = value, fill = order)) +
  geom_bar(stat = "identity", position = "fill", colour = "black") +
  facet_grid(.~group) +
  theme_minimal(base_size = 20) +
  scale_fill_manual(values=taxon_design_modi$colour) +
  scale_y_continuous(labels = percent, expand = c(0,0)) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 15),
        strip.text  = element_text(size=20),
        legend.position = 'none') +
  labs(x= "", y = "Relative abundance")

ggsave(Fig_1D, file=paste(Fig, "taxa_bar_plot_nic_Fig_1D.pdf", sep=""), width=6, height=5, bg="transparent")

##Dual
table_dual <- order_table[order_table$treatment=="dual",]
table_dual$group 　　<- str_replace(table_dual$group , "Mock", "0")
table_dual$group 　　<- str_replace(table_dual$group , "Low", "Dual500")
table_dual$group 　　<- str_replace(table_dual$group , "High", "Dual1000")
table_dual$order    <- factor(table_dual$order,    levels=taxon_design_modi$order)
table_dual$group    <- factor(table_dual$group,    levels=c("0","Dual500","Dual1000"))

Fig_S5A <- ggplot(table_dual, aes(x= rep, y = value, fill = order)) +
  geom_bar(stat = "identity", position = "fill", colour = "black") +
  facet_grid(.~group) +
  theme_minimal(base_size = 20) +
  scale_fill_manual(values=taxon_design_modi$colour) +
  scale_y_continuous(labels = percent, expand = c(0,0)) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 15),
        strip.text  = element_text(size=12),
        legend.position = 'none') +
  labs(x= "", y = "Relative abundance")

ggsave(Fig_S5A, file=paste(Fig, "taxa_bar_plot_dual_Fig_S6A.pdf", sep=""), width=6, height=5, bg="transparent")


## relative aboundance of genus Arthrobacter
ra_genus <- apply(genus, 2, function(x) x/sum(x))
ra_genus100 <-  apply(ra_genus, c(1,2), function(x){
  a <- x*100
  return(a)
})
ra_genus100 <- data.frame(ra_genus100, stringsAsFactors=F, check.names=F)

art <- "k__Bacteria;p__Actinobacteria;c__Actinobacteria;o__Actinomycetales;f__Micrococcaceae;g__Arthrobacter"
Art_table <- data.frame(t(ra_genus100[art,]), stringsAsFactors=F, check.names=F)
idx <- match(design$V1, rownames(Art_table))
design$Arthrobacter <- Art_table$`k__Bacteria;p__Actinobacteria;c__Actinobacteria;o__Actinomycetales;f__Micrococcaceae;g__Arthrobacter`[idx]
Art <- design %>%
  group_by(V2) %>%
  summarise(mean = mean(Arthrobacter),
            sd   = sd(Arthrobacter))
idx <- match(Art$V2, design$V2)
Art$treatment <- design$treatment[idx]

col   <-               data.frame(names = unique(Art$treatment),
                                  col   = c("#3B3B3B", "#EBEBEB", "#636363", "#ABABAB", "#030303"),
                                  stringsAsFactors=F, check.names=F)

## change names
Art$V2 　　<- str_replace(Art$V2, "Sto_High", "Santhopine 1000")
Art$V2 　　<- str_replace(Art$V2, "Sto_Mock", "Santhopine 0")
Art$V2 　　<- str_replace(Art$V2, "Nic_High", "Nicotine 1000")
Art$V2 　　<- str_replace(Art$V2, "Nic_Mock", "Nicotine 0")
Art$V2 　　<- str_replace(Art$V2, "Nt_Root", "Tobacco ES")
Art$V2 　　<- str_replace(Art$V2, "Nt_Bulk", "Tobacco Bulk")
Art$V2 　　<- str_replace(Art$V2, "KUAS_Gm", "Soybean ES")
Art$V2 　　<- str_replace(Art$V2, "KUAS_Sl", "Tomato ES")
Art$V2 　　<- str_replace(Art$V2, "KUAS_Mc", "Bitter melon ES")
Art$V2 　　<- str_replace(Art$V2, "KUAS_Bulk", "Farm Bulk")
Art$V2 　　<- str_replace(Art$V2, "dual_High", "Dual 1000")
Art$V2 　　<- str_replace(Art$V2, "dual_Low", "Dual 500")
Art$V2 　　<- str_replace(Art$V2, "dual_Mock", "Dual 0")

## Fig 2B
list_2B <- c("Santhopine 0","Santhopine 1000", "Nicotine 0", "Nicotine 1000",
                "Tobacco Bulk", "Tobacco ES", "Farm Bulk","Soybean ES", "Tomato ES", "Bitter melon ES")

Art_2B                 <- Art[Art$V2 %in% list_2B,]
Art_2B$V2              <-factor(Art_2B$V2,    levels=list_2B)

Fig_2B <- ggplot(Art_2B, aes(x = V2, y = mean, fill = treatment)) +
  geom_bar(stat = "identity", colour = "black") +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd, width = 0.3)) +
  theme_classic(base_size = 20) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 25)) +
  scale_fill_manual(values = col$col[col$names %in% unique(Art_2B$treatment)])+
  labs( x= "",y = "Relative abundance (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = "black", size = 18),
        axis.text.y = element_text(colour = "black", size = 15),
        axis.text   = element_text(colour = "black"),
        legend.position = 'none') 

ggsave(Fig_2B, file=paste(Fig, "Art_relative_abundance_Fig_2B.pdf", sep=""), width=6, height=5, bg="transparent")

                  
## relative aboundance of genus Arthrobacter in metabolite-treated soils
## Fig S2C (santhopine)                  
Art_S5B$V2              <-factor(Art_S5B$V2,    levels=list_S5B )
list_sto                 <- c("Santhopine 0","Santhopine 50", "Santhopine 250", "Santhopine 1000")
Art_sto                  <- Art[Art$V2 %in% list_sto,]
Art_sto$V2               <-factor(Art_sto$V2,    levels=list_sto)

Fig_sto <- ggplot(Art_sto, aes(x = V2, y = mean, fill = treatment)) +
  geom_bar(stat = "identity", colour = "black") +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd, width = 0.3)) +
  theme_classic(base_size = 20) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 25)) +
  scale_fill_manual(values = col$col[col$names %in% unique(Art_sto$treatment)])+
  labs( x= "",y = "Relative abundance of \nMicrococcaceae (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = "black", size = 15),
        axis.text.y = element_text(colour = "black", size = 15),
        axis.text   = element_text(colour = "black"),
        legend.position = 'none') 
Fig_sto

ggsave(Fig_sto, file=paste(Fig, "abundance_micro_sto.pdf", sep=""), width=6, height=5, bg="transparent")

## Fig S2D (Nicotine)
list_nic                 <- c("Nicotine 0","Nicotine 50", "Nicotine 250", "Nicotine 1000")
Art_nic                  <- Art[Art$V2 %in% list_nic,]
Art_nic$V2               <-factor(Art_nic$V2,    levels=list_nic)

Fig_nic <- ggplot(Art_nic, aes(x = V2, y = mean, fill = treatment)) +
  geom_bar(stat = "identity", colour = "black") +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd, width = 0.3)) +
  theme_classic(base_size = 20) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 25)) +
  scale_fill_manual(values = col$col[col$names %in% unique(Art_nic$treatment)])+
  labs( x= "",y = "Relative abundance of \nMicrococcaceae (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = "black", size = 15),
        axis.text.y = element_text(colour = "black", size = 15),
        axis.text   = element_text(colour = "black"),
        legend.position = 'none') 
Fig_nic

ggsave(Fig_nic, file=paste(Fig, "abundance_micro_nic.pdf", sep=""), width=6, height=5, bg="transparent")

## Fig S5B (Dual)
list_S5B 　　　　　　　　　<- c("Dual 0","Dual 500", "Dual 1000")
Art_S5B                 <- Art[Art$V2 %in% list_S5B,]
Art_S5B$V2              <-factor(Art_S5B$V2,    levels=list_S5B )

Fig_S5B <- ggplot(Art_S5B, aes(x = V2, y = mean, fill = treatment)) +
  geom_bar(stat = "identity", colour = "black") +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd, width = 0.3)) +
  theme_classic(base_size = 20) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 25)) +
  scale_fill_manual(values = col$col[col$names %in% unique(Art_S6B$treatment)])+
  labs( x= "",y = "Relative abundance (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = "black", size = 18),
        axis.text.y = element_text(colour = "black", size = 15),
        axis.text   = element_text(colour = "black"),
        legend.position = 'none') 
                  
ggsave(Fig_S5B, file=paste(Fig, "Art_relative_abundance_Fig_S6B.pdf", sep=""), width=6, height=5, bg="transparent")
