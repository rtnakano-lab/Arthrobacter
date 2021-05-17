# R script for analyzing MiSeq data from Tomohisa Shimasaki, Kyoto Univ
#
# Enrichment test with GLM
#
# originally by Ruben Garrido-Oter
# edited by Ryohei Thomas Nakano and Tomohisa Shimasaki
# nakano@mpipz.mpg.de
# shimasaki.tomohisa.3e@kyoto-u.ac.jp

options(warn=-1)

# clean up
rm(list=ls())

# load packages
library(edgeR,      quietly=T, warn.conflicts=F)
library(stringr,    quietly=T, warn.conflicts=F)
library(gplots,     quietly=T, warn.conflicts=F)
library(RColorBrewer, quietly=T, warn.conflicts=F)

# directories
data           <- "/Volumes/Tomo/github/2021_mBio/glm/data/"
results        <- "/Volumes/Tomo/github/2021_mBio/glm/results/"
list           <- "/Volumes/Tomo/github/2021_mBio/glm/list/"
fig            <- "/Volumes/Tomo/github/2021_mBio/glm/fig/"

# parameters
p.adj.method <- "fdr"
alpha <- 0.05

# data import
design  <-             read.table(paste(list, "metadata.tsv", sep=""), sep="\t", header=F, stringsAsFactors=F)
tab     <- t(as.matrix(read.table(paste(data, "level-6.csv", sep=""),  sep=",", header=T, stringsAsFactors=F, row.names=1, as.is=T)))

# refine metadata
design$rep       <- factor(str_sub(design$V1, -1))

# levels
idx <- match(colnames(tab), design$V1)
treatment_labels <- design$V2[idx]

# create design matrix 
treatment <- treatment_labels
replicate <- factor(design$rep, levels=c(1:4))

model <- model.matrix(~ 0 + treatment + replicate)
colnames(model) <- str_replace(colnames(model), "treatment", "")

# remove missing design columns
idx <- colSums(model) == 0
model <- model[, !idx]

# create contrast matrix
contrasts.mat <- makeContrasts(
     Dual_low   = c(dual_Low  - dual_Mock),
     Dual_High  = c(dual_High - dual_Mock),
     Santhopine = c(Sto_High  - Sto_Mock ),
     Nicotine   = c(Nic_High  - Nic_Mock ),
     root_RP    = c(Nt_Root   - Nt_RP    ),
     KUAS_Gm    = c(KUAS_Gm   - KUAS_Bulk),
     KUAS_Mc    = c(KUAS_Mc   - KUAS_Bulk),
     KUAS_Sl    = c(KUAS_Sl   - KUAS_Bulk),
     levels=model)

contrasts.names <- attr(contrasts.mat, "dimnames")$Contrasts
n <- length(contrasts.names)

# create DGEList object for edgeR
d <- DGEList(counts=tab, group=treatment)
d <- calcNormFactors(d)


# Estimate common and tag-wise disperse
d2 <- estimateGLMCommonDisp(d, model)
d2 <- estimateGLMTagwiseDisp(d2, model)

fit <- glmFit(d2, model)

# LRT for each contrasts
LRT.list <- lapply(contrasts.names, function(x) glmLRT(fit, contrast=contrasts.mat[, which(contrasts.names == x)]))
names(LRT.list) <- contrasts.names

# logFC and PValue tables
logFC_P.list <- lapply(1:n, function(x) {
	table <- LRT.list[[x]]$table[,c(1,4)]
	table$PValue <- p.adjust(table$PValue, method=p.adj.method)
	colnames(table) <- paste(contrasts.names[x], colnames(table), sep="_")
	return(table)
	})
logFC_P <- do.call(data.frame, logFC_P.list)
write.table(logFC_P, file=paste(results, "logFC.P.genus.txt", sep=""), sep="\t", quote=F, col.names=NA, row.names=T)


# Significance picking for each tested model
DE.list <- lapply(contrasts.names, function(x) decideTestsDGE(LRT.list[[which(contrasts.names == x)]], adjust.method=p.adj.method, p.value=alpha))
names(DE.list) <- contrasts.names

# Number of significant differentially abundant OTUs
total    <- sapply(DE.list, function(x) sum(abs(x)))
enriched <- sapply(DE.list, function(x) sum(x ==  1))
depleted <- sapply(DE.list, function(x) sum(x == -1))
count <- data.frame(total, enriched, depleted)
write.table(count, file=paste(results, "number_of_enrichd_depleted_genus.txt", sep=""), quote=F, row.names=T, col.names=NA, sep="\t")


# significance table
DE <- sapply(1:n, function(x) DE.list[[x]][,1])
colnames(DE) <- contrasts.names
write.table(DE, file=paste(results, "enrichment_test.significance_table.genus.txt", sep=""), sep="\t", quote=T, row.names=T, col.names=NA)

idx <- rowSums(abs(DE)) > 0
write.table(DE[idx,], file=paste(results, "enrichment_test.significance_table.sig_only.genus.txt", sep=""), sep="\t", quote=T, row.names=T, col.names=NA)

log2cpm <- cpm(d2, prior.count=2, log=TRUE)
write.table(log2cpm, file=paste(results, "log2cpm.genus.txt", sep=""), sep="\t", quote=F, col.names=NA, row.names=T)


##logFC_data
logFC_P_d <- data.frame(taxnomy=row.names(logFC_P), logFC_P, stringsAsFactors=F)
rownames(logFC_P_d) <- NULL 

##Tobacco enrhich
Tobacco_En <- logFC_P_d[(logFC_P$root_RP_PValue<alpha)&(logFC_P$root_RP_logFC>0),]
write.table(Tobacco_En , file=paste(results, "Tobacco_Enrich.txt", sep=""), quote=F, row.names=T, col.names=NA, sep="\t")

##Santhopine enrhich
Sto_En <- logFC_P_d[(logFC_P$Santhopine_PValue<alpha)&(logFC_P$Santhopine_logFC>0),]
write.table(Sto_En , file=paste(results, "Sto_Enrich.txt", sep=""), quote=F, row.names=T, col.names=NA, sep="\t")

##Nicotine enrhich
Nic_En <- logFC_P_d[(logFC_P$Nicotine_PValue<alpha)&(logFC_P$Nicotine_logFC>0),]
write.table(Nic_En , file=paste(results, "Nic_Enrich.txt", sep=""), quote=F, row.names=T, col.names=NA, sep="\t")

##Dual enrhich
Dual_En <- logFC_P_d[(logFC_P$Dual_low_PValue<alpha)&(logFC_P$Dual_low_logFC>0),]
write.table(Dual_En , file=paste(results, "Dual_Enrich.txt", sep=""), quote=F, row.names=T, col.names=NA, sep="\t")

##plot heatmap 
Enrich <-　rbind(Sto_En, Nic_En, Tobacco_En)

df_en <- data.frame(taxnomy=Enrich[,1],  Sto=Enrich$Santhopine_logFC, Nic=Enrich$Nicotine_logFC, Tobacco=Enrich$root_RP_logFC, stringsAsFactors=F)
df_en <- df_en[!duplicated(df_en$taxnomy), ] 

pdf(paste(fig, "Fig2B.pdf", sep=""))
scale_range = seq(-9,9,0.1)
heatmap.2(as.matrix(df_en[, -1]), labRow = "", labCol = "", col=colorRampPalette(c("#113285","white","#CB4042"))(length(scale_range)-1), scale="none",
          key=TRUE, symkey=FALSE, density.info="none", trace="none", Colv=NA, Rowv=NA,
          cexRow=1.2, cexCol=1.8, margin=c(6,15), main="",srtCol=0, adjCol = c(NA,1),
          key.xlab="log2FC", key.title="")
dev.off()
