# an R script to reconstruct ancestral status of the most recent common ancestor 
# https://www.rdocumentation.org/packages/ape/versions/5.4/topics/ace
# originally by Ryohei Thomas Nakano and Tomohisa Shimasaki
# edited by Tomohisa shimasaki  
# shimasaki.tomohisa.3e@kyoto-u.ac.jp

rm(list=ls())
options(warn=-1)

library(ape,    quietly=T, warn.conflicts=F)
library(furrr,  quietly=T, warn.conflicts=F)
library(gplots, quietly=T, warn.conflicts=F)

# directories
data    <- "/Volumes/Tomo/github/2021_mBio/ACE/data/"
result  <- "/Volumes/Tomo/github/2021_mBio/ACE/result/"
fig     <- "/Volumes/Tomo/github/2021_mBio/ACE/fig/"

# load data
Go     <- as.matrix(read.table(file=paste(data, "Orthogroups.GeneCount.txt", sep=""), header=T, row.names=1, sep="\t", check.names=F,stringsAsFactors=F))
Ko     <- as.matrix(read.table(file=paste(data, "KEGG_matrix.txt", sep=""), header=T, row.names=1, sep="\t", check.names=F,stringsAsFactors=F))
tree   <- read.tree(file=paste(data, "amphora_tree.txt", sep=""))

Go[Go >1] <-1
tree <- multi2di(tree)

# replace zero length
idx <- tree$edge.length == 0
tree$edge.length[idx] <- min(tree$edge.length[!idx])*1e-6

# out tree with node labels
pdf(paste(result, "tree_nodelabel.pdf", sep=""), width=18, height=18)
plot(tree, type="phylogram", use.edge.length=FALSE, show.tip.label=TRUE, label.offset=5)
nodelabels(paste("node", 1:tree$Nnode, sep=""))
dev.off()

# remove common absence/presence
# Go
idx <- rowSums(Go == 1) %in% c(0, ncol(Go))
Go <- Go[!idx,]
idx <- match(tree$tip.label, colnames(Go))
Go <- Go[, idx]

# Ko
idx <- rowSums(Ko == 1) %in% c(0, ncol(Ko))
Ko <- Ko[!idx,]
idx <- match(tree$tip.label, colnames(Ko))
Ko <- Ko[, idx]

col <- c("0"="black", "1"="green4")

# ACE 
plan(multiprocess(workers=40))

# Go
out_list_Go <- future_map(1:nrow(Go), function(i){
  x <- Go[i,]

ace <- ace(x=x, phy=tree, type="discrete")

lik_root <- ace$lik.anc[1, 2]
lik_A <- ace$lik.anc[69,2]
lik_B <- ace$lik.anc[27,2]
lik_C <- ace$lik.anc[2,2]

out <- data.frame(Go=rownames(Go)[i], root=lik_root, A=lik_A, B=lik_B, C=lik_C)

pdf(paste(result, "ind_tree/tree_ace.", rownames(Go)[i], ".pdf", sep=""), width=18, height=18)
plot(tree, type="phylogram", use.edge.length=FALSE, show.tip.label=TRUE, label.offset=5)
nodelabels(round(ace$lik.anc[, "1"], 3), col="white", pie=ace$lik.anc, piecol=col, cex=0.5)
tiplabels(pie=as.matrix(data.frame("0"=1-x, "1"=x, check.names=F)), piecol=col, cex=0.3, offset=2.5)
dev.off()

return(out)

})	

ace_lik_Go <- do.call(rbind, out_list_Go)

# listup ACE at each sublinage
Ace_Go   <- ace_lik_Go$Go[ace_lik_Go$root  > .8]
Ace_Go_A <- ace_lik_Go$Go[ace_lik_Go$A > .8]
Ace_Go_B <- ace_lik_Go$Go[ace_lik_Go$B > .8]
Ace_Go_C <- ace_lik_Go$Go[ace_lik_Go$C > .8]

# venn diagram analysis
venn_Go <- list(Ace_Go_A, Ace_Go_B, Ace_Go_C)
pdf(paste(fig, "venn_Go.pdf", sep=""))
venn(venn_Go)
dev.off()

Uniq_A <- setdiff(setdiff(Ace_Go_A, Ace_Go_B), Ace_Go_C)
Uniq_B <- setdiff(setdiff(Ace_Go_B, Ace_Go_A), Ace_Go_C)
Uniq_C <- setdiff(setdiff(Ace_Go_C, Ace_Go_A), Ace_Go_B)
Uniq_AB <- setdiff(intersect(Ace_Go_A, Ace_Go_B), Ace_Go_C)

write.table(Uniq_A, file=paste(result,"lineage_A_Uniq-ACE_Go.txt", sep=""), sep="\t", col.names=T, row.names=F, quote=F)
write.table(Uniq_B, file=paste(result,"lineage_B_Uniq-ACE_Go.txt", sep=""), sep="\t", col.names=T, row.names=F, quote=F)
write.table(Uniq_C, file=paste(result,"lineage_C_Uniq-ACE_Go.txt", sep=""), sep="\t", col.names=T, row.names=F, quote=F)
write.table(Uniq_AB, file=paste(result,"lineage_AB_Uniq-ACE_Go.txt", sep=""), sep="\t", col.names=T, row.names=F, quote=F)
write.table(ace_lik,  file=paste(result,"ACE_likelihood_list.txt", sep=""), sep="\t", col.names=T, row.names=F, quote=F)

# Ko
out_list_Ko <- future_map(1:nrow(Ko), function(i){
  x <- Ko[i,]
  
  ace <- ace(x=x, phy=tree, type="discrete")
  
  lik_root <- ace$lik.anc[1, 2]
  lik_A <- ace$lik.anc[69,2]
  lik_B <- ace$lik.anc[27,2]
  lik_C <- ace$lik.anc[2,2]
  
  out <- data.frame(Ko=rownames(Ko)[i], root=lik_root, A=lik_A, B=lik_B, C=lik_C)
  
  pdf(paste(result, "ind_tree/tree_ace.", rownames(Ko)[i], ".pdf", sep=""), width=18, height=18)
  plot(tree, type="phylogram", use.edge.length=FALSE, show.tip.label=TRUE, label.offset=5)
  nodelabels(round(ace$lik.anc[, "1"], 3), col="white", pie=ace$lik.anc, piecol=col, cex=0.5)
  tiplabels(pie=as.matrix(data.frame("0"=1-x, "1"=x, check.names=F)), piecol=col, cex=0.3, offset=2.5)
  dev.off()
  
  return(out)
  
})	

ace_lik_Ko <- do.call(rbind, out_list_Ko)

# listup ACE at each sublinage
Ace_Ko   <- ace_lik_Ko$Ko[ace_lik_Ko$root  > .8]
Ace_Ko_A <- ace_lik_Ko$Ko[ace_lik_Ko$A > .8]
Ace_Ko_B <- ace_lik_Ko$Ko[ace_lik_Ko$B > .8]
Ace_Ko_C <- ace_lik_Ko$Ko[ace_lik_Ko$C > .8]

# venn diagram analysis
venn_Ko <- list(Ace_Ko_A, Ace_Ko_B, Ace_Ko_C)
pdf(paste(fig, "venn_Ko.pdf", sep=""))
venn(venn_Ko)
dev.off()

Uniq_A <- setdiff(setdiff(Ace_Ko_A, Ace_Ko_B), Ace_Ko_C)
Uniq_B <- setdiff(setdiff(Ace_Ko_B, Ace_Ko_A), Ace_Ko_C)
Uniq_C <- setdiff(setdiff(Ace_Ko_C, Ace_Ko_A), Ace_Ko_B)
Uniq_AB <- setdiff(intersect(Ace_Ko_A, Ace_Ko_B), Ace_Ko_C)

write.table(Uniq_A, file=paste(result,"lineage_A_Uniq-ACE_Ko.txt", sep=""), sep="\t", col.names=T, row.names=F, quote=F)
write.table(Uniq_B, file=paste(result,"lineage_B_Uniq-ACE_Ko.txt", sep=""), sep="\t", col.names=T, row.names=F, quote=F)
write.table(Uniq_C, file=paste(result,"lineage_C_Uniq-ACE_Ko.txt", sep=""), sep="\t", col.names=T, row.names=F, quote=F)
write.table(Uniq_AB, file=paste(result,"lineage_AB_Uniq-ACE_Ko.txt", sep=""), sep="\t", col.names=T, row.names=F, quote=F)
write.table(ace_lik_Ko,  file=paste(result,"ACE_likelihood_list.txt", sep=""), sep="\t", col.names=T, row.names=F, quote=F)
