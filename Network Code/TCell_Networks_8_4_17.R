tcell.networks.only <- load("~/Downloads/geneset.data_playground.p900_consensus.radius_RNA_explicit_dirnet_redundant.RData")
tdiff <- path.all
network.data.tcells <- load("~/Desktop/TCell_DiffData/geneset.data_playground.p900_consensus.radius_RNA_explicit_fullnet_redundant.RData")
# path.all <- path.all[-which(path.all$schema == "tcelldiff"), ]
path.all <- rbind(path.all, tdiff)
path.all.negative <- path.all[path.all$HY_HO.mean_deltaProp < 0, ]
path.all.negative <- path.all.negative[abs(path.all.negative$HY_HO.mean_deltaProp) > quantile(abs(path.all.negative$HY_HO.mean_deltaProp),prob=1-75/100),]

path.all.positive <- path.all[path.all$HY_HO.mean_deltaProp > 0, ]
path.all.positive <- path.all.positive[abs(path.all.positive$HY_HO.mean_deltaProp) > quantile(abs(path.all.positive$HY_HO.mean_deltaProp),prob=1-75/100),]

# path.all.negative <- path.all[path.all$logFC.rna < 0, ]
# path.all.negative <- path.all.negative[abs(path.all.negative$logFC.rna) > quantile(abs(path.all.negative$logFC.rna),prob=1-75/100),]
# 
# path.all.positive <- path.all[path.all$logFC.rna > 0, ]
# path.all.positive <- path.all.positive[abs(path.all.positive$logFC.rna) > quantile(abs(path.all.positive$logFC.rna),prob=1-75/100),]


path.all <- rbind(path.all.positive, path.all.negative)
# path.all <- path.all[which(path.all$PValue.rna < 0.01), ]

atac_do.tcelldiff_enrichment <- read.csv(file="~/Downloads/atac_do.tcelldiff_enrichment.txt", sep = "\t", header = T, stringsAsFactors =F)

path.all <- path.all[path.all$schema %in% c("imm_coarse", "imm_fine"), ]
path.imm_fine <- path.all[path.all$schema %in% "imm_fine", ]
path.imm_fine <- path.imm_fine[path.imm_fine$pathway %in% atac_do.imm_fine_enrichment[which(atac_do.imm_fine_enrichment$hypergeom.fdr <=  quantile(abs(atac_do.imm_fine_enrichment$hypergeom.fdr),prob=1-95/100)), "Module.Name"], ]
path.imm_coarse <- path.all[path.all$schema %in% "imm_coarse", ]
path.imm_coarse <- path.imm_coarse[path.imm_coarse$pathway %in% atac_do.imm_coarse_enrichment[which(atac_do.imm_coarse_enrichment$hypergeom.fdr <= quantile(abs(atac_do.imm_coarse_enrichment$hypergeom.fdr),prob=1-95/100)), "Module.Name"], ]
path.all <- rbind(path.imm_coarse, path.imm_fine)

# path.all.five <- unique(path.all$pathway)
path.all.one <- unique(path.all$pathway)

clipboard_genes <- function(x){
  current.module <- unique(path.all$pathway)[x]
  fp.data.example <- path.all[path.all$pathway %in% current.module, ]
  
  fp.data.example <- fp.data.example[abs(fp.data.example$signed_logOR) > quantile(abs(fp.data.example$signed_logOR), prob=1-50/100), ]
    data <- unique(c(fp.data.example$tf.gene, fp.data.example$GeneName))
    clip <- pipe("pbcopy", "w")                       
    write.table(data, file=clip, quote = F, row.names = F)                               
    close(clip)
}


for(geneset.name in unique(path.all$schema)){
  # geneset.name <- unique(path.all$schema)[1]
  # since data is presented differently for some of the gene sets, individual debugging may be necessary
  # geneset.name <- "reactome"
  # geneset.name <- "tcelldiff"
  all.info <- path.all[path.all$schema %in% geneset.name, ]
  
  # atac.enrichment <- get(paste("atac_do.", geneset.name, "_enrichment", sep = ""))
  # rna.enrichment <- get(paste("rna_do.", geneset.name, "_enrichment", sep = ""))
  
  # all.info <- all.info[all.info$pathway %in% atac.enrichment[which(atac.enrichment$hypergeom.fdr <= 0.05), "Module.Name"], ]
  # 
  # print(paste("Working on geneset", geneset.name))
  # current_geneset <- test.table.tcell$GeneName
  
  for(current.module in unique(all.info$pathway)){
    current.module <- unique(all.info$pathway)[4]
    # current.module <- "Cellular Senescence"
    fp.data.example <- all.info[all.info$pathway %in% current.module, ]

    fp.data.example <- fp.data.example[abs(fp.data.example$signed_logOR) > quantile(abs(fp.data.example$signed_logOR), prob=1-80/100), ]

    if(nrow(fp.data.example) > 2){
      data <- unique(c(fp.data.example$tf.gene, fp.data.example$GeneName))
      clip <- pipe("pbcopy", "w")                       
      write.table(data, file=clip, quote = F, row.names = F)                               
      close(clip)
      
      df <- data.frame(fp.data.example %>% select(tf.name,GeneName), stringsAsFactors = FALSE)
      # ***If you want to differentiate between proximal and distal, use:***
      # df <- data.frame(fp.data.example %>% select(tf_distype, GeneName), stringsAsFactors = FALSE)
      df <- df[!duplicated(df),]
      fp.data.example <- fp.data.example[row.names(df), ]
      nw <- graph_from_data_frame(df,directed=T)

      edge_attr(nw,"deltaProp") <- fp.data.example$HY_HO.mean_deltaProp
      edge_attr(nw,"deltaProp_color")[which(edge_attr(nw,"deltaProp") <= 0)] <- "blue"
      edge_attr(nw,"deltaProp_color")[which(edge_attr(nw,"deltaProp") > 0)] <- "red"

      gene.names <- unique(df$GeneName)
      tf.names <- unique(df$tf.name)

      vertex_attr(nw, "logFC")[which(vertex_attr(nw, "name") %in% gene.names)] <- unique(fp.data.example[fp.data.example$GeneName %in% gene.names, "logFC.rna"])
      vertex_attr(nw, "logFC")[which(vertex_attr(nw, "name") %in% tf.names)] <- rep(0, length(vertex_attr(nw, "logFC")[which(vertex_attr(nw, "name") %in% tf.names)]))

      vertex_attr(nw,"logFC_color")[which(vertex_attr(nw,"logFC") < 0)] <- "blue"
      vertex_attr(nw,"logFC_color")[intersect(which(vertex_attr(nw,"logFC") < 0), which(vertex_attr(nw, "name") %in% test.table.tcell$GeneName))] <- "cyan"
      vertex_attr(nw,"logFC_color")[which(vertex_attr(nw,"logFC") > 0)] <- "red"
      vertex_attr(nw,"logFC_color")[intersect(which(vertex_attr(nw,"logFC") > 0), which(vertex_attr(nw, "name") %in% test.table.tcell$GeneName))] <- "magenta"
      # vertex_attr(nw, "color")[!vertex_attr(nw, "name") %in% df$tf.name] <- "white"
      vertex_attr(nw, "is.tf") <- vertex_attr(nw,"name") %in% tf.names
      vertex_attr(nw,"logFC_color")[vertex_attr(nw,"is.tf")] <- "white"

      nw.gplot <- ggnet2(nw,
                         node.size = rescale(abs(vertex_attr(nw,"logFC")), c(5, 30)),
                         # node.shape = vertex_attr(nw,"shape"),
                         label = F, # I use ggplot (below) to add labels
                         mode = "fruchtermanreingold",
                         arrow.size = 5,
                         arrow.type = "open",
                         arrow.gap = 0.025,
                         edge.color = edge_attr(nw,"deltaProp_color"),
                         edge.size = rescale(abs(edge_attr(nw, "deltaProp")), c(0.1, 2)),
                         node.alpha = 0.8,
                         node.color = unlist(vertex_attr(nw,"logFC_color")),
                         legend.position = "none")

      # Then just plot with p
      pdf(paste("~/Desktop/relevant_immgen/", geneset.name, "_", gsub("[[:punct:]]", " ", gsub(",.*$", "", current.module)), ".pdf", sep = ""))
      p <- nw.gplot +
        # This replots the nodes, adding edges to the nodes drawn by ggnet2
        # I have included an attribute “is.tf” that informs whether a node is a TF or a target, and use this to plot TFs and targets with different shapes and sizes, and to format the labels differently
        geom_point(shape=as.numeric(sub("TRUE",22,sub("FALSE",21,vertex_attr(nw,"is.tf")))),
                   size=as.numeric(sub("TRUE",8,sub("FALSE",9,vertex_attr(nw,"is.tf")))),
                   color=sub("TRUE","maroon2",sub("FALSE","white",vertex_attr(nw,"is.tf"))),alpha=0.6) +
        # This adds labels to the nodes
        geom_text(aes(label=label),
                  size=as.numeric(sub("FALSE",2.5,sub("TRUE",2,vertex_attr(nw,"is.tf")))),
                  color=sub("FALSE","gray40",sub("TRUE","maroon4",vertex_attr(nw,"is.tf"))),hjust=-0.2)
      print(p)
      dev.off()

      print(paste("Just finished network for ", current.module))
    }
  }
}

#######################
#
# Comparison of expression results 
#
#######################

immgen_geneset <- read.csv(file='~/Desktop/immgen_geneset.csv')

clipboard_genes <- function(x){
  current.module <- unique(tdiff$pathway)[3]
  fp.data.example <- tdiff[tdiff$pathway %in% current.module, ]
  
  # fp.data.example <- fp.data.example[abs(fp.data.example$signed_logOR) > quantile(abs(fp.data.example$signed_logOR), prob=1-50/100), ]
  tf.genes <- as.character(fp.data.example$tf.gene)
  tf.genes <- gsub("_.*", "", tf.genes)
  data <- unique(c(tf.genes, as.character(fp.data.example$GeneName)))
  clip <- pipe("pbcopy", "w")                       
  write.table(data, file=clip, quote = F, row.names = F)                               
  close(clip)
}

tdiff_geneset <- read.csv(file = "~/Desktop/tdiff_geneset.csv")
overlapping_genes <- as.character(tdiff_geneset[which(tdiff_geneset$Gene %in% immgen_geneset$Gene), "Gene"])

immgen_geneset <- immgen_geneset[!duplicated(immgen_geneset$Gene),]
rownames(tdiff_geneset) <- tdiff_geneset$Gene
rownames(immgen_geneset) <- immgen_geneset$Gene

futher_set <- tdiff_geneset[overlapping_genes, ]
further_tfs <- tdiff_geneset[overlapping_genes %in% tf.genes, ]

tdiff <- tdiff[!duplicated(tdiff), ]
tdiff$tf.genename <- gsub("_.*", "", tdiff$tf.gene)

further_tfs <- merge(further_tfs, tdiff, by.x="Gene", by.y="tf.genename" )
colnames(further_tfs)[1] <- "tf.GeneName"
#################
#
# Figure illustrating 
# correlation between 
# diff + age
#
#################


diff.mx <- read.table(file="~/Desktop/TCell_DiffData/tcell_diff.txt",header = TRUE,sep = "\t",stringsAsFactors = F)
reference <- read.table(file="~/Downloads/aging_merged_consensus_whitelisted_annotated.txt",header = TRUE,sep = "\t",stringsAsFactors = F)

rna.ref <- read.table(file="~/Desktop/TCell_DiffData/merged_consensus_rna_results.txt",sep = "\t")
colnames(rna.ref) <- c("GeneName", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")
rna.ref <- rna.ref[-1,]

diff.mx <- merge(diff.mx, rna.ref[, c("GeneName", "logFC")], by="GeneName")
colnames(diff.mx)[which(colnames(diff.mx) == "logFC")] <- "logFC.rna"
# atac.info <- read.table(file="~/Downloads/atac_do.tcelldiff_enrichment.txt",header = TRUE,sep = "\t",stringsAsFactors = F)
# atac.ref <- read.table(file="~/Downloads/aging_merged_consensus_atac_results.txt",sep = "\t")
# colnames(atac.ref) <- c("chr", "start", "end", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")
# atac.ref <- atac.ref[-1,]
# atac.ref$logFC <- as.numeric(as.character((atac.ref$logFC)))
# atac.ref <- aggregate(logFC~chr + start + end,data=atac.ref,FUN=mean)
# atac.ref <- merge(reference, atac.ref, by=c("chr", "start", "end"))
# atac.ref <- atac.ref[atac.ref$GeneName %in% diff.mx$GeneName, ]
# atac.ref <- do.call(rbind,lapply(split(atac.ref,atac.ref$GeneName),function(chunk) chunk[which.max(abs(chunk$logFC)),]))
# diff.mx <- merge(diff.mx[diff.mx$GeneName %in% atac.ref$GeneName, ], atac.ref, by="GeneName",suffixes = c(".rna",".atac"))

row.names(diff.mx) <- diff.mx$GeneName
diff.mx <- diff.mx[order(diff.mx$logFC.atac), ]
diff.mx[is.na(diff.mx)] <- 0

diff.mx <- diff.mx[order(diff.mx$logFC.atac), ]
newCols <- colorRampPalette(colorspace::diverge_hsv(length(unique(diff.mx$logFC.atac))))
mycolors <- newCols(length(unique(diff.mx$logFC.atac)))
names(mycolors) <- unique(diff.mx$logFC.atac)
mycolors <- list(logFC.atac = mycolors)

diff.mx <- diff.mx[order(diff.mx$logFC.rna), ]
newCols2 <- colorRampPalette(colorspace::diverge_hsv(length(unique(diff.mx$logFC.rna))))
mycolors2 <- newCols2(length(unique(diff.mx$logFC.rna)))
names(mycolors2) <- unique(diff.mx$logFC.rna)
mycolors2 <- list(logFC.rna = mycolors2)

diff.mx <- diff.mx[order(diff.mx$logFC.rna), ]
pdf("~/Desktop/ExpressionMatrix.rna.pdf")
pheatmap(data.frame(diff.mx[, c(3:8)], row.names = row.names(diff.mx)),
         show_rownames = T, scale = "none",
         color = colorRampPalette(brewer.pal(6, "BuGn"))(100), 
         annotation_row = x <- diff.mx[ ,"logFC.rna", drop=FALSE], 
         annotation_colors = mycolors2, 
         cluster_rows=FALSE, cluster_cols=FALSE, 
         annotation_legend = F)
dev.off()

diff.mx <- diff.mx[order(diff.mx$logFC.atac), ]
pdf("~/Desktop/ExpressionMatrix.atac.pdf")
pheatmap(data.frame(diff.mx[, c(3:8)], row.names = row.names(diff.mx)),show_rownames = T,scale = "row",color = colorRampPalette(c("white", "lightcyan", "lightcyan1", "lightcyan2", "lightcyan3", "lightcyan4"))(100), annotation_row = diff.mx[, c("logFC.atac", "logFC.rna")], annotation_colors = c(mycolors, mycolors2), cluster_rows=FALSE, cluster_cols=FALSE, annotation_legend = F)
dev.off()


##################
#
# ATAC-seq RNA-seq agreement
#
##################

rm(list = ls())

library(ggplot2)
library(reshape2)
library(MKmisc)
library(pheatmap)
library(ggExtra)
library(gridExtra)
library(gtable)
library(GenomicRanges)
library(dplyr)
library(broom)

sse <- function(x) sum((x-mean(x))^2)
andred <- function(x) {Reduce("&",x)}
orred <- function(x) {Reduce("|",x)}
minmaxnorm <- function(x) {(x - min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE))}
minmaxnormc <- function(x) {minmaxnorm(x)-mean(minmaxnorm(x))}
znorm <- function(x) {(x - mean(x)) / sd(x)}

colors_hmm6=c("gold", "mediumaquamarine", "bisque", "slategray2", "firebrick1", "green4")
colors_redblu10 = c("#3360a4","#1a8af9","#06b8ff","#b3defb","#ffffff","#fdc82f","#f67028","#e9302a","firebrick3","firebrick4")

# setwd("/Volumes/emext/Dropbox/Jax/PBMC/DifferentialAnalysis/manuscripts/aging-01/scripts/figures.regen.stm/")
# setwd("~/Dropbox/Jax/PBMC/DifferentialAnalysis/manuscripts/aging-01/scripts/figures.regen.stm/")
# setwd("H:/Dropbox/Jax/PBMC/DifferentialAnalysis/manuscripts/aging-01/scripts/figures.regen.stm/")

# function to relabel samples on-the-fly wrt randomIDs
sampinfo <- read.table("~/Downloads/sample.info_stm.txt",header = TRUE,row.names = 1,sep = "\t")
random.codes <- read.table("~/Downloads/RandomRecodeAgingStudy.txt", header = T, sep = "\t", quote = "")
random.codes <- data.frame(merge(data.frame(random.codes,sampid=sub("VHY402","VHY402D0",paste(random.codes$Cohort,random.codes$PID,sep = "")),row.names = "sampid"),sampinfo,by="row.names"),row.names = 1)
rand.recode <- function(x) with(random.codes[x,], paste(group,RandomID,sep=""))

# Previoulsy analyzed ATACseq data
da.pbmc.out <- load("~/Downloads/da.pbmc.RData")
full.results <- read.table("~/Downloads/pbmc_whitelisted_filtered_glm.results.txt", sep = "\t", quote = "", header = T,row.names = "peakco")
full.results$peakco <- rownames(full.results)
adj.atac <- data.frame(merge(full.results[,c("chr","start","end")],da.pbmc$adj,by = "row.names"),row.names = 1)
# atacsamps <- colnames(adj.atac[,-c(1:3)])
rm(list = da.pbmc.out)

# Previoulsy analyzed RNAseq data (samps common in RNA and filtered ATAC)
da.rna.pbmc.out <- load("~/Downloads/da.rna.pbmc.RData")
rna.glmtop <- read.table("~/Downloads/pbmc_rna.atacpaired_glm.results.txt", sep = "\t", quote = "", header = T)
rownames(rna.glmtop) <- rna.glmtop$GeneName
# rnasamps <- colnames(da.pbmc.atacpairs$adj)
adj.rna <- da.pbmc.atacpairs$adj
rpkm <- da.pbmc$raw
rm(list = da.rna.pbmc.out)

samps <- intersect(rnasamps,atacsamps) # Samples that are those present in both datasets
expressed.genes <- rownames(rpkm)

# Splinters TSS peaks out to combine with RNAseq data
# atac.tss <- read.table("datasources.201702/raw_data_tss.peaks_minOver0_hard_stm.txt",quote = "",header = T, sep = "\t")
flankd = 1000
tss.bed <- GRanges(read.table(paste("../../data/Promoters/tssflank",flankd,"_prom_complete.bed",sep = ""),header = F,sep = "\t",quote = "",col.names = c("chr","start","end","GeneName","score","strand")))
tss.bed.expressed <- subset(tss.bed,GeneName %in% expressed.genes)
tss_peaks <- findOverlaps(GRanges(full.results[,c("chr","start","end")]),tss.bed.expressed)
# The following ensures gene names are derived from TSS definitions instead of nearest-TSS annotations
full.results.tss_red <- data.frame(GeneName=as.character(tss.bed.expressed[tss_peaks@to]$GeneName),full.results[tss_peaks@from,-which(colnames(full.results)=="GeneName")],row.names = NULL)
# full.results.tss_red <- merge(full.results.tss_red,data.frame(meanExpression=rowMeans(rpkm)),by.x = "GeneName", by.y = "row.names") # to keep peak associated to highest-expressed gene
full.results.tss_red <- unique(full.results.tss_red[order(full.results.tss_red$GeneName,abs(full.results.tss_red$DistancetoTSS)),]) # to keep one peak per gene (multiple genes per peak allowed), kept peak closest to nearest TSS
full.results.tss <- full.results.tss_red[!duplicated(full.results.tss_red$GeneName),]
rownames(full.results.tss) <- full.results.tss$GeneName
rm("full.results.tss_red")

adj.atac.tss <- merge(full.results.tss[,c("chr","start","end","GeneName")], adj.atac, by = c("chr","start","end"))
# adj.atac.tss.nodup <- adj.atac.tss[order(adj.atac.tss$GeneName,-rowSums(2^adj.atac.tss[,atacsamps])),c("GeneName",atacsamps)]
# adj.atac.tss.nodup <- adj.atac.tss.nodup[!duplicated(adj.atac.tss.nodup$GeneName),]

atacgroup <- relevel(factor(sampinfo[atacsamps,]$group),ref = "HY")
atacsex <- relevel(factor(sampinfo[atacsamps,]$sex),ref = "F")
rnagroup <- relevel(factor(sampinfo[rnasamps,]$group),ref = "HY")
rnasex <- relevel(factor(sampinfo[rnasamps,]$sex),ref = "F")
group <- relevel(factor(sampinfo[samps,]$group),ref = "HY")
sex <- relevel(factor(sampinfo[samps,]$sex),ref = "F")

##### Permutation tests to derive empirical Pvalues for each gene ATAC/RNA fold change
nperm = 1000
#### RNAseq FC
### Code used to generate perms:
# rna.permlist <- replicate(nperm,adj.rna,simplify = F)
# rna.permfc <- mclapply(rna.permlist, function(Y) apply(Y,1,function(x) diff(aggregate(x, list(group=sample(group)), mean)[,2])),mc.cores = 6)
# save(list = c("rna.permlist","rna.permfc"),file = "datasources.201702/RNAseq_permutation.test_stm.RData")
rna.permout <- load("~/Downloads/RNAseq_permutation.test_stm.RData")
# Observed FC
adj.rna <- adj.rna[rownames(adj.rna) %in% test.table.tcell$GeneName, ]
rna.obsfc <- data.frame(logFC.rna.obs=apply(adj.rna[,samps],1,function(x) diff(aggregate(x, list(group=group), mean)[,2])))
# Permuted FC
rna.permfc.mx <- do.call(cbind, rna.permfc)
rna.permfc.mx <- data.frame(rna.permfc.mx)[rownames(rna.obsfc),]
rna.obsfc$logFC.rna.pcounts[rna.obsfc$logFC.rna.obs >= 0] <- rowSums(rna.obsfc$logFC.rna.obs[rna.obsfc$logFC.rna.obs >= 0] >= rna.permfc.mx[rna.obsfc$logFC.rna.obs >= 0,])
rna.obsfc$logFC.rna.pcounts[rna.obsfc$logFC.rna.obs < 0] <- rowSums(rna.obsfc$logFC.rna.obs[rna.obsfc$logFC.rna.obs < 0] < rna.permfc.mx[rna.obsfc$logFC.rna.obs < 0,])
rna.obsfc$logFC.rna.pvalue <- 1-(rna.obsfc$logFC.rna.pcounts/(1+ncol(rna.permfc.mx)))
rna.obsfc$logFC.rna.fdr <- p.adjust(rna.obsfc$logFC.rna.pvalue, method = "fdr")
rm(list = c(rna.permout,"rna.permfc.mx"))
#### ATACseq TSS peaks FC
### Code used to generate perms:
# atac.permlist <- replicate(nperm,adj.atac.tss[,samps],simplify = F)
# atac.permfc <- mclapply(atac.permlist, function(Y) apply(Y,1,function(x) diff(aggregate(x, list(group=sample(group)), mean)[,2])),mc.cores = 6)
###atac.permfc <- lapply(atac.permlist, function(Y) apply(Y,1,function(x) diff(aggregate(x, list(group=sample(group)), mean)[,2])))
# save(list = c("atac.permlist","atac.permfc"),file = "datasources.201702/ATACseq_permutation.test_stm.RData")
atac.permout <- load("~/Downloads/ATACseq_permutation.test_stm.RData")
# Observed FC
atac.obsfc <- data.frame(adj.atac.tss[,c("chr","start","end")],logFC.atac.obs=apply(adj.atac.tss[,samps],1,function(x) diff(aggregate(x, list(group=group), mean)[,2])))
# Permuted FC
atac.permfc.mx <- do.call(cbind, atac.permfc)
atac.permfc.mx <- data.frame(atac.permfc.mx)[rownames(atac.obsfc),]
atac.obsfc$logFC.atac.pcounts[atac.obsfc$logFC.atac.obs >= 0] <- rowSums(atac.obsfc$logFC.atac.obs[atac.obsfc$logFC.atac.obs >= 0] >= atac.permfc.mx[atac.obsfc$logFC.atac.obs >= 0,])
atac.obsfc$logFC.atac.pcounts[atac.obsfc$logFC.atac.obs < 0] <- rowSums(atac.obsfc$logFC.atac.obs[atac.obsfc$logFC.atac.obs < 0] < atac.permfc.mx[atac.obsfc$logFC.atac.obs < 0,])
atac.obsfc$logFC.atac.pvalue <- 1-(atac.obsfc$logFC.atac.pcounts/(1+ncol(atac.permfc.mx)))
atac.obsfc$logFC.atac.fdr <- p.adjust(atac.obsfc$logFC.atac.pvalue, method = "fdr")
atac.obsfc$peakco <- paste(atac.obsfc$chr,atac.obsfc$start,atac.obsfc$end,sep = "_")
atac.obsfc <- atac.obsfc[order(atac.obsfc$peakco,atac.obsfc$logFC.atac.pvalue),]
atac.obsfc <- atac.obsfc[!duplicated(atac.obsfc$peakco),]
rm(list = c(atac.permout,"atac.permfc.mx"))
#####

# Combines ATAC and RNAseq data/results, assigns Pvalues based on permutation tests
atacrna.glm <- merge(full.results.tss,rna.glmtop, by = "GeneName",suffixes = c(".atac",".rna"))
atacrna.glm <- merge(merge(atacrna.glm, atac.obsfc, by = c("chr","start","end")),
                     rna.obsfc, by.x = "GeneName", by.y = "row.names")

diff.mx$signed_exp <- -(diff.mx$TN + diff.mx$TSCM + diff.mx$TCM +diff.mx$TTM) + (diff.mx$TTE+diff.mx$TEM)
diff.mx$color[diff.mx$signed_exp < 0] <- colorRampPalette(c('blue', 'lightblue'))(length(which(diff.mx$signed_exp < 0)))[rank(diff.mx$color[diff.mx$signed_exp < 0])]
diff.mx$color[diff.mx$signed_exp > 0] <- colorRampPalette(c('lightpink', 'red'))(length(which(diff.mx$signed_exp > 0)))[rank(diff.mx$color[diff.mx$signed_exp > 0])]
# diff.mx <- diff.mx[diff.mx$GeneName %in% atacrna.glm$GeneName,]

atacrna.glm.color <- merge(atacrna.glm, diff.mx, by="GeneName")

pdf("~/Fig3A_pbmc.arc.fcfcplot.pdf",paper = "USr")
fcfcplot <- atacrna.glm # atacrna.glm|atacrna.glm.nodup|atacrna.comp.glm|atacrna.comp.glm.nodup :: if use "comp", variables below have to be modified
fcfc.corr <- with(fcfcplot,cor.test(logFC.atac.obs,logFC.rna.obs, method = "pearson"))
atac <- fcfcplot$logFC.atac.obs
rna <- fcfcplot$logFC.rna.obs
fcfc.pcorr_txt <- ifelse(fcfc.corr$p.value<0.001,"P < 0.001",ifelse(fcfc.corr$p.value<0.01,"P < 0.01",ifelse(fcfc.corr$p.value<0.05,"P < 0.05","P > 0.05")))
alpha = 0.01;
xlims <- c(min(fcfcplot$logFC.atac.obs),
           min(quantile(fcfcplot$logFC.atac.obs,0.05),max(fcfcplot$logFC.atac.obs[apply(cbind(fcfcplot$logFC.atac.pvalue<alpha, fcfcplot$logFC.atac.obs<0),1,andred)])),
           max(fcfcplot$logFC.atac.obs[apply(cbind(fcfcplot$logFC.atac.pvalue<alpha, fcfcplot$logFC.atac.obs<0),1,andred)]),
           min(fcfcplot$logFC.atac.obs[apply(cbind(fcfcplot$logFC.atac.pvalue<alpha, fcfcplot$logFC.atac.obs>0),1,andred)]),
           max(quantile(fcfcplot$logFC.atac.obs,0.95),min(fcfcplot$logFC.atac.obs[apply(cbind(fcfcplot$logFC.atac.pvalue<alpha, fcfcplot$logFC.atac.obs>0),1,andred)])),
           max(fcfcplot$logFC.atac.obs))
ylims <- c(min(fcfcplot$logFC.rna.obs),
           min(quantile(fcfcplot$logFC.rna.obs,0.05),max(fcfcplot$logFC.rna.obs[apply(cbind(fcfcplot$logFC.rna.pvalue<alpha, fcfcplot$logFC.rna.obs<0),1,andred)])),
           max(fcfcplot$logFC.rna.obs[apply(cbind(fcfcplot$logFC.rna.pvalue<alpha, fcfcplot$logFC.rna.obs<0),1,andred)]),
           min(fcfcplot$logFC.rna.obs[apply(cbind(fcfcplot$logFC.rna.pvalue<alpha, fcfcplot$logFC.rna.obs>0),1,andred)]),
           max(quantile(fcfcplot$logFC.rna.obs,0.95),min(fcfcplot$logFC.rna.obs[apply(cbind(fcfcplot$logFC.rna.pvalue<alpha, fcfcplot$logFC.rna.obs>0),1,andred)])),
           max(fcfcplot$logFC.rna.obs))
# genesel <- c("LCN2","PDGFRB","KIR3DL1","PDGFD","GNLY","GZMB","KLRF1","CD248","NOG","CCR7","IL7R","LEF1","BACH2","CD8A","LRRN3","GSTT1","NRCAM","SLC16A10","NBEA","GSTM1","MTUS1","B3GAT1","IGFBP3","NCR1","FGFBP2","GZMH","PRSS23","NT5E","S100B") # highlighted in ms.
genesel <- c("LCN2","PDGFRB","KIR3DL1","PDGFD","GNLY","GZMB","KLRF1","CD248","NOG","CCR7","IL7R","LEF1","BACH2","CD8A","LRRN3","GSTT1","NRCAM","SLC16A10","NBEA","GSTM1","MTUS1","B3GAT1","FGFBP2","GZMH","PRSS23","S100B","PRF1","NT5E") # highlighted in ms.
fcfcplot.sel <- subset(fcfcplot,with(fcfcplot,GeneName %in% genesel & apply(cbind(logFC.atac.obs<xlims[3]&logFC.rna.obs<ylims[3],logFC.atac.obs>xlims[4]&logFC.rna.obs>ylims[4]),1,any)))
fcfcplot.sel <- fcfcplot.sel[order(fcfcplot.sel$GeneName,fcfcplot.sel$logFC.atac.obs),]
fcfcplot.sel <- fcfcplot.sel[!duplicated(fcfcplot.sel$GeneName),]
p <- ggplot(data=fcfcplot[order(fcfcplot$logFC.atac.obs),], aes(logFC.atac.obs, logFC.rna.obs)) +
  #annotate("rect", xmin = -max(abs(xlims*1.1)), xmax = xlims[3], ymin = -max(abs(ylims*1.1)), ymax = ylims[3], fill = "darkcyan", alpha = 0.1) + # alt version with differently-colored quadrants. Use "darkorange1" to match both quadrants to the original version
  #annotate("rect", xmin = xlims[4], xmax = max(abs(xlims)), ymin = ylims[4], ymax = max(abs(ylims)), fill = "darkorange1", alpha = 0.1) +
  geom_point(aes(size=factor(apply(cbind(logFC.atac.obs<xlims[3]&logFC.rna.obs<ylims[3],logFC.atac.obs>xlims[4]&logFC.rna.obs>ylims[4]),1,orred)),
                 alpha=factor(apply(cbind(logFC.atac.obs<xlims[3]&logFC.rna.obs<ylims[3],logFC.atac.obs>xlims[4]&logFC.rna.obs>ylims[4]),1,orred)),
                 color=factor(apply(cbind(logFC.atac.obs<xlims[3]&logFC.rna.obs<ylims[3],logFC.atac.obs>xlims[4]&logFC.rna.obs>ylims[4]),1,orred)))) +
  # geom_text(data = fcfcplot[fcfcplot$logFC.atac.obs*fcfcplot$logFC.rna.obs > 0.2,], aes(logFC.atac.obs,logFC.rna.obs,label=GeneName), size=2.5) + # uses a threshhold to apply labels to extreme points
  geom_text_repel(data = fcfcplot,aes(logFC.atac.obs,logFC.rna.obs,label=GeneName), size=3.5, color = atacrna.glm.color$color) + # preselect labels used in the paper (negative)
  # geom_text(data = subset(fcfcplot.sel,logFC.atac.obs>0),aes(logFC.atac.obs,logFC.rna.obs,label=GeneName), size=5, vjust = -0.3, hjust = 0, color = "green4") + # preselect labels used in the paper (positive)
  geom_point(data = fcfcplot.sel,aes(logFC.atac.obs,logFC.rna.obs), color = "green3", size = 0.8) +
  scale_color_manual(values = c("FALSE"="dimgray","TRUE"="navy"), guide=F) +
  scale_size_manual(values = c("FALSE"=0.5,"TRUE"=0.7), guide=F) +
  scale_alpha_manual(values = c("FALSE"=0.25,"TRUE"=0.75), guide=F) +
  scale_x_continuous(limits = c(-max(abs(xlims*1.1)),max(abs(xlims)))) +
  scale_y_continuous(limits = c(-max(abs(ylims*1.1)),max(abs(ylims)))) +
  xlab("ATAC-seq logFC (Promoter accessibility)") +
  ylab("RNA-seq logFC (Expression)") +
  ggtitle(paste("logFC RNAseq vs ATACseq (r = ",sprintf("%0.2f",fcfc.corr$estimate),", ",fcfc.pcorr_txt,")", sep = "")) +
  theme_bw(base_family = "Helvetica",base_size = 10) +
  geom_smooth(method=lm, size=., se = FALSE, color="black") +
  theme(panel.background = element_rect(fill = "white", color = "black"),
        # panel.grid.major.x = element_line(color = "ivory3", size = 0.2, linetype = 2),
        # panel.grid.major.y = element_line(color = "ivory3", size = 0.2, linetype = 2),
        panel.grid.major = element_blank(),
        axis.line = element_line(color = "black", size = 0.25),
        axis.text.x = element_text(size=10))
  # geom_hline(yintercept = ylims[3:4], color = "palegreen3", linetype = "22", size=0.5, alpha=0.5) +
  # geom_vline(xintercept = xlims[3:4], color = "palegreen3", linetype = "22", size=0.5, alpha=0.5)
# p <- ggMarginal(p,type = "histogram", size = 15, bins = 150, fill="palegreen3", color="seagreen3")
print(p)
dev.off()
