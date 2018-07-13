## Code used to whitelist (vet) ATACseq peaks. 
## Here, peaks are vetted by qvalue (a measure of how much the peak is different 
## from its local background) and 
## to exclude peaks overlapping ENCODE-blacklisted regions. 
## These regions often suffer from poor mappability and reads there cannot be trusted.

rm(list = ls())

# library(GenomicRanges)
library(dplyr)

blacklst <- GRanges(read.table("~/Desktop/SSP/Data/support/hg19.blacklists_merged.bed",sep = "\t",header = F,quote = "",col.names = c("chr","start","end","region.name"))) # need to edit. check data/support directory
qvals <- read.table("~/Desktop/SSP/Data/jax/atac/atac_summerXwinter_pbmc_sorted_merged_consensus_peakstats_qvalues.txt",header = T,quote = "",sep = "\t") # need to edit also, if using. check data/jax/atac
qval.thresh = 1e-2

listFunc <- function(a, b) list(a, b) # basic function which creates list from two objects
qvals$peakco <- apply(qvals[, c('start', 'end')], 1, function(x) listFunc(x[1], x[2])) # create peakco column using listFunc
qvals.filtr <- GRanges(qvals[apply((qvals %>% select(-chr,-start,-end,-peakco)) >= -log10(qval.thresh),1,function(x) sum(x)>1),] %>% select(chr,start,end))

atac.raw <- read.table("~/Desktop/SSP/Data/jax/atac/atac_summerXwinter_pbmc_sorted_merged_peaks_rawcounts.txt",sep = "\t",header = T,quote = "") # check data/jax/atac... you get the idea, so I won't further point this out

vet_in <- findOverlaps(GRanges(atac.raw),qvals.filtr)
peaks <- GRanges(atac.raw[vet_in@from,c("chr","start","end")])
vet_out <- findOverlaps(peaks,blacklst)
peaks_vetted <- merge(setNames(as.data.frame(sort(peaks[-unique(vet_out@from),]))[,1:3],c("chr","start","end")),atac.raw,by = c("chr","start","end"))

## In addition to qvalues, for PBMC we applied other filters, which I'm adding in the following code. You may want to implement those too, in any way you see fit:
# This filters out peaks in sex chromosomes:
chromfilter.pbmc <- data.frame(chromfilter=!is.element(atac.basepeaks.pbmc$chr,c("chrX","chrY")),row.names = rownames(atac.basepeaks.pbmc))
peaks_vetted <- peaks_vetted[which(!peaks_vetted$chr %in% c("chrX", "chrY")),]
# This computes the maximum number of reads in each peak across all samples, and then filters out any peak with a maximum of 20 or fewer reads
maxcts.pbmc <- data.frame(maxcts=apply(peaks_vetted[, colnames(peaks_vetted) %in% rownames(sampinfo)],1,max))
sizefilter.pbmc <- data.frame(maxcts=maxcts.pbmc>20,row.names = rownames(maxcts.pbmc))
# This computes the maximum reads per million per peak, and filters out peaks with 500 or more reads (probably a mapping artifact)
maxcpm.pbmc <- data.frame(maxcpm=apply(cpm(atac.basepeaks.pbmc[,sampid.pbmc],lib.size = libsize.atac.pbmc),1,max))
cpmfilter.pbmc <- data.frame(cpmfilter=maxcpm.pbmc<500,row.names = rownames(maxcpm.pbmc))
# This collects all of the previous filters and actually filters the peaks
peakfilter.pbmc <- apply(data.frame(merge(data.frame(merge(merge(callfilter.pbmc,qfilter.pbmc,by="row.names"),merge(sizefilter.pbmc,chromfilter.pbmc,by="row.names"),by = "Row.names"),row.names = 1),cpmfilter.pbmc,by = "row.names"),row.names = 1),1,andred)
##

write.table(peaks_vetted,"~/Desktop/SSP/monot_merged_consensus_raw_atac_whitelisted.txt",sep = "\t",col.names = T,row.names = F,quote = F)
write.table(data.frame(sort(GRanges(peaks_vetted[,c("chr","start","end")])))[,c("seqnames","start","end")],"~/Desktop/SSP/monot_merged_consensus_whitelisted.bed",sep = "\t",col.names = F,row.names = F,quote = F)
# the above saves two files, one with the vetted raw counts and another, referred to simply as "consensus peask" which contains only the peak locations

## Code used to select "open genes" from ATACseq data
## Simply matches consensus peaks to a list of known TSS and keeps only those peaks whose center is within 1kb from a TSS

library(biomaRt)
# library(GenomicRanges)

tss_flank = 1000
tss <- keepStandardChromosomes(GRanges(read.table("~/Desktop/SSP/Data/support/refSeq_TSS_clean.bed",header = F,sep = "\t",quote = "",col.names = c("chr","start","end","strand","GeneName"))), pruning.mode = "coarse") # Loads TSS
peaks <- dropSeqlevels(GRanges(read.table("~/Desktop/SSP/RNA-seq_6_23_17/monot_merged_consensus_raw_atac_whitelisted.txt",sep = "\t",header = T,quote = "")[,c("chr","start","end")]),'chrM')
peaks2tss <- distanceToNearest(peaks,tss)

# all peaks without flank
large_tss_flank = 10^10 # impossibly high to still change data type, but still keep the peaks
tsspeaks <- peaks[subset(peaks2tss,abs(distance)<=large_tss_flank)@from,]
tsspeaks$GeneName <- tss[subset(peaks2tss,abs(distance)<=large_tss_flank)@to,]$GeneName
tsspeaks.df <- data.frame(sort(tsspeaks))

# with flank
tsspeaks <- peaks[subset(peaks2tss,abs(distance)<=tss_flank)@from,]
tsspeaks$GeneName <- tss[subset(peaks2tss,abs(distance)<=tss_flank)@to,]$GeneName
tsspeaks.df <- data.frame(sort(tsspeaks))

open.genes <- data.frame(GeneName=sort(unique(tsspeaks.df$GeneName)))
write.table(open.genes,"~/Desktop/SSP/geneset.open.monot.txt",quote = F,sep = "\t",col.names = T,row.names = F)
Hide
## Code used to build expression matrix

library(edgeR)
library(biomaRt)
library(dplyr)

expressionSet.1 <- read.csv("./data_monot/RNAseq/PBMC+B+T_counts.csv",stringsAsFactors = F)
expressionSet.2 <- read.csv("./data_monot/RNAseq/CD14_samples_counts.csv",stringsAsFactors = F)
expressionSet <- merge(expressionSet.1,expressionSet.2,by = "EnsemblID")
# Aging PBMC data is also split in two files, named counts_protein_coding_genes_library_size_normalized.csv and expected_counts.txt.
# Keep this in mind when processing: 
# (1) The samples corresponding to the ATAC samples are distributed
#     between the two files (not all ATAC samples have RNA samples); 
# (2) don't use samples whose name contains the name Nugen; 
# (3) the first file is normalized so that each column adds up to 1e6, 
#     and only include protein-coding genes. You should match the rows of the two files 
#     AND THEN normalize the second file also so that columns add up to 1e6 
#     (you can use the cpm function in the edgeR package, for example). 
#     Otherwise, the two sets of samples cannot be compared.

expressionSet.1 <- read.csv(file = "~/Desktop/SSP/Data/jax/rna/counts_protein_coding_genes_library_size_normalized.csv", stringsAsFactors = F)
colnames(expressionSet.1)[1] <- "EnsemblID"

expressionSet.2 <- read.csv(file = "~/Desktop/SSP/Data/jax/rna/expected_counts.txt", stringsAsFactors = F, sep = "\t")
colnames(expressionSet.2)[1] <- "EnsemblID"

for(i in 2:31){
  expressionSet.2[,i] <- cpm(expressionSet.2[,i])
  print(sum(expressionSet.2[,i]))
}

expressionSet <- merge(expressionSet.1, expressionSet.2, by = "EnsemblID")

# This converts Ensembl names to official gene symbols. Even if the files above include gene symbols, I would replace them by their annotatios derived from fresh conversion from Ensembl ids.
expressionSet.genes <- getBM(attributes = c("ensembl_gene_id","external_gene_name","gene_biotype","chromosome_name"),
                             filters = "ensembl_gene_id",
                             values = expressionSet$EnsemblID,
                             mart=useMart(biomart="ENSEMBL_MART_ENSEMBL",
                                          dataset="hsapiens_gene_ensembl",
                                          host="grch37.ensembl.org"))

atac.monot.peaks <- merge(expressionSet, open.genes, by.x="Associated.Gene.Name", by.y = "GeneName")

monot.expressionSet <- merge(subset(expressionSet.genes,gene_biotype=="protein_coding" & (!is.na(as.numeric(expressionSet.genes$chromosome_name)) | expressionSet.genes$chromosome_name %in% c("X","Y"))),
                             expressionSet[,c(intersect(colnames(expressionSet),colnames(atac.monot.peaks)))],  # you don't have, or need, a "peakcalls" file. This intersect simply selects the RNAseq columns for which there is ATACseq data. You can use any of the data files headers to do this. 
                             by.x="ensembl_gene_id",
                             by.y="EnsemblID") 
monot.expressionSet[, c("gene_biotype","chromosome_name", "Associated.Gene.Name")] <- NULL
names(monot.expressionSet) <- sub("ensembl_gene_id","EnsemblID",sub("external_gene_name","GeneName",colnames(monot.expressionSet)))

monot.expressionSet <- cbind(monot.expressionSet %>% select(EnsemblID,GeneName),data.frame(cpm(monot.expressionSet %>% select(-EnsemblID,-GeneName, -Description))))
write.table(monot.expressionSet,"~/Desktop/SSP/monot_rnaseq.txt",sep = "\t",col.names = T,quote = F)
Hide
## Code used to select "expressed genes" from RNAseq data. This has already been carried out in this dataset

library(dplyr)

expression.data <- read.table("~/Desktop/SSP/monot_rnaseq.txt",sep = "\t",header = T,quote = "",stringsAsFactors = F)
expressed.genes <- data.frame(GeneName=expression.data$GeneName[expression.data %>% select(-GeneName,-EnsemblID) %>% cpm(.,log=F) %>% apply(.,1,function(x) sum(x>3)>1)],stringsAsFactors = F)
write.table(expressed.genes,"~/Desktop/SSP/geneset.expressed.monot.txt",quote = F,sep = "\t",col.names = T,row.names = F)
Hide
## Code used for differential analyses on ATACseq data
## This uses a different method (voom) from what we used in the paper (edgeR)

library(pheatmap)
library(ggplot2)
library(sva)

# For this I'm including code from our original PBMC analysis and the newer differential analyses, used for t cells. Let me know if you need help figuring out what goes where.

sampinfo <- read.table("~/Desktop/SSP/data/jax/sample.info.comprehensive.txt",header = TRUE,sep = "\t",stringsAsFactors = F) # check data/jax
rownames(sampinfo) <- with(sampinfo,paste(sampid,type,sep = "_"))
sampinfo[sampinfo$age>40 & sampinfo$group=="HY",]$group <- "HM"

# The following factors should be filtered to include only HY-HO samples for which there are data, including only samples used in the paper (C,NH)
samp.filter <- apply(cbind(sampinfo$group %in% c("HY","HO"),!sampinfo$race %in% c("A","B"),!sampinfo$ethnicity=="H",!is.na(sampinfo$depth),!grepl("VHY402D[1,2,6,7].*",sampinfo$sampid),!sub("_.*","",sampinfo$sampid) %in% c("HO207","ND1","ND2")),1, all)
sampid <- rownames(sampinfo[samp.filter,])
samps <- unique(sampinfo[samp.filter,]$sampid)
age.group <- relevel(factor(setNames(sampinfo[samp.filter,]$group,sampid)),ref = "HO") # this relevel/ref statement is actually what makes HY-specific peaks to have negative logs. By default, HY peaks would have positive logs.
sex <- relevel(factor(setNames(sampinfo[samp.filter,]$sex,sampid)),ref = "F")
age <- setNames(sampinfo[samp.filter,]$age,sampid)
season <- setNames(sampinfo[samp.filter,]$season,sampid)
batch <- setNames(sampinfo[samp.filter,]$batch,sampid)
cmv.status <- setNames(sampinfo[samp.filter,]$CMV.status,sampid)
cmv.titer <- setNames(sampinfo[samp.filter,]$CMV.titer,sampid)
cell.type <- setNames(sampinfo[samp.filter,]$type,sampid)
depth <- setNames(sampinfo[samp.filter,]$depth,sampid)


atac.tcells.raw <- read.table("~/Desktop/SSP/monot_merged_consensus_raw_atac_whitelisted.txt",sep = "\t",header = T,quote = "") # this corresponds to the vetted peak file you generated above
atac.samples <- rownames(group.names)
atac.tcells.raw$nhits <- rowSums(atac.tcells.raw[,which(!colnames(atac.tcells.raw) %in% c("chr","start","end"))])
atac.tcells.raw$peakco <- paste(atac.tcells.raw$chr,":",atac.tcells.raw$start,"-",atac.tcells.raw$end,sep = "")

atac.tcells.raw.mx <- as.matrix(data.frame(atac.tcells.raw[,which(!colnames(atac.tcells.raw) %in% c("chr","start","end","nhits"))],row.names = "peakco"))
colnames(atac.tcells.raw.mx)[58] <- "VHY402D0_PBMC"
atac.tcells.raw.mx <- atac.tcells.raw.mx[, which(colnames(atac.tcells.raw.mx) %in% rownames(group.names))]
atac.tcells.libsize <- colSums(atac.tcells.raw.mx)

group.names.2 <- sampinfo[names(atac.tcells.libsize),] %>% select(sex,group,season) # Use this; accounts for discrepancies between samples seen in atac peaks
                                                                                    # and data that was filtered out of samples due to race, age, health, etc.
group.names.2$group <- relevel(factor(group.names.2$group),ref="HY")
group.names.2$sex <- relevel(factor(group.names.2$sex),ref="F")

group.names <- data.frame(sex, age.group, season) # All filtered values

atac.tcells.design <- with(group.names.2,model.matrix(~sex+season+group)) # in our paper, our design was sex+season+group
atac.tcells.design.null <- with(group.names.2,model.matrix(~sex+season,data = data.frame(atac.tcells.libsize))) # the null design is sex+season (i.e. everything but what we are interested in)

# This is for SVA, surrogate variable analysis. You can read a bit about it and we can chat later on its purpose and benefits
n.sv <- num.sv(cpm(atac.tcells.raw.mx,lib.size = atac.tcells.libsize)*1e6,atac.tcells.design)
sv <- svaseq(cpm(atac.tcells.raw.mx,lib.size = atac.tcells.libsize,log = F)*1e6,atac.tcells.design,atac.tcells.design.null,n.sv = n.sv)$sv
atac.tcells.design <- cbind(atac.tcells.design,sv)
atac.tcells.design.null <- cbind(atac.tcells.design.null,sv)

## VOOM-based DA
atac.tcells.dge <- voom(atac.tcells.raw.mx,design = atac.tcells.design,plot = F)
atac.tcells.norm <- atac.tcells.dge$E
atac.tcells.lm <- lmFit(atac.tcells.dge,design = atac.tcells.design)
atac.tcells.lmfit <- eBayes(contrasts.fit(atac.tcells.lm,coefficients = "groupHO"),trend = T,robust = T)
atac.tcells.lmtop <- topTable(atac.tcells.lmfit,number = nrow(atac.tcells.raw.mx),sort.by = "none",adjust.method = "BH")
atac.tcells.dge.null <- voom(atac.tcells.raw.mx,design = atac.tcells.design.null,plot = F)
atac.tcells.lm.null <- lmFit(atac.tcells.dge.null,design = atac.tcells.design.null)
sum(atac.tcells.lmtop$adj.P.Val<.05) # 9240
hist(atac.tcells.lmtop$P.Value,100)

atac.tcells.norm.data <- cbind(atac.tcells.raw[,c("chr","start","end")],atac.tcells.norm)
atac.tcells.adj.data <-  cbind(atac.tcells.raw[,c("chr","start","end")],atac.tcells.norm - fitted.MArrayLM(atac.tcells.lm.null))

write.table(atac.tcells.norm.data,"~/Desktop/SSP/tcell-peaks6_23_17/tcells_merged_consensus_voom_atac_whitelisted.txt",quote = F,sep = "\t",row.names = F,col.names = T)
write.table(atac.tcells.adj.data,"~/Desktop/SSP/tcell-peaks6_23_17/tcells_merged_consensus_adj_atac_whitelisted.txt",quote = F,sep = "\t",row.names = F,col.names = T)
write.table(cbind(atac.tcells.raw[,c("chr","start","end")],atac.tcells.lmtop),"~/Desktop/SSP/tcell-peaks6_23_17/tcells_merged_consensus_atac_results.txt",quote = F,sep = "\t",row.names = F,col.names = T)

# Inputs for visualization scripts

annotated.filter.bed$peakco <- paste(annotated.filter.bed$chr,"_",annotated.filter.bed$start,"_",annotated.filter.bed$end,sep = "")


atac.tcells.adj.data$peakco <- paste(atac.tcells.adj.data$chr,":",atac.tcells.adj.data$start,"-",atac.tcells.adj.data$end,sep = "")
full.results.age <- merge(annotated.filter.bed, atac.tcells.adj.data, by="peakco")
write.csv(full.results.age, file="~/Desktop/full_results_age.csv")

# Some visualizations.
# Make sure atac.tcells.lmtop$adj.P.Val<0.05 is not too big, or this may crash your computer. Start low (only plot, say, 2000 with lowest adj.P.Val) and increase from there
pheatmap(atac.tcells.adj.data[atac.tcells.lmtop$adj.P.Val<0.05,names(atac.tcells.libsize)],show_rownames = F,scale = "row",color = colorRampPalette(c("#3360a4","#1a8af9","cornsilk","firebrick3","darkred"))(100),cutree_cols = 2,clustering_distance_rows = "correlation",clustering_distance_cols = "correlation",clustering_method = "complete")
# MA-plot:
ggplot(aes(AveExpr,logFC,color=interaction(adj.P.Val<0.05,logFC<0)),
       data=atac.tcells.lmtop[order(-atac.tcells.lmtop$P.Value),]) +
  geom_point(size=0.5) +
  scale_color_manual(values = c("FALSE.FALSE"="black","TRUE.FALSE"="violetred2","FALSE.TRUE"="black","TRUE.TRUE"="darkorange"),guide=F) +
  geom_hline(yintercept = 0, size=0.25) +
  theme_bw(base_family = "Helvetica",base_size = 10) +
  theme(panel.border = element_rect(fill = NA,color = "black",size=0.25))
Hide

## RNAseq differential analysis is essentially the same as ATACseq's

rna.raw <- read.table("~/Desktop/SSP/monot_rnaseq.txt",sep = "\t",header = T,quote = "")
indices <- which(grepl("CD8", colnames(rna.raw)) | grepl("CD4", colnames(rna.raw)) | grepl("Memory", colnames(rna.raw)))
rna.raw <- rna.raw[, which(!colnames(rna.raw) %in% colnames(rna.raw)[indices])]
rna.samples <- rownames(group.names)
rna.raw$nhits <- rowSums(rna.raw[,which(!colnames(rna.raw) %in% c("EnsemblID","GeneName"))])

rna.raw <- rna.raw[which(rna.raw$GeneName %in% tsspeaks.df$GeneName), ] # might need to account for unincluded gene segments

# rna.raw$peakco <- paste(atac.tcells.raw$chr,":",atac.tcells.raw$start,"-",atac.tcells.raw$end,sep = "")
rna.raw$rna.max.val <- apply(rna.raw[, 3:101], 1, max)
test.rna.raw <- rna.raw
test.final <- do.call(rbind,lapply(split(test.rna.raw,test.rna.raw$GeneName),function(chunk) chunk[which.max(chunk$rna.max.val),]))
rna.raw <- test.final
rna.raw.mx <- as.matrix(data.frame(rna.raw[,which(!colnames(rna.raw) %in% c("EnsemblID", "nhits", "rna.max.val"))], row.names = "GeneName"))
in.group <- c()

for(i in 1:ncol(rna.raw.mx)){
  for(j in 1:nrow(group.names)){
    if(grepl(rownames(group.names)[j], colnames(rna.raw.mx)[i])){
      in.group <- append(in.group, i)
      colnames(rna.raw.mx)[i] <- rownames(group.names)[j]
    }
  }
}

rna.raw.mx <- rna.raw.mx[, in.group]
rna.libsize <- colSums(rna.raw.mx)

group.names.rna <- sampinfo[names(rna.libsize),] %>% select(sex,group,season) # Use this; accounts for discrepancies between samples seen in atac peaks                                                                            # and data that was filtered out of samples due to race, age, health, etc.
group.names.rna$group <- relevel(factor(group.names.rna$group),ref="HY")
group.names.rna$sex <- relevel(factor(group.names.rna$sex),ref="F")
group.names <- data.frame(sex, age.group, season) # All filtered values

rna.design <- with(group.names.rna ,model.matrix(~sex+season+group)) # in our paper, our design was sex+season+group
rna.design.null <- with(group.names.rna ,model.matrix(~sex+season,data = data.frame(rna.libsize))) # the null design is sex+season (i.e. everything but what we are interested in)

# This is for SVA, surrogate variable analysis. You can read a bit about it and we can chat later on its purpose and benefits
rna.cpm <- cpm(rna.raw.mx,lib.size = rna.libsize)
rna.cpm <- rna.cpm[rowSums(1*(rna.cpm>3))>1,]

rna.cpm2 <- cpm(rna.raw.mx,lib.size = rna.libsize,log = F)
rna.cpm2 <- rna.cpm2[rowSums(1*(rna.cpm2>3))>1,]

n.sv <- num.sv(rna.cpm*1e6,rna.design)
sv <- svaseq(rna.cpm2*1e6,rna.design,rna.design.null,n.sv = n.sv)$sv
rna.design <- cbind(rna.design,sv)
rna.design.null <- cbind(rna.design.null,sv)

## VOOM-based DA
rna.dge <- voom(rna.raw.mx,design = rna.design,plot = F)
rna.norm <- rna.dge$E
rna.lm <- lmFit(rna.dge,design = rna.design)
rna.lmfit <- eBayes(contrasts.fit(rna.lm,coefficients = "groupHO"),trend = T,robust = T)
rna.lmtop <- topTable(rna.lmfit,number = nrow(rna.raw.mx),sort.by = "none",adjust.method = "BH")
rna.dge.null <- voom(rna.raw.mx,design = rna.design.null,plot = F)
rna.lm.null <- lmFit(rna.dge.null,design = rna.design.null)
sum(rna.lmtop$adj.P.Val<.05) # 526
hist(rna.lmtop$P.Value,100)

rna.norm.data <- cbind(rna.raw[,c("GeneName")],rna.norm)
rna.adj.data <-  cbind(rna.raw[,c("GeneName")],rna.norm - fitted.MArrayLM(rna.lm.null))

write.table(rna.norm.data,"~/Desktop/SSP/RNA-peaks/merged_consensus_voom_rna_whitelisted.txt",quote = F,sep = "\t",row.names = F,col.names = T)
write.table(rna.adj.data,"~/Desktop/SSP/RNA-peaks/merged_consensus_adj_rna_whitelisted.txt",quote = F,sep = "\t",row.names = F,col.names = T)
write.table(cbind(rna.raw[,c("GeneName")],rna.lmtop),"~/Desktop/SSP/RNA-peaks/merged_consensus_rna_results.txt",quote = F,sep = "\t",row.names = F,col.names = T)

genes.w.info <- cbind(rna.raw[,c("GeneName")],rna.lmtop)
genes.w.info <- genes.w.info[which(genes.w.info$adj.P.Val < 0.05), ]
colnames(genes.w.info)[1] <- "GeneName"
gene.ensembl <- merge(genes.w.info, monot.expressionSet[, 1:2], by= "GeneName")
gene.ensembl <- gene.ensembl[!duplicated(gene.ensembl$GeneName), ]

write.table(gene.ensembl[which(gene.ensembl$logFC <0), "EnsemblID"], "~/Desktop/rna_diff_neg.txt", quote = F,sep = "\t",row.names = F,col.names = T)
write.table(gene.ensembl[which(gene.ensembl$logFC >0), "EnsemblID"], "~/Desktop/rna_diff_pos.txt", quote = F,sep = "\t",row.names = F,col.names = T)

# RNA-seq visualizations
pheatmap(rna.adj.data[rna.lmtop$adj.P.Val<0.05,names(rna.libsize)],show_rownames = F,scale = "row",color = colorRampPalette(c("#3360a4","#1a8af9","cornsilk","firebrick3","darkred"))(100),cutree_cols = 2,clustering_distance_rows = "correlation",clustering_distance_cols = "correlation",clustering_method = "complete")
# MA-plot:
ggplot(aes(AveExpr,logFC,color=interaction(adj.P.Val<0.001,logFC<0)),
       data=rna.lmtop[order(-rna.lmtop$P.Value),]) +
  geom_point(size=0.5) +
  scale_color_manual(values = c("FALSE.FALSE"="black","TRUE.FALSE"="violetred2","FALSE.TRUE"="black","TRUE.TRUE"="darkorange"),guide=F) +
  geom_hline(yintercept = 0, size=0.25) +
  theme_bw(base_family = "Helvetica",base_size = 10) +
  theme(panel.border = element_rect(fill = NA,color = "black",size=0.25))
Hide

``
# This is an example of what I've used to annotate peaks with chromHMM states. The thing to keep in mind is that a single peak can be annotated with multiple chromHMM states, and so we use priority rules to decide on one annotation per peak. 

library(parallel)
library(GenomicRanges)

# The annotated peak files that are the input for the code below are obtained from running the consensus vetted peaks, obtained above, through the add_chromhmm.sh Unix script (@ scripts/annotation). You'll need to edit this script too, which is set up for PBMC and various other cell types. Note that this not only adds chromHMM annotations but also HOMER annotations (HOMER=annotatePeaks.pl). In addition to HOMER, you'll need to install BEDTools to run that.
atac.homer.tcells <- read.table("~/Desktop/SSP/annotations_6_23_17/monot_merged_consensus_raw_atac_whitelisted_annotated.txt",sep = "\t",header = T,quote = "")
# Assigns priority to retain a unique chromHMM state per peak
chmm.prior <- data.frame(
  Priority=1:18,
  chromHMMstate=c("TssA","EnhA1","EnhA2","EnhG1","EnhG2","EnhWk","Tx","TssFlnk","TssFlnkU","TssFlnkD",
                  "TxWk","TssBiv","EnhBiv","ReprPCWk","ReprPC","ZNF_Rpts","Het","Quies")
)
prior.assign <- function(h) {chmm.prior[grep(paste("^",h["chromHMMstate"],"$",sep=""),chmm.prior$chromHMMstate),"Priority"]}
atac.chromhmm.tcells <- mclapply(dir("~/Desktop/SSP/annotations_6_23_17/",pattern = "*chromHMM_*",full.names = T),function(f) {
  anno <- read.table(f,sep="\t",header = T,quote = "",col.names = c("chr","start","end","chromHMMstate"))
  anno$peakco <- paste(anno$chr,":",anno$start,"-",anno$end,sep = "")
  anno$Priority <- apply(anno,1,prior.assign)
  anno <- anno[order(anno$peakco,anno$Priority),]
  anno <- anno[!duplicated(anno$peakco),]
  rownames(anno) <- anno$peakco
  anno <- anno[order(anno$chr,anno$start),c("chr","start","end","chromHMMstate")]
  colnames(anno) <- sub("chromHMMstate$",paste("chromHMMstate",sub("\\.txt","",sub(".*chromHMM_","",f)),sep = "_"),colnames(anno))
  return(anno)
},mc.cores = 4) # you want to set n.cores to a number that doesn't crash your machine. Use detectCores() to figure out how many cores you have

atac.glm.tcells <- read.table("~/Desktop/SSP/tcell-peaks6_23_17/tcells_merged_consensus_atac_results.txt",sep = "\t",header = T,quote = "")

# Peak filter
annotated.filter.bed <- data.frame(sort(GRanges(merge(atac.glm.tcells,
                                                      merge(atac.homer.tcells,
                                                            Reduce(function(x,y) merge(x,y,by = c("chr","start","end")),atac.chromhmm.tcells),
                                                            by = c("chr","start","end")),
                                                      by = c("chr","start","end")))))
colnames(annotated.filter.bed) <- sub("AveExpr","logCPM.atac",sub("^logFC$","logFC.atac",sub("^P.Value$","PValue.atac",sub("^adj.P.Val$","FDR.atac",sub("^GeneName$","GeneName.homer",sub("^seqnames$","chr",colnames(annotated.filter.bed))))))) # you may or may not include this step. I rename the columns to an old standard I use
colnames(annotated.filter.bed)[16] <- "chromHMM"


annotated.filter.pos.bed <- annotated.filter.bed[which(annotated.filter.bed$logFC.atac > 0 & annotated.filter.bed$FDR.atac < 0.05), 1:3]
annotated.filter.neg.bed <- annotated.filter.bed[which(annotated.filter.bed$logFC.atac < 0 & annotated.filter.bed$FDR.atac < 0.05), 1:3]

annotated.filter.final <- annotated.filter.bed[which(annotated.filter.bed$FDR.atac < 0.05), ]

annotated.logCPM <- annotated.filter.final[order(-annotated.filter.final$logCPM.atac), ]
annotated.logCPM.top <- annotated.logCPM[1:nrow(annotated.logCPM)/4, 1:3]
annotated.logCPM.bottom <- annotated.logCPM[(3*nrow(annotated.logCPM)/4):nrow(annotated.logCPM), 1:3]

write.table(annotated.filter.pos.bed, file="~/annotated.filter.pos.bed", quote=F, sep="\t", row.names=F, col.names=F)
write.table(annotated.filter.neg.bed, file="~/annotated.filter.neg.bed", quote=F, sep="\t", row.names=F, col.names=F)

write.table(annotated.logCPM.top, file = "~/annotated.logCPM.top.bed", quote=F, sep="\t", row.names = F, col.names = F)
write.table(annotated.logCPM.bottom, file= "~/annotated.logCPM.bottom.bed", quote=F, sep="\t", row.names = F, col.names = F)

for(i in 1:length(unique(annotated.filter.final$chromHMM))){
  type.chromHMM <- unique(annotated.filter.final$chromHMM)[i]
  fn <- paste("~/annotated.chromHMM.", type.chromHMM, ".bed", sep = "")
  write.table(annotated.filter.final[which(annotated.filter.final$chromHMM %in% type.chromHMM), 1:3], file = fn, quote=F, sep="\t", row.names = F, col.names = F)
}

# Mapping ATAC-seq and RNA-seq peaks to each other by proximity

atac2rna <- annotated.filter.final
atac2rna <- atac2rna[atac2rna$GeneName.homer %in% rownames(rna.peaks), ]
atac2rna.1000 <- atac2rna[which(abs(as.numeric(atac2rna$DistancetoTSS)) <= 1000), ]

rna.peaks <- rna.lmtop[rna.lmtop$adj.P.Val<.05, ] #510

# Combine enhancers (except weak) for logFC < 0
write.table(annotated.filter.final[which(grepl("Enh", as.character(annotated.filter.final$chromHMM)) & annotated.filter.final$logFC.atac<0), 1:3], file = "~/closing.enhancers.bed", quote=F, sep="\t", row.names = F, col.names = F)

# Combine flanking regions
regions <- c("TssFlnk", "TssFlnkU", "TssFlnkD")
write.table(annotated.filter.final[which(as.character(annotated.filter.final$chromHMM) %in% regions), 1:3], file = "~/flanking.tss.bed", quote=F, sep="\t", row.names = F, col.names = F)

# Boxplots to compare logFC of each enhancer -- opening/closing?
annotated.enhancers <- annotated.filter.final[which(grepl("Enh", as.character(annotated.filter.final$chromHMM))),]
with(annotated.enhancers, boxplot(logFC.atac~chromHMM, las=2, ylab = "logFC", col = c("palegreen", "green", "red", "lemonchiffon", "yellow", "purple")))

# Check for six groups in the paper and implement here

# Import annotations for ATAC-seq peaks

temp = list.files(pattern="*.tsv", path = "~/Desktop/SSP/GREAT_annotations/")
addresses <- paste("~/Desktop/SSP/GREAT_annotations/", temp, sep = "")
list2env(
  lapply(setNames(addresses, make.names(gsub("*.tsv$", "", temp))), 
         read.csv, sep = "\t", header = F, stringsAsFactors =F), envir = .GlobalEnv)

varnames <- ls(pattern="great.*")
for(i in varnames){
  mx <- get(i)
  mx[, ncol(mx)] <- NULL
  colnames(mx) <- mx[1,]
  mx = mx[-1, ]  
  assign(as.character(i), mx, envir= .GlobalEnv)
}

# Import annotations for two gene sets

temp = list.files(pattern="*.txt", path = "~/Desktop/SSP/DAVID_annotations/")

addresses <- paste("~/Desktop/SSP/DAVID_annotations/", temp, sep = "")
list2env(
  lapply(setNames(addresses, make.names(gsub("*.txt$", "", temp))), 
         read.csv, sep = "\t", header = T, stringsAsFactors =F), envir = .GlobalEnv)

varnames = ls(pattern = "Functional*")
for(i in varnames){
  mx <- get(i)
  mx <- mx[which(!mx$Category %in% "Category"), ]
  colnames(mx)[which(colnames(mx) %in% "X.")] <- "%"
  assign(as.character(i), mx, envir= .GlobalEnv)
}

rm(list = ls(pattern = "*.txt"))
rm(list = ls(pattern = "*.tsv"))

# Process data in negative clustering file
neg_cluster <- read.csv(file="~/Desktop/SSP/DAVID_annotations/Functional_Annotation_Clustering_neg.txt", sep="\t", header = F, stringsAsFactors = F)
cluster_info <- neg_cluster[grep("Annotation Cluster", neg_cluster$V1),]
rownames(cluster_info) <- as.numeric(rownames(cluster_info))
neg_cluster <- neg_cluster[!rownames(neg_cluster) %in% grep("Annotation Cluster", neg_cluster$V1),]

annotation_cluster <- rep("", nrow(neg_cluster))
enrichment_val <- rep("", nrow(neg_cluster))
for(i in 1:nrow(neg_cluster)){
  row <- as.numeric(rownames(neg_cluster)[i])
  annotation_cluster[i] = cluster_info[tail(rownames(cluster_info)[as.numeric(rownames(cluster_info)) <= row], 1), 1]
  enrichment_val[i] = cluster_info[tail(rownames(cluster_info)[as.numeric(rownames(cluster_info)) <= row], 1), 2]
  enrichment_val[i] <- substr(enrichment_val[i], 18, nchar(enrichment_val[i]))
}
neg_cluster$annotation_cluster <- annotation_cluster
neg_cluster$enrichment_val <- as.numeric(enrichment_val)

colnames(neg_cluster) <- neg_cluster[1,]
neg_cluster <- neg_cluster[-1,]
neg_cluster <- neg_cluster[-grep("Term", neg_cluster$Term),]
rownames(neg_cluster) <- 1:nrow(neg_cluster)
colnames(neg_cluster)[14:15] <- c("Annotation Cluster", "Enrichment Score")

library(dplyr)
library(tidyr)
neg_cluster <- separate_rows(neg_cluster, Genes, sep = ",")
neg_cluster_gene_list <- getBM(attributes = c("ensembl_gene_id","external_gene_name","gene_biotype","chromosome_name"),
                             values = neg_cluster$Genes,
                             mart=useMart(biomart="ENSEMBL_MART_ENSEMBL",
                                          dataset="hsapiens_gene_ensembl",
                                          host="grch37.ensembl.org"))
for(i in 1:nrow(neg_cluster)){
  neg_cluster[i, "Genes"] <- gsub(" ", "", neg_cluster[i, "Genes"], fixed = TRUE)
  neg_cluster[i, "Genes"] <- neg_cluster_gene_list[which(neg_cluster_gene_list$ensembl_gene_id %in% neg_cluster[i, "Genes"]), "external_gene_name"]
  neg_cluster[i, "Term"] <- gsub("^.*?~","",neg_cluster[i, "Term"])
  neg_cluster[i, "Term"] <- gsub("^.*?:","",neg_cluster[i, "Term"])
  
  # neg_cluster[i, "Bonferroni"] <- as.numeric(strsplit(neg_cluster[i,"Bonferroni"], "E")[[1]][1]) * 10^as.numeric(strsplit(neg_cluster[i,"Bonferroni"], "E")[[1]][2])
  # neg_cluster[i, "Benjamini"] <- as.numeric(strsplit(neg_cluster[i,"Benjamini"], "E")[[1]][1]) * 10^as.numeric(strsplit(neg_cluster[i,"Benjamini"], "E")[[1]][2])
  # neg_cluster[i, "FDR"] <- as.numeric(strsplit(neg_cluster[i,"FDR"], "E")[[1]][1]) * 10^as.numeric(strsplit(neg_cluster[i,"FDR"], "E")[[1]][2])
  # neg_cluster[i, "PValue"] <- as.numeric(strsplit(neg_cluster[i,"PValue"], "E")[[1]][1]) * 10^as.numeric(strsplit(neg_cluster[i,"PValue"], "E")[[1]][2])
  # 
  # neg_cluster[i, c("Count", "%", "List Total", "Pop Hits", "Pop Total", "Fold Enrichment")] <- as.numeric(neg_cluster[i, c("Count", "%", "List Total", "Pop Hits", "Pop Total", "Fold Enrichment")])
}

for(i in c("Count", "%", "List Total", "Pop Hits", "Pop Total", "Fold Enrichment", "Bonferroni", "Benjamini", "FDR", "PValue")){
  neg_cluster[, i] <- as.numeric(neg_cluster[, i])
}

options(scipen=999)
neg_cluster <-format(neg_cluster, scientific = FALSE)
write.csv(neg_cluster, file = "~/Desktop/neg_cluster")

install.packages("igraph")
library("igraph")
