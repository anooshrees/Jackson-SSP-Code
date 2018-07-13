################################
#
# Processing for TF networks
#
################################

geneset.info <- load("~/Downloads/geneset.info.RData")
network.data <- load("~/Downloads/aging_fpanalysis_gene.target.mode.radius_p900_consensus.filtered_mbf1_peaks.PD_explicit (1).RData")

temp = list.files(pattern="*.txt", path = "~/Downloads/genesetAnalysis/")
addresses <- paste("~/Downloads/genesetAnalysis/", temp, sep = "")
list2env(
  lapply(setNames(addresses, make.names(gsub("*.txt$", "", temp))), 
         read.csv, sep = "\t", header = T, stringsAsFactors =F), envir = .GlobalEnv)

atac.samples <- sample.info$samp
rna.samples <- colnames(source.data$RNA)[3:ncol(source.data$RNA)]
gene.stats <- summary.data$RNA

################################
#
# Filtering 
#
################################

# Some general filtering

# by ATAC-seq FDR and P-value
fp.data <- fp.data[intersect(which(fp.data$FDR.atac < 0.05), which(fp.data$PValue.atac < 0.01)), ]

# by delta binding and delta expression
fp.data.negative <- fp.data[fp.data$HY_HO.deltaProp < 0, ]
fp.data.negative <- fp.data.negative[abs(fp.data.negative$HY_HO.deltaProp) > quantile(abs(fp.data.negative$HY_HO.deltaProp),prob=1-1/100),]

fp.data.positive <- fp.data[fp.data$HY_HO.deltaProp > 0, ]
fp.data.positive <- fp.data.positive[abs(fp.data.positive$HY_HO.deltaProp) > quantile(abs(fp.data.positive$HY_HO.deltaProp),prob=1-1/100),]

fp.data <- rbind(fp.data.negative, fp.data.positive)

# When selecting genes for networks, there are three options for ensuring each TF connnects to only one gene
# 1) the closest gene:
fp.data.single <- do.call(rbind,lapply(split(fp.data,fp.data$fpco),function(chunk) chunk[which.min(chunk$DistancetoTSS),]))
# 2) the gene with highest logFC expression in relation to TF binding
fp.data.single <- do.call(rbind,lapply(split(fp.data,fp.data$tf.name),function(chunk) chunk[which.max(abs(gene.stats[gene.stats$GeneName %in% chunk$GeneName, "logFC.rna"])),]))
# 3) the gene whose accessibility matches expression the closest (maximizing logFC product)
fp.data.single <- do.call(rbind,lapply(split(fp.data,fp.data$tf.name),function(chunk) chunk[which.max(abs(gene.stats[gene.stats$GeneName %in% chunk$GeneName, "logFC.rna"]*chunk$logFC.atac)),]))

#########################
#
# Filtering by gene sets
#
#########################

# GOBP
all.names <- geneset.names.gobp
atac.enrichment <- atac_do.gobp_enrichment
rna.enrichment <- rna_do.gobp_enrichment
genes <- geneset.genes.gobp
gene.color <- data.frame(unlist(geneset.colors.gobp), stringsAsFactors = FALSE)

# GOCC
all.names <- geneset.names.gocc
atac.enrichment <- atac_do.gocc_enrichment
rna.enrichment <- rna_do.gocc_enrichment
genes <- geneset.genes.gocc
gene.color <- data.frame(unlist(geneset.colors.gocc), stringsAsFactors = FALSE)

# GOMF
all.names <- geneset.names.gomf
atac.enrichment <- atac_do.gomf_enrichment
rna.enrichment <- rna_do.gomf_enrichment
genes <- geneset.genes.gomf
gene.color <- data.frame(unlist(geneset.colors.gomf), stringsAsFactors = FALSE)

# GWAS
all.names <- geneset.names.gwas
atac.enrichment <- atac_do.gwas_enrichment
rna.enrichment <- rna_do.gwas_enrichment
genes <- geneset.genes.gwas
gene.color <- data.frame(unlist(geneset.colors.gwas), stringsAsFactors = FALSE)

# IMM_COARSE
all.names <- geneset.names.imm_coarse
atac.enrichment <- atac_do.imm_coarse_enrichment
rna.enrichment <- rna_do.imm_coarse_enrichment
genes <- geneset.genes.imm_coarse
gene.color <- data.frame(unlist(geneset.colors.imm_coarse), stringsAsFactors = FALSE)

# IMM_FINE
all.names <- geneset.names.imm_fine
atac.enrichment <- atac_do.imm_fine_enrichment
rna.enrichment <- rna_do.imm_fine_enrichment
genes <- geneset.genes.imm_fine
gene.color <- data.frame(unlist(geneset.colors.imm_fine), stringsAsFactors = FALSE)

# INTERPRO
all.names <- geneset.names.interpro
atac.enrichment <- atac_do.interpro_enrichment
rna.enrichment <- rna_do.interpro_enrichment
genes <- geneset.genes.interpro
gene.color <- data.frame(unlist(geneset.colors.interpro), stringsAsFactors = FALSE)

# KEGG
all.names <- geneset.names.kegg
atac.enrichment <- atac_do.kegg_enrichment
rna.enrichment <- rna_do.kegg_enrichment
genes <- geneset.genes.kegg
gene.color <- data.frame(unlist(geneset.colors.kegg), stringsAsFactors = FALSE)

# MSIGDB
all.names <- geneset.names.msigdb
atac.enrichment <- atac_do.msigdb_enrichment
rna.enrichment <- rna_do.msigdb_enrichment
genes <- geneset.genes.msigdb
gene.color <- data.frame(unlist(geneset.colors.msigdb), stringsAsFactors = FALSE)

# PID
all.names <- geneset.names.pid
atac.enrichment <- atac_do.pid_enrichment
rna.enrichment <- rna_do.pid_enrichment
genes <- geneset.genes.pid
gene.color <- data.frame(unlist(geneset.colors.pid), stringsAsFactors = FALSE)

# REACTOME
all.names <- geneset.names.reactome
atac.enrichment <- atac_do.reactome_enrichment
rna.enrichment <- rna_do.reactome_enrichment
genes <- geneset.genes.reactome
gene.color <- data.frame(unlist(geneset.colors.reactome), stringsAsFactors = FALSE)

# VP2008
all.names <- geneset.names.vp2008
atac.enrichment <- atac_do.vp2008_enrichment
rna.enrichment <- rna_do.vp2008_enrichment
genes <- geneset.genes.vp2008
gene.color <- data.frame(unlist(geneset.colors.vp2008), stringsAsFactors = FALSE)

# VP2015
all.names <- geneset.names.vp2015
atac.enrichment <- atac_do.vp2015_enrichment
rna.enrichment <- rna_do.vp2015_enrichment
genes <- geneset.genes.vp2015
gene.color <- data.frame(unlist(geneset.colors.vp2015), stringsAsFactors = FALSE)

# WP
all.names <- geneset.names.wp
atac.enrichment <- atac_do.wp_enrichment
rna.enrichment <- rna_do.wp_enrichment
genes <- geneset.genes.wp
gene.color <- data.frame(unlist(geneset.colors.wp), stringsAsFactors = FALSE)

# General code to reformat color matrix
gene.color$Module.ID <- (all.names %>% slice(match(rownames(gene.color), Module.Name)) %>% select(Module.ID))$Module.ID
colnames(gene.color) <- c("color", "Module.ID")
gene.color$GeneName <- (genes %>% slice(match(gene.color$Module.ID, Module.ID)) %>% select(GeneName))$GeneName

# Refilter fp.data -- may not be necessary...
fp.data.example <- fp.data.single[fp.data.single$GeneName %in% genes[genes$Module.ID %in% all.names$Module.ID, "GeneName"], ]
fp.data.example$Module.ID <- (genes %>% slice(match(fp.data.example$GeneName,GeneName)) %>% select(Module.ID))$Module.ID


################################
#
# Function that creates all networks for all gene sets
#
################################

# genesets <- c("gobp", "gocc", "gomf", "gwas", "imm_coarse", "imm_fine", "interpro", "kegg", "msigdb", "pid", "reactome", "vp2008", "vp2015", "wp")
# genesets <- c("kegg", "msigdb", "pid", "reactome", "vp2008", "vp2015", "wp")
genesets <- c("kegg", "vp2008", "gobp", "reactome")

generate_networks <- function(x){
  for(geneset.name in genesets){
    # since data is presented differently for some of the gene sets, individual debugging may be necessary
    # geneset.name <- "gobp"
    all.names <- get(paste("geneset.names.", geneset.name, sep = ""))
    atac.enrichment <- get(paste("atac_do.", geneset.name, "_enrichment", sep = ""))
    rna.enrichment <- get(paste("rna_do.", geneset.name, "_enrichment", sep = ""))

    all.names <- all.names[all.names$Module.Name %in% rna.enrichment[rna.enrichment$hypergeom.fdr < 0.05, "Module.Name"], ]
    genes <- get(paste("geneset.genes.", geneset.name, sep= ""))
    genes <- genes[genes$Module.ID %in% all.names$Module.ID, ]
    
    gene.color <- data.frame(unlist(get(paste("geneset.colors.", geneset.name, sep = ""))), stringsAsFactors = FALSE)
    gene.color$Module.Name <- names(get(paste("geneset.colors.", geneset.name, sep = ""))) #kegg
    gene.color <- gene.color[gene.color$Module.Name %in% rna.enrichment[rna.enrichment$hypergeom.fdr < 0.05, "Module.Name"], ]
    # gene.color$Module.ID <- (all.names %>% slice(match(rownames(gene.color), Module.Name)) %>% select(Module.ID))$Module.ID
    gene.color$Module.ID <- (all.names %>% slice(match(gene.color$Module.Name, Module.Name)) %>% select(Module.ID))$Module.ID #kegg
    colnames(gene.color)[1] <- c("color")
    # gene.color$GeneName[gene.color$Module.ID %in% genes$Module.ID] <- (genes %>% slice(match(gene.color$Module.ID, Module.ID)) %>% select(GeneName))$GeneName #kegg1
    # gene.color <- gene.color[(which(!is.na(gene.color$GeneName))),]
    
    print(paste("Working on geneset", geneset.name))
    
    for(current.module in all.names$Module.ID){
      fp.data.example <- fp.data[fp.data$GeneName %in% current_geneset, ]
      
      if(nrow(fp.data.example) > 2){
        fp.data.example$Module.ID <- (genes %>% slice(match(fp.data.example$GeneName,GeneName)) %>% select(Module.ID))$Module.ID
        fp.data.example <-  do.call(rbind,lapply(split(fp.data.example,fp.data.example$fpco),function(chunk) chunk[which.min(chunk$DistancetoTSS),]))
        
        df <- data.frame(fp.data.example %>% select(tf_distype,GeneName), stringsAsFactors = FALSE)
        # ***If you want to differentiate between proximal and distal, use:***
        # df <- data.frame(fp.data.example %>% select(tf_distype, GeneName), stringsAsFactors = FALSE)
        df <- df[!duplicated(df), ]
        nw <- graph_from_data_frame(df,directed=T)
        
        edge_attr(nw,"deltaProp") <- fp.data.example$HY_HO.deltaProp[intersect(which(fp.data.example$GeneName %in% df$GeneName), which(fp.data.example$tf_distype %in% df$tf_distype))]
        edge_attr(nw,"deltaProp_color")[which(edge_attr(nw,"deltaProp") < 0)] <- "blue"
        edge_attr(nw,"deltaProp_color")[which(edge_attr(nw,"deltaProp") > 0)] <- "red"
        
        vertex_attr(nw,"logFC")[which(vertex_attr(nw, "name") %in% fp.data.example$GeneName)] <- (gene.stats %>% slice(match(vertex_attr(nw, "name"),GeneName)) %>% select(logFC.rna))$logFC.rna
        # vertex_attr(nw, "logFC") <- 
        vertex_attr(nw,"logFC_color")[which(vertex_attr(nw,"logFC") < 0)] <- "blue"
        vertex_attr(nw,"logFC_color")[which(vertex_attr(nw,"logFC") > 0)] <- "red"
        # vertex_attr(nw, "color")[!vertex_attr(nw, "name") %in% df$tf.name] <- "white"
        vertex_attr(nw, "is.tf") <- vertex_attr(nw,"name") %in% fp.data.example$tf_distype
        vertex_attr(nw,"logFC_color")[vertex_attr(nw,"is.tf")] <- NULL
        vertex_attr(nw,"transparency")[vertex_attr(nw,"is.tf")] <- "transparent"
        vertex_attr(nw,"transparency")[!vertex_attr(nw,"is.tf")] <- 0.2
        
        vertex_attr(nw,"logFC") <- rescale(abs(vertex_attr(nw,"logFC")), c(5, 30))
        vertex_attr(nw,"logFC")[vertex_attr(nw, "is.tf")] <- 0
        nw.gplot <- ggnet2(nw,
                           node.size = vertex_attr(nw,"logFC"),
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
                           alpha = vertex_attr(nw, "transparency"),
                           legend.position = "none")
        
        # Then just plot with p
        current.module.name <- all.names[all.names$Module.ID %in% current.module, "Module.Name"]
        pdf(paste("~/Desktop/Differentiation_General_Network", ".pdf", sep = ""))
        p <- nw.gplot + 
          # This replots the nodes, adding edges to the nodes drawn by ggnet2
          # I have included an attribute “is.tf” that informs whether a node is a TF or a target, and use this to plot TFs and targets with different shapes and sizes, and to format the labels differently
          geom_point(shape=as.numeric(sub("TRUE",22,sub("FALSE",21,vertex_attr(nw,"is.tf")))),
                     size=as.numeric(sub("TRUE",8,sub("FALSE",9,vertex_attr(nw,"is.tf")))),
                     color=sub("TRUE","maroon2",sub("FALSE","transparent",vertex_attr(nw,"is.tf"))), fill=NA) +
          # This adds labels to the nodes
          geom_text(aes(label=label),
                    size=as.numeric(sub("FALSE",2.5,sub("TRUE",2,vertex_attr(nw,"is.tf")))),
                    color=sub("FALSE","gray40",sub("TRUE","maroon4",vertex_attr(nw,"is.tf"))),hjust=-0.2)
        print(p)
        dev.off()
        
        print(paste("Just finished network for ", current.module, "with name ", current.module.name))
      }
    }
  }
}

generate_networks(genesets)

################################
#
# Filtering by metadata
#
################################

# upper or lower half of logCPM
fp.data.example <- fp.data.example[fp.data.example$logCPM.atac > median(fp.data.example$logCPM.atac), ]
fp.data.example <- fp.data.example[fp.data.example$logCPM.atac < median(fp.data.example$logCPM.atac), ]

# only positive logFC or only negative
fp.data.example <- fp.data.example[fp.data.example$logFC.atac > 0, ]
fp.data.example <- fp.data.example[fp.data.example$logFC.atac < 0, ]

# Distal or proximal location
fp.data.example <- fp.data.example[fp.data.example$DisType %in% "Proximal", ]
fp.data.example <- fp.data.example[fp.data.example$DisType %in% "Distal", ]

# Low FDR and P Val for gene expression
fp.data.example <- fp.data.example[fp.data.example$GeneName %in% gene.stats[intersect(which(gene.stats$FDR.rna < 0.05), which(gene.stats$PValue.rna < 0.01)), "GeneName"], ]

# Only negative logFC or only positive (for RNA)
fp.data.example <- fp.data.example[fp.data.example$GeneName %in% gene.stats[gene.stats$logFC.rna > 0, "GeneName"], ]
fp.data.example <- fp.data.example[fp.data.example$GeneName %in% gene.stats[gene.stats$logFC.rna < 0, "GeneName"], ]

################################
#
# Plotting
#
################################
library(igraph)
library(ggplot2)
library(GGally)
library(dplyr)
library(sna)
library(network)
library(igraph) # again to preserve certain functions

df <- data.frame(fp.data.example %>% select(tf.name,GeneName), stringsAsFactors = FALSE)
nw <- graph_from_data_frame(df,directed=T)

vertex_attr(nw,"logFC") <- (gene.stats %>% slice(match(vertex_attr(nw,"name"),GeneName)) %>% select(logFC.rna))$logFC.rna
vertex_attr(nw, "color")[vertex_attr(nw, "name") %in% df$GeneName] <- (gene.color %>% slice(match(vertex_attr(nw,"name"), GeneName)) %>% select(color))$color
vertex_attr(nw, "color")[!vertex_attr(nw, "name") %in% df$GeneName] <- rep("white", length(which(!vertex_attr(nw, "name") %in% df$GeneName)))
vertex_attr(nw, "is.tf") <- vertex_attr(nw,"name") %in% df$tf.name

ggnet2(nw)

nw.gplot <- ggnet2(nw,
                   # node.size = vertex_attr(nw,"logFC"),
                   # node.shape = vertex_attr(nw,"shape"),
                   label = F, # I use ggplot (below) to add labels
                   mode = "fruchtermanreingold",
                   arrow.size = 5,
                   arrow.type = "open",
                   arrow.gap = 0.025,
                   # edge.color = edge_attr(nw,"color"),
                   edge.size = 0.5,
                   node.alpha = 0.8,
                   node.color = unlist(vertex_attr(nw,"color")),
                   legend.position = "none")

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

# Then just plot with p
p

#######################
#
# Optimized Networks
#
#######################

# Networks with proximal/distal annotation

network.data.pd <- load("~/Downloads/geneset.data_playground.p900_consensus.radius_RNA_explicit_fullnet.RData")
path.all.negative <- path.all[path.all$HY_HO.mean_deltaProp < 0, ]
path.all.negative <- path.all.negative[abs(path.all.negative$HY_HO.mean_deltaProp) > quantile(abs(path.all.negative$HY_HO.mean_deltaProp),prob=1-50/100),]

path.all.positive <- path.all[path.all$HY_HO.mean_deltaProp > 0, ]
path.all.positive <- path.all.positive[abs(path.all.positive$HY_HO.mean_deltaProp) > quantile(abs(path.all.positive$HY_HO.mean_deltaProp),prob=1-50/100),]

path.all <- rbind(path.all.positive, path.all.negative)

for(geneset.name in unique(path.all$schema)){
  # since data is presented differently for some of the gene sets, individual debugging may be necessary
  # geneset.name <- "reactome"
  # all.info <- path.all[path.all$schema %in% geneset.name, ]
  # 
  # print(paste("Working on geneset", geneset.name))
  current_geneset <- test.table.tcell$GeneName
  
  for(current.module in unique(all.info$pathway)){
    current.module <- "Cellular Senescence"
    fp.data.example <- all.info[all.info$pathway %in% current.module, ]
    
    fp.data.example <- fp.data.example[abs(fp.data.example$signed_logOR) > quantile(abs(fp.data.example$signed_logOR), prob=1-75/100), ]
    
    if(nrow(fp.data.example) > 2){
      df <- data.frame(fp.data.example %>% select(tf.name,GeneName), stringsAsFactors = FALSE)
      # ***If you want to differentiate between proximal and distal, use:***
      # df <- data.frame(fp.data.example %>% select(tf_distype, GeneName), stringsAsFactors = FALSE)
      df <- df[!duplicated(df), ]
      nw <- graph_from_data_frame(df,directed=T)
      
      edge_attr(nw,"deltaProp") <- fp.data.example$HY_HO.mean_deltaProp
      edge_attr(nw,"deltaProp_color")[which(edge_attr(nw,"deltaProp") < 0)] <- "blue"
      edge_attr(nw,"deltaProp_color")[which(edge_attr(nw,"deltaProp") > 0)] <- "red"
      
      gene.names <- unique(df$GeneName)
      tf.names <- unique(df$tf.name)
      
      vertex_attr(nw, "logFC")[which(vertex_attr(nw, "name") %in% gene.names)] <- unique(fp.data.example[fp.data.example$GeneName %in% gene.names, "logFC.rna"])
      vertex_attr(nw, "logFC")[which(vertex_attr(nw, "name") %in% tf.names)] <- rep(0, length(vertex_attr(nw, "logFC")[which(vertex_attr(nw, "name") %in% tf.names)]))
      
      vertex_attr(nw,"logFC_color")[which(vertex_attr(nw,"logFC") < 0)] <- "blue"
      vertex_attr(nw,"logFC_color")[which(vertex_attr(nw,"logFC") > 0)] <- "red"
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
      pdf(paste("~/Desktop/", geneset.name, "_", gsub("[[:punct:]]", " ", gsub(",.*$", "", current.module)), ".pdf", sep = ""))
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
      
      print(paste("Just finished network for ", current.module, "with name ", current.module.name))
    }
  }
}

# Networks without annotation (just TF names)

network.data.full <- load("~/Downloads/geneset.data_playground.p900_consensus.radius_RNA__fullnet.RData")
path.all.negative <- path.all[path.all$HY_HO.mean_deltaProp < 0, ]
path.all.negative <- path.all.negative[abs(path.all.negative$HY_HO.mean_deltaProp) > quantile(abs(path.all.negative$HY_HO.mean_deltaProp),prob=1-25/100),]

path.all.positive <- path.all[path.all$HY_HO.mean_deltaProp > 0, ]
path.all.positive <- path.all.positive[abs(path.all.positive$HY_HO.mean_deltaProp) > quantile(abs(path.all.positive$HY_HO.mean_deltaProp),prob=1-25/100),]

path.all <- rbind(path.all.positive, path.all.negative)

for(geneset.name in unique(path.all$schema)){
  # since data is presented differently for some of the gene sets, individual debugging may be necessary
  # geneset.name <- "gobp"
  all.info <- path.all[path.all$schema %in% geneset.name, ]
  
  print(paste("Working on geneset", geneset.name))
  
  for(current.module in unique(all.info$pathway)){
    fp.data.example <- all.info[all.info$pathway %in% current.module, ]
    
    if(nrow(fp.data.example) > 2){
      df <- data.frame(fp.data.example %>% select(tf.name,GeneName), stringsAsFactors = FALSE)
      # ***If you want to differentiate between proximal and distal, use:***
      # df <- data.frame(fp.data.example %>% select(tf_distype, GeneName), stringsAsFactors = FALSE)
      df <- df[!duplicated(df), ]
      nw <- graph_from_data_frame(df,directed=T)
      
      edge_attr(nw,"deltaProp") <- fp.data.example$HY_HO.mean_deltaProp
      edge_attr(nw,"deltaProp_color")[which(edge_attr(nw,"deltaProp") < 0)] <- "blue"
      edge_attr(nw,"deltaProp_color")[which(edge_attr(nw,"deltaProp") > 0)] <- "red"
      
      gene.names <- unique(df$GeneName)
      tf.names <- unique(df$tf.name)
      
      vertex_attr(nw, "logFC")[which(vertex_attr(nw, "name") %in% gene.names)] <- unique(fp.data.example[fp.data.example$GeneName %in% gene.names, "logFC.rna"])
      vertex_attr(nw, "logFC")[which(vertex_attr(nw, "name") %in% tf.names)] <- rep(0, length(vertex_attr(nw, "logFC")[which(vertex_attr(nw, "name") %in% tf.names)]))
      
      vertex_attr(nw,"logFC_color")[which(vertex_attr(nw,"logFC") < 0)] <- "blue"
      vertex_attr(nw,"logFC_color")[which(vertex_attr(nw,"logFC") > 0)] <- "red"
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
      pdf(paste("~/Desktop/SSP/Networks_7_22_17/", "less_detail_", geneset.name, "_", gsub("[[:punct:]]", " ", gsub(",.*$", "", current.module)), ".pdf", sep = ""))
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
      
      print(paste("Just finished network for ", current.module, "with name ", current.module.name))
    }
  }
}

#######################
#
# Export to Gephi
#
#######################
library(rgexf)

find.edgeid <- function(edgelabel,nodes) grep(paste("^",edgelabel,"$",sep = ""),nodes$label)

# network
nw <- topnw.all

# nodes
nodes <- data.frame(id=seq_len(vcount(nw)),
                    label=vertex_attr(nw,"name"),
                    stringsAsFactors = F)
nodeAtt <- data.frame(logFC_color=vertex_attr(nw,"color"),
                      logFC.rna=vertex_attr(nw,"logFC"))
                      # specificity=vertex_attr(nw,"specificity"),
                      # target.scaled_freq=(data.frame(tf.gene=vertex_attr(nw,"name"),stringsAsFactors = F) %>% full_join(ylong_bytargets_freqs %>% select(-target.specificity) %>% mutate(tf.gene=as.character(tf.gene))) %>% slice(match(vertex_attr(nw,"name"),tf.gene)) %>% select(freq) %>% mutate(freq=2*(freq-0.5)) %>% mutate(freq=ifelse(is.na(freq),0,freq)))$freq)
nodeVizAtt <- list(size=as.numeric(sub("TRUE",5.5,sub("FALSE",6,vertex_attr(nw,"is.tf")))),
                   shape=sub("TRUE","square",sub("FALSE","circle",vertex_attr(nw,"is.tf"))),
                   color=as.matrix(setNames(data.frame(t(col2rgb(vertex_attr(nw,"color"))),a=0.8),c("r","g","b","a"))))  # color nodes by logFC
# color=as.matrix(setNames(data.frame(t(col2rgb(vec2color.symmetrized(nodeAtt$target.scaled_freq,colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu")[-c(4,7)])(nrow(nodeAtt))))),a=0.8),c("r","g","b","a"))))  # color nodes by freq of cell specific targets (-1,1, positive for bcells, negative for tcells, zero for unbiased and for non-tf genes)

# edges 
edges.df <-as_data_frame(nw)
colnames(edges.df) <- sub("color","slogOR_color",sub("from","source",sub("to","target",colnames(edges.df))))
edges <- sapply(edges.df[,c("source","target")],function(e) as.vector(sapply(e,find.edgeid,nodes)))
edgeAtt <- edges.df %>% select(source,slogOR_color,slogOR,affinity) %>% left_join(ylong_bytargets_freqs %>% mutate(source=as.character(tf.gene)) %>% select(source,freq),by = "source") %>% mutate(target.cum_specificity=sub("TRUE","TCells",sub("FALSE","BCells",(freq-0.5)<0))) %>% select(-source,-freq) # target.cum_specificity labels an edge by the overall proportion of cell-specifity of the targets. For example, if 55% of the targest of TF x are B-cell specific (ie. more expressed in B cell) then all edges from this TF are labeled "BCells"
edgeVizAtt <- list(size=rep(2,ecount(nw))
                   # color=as.matrix(setNames(data.frame(t(col2rgb(vec2color.symmetrized((edges.df %>% select(source,slogOR_color,slogOR,affinity) %>% left_join(ylong_bytargets_freqs %>% mutate(source=as.character(tf.gene)) %>% select(source,freq),by = "source") %>% mutate(freq=2*(freq-0.5)))$freq,colorRampPalette(colors_redblu10[-9])(nrow(edges.df))))),a=1),c("r","g","b","a")))
                   # color=setNames(data.frame(t(col2rgb(edges.df$slogOR_color)),a=1),c("r","g","b","a"))
)

# export
grfx <- write.gexf(nodes = nodes,
                   edges = edges,
#                   edgesAtt = edgeAtt,
                   nodesAtt = nodeAtt,
#                   edgesVizAtt = edgeVizAtt,
                   nodesVizAtt = nodeVizAtt,
                   output = "test1.gexf",
                   defaultedgetype = "directed")

###############################
#
# T Cell differentiation geneset
#
###############################

test.table.tcell <- read.csv(file = "~/Desktop/SSP/Table1.txt", sep=" ", header = TRUE)

# library(biomaRt)
# tcells.genes <- getBM(attributes = c("ensembl_gene_id","external_gene_name","gene_biotype","chromosome_name"),
#                                              filters = "family",
#                                              values = test.table.tcell$Antigen,
#                                              mart=useMart(biomart="ENSEMBL_MART_ENSEMBL",
#                                                           dataset="hsapiens_gene_ensembl",
#                                                           host="grch37.ensembl.org"))

# write.table(paste(gsub("-", "", test.table.tcell$Antigen), "_HUMAN", sep=""), file="~/Desktop/protein_names.txt", quote = F,row.names = F,col.names = T)

#########################################
#
# Create ATAC logFC vs. RNA logFC plot to 
# prove that T cells further along in 
# linear differentiation are more present
# in older individuals
#
#########################################

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
colors_tcells= c("violet", "violetred", "violetred1", "violetred2", "violetred3", "violetred4")

# setwd("/Volumes/emext/Dropbox/Jax/PBMC/DifferentialAnalysis/manuscripts/aging-01/scripts/figures.regen.stm/")
# setwd("~/Dropbox/Jax/PBMC/DifferentialAnalysis/manuscripts/aging-01/scripts/figures.regen.stm/")
# setwd("H:/Dropbox/Jax/PBMC/DifferentialAnalysis/manuscripts/aging-01/scripts/figures.regen.stm/")

# function to relabel samples on-the-fly wrt randomIDs
sampinfo <- read.table("~/Desktop/SSP/data/jax/sample.info.comprehensive.txt",header = TRUE,sep = "\t",stringsAsFactors = F) # check data/jax
rownames(sampinfo) <- with(sampinfo,paste(sampid,type,sep = "_"))
sampinfo[sampinfo$age>40 & sampinfo$group=="HY",]$group <- "HM"

random.codes <- read.table("~/Downloads/RandomRecodeAgingStudy.txt", header = T, sep = "\t", quote = "")
random.codes <- data.frame(merge(data.frame(random.codes,sampid=sub("VHY402","VHY402D0",paste(random.codes$Cohort,random.codes$PID,sep = "")),row.names = "sampid"),sampinfo,by.x="row.names", by.y="sampid"))
rownames(random.codes) <- with(random.codes,paste(random.codes[,1],random.codes[,8],sep = "_"))
random.codes <- random.codes[, -1]
rand.recode <- function(x) with(random.codes[x,], paste(group,RandomID,sep=""))

# Previoulsy analyzed ATACseq data
# da.pbmc.out <- load("~/Downloads/da.pbmc.RData")
#full.results <- read.table("datasources.201702/pbmc_whitelisted_filtered_glm.results.txt", sep = "\t", quote = "", header = T,row.names = "peakco")
full.results <- read.csv("~/Desktop/SSP/DA_7_12_17/full_results_age.csv")
full.results$peakco <- paste(full.results$chr,"_",full.results$start,"_",full.results$end,sep = "")
rownames(full.results) <- full.results$peakco 

adj.info <- read.table("~/Desktop/SSP/DA_7_12_17/tcells_merged_consensus_adj_atac_whitelisted.txt",header = T, sep = "\t", quote = "")
adj.info$peakco <- paste(adj.info$chr, "_", adj.info$start, "_", adj.info$end, sep="")
adj.atac <- data.frame(merge(full.results[,c("chr","start","end")],adj.info, by.x = "row.names", by.y = "peakco"),row.names = 1)
atacsamps <- gsub("*_PBMC$", "", colnames(adj.atac[,-c(1:6)]))
rm(list = da.pbmc.out)

# Previoulsy analyzed RNAseq data (samps common in RNA and filtered ATAC)
da.rna.pbmc.out <- load("~/Downloads/da.rna.pbmc.RData")
rna.glmtop <- read.table("~/Downloads/pbmc_rna.atacpaired_glm.results.txt", sep = "\t", quote = "", header = T)

# Collect data from DA script
atac.glm.tcells$peakco <- paste(atac.glm.tcells$chr, "_", atac.glm.tcells$start, "_", atac.glm.tcells$end, sep="")
rna.glmtop2 <- merge(atac.glm.tcells, adj.info, by = "peakco")
tsspeaks.df$peakco <- paste(tsspeaks.df$chr, "_", tsspeaks.df$start, "_", tsspeaks.df$end, sep="")
rna.glmtop2 <- merge(tsspeaks.df, rna.glmtop2, by = "peakco")
rna.glmtop2 <- merge(annotated.filter.bed[, c("logCPM.atac", "chr", "start", "end")], rna.glmtop2, by = c("chr", "start", "end"))

# Reorganize and rename
colnames(rna.glmtop2)[21:64] <- gsub("*_PBMC$", "", colnames(rna.glmtop2)[21:64])
rna.glmtop2[, which(colnames(rna.glmtop2) %in% c("chr.x", "start.x", "end.x", "peakco", "seqnames", "chr.y", "start.y", "end.y", "start", "end", "width", "chr"))] <- NULL
colnames(rna.glmtop2)[1] <- c("logCPM")
colnames(rna.glmtop2)[7:8] <- c("PValue", "FDR")
rna.glmtop <- rna.glmtop2

rna.glmtop <- rna.glmtop[rna.glmtop$GeneName %in% test.table.tcell$GeneName, ]
rna.glmtop <- do.call(rbind,lapply(split(rna.glmtop,rna.glmtop$GeneName),function(chunk) chunk[which.min(chunk$FDR),]))

rownames(rna.glmtop) <- rna.glmtop$GeneName
# rnasamps <- gsub("_PBMC_subsampled", "", rna.samples)
# atacsamps <- gsub("_PBMC_subsampled", "", atac.samples)

adj.rna <- read.table(file = "~/Desktop/SSP/DA_7_12_17/merged_consensus_adj_rna_whitelisted.txt", sep = "\t", quote = "", header = T)

adj.rna$GeneName <- unique(aging.rna$GeneName)
rpkm <- da.pbmc$raw
rm(list = da.rna.pbmc.out)

samps <- intersect(rnasamps,atacsamps) # Samples that are those present in both datasets
expressed.genes <- rownames(rpkm)

samps <- samps[samps %in% gsub("_PBMC", "", colnames(adj.rna))]
rnasamps <- rnasamps[rnasamps %in% gsub("_PBMC", "", colnames(adj.rna))]
atacsamps <- atacsamps[atacsamps %in% gsub("_PBMC", "", colnames(adj.rna))]

# Splinters TSS peaks out to combine with RNAseq data
# atac.tss <- read.table("datasources.201702/raw_data_tss.peaks_minOver0_hard_stm.txt",quote = "",header = T, sep = "\t")
flankd = 1000
# tss.bed <- GRanges(read.table(paste("../../data/Promoters/tssflank",flankd,"_prom_complete.bed",sep = ""),header = F,sep = "\t",quote = "",col.names = c("chr","start","end","GeneName","score","strand")))
tss.bed <- GRanges(tsspeaks.df)
tss.bed.expressed <- subset(tss.bed,GeneName %in% expressed.genes)
tss_peaks <- findOverlaps(GRanges(full.results[,c("chr","start","end")]),tss.bed.expressed)
# The following ensures gene names are derived from TSS definitions instead of nearest-TSS annotations
full.results.tss_red <- data.frame(GeneName=as.character(tss.bed.expressed[tss_peaks@to]$GeneName),full.results[tss.bed.expressed[tss_peaks@to]$peakco,which(colnames(full.results)!="GeneName")],row.names = NULL)

# full.results.tss_red <- full.results.tss_red[full.results.tss_red$GeneName %in% test.table.tcell$GeneName, ]
# full.results.tss_red <- merge(full.results.tss_red,data.frame(meanExpression=rowMeans(rpkm)),by.x = "GeneName", by.y = "row.names") # to keep peak associated to highest-expressed gene
full.results.tss_red <- unique(full.results.tss_red[order(full.results.tss_red$GeneName,abs(full.results.tss_red$DistancetoTSS)),]) # to keep one peak per gene (multiple genes per peak allowed), kept peak closest to nearest TSS
full.results.tss <- full.results.tss_red[!duplicated(full.results.tss_red$GeneName),]
rownames(full.results.tss) <- full.results.tss$GeneName
rm("full.results.tss_red")

colnames(tsspeaks.df)[1] <- "chr"
colnames(adj.atac)[1:3] <- c("chr", "start", "end")
adj.atac <- adj.atac[, -c(4:6)]
adj.atac.tss <- merge(tsspeaks.df[,c("chr","start","end","GeneName")], adj.atac, by = c("chr","start","end"))
# adj.atac.tss.nodup <- adj.atac.tss[order(adj.atac.tss$GeneName,-rowSums(2^adj.atac.tss[,atacsamps])),c("GeneName",atacsamps)]
# adj.atac.tss.nodup <- adj.atac.tss.nodup[!duplicated(adj.atac.tss.nodup$GeneName),]

samps_pbmc <- paste(samps, "_", "PBMC", sep="")
rnasamps_pbmc <- paste(rnasamps, "_", "PBMC", sep="")
atacsamps_pbmc <- paste(atacsamps, "_", "PBMC", sep="")

atacgroup <- relevel(factor(sampinfo[atacsamps_pbmc[atacsamps_pbmc %in% colnames(adj.atac.tss)],]$group),ref = "HY")
atacsex <- relevel(factor(sampinfo[atacsamps_pbmc,]$sex),ref = "F")
rnagroup <- relevel(factor(sampinfo[rnasamps_pbmc,]$group),ref = "HY")
rnasex <- relevel(factor(sampinfo[rnasamps_pbmc,]$sex),ref = "F")
group <- relevel(factor(sampinfo[samps_pbmc %in% colnames(adj.atac.tss),]$group),ref = "HY")
sex <- relevel(factor(sampinfo[samps_pbmc,]$sex),ref = "F")

##### Permutation tests to derive empirical Pvalues for each gene ATAC/RNA fold change
nperm = 1000
#### RNAseq FC
### Code used to generate perms:
# rna.permlist <- replicate(nperm,adj.rna,simplify = F)
# rna.permfc <- mclapply(rna.permlist, function(Y) apply(Y,1,function(x) diff(aggregate(x, list(group=sample(group)), mean)[,2])),mc.cores = 6)
# save(list = c("rna.permlist","rna.permfc"),file = "datasources.201702/RNAseq_permutation.test_stm.RData")
rna.permout <- load("~/Downloads/RNAseq_permutation.test_stm.RData")
# Observed FC
row.names(adj.rna) <- adj.rna$GeneName
adj.rna <- adj.rna[which(row.names(adj.rna) %in% test.table.tcell$GeneName), ]
rna.obsfc <- data.frame(logFC.rna.obs=apply(adj.rna[,samps_pbmc[samps %in% gsub("_PBMC", "", colnames(adj.rna))]],1,function(x) diff(aggregate(x, list(group=group), mean)[,2])), row.names = row.names(adj.rna[,samps_pbmc[samps %in% gsub("_PBMC", "", colnames(adj.rna))]]))

# Permuted FC
row.names(rna.obsfc) <- row.names(adj.rna)
rna.permfc.mx <- do.call(cbind, rna.permfc)
rna.permfc.mx <- data.frame(rna.permfc.mx)[row.names(rna.obsfc) %in% row.names(data.frame(rna.permfc.mx)),]
rna.permfc.mx <- rna.permfc.mx[row.names(rna.permfc.mx) %in% test.table.tcell$GeneName, ]
temp <- rna.obsfc[row.names(rna.obsfc) %in% row.names(data.frame(rna.permfc.mx)), ]

rna.obsfc <- data.frame(logFC.rna.obs = temp, row.names = row.names(rna.obsfc)[rna.obsfc$logFC.rna.obs %in% temp])

rna.obsfc$logFC.rna.pcounts[rna.obsfc$logFC.rna.obs >= 0] <- rowSums(rna.obsfc$logFC.rna.obs[rna.obsfc$logFC.rna.obs >= 0] >= rna.permfc.mx[which(rna.obsfc$logFC.rna.obs >= 0),])
rna.obsfc$logFC.rna.pcounts[rna.obsfc$logFC.rna.obs < 0] <- rowSums(rna.obsfc$logFC.rna.obs[rna.obsfc$logFC.rna.obs < 0] < rna.permfc.mx[rna.obsfc$logFC.rna.obs < 0,])
rna.obsfc$logFC.rna.pvalue <- 1-(rna.obsfc$logFC.rna.pcounts/(1+ncol(rna.permfc.mx)))
rna.obsfc$logFC.rna.fdr <- p.adjust(rna.obsfc$logFC.rna.pvalue, method = "fdr")

rm(list = c(rna.permout,"rna.permfc.mx"))
#### ATACseq TSS peaks FC
### Code used to generate perms:
# atac.permlist <- replicate(nperm,adj.atac.tss[,samps_pbmc],simplify = F)
# atac.permfc <- mclapply(atac.permlist, function(Y) apply(Y,1,function(x) diff(aggregate(x, list(group=sample(group)), mean)[,2])),mc.cores = 6)
###atac.permfc <- lapply(atac.permlist, function(Y) apply(Y,1,function(x) diff(aggregate(x, list(group=sample(group)), mean)[,2])))
# save(list = c("atac.permlist","atac.permfc"),file = "datasources.201702/ATACseq_permutation.test_stm.RData")
atac.permout <- load("~/Downloads/ATACseq_permutation.test_stm.RData")
# Observed FC

group <- relevel(factor(sampinfo[samps_pbmc[samps_pbmc %in% colnames(adj.atac.tss)],]$group),ref = "HY")

atac.obsfc <- data.frame(adj.atac.tss[,c("chr","start","end")],logFC.atac.obs=apply(adj.atac.tss[,colnames(adj.atac.tss) %in% samps_pbmc],1,function(x) diff(aggregate(x, list(group=group), mean)[,2])))
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
atac.obsfc <- atac.obsfc[which(!is.na(atac.obsfc$logFC.atac.pvalue)), ]
rm(list = c(atac.permout,"atac.permfc.mx"))
#####

# Combines ATAC and RNAseq data/results, assigns Pvalues based on permutation tests
colnames(full.results.tss) <- gsub("_PBMC", "", colnames(full.results.tss))
atacrna.glm <- merge(full.results.tss,rna.glmtop, by = "GeneName",suffixes = c(".atac",".rna"))
atacrna.glm <- merge(merge(atacrna.glm, atac.obsfc, by = c("chr","start","end")),
                     rna.obsfc[row.names(rna.obsfc) %in% atacrna.glm$GeneName, ], by.x = "GeneName", by.y = "row.names")

pdf("~/Desktop/Fig3A_pbmc.arc.fcfcplot.pdf",paper = "USr")
fcfcplot <- atacrna.glm # atacrna.glm|atacrna.glm.nodup|atacrna.comp.glm|atacrna.comp.glm.nodup :: if use "comp", variables below have to be modified
fcfc.corr <- with(fcfcplot,cor.test(logFC.atac.obs,logFC.rna.obs, method = "pearson"))
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
ylims[2:3] <- c(0,0)
# xlims[2:3] <- c(0,0)
# genesel <- c("LCN2","PDGFRB","KIR3DL1","PDGFD","GNLY","GZMB","KLRF1","CD248","NOG","CCR7","IL7R","LEF1","BACH2","CD8A","LRRN3","GSTT1","NRCAM","SLC16A10","NBEA","GSTM1","MTUS1","B3GAT1","IGFBP3","NCR1","FGFBP2","GZMH","PRSS23","NT5E","S100B") # highlighted in ms.
genesel <- test.table.tcell$GeneName
fcfcplot.sel <- subset(fcfcplot,with(fcfcplot,GeneName %in% genesel))
fcfcplot.sel <- fcfcplot.sel[order(fcfcplot.sel$GeneName,fcfcplot.sel$logFC.atac.obs),]
fcfcplot.sel <- fcfcplot.sel[!duplicated(fcfcplot.sel$GeneName),]
p <- ggplot(data=fcfcplot[order(fcfcplot$logFC.atac.obs),], aes(logFC.atac.obs, logFC.rna.obs)) +
  annotate("rect", xmin = -max(abs(xlims*1.1)), xmax = xlims[3], ymin = -max(abs(ylims*1.1)), ymax = ylims[3], fill = "darkcyan", alpha = 0.1) + # alt version with differently-colored quadrants. Use "darkorange1" to match both quadrants to the original version
  annotate("rect", xmin = xlims[3], xmax = max(abs(xlims)), ymin = 0, ymax = max(abs(ylims)), fill = "darkorange1", alpha = 0.1) +
  geom_point(aes(size=factor(apply(cbind(logFC.atac.obs<xlims[3]&logFC.rna.obs<ylims[3],logFC.atac.obs>0&logFC.rna.obs>0),1,orred)),
                 alpha=factor(apply(cbind(logFC.atac.obs<xlims[3]&logFC.rna.obs<ylims[3],logFC.atac.obs>0&logFC.rna.obs>0),1,orred)),
                 color=factor(apply(cbind(logFC.atac.obs<xlims[3]&logFC.rna.obs<ylims[3],logFC.atac.obs>0&logFC.rna.obs>0),1,orred)))) +
  geom_text(data = fcfcplot[fcfcplot$logFC.atac.obs*fcfcplot$logFC.rna.obs > 0,], aes(logFC.atac.obs,logFC.rna.obs,label=GeneName), size=5) + # uses a threshhold to apply labels to extreme points
  geom_text(data = subset(fcfcplot.sel,logFC.atac.obs<0),aes(logFC.atac.obs,logFC.rna.obs,label=GeneName), size=5, vjust = 1.3, hjust = 1, color = "green4") + # preselect labels used in the paper (negative)
  geom_text(data = subset(fcfcplot.sel,logFC.atac.obs>0),aes(logFC.atac.obs,logFC.rna.obs,label=GeneName), size=5, vjust = -0.3, hjust = 0, color = "green4") + # preselect labels used in the paper (positive)
  geom_point(data = fcfcplot.sel,aes(logFC.atac.obs,logFC.rna.obs), color = "green3", size = 0.8) +
  scale_color_manual(values = c("FALSE"="dimgray","TRUE"="navy"), guide=F) +
  scale_size_manual(values = c("FALSE"=0.5,"TRUE"=0.7), guide=F) +
  scale_alpha_manual(values = c("FALSE"=0.25,"TRUE"=0.75), guide=F) +
  scale_x_continuous(limits = c(-max(abs(xlims*1.1)),max(abs(xlims)))) +
  scale_y_continuous(limits = c(-max(abs(ylims*1.1)),max(abs(ylims)))) +
  xlab("ATAC-seq logFC (Promoter accessibility)") +
  ylab("RNA-seq logFC (Expression)") +
  ggtitle(paste("logFC RNAseq vs ATACseq")) +
  theme_bw(base_family = "Helvetica",base_size = 10) +
  theme(panel.background = element_rect(fill = "white", color = "black"),
        # panel.grid.major.x = element_line(color = "ivory3", size = 0.2, linetype = 2),
        # panel.grid.major.y = element_line(color = "ivory3", size = 0.2, linetype = 2),
        panel.grid.major = element_blank(),
        axis.line = element_line(color = "black", size = 0.25),
        axis.text.x = element_text(size=10)) +
  geom_hline(yintercept = ylims[3]:0, color = "palegreen3", linetype = "22", size=0.5, alpha=0.5) +
  geom_vline(xintercept = xlims[3]:0, color = "palegreen3", linetype = "22", size=0.5, alpha=0.5) +
# p <- ggMarginal(p,type = "histogram", size = 15, bins = 150, fill="palegreen3", color="seagreen3")
print(p)


fcfcplot <- atacrna.glm # atacrna.glm|atacrna.glm.nodup|atacrna.comp.glm|atacrna.comp.glm.nodup :: if use "comp", variables below have to be modified
fcfc.corr <- with(fcfcplot,cor.test(logFC.atac.obs,logFC.rna.obs, method = "pearson"))
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
ylims[2:3] <- c(0,0)
# genesel <- c("LCN2","PDGFRB","KIR3DL1","PDGFD","GNLY","GZMB","KLRF1","CD248","NOG","CCR7","IL7R","LEF1","BACH2","CD8A","LRRN3","GSTT1","NRCAM","SLC16A10","NBEA","GSTM1","MTUS1","B3GAT1","IGFBP3","NCR1","FGFBP2","GZMH","PRSS23","NT5E","S100B") # highlighted in ms.
genesel <- c("LCN2","PDGFRB","KIR3DL1","PDGFD","GNLY","GZMB","KLRF1","CD248","NOG","CCR7","IL7R","LEF1","BACH2","CD8A","LRRN3","GSTT1","NRCAM","SLC16A10","NBEA","GSTM1","MTUS1","B3GAT1","FGFBP2","GZMH","PRSS23","S100B","PRF1","NT5E") # highlighted in ms.
fcfcplot.sel <- subset(fcfcplot,with(fcfcplot,GeneName %in% genesel & apply(cbind(logFC.atac.obs<xlims[3]&logFC.rna.obs<ylims[3],logFC.atac.obs>xlims[4]&logFC.rna.obs>ylims[4]),1,any)))
fcfcplot.sel <- fcfcplot.sel[order(fcfcplot.sel$GeneName,fcfcplot.sel$logFC.atac.obs),]
fcfcplot.sel <- fcfcplot.sel[!duplicated(fcfcplot.sel$GeneName),]
p <- ggplot(data=fcfcplot[order(fcfcplot$logFC.atac.obs),], aes(logFC.atac.obs, logFC.rna.obs)) +
  annotate("rect", xmin = -max(abs(xlims*1.1)), xmax = xlims[3], ymin = -max(abs(ylims*1.1)), ymax = ylims[3], fill = "darkcyan", alpha = 0.1) + # alt version with differently-colored quadrants. Use "darkorange1" to match both quadrants to the original version
  annotate("rect", xmin = xlims[4], xmax = max(abs(xlims)), ymin = ylims[4], ymax = max(abs(ylims)), fill = "darkorange1", alpha = 0.1) +
  geom_point(aes(size=factor(apply(cbind(logFC.atac.obs<xlims[3]&logFC.rna.obs<ylims[3],logFC.atac.obs>xlims[4]&logFC.rna.obs>ylims[4]),1,orred)),
                 alpha=factor(apply(cbind(logFC.atac.obs<xlims[3]&logFC.rna.obs<ylims[3],logFC.atac.obs>xlims[4]&logFC.rna.obs>ylims[4]),1,orred)),
                 color=factor(apply(cbind(logFC.atac.obs<xlims[3]&logFC.rna.obs<ylims[3],logFC.atac.obs>xlims[4]&logFC.rna.obs>ylims[4]),1,orred)))) +
  # geom_text(data = fcfcplot[fcfcplot$logFC.atac.obs*fcfcplot$logFC.rna.obs > 0.2,], aes(logFC.atac.obs,logFC.rna.obs,label=GeneName), size=2.5) + # uses a threshhold to apply labels to extreme points
  geom_text(data = subset(fcfcplot.sel,logFC.atac.obs<0),aes(logFC.atac.obs,logFC.rna.obs,label=GeneName), size=2.5, vjust = 1.3, hjust = 1, color = "green4") + # preselect labels used in the paper (negative)
  geom_text(data = subset(fcfcplot.sel,logFC.atac.obs>0),aes(logFC.atac.obs,logFC.rna.obs,label=GeneName), size=2.5, vjust = -0.3, hjust = 0, color = "green4") + # preselect labels used in the paper (positive)
  geom_point(data = fcfcplot.sel,aes(logFC.atac.obs,logFC.rna.obs), color = "green3", size = 0.8) +
  scale_color_manual(values = c("FALSE"="dimgray","TRUE"="navy"), guide=F) +
  scale_size_manual(values = c("FALSE"=0.5,"TRUE"=0.7), guide=F) +
  scale_alpha_manual(values = c("FALSE"=0.25,"TRUE"=0.75), guide=F) +
  scale_x_continuous(limits = c(-max(abs(xlims*1.1)),max(abs(xlims)))) +
  scale_y_continuous(limits = c(-max(abs(ylims*1.1)),max(abs(ylims)))) +
  xlab("ATAC-seq logFC (Promoter accessibility)") +
  ylab("RNA-seq logFC (Expression)") +
  ggtitle(paste("logFC RNAseq vs ATACseq")) +
  theme_bw(base_family = "Helvetica",base_size = 10) +
  theme(panel.background = element_rect(fill = "white", color = "black"),
        # panel.grid.major.x = element_line(color = "ivory3", size = 0.2, linetype = 2),
        # panel.grid.major.y = element_line(color = "ivory3", size = 0.2, linetype = 2),
        panel.grid.major = element_blank(),
        axis.line = element_line(color = "black", size = 0.25),
        axis.text.x = element_text(size=10)) +
  geom_hline(yintercept = ylims[3:4], color = "palegreen3", linetype = "22", size=0.5, alpha=0.5) +
  geom_vline(xintercept = xlims[3:4], color = "palegreen3", linetype = "22", size=0.5, alpha=0.5)
# p <- ggMarginal(p,type = "histogram", size = 15, bins = 150, fill="palegreen3", color="seagreen3")
print(p)

dev.off()
