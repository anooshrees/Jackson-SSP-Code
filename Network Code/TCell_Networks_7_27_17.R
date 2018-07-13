network.data.tcells <- load("~/Desktop/TCell_DiffData/geneset.data_playground.p900_consensus.radius_RNA_explicit_fullnet_redundant.RData")
# path.all.negative <- path.all[path.all$HY_HO.mean_deltaProp < 0, ]
# path.all.negative <- path.all.negative[abs(path.all.negative$HY_HO.mean_deltaProp) > quantile(abs(path.all.negative$HY_HO.mean_deltaProp),prob=1-50/100),]
# 
# path.all.positive <- path.all[path.all$HY_HO.mean_deltaProp > 0, ]
# path.all.positive <- path.all.positive[abs(path.all.positive$HY_HO.mean_deltaProp) > quantile(abs(path.all.positive$HY_HO.mean_deltaProp),prob=1-50/100),]
# 
# path.all <- rbind(path.all.positive, path.all.negative)

for(geneset.name in unique(path.all$schema)){
  # geneset.name <- unique(path.all$schema)[1]
  # since data is presented differently for some of the gene sets, individual debugging may be necessary
  # geneset.name <- "reactome"
  all.info <- path.all[path.all$schema %in% geneset.name, ]
  # 
  # print(paste("Working on geneset", geneset.name))
  # current_geneset <- test.table.tcell$GeneName
  
  for(current.module in unique(all.info$pathway)){
    # current.module <- unique(all.info$pathway)[1]
    # current.module <- "Cellular Senescence"
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
      pdf(paste("~/Desktop/TF_Diff_7_27_17/", geneset.name, "_", gsub("[[:punct:]]", " ", gsub(",.*$", "", current.module)), ".pdf", sep = ""))
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

atac.ref <- read.table(file="~/Desktop/SSP/DA_7_12_17/tcells_merged_consensus_atac_results.txt",sep = "\t")
colnames(atac.ref) <- c("chr", "start", "end", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")
atac.ref <- atac.ref[-1,]

atac.ref <- merge(reference, atac.ref, by=c("chr", "start", "end"))
atac.ref <- atac.ref[atac.ref$GeneName %in% diff.mx$GeneName, ]
atac.ref <- do.call(rbind,lapply(split(atac.ref,atac.ref$GeneName),function(chunk) chunk[which.min(chunk$DistancetoTSS),]))

diff.mx <- merge(diff.mx[diff.mx$GeneName %in% atac.ref$GeneName, ], atac.ref, by="GeneName",suffixes = c(".rna",".atac"))

row.names(diff.mx) <- diff.mx$GeneName
diff.mx <- diff.mx[order(diff.mx$logFC.atac), ]
diff.mx[is.na(diff.mx)] <- 0

newCols <- colorRampPalette(colorspace::diverge_hsv(length(unique(diff.mx$logFC.atac))))
mycolors <- newCols(length(unique(diff.mx$logFC.atac)))
names(mycolors) <- unique(diff.mx$logFC.atac)
mycolors <- list(logFC.atac = mycolors)

newCols2 <- colorRampPalette(colorspace::diverge_hsv(length(unique(diff.mx$logFC.rna))))
mycolors2 <- newCols2(length(unique(diff.mx$logFC.rna)))
names(mycolors2) <- unique(diff.mx$logFC.rna)
mycolors2 <- list(logFC.rna = mycolors2)

pdf("~/Desktop/WTF.pdf")
pheatmap(data.frame(diff.mx[, c(3:8)], row.names = row.names(diff.mx)),show_rownames = T,scale = "row",color = colorRampPalette(cm.colors(5))(100), annotation_row = diff.mx[c("logFC.atac", "logFC.rna")], annotation_colors = c(mycolors, mycolors2), cluster_rows=FALSE, cluster_cols=FALSE, annotation_legend = F)
dev.off()
