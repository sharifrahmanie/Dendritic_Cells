require(tidyverse)
require(gridExtra)
require(caret)
require(umap)
require(DESeq2)
require(Seurat)
require(org.Mm.eg.db)
require(gage)
require(gageData)


data_importing <- function(path) {
  path1 <- paste0(path, "/droplets/")
  anno_dropseq <- read_csv(paste0(path1, "annotations_droplets.csv"))
  anno_dropseq <- anno_dropseq[which(anno_dropseq$cell_ontology_class == "dendritic cell"),1:3]
  anno_dropseq$cell <- gsub("10X_P[0-9]_[0-9]*_", "", anno_dropseq$cell)
  anno_dropseq_lung <- anno_dropseq[which(anno_dropseq$tissue == "Lung"),]
  anno_dropseq_spleen <- anno_dropseq[which(anno_dropseq$tissue == "Spleen"),]
  lung_dropseq <- Read10X(paste0(path1, "Lung-10X_P7_9/"))
  colnames(lung_dropseq) <- gsub("-1", "", colnames(lung_dropseq))
  spleen_dropseq <- Read10X(paste0(path1, "Spleen-10X_P7_6/"))
  colnames(spleen_dropseq) <- gsub("-1", "", colnames(spleen_dropseq))
  lung_cells <- which(anno_dropseq_lung$cell %in% colnames(lung_dropseq))
  lung_dropseq <- lung_dropseq[, lung_cells]
  spleen_cells <- which(anno_dropseq_spleen$cell %in% colnames(spleen_dropseq))
  spleen_dropseq <- spleen_dropseq[,spleen_cells]
  
  path2 <- paste0(path, "/FACS/")
  anno_facs <- read_csv(paste0(path2,"annotations_FACS.csv"))
  anno_facs <- anno_facs[which(anno_facs$cell_ontology_class == "dendritic cell"),1:3]
  lung_facs <- data.frame(read_csv(paste0(path2,"Lung-counts.csv")))
  dend_cols <- which(colnames(lung_facs) %in% anno_facs$cell)
  lung_facs <- lung_facs[,c(1, dend_cols)]
  rownames(lung_facs) <- lung_facs[,1]
  lung_facs <- lung_facs[,-1]
  counts <- cbind(lung_dropseq, spleen_dropseq, lung_facs)
  groups <- data.frame(groups = rep(c("lung_drop", "spleen_drop", "lung_facs"),
                                    times = c(ncol(lung_dropseq),
                                              ncol(spleen_dropseq),
                                              ncol(lung_facs))))
  return(list(counts = counts, groups = groups))
}

qc_check <- function(counts) {
  dir.create("plots", showWarnings = F)
  qc_data <- cbind(rownames(counts), counts)
  colnames(qc_data)[1] <- "Gene"
  long_data <- gather(qc_data,
                      key = "sample",
                      value = "expression",
                      -Gene)
  gene_counts <- long_data %>% 
    group_by(sample) %>% 
    summarize(num_genes = sum(expression > 0))
  cell_counts <- long_data %>% 
    group_by(Gene) %>% 
    summarize(num_cells = sum(expression > 0))
  g1 <- ggplot(gene_counts, aes(x = "", y = num_genes)) +
    geom_violin(fill= "#873600") +
    geom_jitter()+
    xlab("") +
    ylab("Number of Genes") +
    mytheme() +
    theme_classic()
  g2 <- ggplot(cell_counts, aes(x = "", y = num_cells)) +
    geom_violin(fill = "#873600") +
    geom_jitter()+
    xlab("") +
    ylab("Number of Cells") +
    mytheme() +
    theme_classic()
  pdf("plots/QC.pdf", width = 13, height = 8)
  grid.arrange(g1, g2,nrow = 1)
  dev.off()
}

filtering <- function(counts, min.cell, min.feature) {
  counts <- counts[rowSums(counts >= 1) >= min.cell, ]
  counts <- counts[, colSums(counts >= 1) >= min.feature]
  return(counts)
}

normalize_log <- function(counts, scale.factor) {
  library_sizes <- colSums(counts)
  counts <- t(t(counts) / library_sizes)
  counts <- counts * scale.factor
  counts <- counts + 1
  counts <- log2(counts)
}

mytheme <- function() {
  mytheme <- theme(axis.text = element_text(family = "Times",size = 24 , colour = "black"),
                   axis.text.x = element_text(family = "Times",colour = "black", size = 14),
                   axis.text.y = element_text(family = "Times",colour = "black", size = 14),
                   plot.subtitle = element_text(family = "Times",size = 24, colour = "black", hjust = 0.5),
                   axis.title.y = element_text(family = "Times", size = 16, angle = 90),
                   axis.title.x = element_text(family = "Times", size = 16, angle = 00),
                   legend.text = element_text(size = 16, family = "Times"), 
                   legend.title = element_text(size = 16, family = "Times"))
  return(mytheme)
}

the_elbow <- function(counts) {
  dir.create("plots", showWarnings = F)
  pca_df <- scale(data.frame(t(counts)))
  pca <- preProcess(x = pca_df,
                    method = 'pca',
                    pcaComp = 2)
  pca_df <- data.frame(predict(pca, pca_df))
  wcss <- vector()
  for(i in 1:10){
    wcss[i] <- sum(kmeans(pca_df, i)$withinss)
  }
  g <- ggplot() + 
    geom_line(aes(x = 1:10, y = wcss),
              color = "#DC7633", linewidth = 1) +
    geom_point(aes(x = 1:10, y = wcss), color = "#7E5109") +
    xlab('Number of clusters') + ylab('WCSS') + 
    scale_x_continuous(breaks = seq(1, 10, 1))+
    theme_classic() + mytheme() +
    labs(subtitle = 'The Elbow Method')
  pdf("plots/TheElbow.pdf", width = 10, height = 8)
  print(g)
  dev.off()
}

umap_plot <- function(counts,
                      n.clusters,
                      n.neighbors) {
  clusters <- kmeans(t(counts),
                     n.clusters)$cluster %>%
    as.factor()
  umap_result <- umap(t(counts),
                      n_neighbors = n.neighbors,
                      labels = clusters)
  umap_df <- data.frame(umap_result$layout,
                        sample = colnames(counts),
                        Cluster = clusters)
  g <- ggplot(umap_df, aes(x = X1, y = X2, color = Cluster)) +
    geom_point(size = 2) +
    xlab('UMAP1') + ylab('UMAP2') +
    theme_classic() + mytheme() +
    labs(subtitle = 'UMAP plot')
  pdf("plots/umap.pdf", width = 10, height = 8)
  print(g)
  dev.off()
}

deseq2_results <- function(counts,
                           groups,
                           reference,
                           contrasts,
                           FC,
                           padj) {
  dir.create("DESeq2", showWarnings = F)
  groups$groups <- factor(groups$groups)
  counts1 <- counts + 1
  deseq <- DESeqDataSetFromMatrix(countData = counts1, 
                                  colData = groups, 
                                  design = ~ groups)
  deseq$groups <- relevel(deseq$groups, ref = reference)
  deseq <- DESeq(deseq,
                 fitType='local',
                 parallel = TRUE)
  res <- results(deseq, contrast = paste0(c("groups", contrasts)))
  all <- data.frame(res)
  res <- res[which(abs(res$log2FoldChange) > FC & res$padj < padj), ]
  up <- data.frame(res[which(res$log2FoldChange > FC ),])
  up <- up[order(up$log2FoldChange, decreasing = T),]
  down <- data.frame(res[which(!res$log2FoldChange > FC ),])
  down <- down[order(down$log2FoldChange, decreasing = F),]
  both <- rbind(up, down)
  message(cat("Up-regulated: ", nrow(up), ", Down-regulated: ", nrow(down)))
  message(cat("Top 10 Up-regulated and Down-regulated genes"))
  write.csv(all, paste0("DESeq2/all_", paste0(contrasts,  collapse = "_"), ".csv"))
  write.csv(both, paste0("DESeq2/UPDown_", paste0(contrasts,  collapse = "_"), ".csv"))
  return(list(UP= head(up, 10), DOWN= head(down, 10)))
}

volcano_plot <- function(compare,
                         padjlevel,
                         Up_FC,
                         Down_FC) {
  dir.create("plots", showWarnings = F)
  volcano <- data.frame(read_csv(paste0("DESeq2/all_", 
                                        paste0(compare, collapse = "_") , ".csv")))
  volcano <- na.omit(volcano)
  colnames(volcano)[1] <- "ID"
  volcano <- volcano %>%
    mutate(log10padj = -log10(padj)) %>%
    mutate(DEG = "NotSignificant") %>%
    mutate(DEG = ifelse(log10padj > -log10(padjlevel) & log2FoldChange > Up_FC, "Upregulated", DEG)) %>%
    mutate(DEG = ifelse(log10padj > -log10(padjlevel) & log2FoldChange < Down_FC, "Downregulated", DEG))
  my_pal <- c("#1B9E77","#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666", "#9A7D0A")
  g <- ggplot(data = volcano, aes(x = log2FoldChange, y= log10padj, color= DEG, fill = DEG, label= ID)) + 
    labs(x= 'log2(Fold Change)', y= "-log10(P.adj)") +
    geom_point(size = 2.5, shape = 21) +
    geom_text(check_overlap = T,vjust = -0.3, nudge_y = 0.1) +
    scale_color_manual(values=c(my_pal)) +
    scale_fill_manual(values=my_pal) +
    theme_classic() + mytheme() +
    labs(subtitle = paste0(compare, collapse = "-"))
  pdf(paste0("plots/Volcano_", paste0(compare, collapse = "_"), ".pdf"), width = 10, height = 8)
  print(g)
  dev.off() 
}

pathway_plot <- function(compare,
                         FDR) {
  dir.create("plots", showWarnings = F)
  up_down <- data.frame(read_csv(paste0("DESeq2/UPDown_", 
                                        paste0(compare, collapse = "_") , ".csv")))
  colnames(up_down)[1] <- "genes"
  up_down$entrez <- mapIds(org.Mm.eg.db,
                           keys=up_down$genes,
                           column="ENTREZID",
                           keytype="SYMBOL",
                           multiVals="first")
  data(kegg.sets.mm)
  foldchanges <- up_down$log2FoldChange
  names(foldchanges) <- up_down$entrez
  kegg_up_down <- gage(foldchanges, gsets = kegg.sets.mm, same.dir = TRUE)
  pathway_up <- data.frame(kegg_up_down$greater)
  pathway_up$ID <- str_extract(rownames(pathway_up), "mmu[0-9]+(?= )")
  pathway_up$Pathway <- gsub("mmu[0-9]+ ", "", rownames(pathway_up))
  pathway_up <- pathway_up[order(pathway_up$q.val, decreasing = FALSE), c(7, 8, 5, 4)]
  colnames(pathway_up) <- c("ID", "Pathway", "n_Gene", "FDR")
  rownames(pathway_up) <- NULL
  pathway_up <- pathway_up[which(pathway_up$FDR < FDR), ]

  pathway_down <- data.frame(kegg_up_down$less)
  pathway_down$ID <- str_extract(rownames(pathway_down), "mmu[0-9]+(?= )")
  pathway_down$Pathway <- gsub("mmu[0-9]+ ", "", rownames(pathway_down))
  pathway_down <- pathway_down[order(pathway_down$q.val, decreasing = FALSE), c(7, 8, 5, 4)]
  colnames(pathway_down) <- c("ID", "Pathway", "n_Gene", "FDR")
  rownames(pathway_down) <- NULL
  pathway_down <- pathway_down[which(pathway_down$FDR < FDR), ]
  
  pathways <- rbind(pathway_up, pathway_down)
  pathways <- na.omit(pathways)
  g <- pathways %>%
    ggplot(aes(reorder(Pathway, n_Gene, sum), n_Gene)) +
    geom_col(aes(color = FDR, fill = FDR), width = 0.4) +
    coord_flip() + xlab("Pathway") + ylab("Enriched genes") +
    theme_classic() +
    scale_color_gradient(low = "#D7DBDD", high = "#186A3B") +
    scale_fill_gradient(low = "#D7DBDD", high = "#186A3B") +
    mytheme() +
    labs(subtitle = "")
  pdf(paste0("plots/Pathway_", paste0(compare, collapse = "_"), ".pdf"), width = 10, height = 8)
  print(g)
  dev.off() 
}


set.seed(1234)
# Importing data
mydata <- data_importing(path = "data")
counts <- mydata$counts
groups <- mydata$groups

# QC check 
qc_check(counts = counts)

# Filtering based on a minimum genes detected in all cells
# and minimum number of genes in a cell
counts <- filtering(counts = counts,
                    min.cell = 10,
                    min.feature = 200)

# Library size normalization and log2 transformation
norm_counts <- normalize_log(counts = counts,
                             scale.factor = 10000)

# The Elbow method for detecting the number of clusters
the_elbow(counts = norm_counts)

# UMAP plot 
umap_plot(counts = norm_counts,
          n.clusters = 3, 
          n.neighbors = 7)

# DEG result of different comparisons
deseq2_results(counts = counts,
              groups = groups,
              reference = "lung_drop", 
              contrasts = c("lung_facs", "spleen_drop"),
              FC = 1, 
              padj = 0.05)
deseq2_results(counts = counts,
              groups = groups,
              reference = "lung_drop", 
              contrasts = c("lung_drop", "spleen_drop"),
              FC = 1, 
              padj = 0.05)

deseq2_results(counts = counts,
              groups = groups,
              reference = "lung_drop", 
              contrasts = c("lung_facs", "lung_drop"),
              FC = 1, 
              padj = 0.05)

# Volcano of different comparisons
volcano_plot(compare = c("lung_facs", "spleen_drop"), 
             padjlevel = 0.05,
             Up_FC = 1, 
             Down_FC = -1)

volcano_plot(compare = c("lung_drop", "spleen_drop"), 
             padjlevel = 0.05,
             Up_FC = 1, 
             Down_FC = -1)

volcano_plot(compare = c("lung_facs", "lung_drop"), 
             padjlevel = 0.05,
             Up_FC = 1, 
             Down_FC = -1)
# Pathway plot of different comparisons

pathway_plot(compare = c("lung_facs", "spleen_drop"),
             FDR = 0.05)
pathway_plot(compare = c("lung_drop", "spleen_drop"),
             FDR = 0.05)
pathway_plot(compare = c("lung_facs", "lung_drop"),
             FDR = 0.05)

