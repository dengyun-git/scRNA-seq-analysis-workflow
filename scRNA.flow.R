### This is a TOY script for scRNA seq analysis, which is NOT my original work. I curated approaches both of QC and downstream analysis from different online resources, which could be a good reference for my own potential project.
### Oct 2023

library(Seurat)
library(SoupX)
library(ggplot2)
library(SingleR)
library(dplyr)
library(celldex)
library(RColorBrewer)
library(SingleCellExperiment)
library(tidyverse)
library(magrittr)
library(monocle)
library(data.table)
library(celldex)
library(SingleR)
library(SCENIC)
library(SCopeLoomR)
library(network)
library(igraph)
library(ArchR)

### pay attention to the library attaching order, to avoid package in need to be masked

inputFile <- "/Users/ydeng/Desktop/工作总结/scRNA-seq analysis/input source/"

#------------------
#------------------
# QUALITY CONTROL
#------------------
#------------------
######################################
# ambient RNA correction using soupX
######################################
### download the files
download.file("https://cf.10xgenomics.com/samples/cell-exp/4.0.0/Parent_NGSC3_DI_PBMC/Parent_NGSC3_DI_PBMC_filtered_feature_bc_matrix.h5",
              destfile = paste0(inputFile,"pbmc10k_filt.h5"))
download.file("https://cf.10xgenomics.com/samples/cell-exp/4.0.0/Parent_NGSC3_DI_PBMC/Parent_NGSC3_DI_PBMC_raw_feature_bc_matrix.h5",
              destfile = paste0(inputFile,"pbmc10k_raw.h5"))
### read in HDF5 data 
filt.matrix <- Read10X_h5(paste0(inputFile,"pbmc10k_filt.h5"),use.names = T) ### use.names=F,Ensemble IDs as row identifiers
raw.matrix  <- Read10X_h5(paste0(inputFile,"pbmc10k_raw.h5"),use.names = T)
# str(raw.matrix)

### creat a seruat object
srat  <- CreateSeuratObject(counts = filt.matrix)
View(srat)
soup.channel  <- SoupChannel(raw.matrix, filt.matrix)
soup.channel

### define marker genes
srat    <- SCTransform(srat, verbose = F)
srat    <- RunPCA(srat, verbose = F)
srat    <- RunUMAP(srat, dims = 1:30, verbose = F)
srat    <- FindNeighbors(srat, dims = 1:30, verbose = F)

meta    <- srat@meta.data
umap    <- srat@reductions$umap@cell.embeddings
soup.channel  <- setClusters(soup.channel, setNames(meta$seurat_clusters, rownames(meta)))
soup.channel  <- setDR(soup.channel, umap)

###  calculating ambient RNA profile.
soup.channel  <- autoEstCont(soup.channel)

### This outputs the matrix of soup-corrected reads which we will use for all future analysis.
combined  <- adjustCounts(soup.channel, roundToInt = T)

##################################
# doublet detection using scrublet
##################################
### dimension reduction plot: whether the double cells together
FeaturePlot(combined, features = "DoubletScores", pt.size = 0.01)
DimPlot(combined, group.by = "DoubletPrediction",pt.size = 0.01,cols = c("darkgreen","firebrick"))

### nUMI distribution: whether the double cells together
VlnPlot(combined,features = "nCount_RNA",pt.size = 0,group.by = "DoubletPrediction") + NoLegend()

### Visualise factions of doublet per cluster in stacked bar chart
df <- data.table(combined@meta.data)
sel.meta <- c("DoubletPrediction", "cluster", "Individual")
df <- df[, sel.meta, with = FALSE]

df[, 2:3] %>% map( ~ {
  freq1 <- df[, .N, keyby = .(.x, DoubletPrediction)]
  freq1[, total := sum(N), by = .(.x)]
  freq1[, ratio := N / total]
  
  linesize = .35
  fontsize = 8
  
  ggplot(freq1, aes(fill=DoubletPrediction, y=ratio, x= .x)) + 
    geom_bar(position="stack", stat="identity")+
    scale_fill_manual(values = c("Doublet" = 'red', "Singlet" = "grey")) +
    xlab('Clsuter') +
    scale_y_continuous(breaks = seq(0,1,0.1), expand = c(0,0), name = 'Percentage')+
    theme_bw()+
    theme( panel.grid.major.x = element_blank(), 
           panel.grid.major.y = element_blank(),
           panel.grid.minor = element_blank(),
           strip.background = element_blank(),panel.border = element_rect(size = linesize),
           axis.ticks = element_blank(), 
           axis.text.x = element_text(size = 5))
})

### component clusters for doublets 
### Differential Expressed Genes
cluster.markers <- FindMarkers(combined, ident.1 = c("InMGE"), ident.2 = "InCGE", min.pct = 0.25)
sel.idents <- c("InMGE", "InCGE",  "D33")
combined.small <- subset(combined, cells = WhichCells(combined, idents = sel.idents))
DoHeatmap(combined.small, features = rownames(cluster.markers)[1:40], raster = F)

### canonical gene
sel.feature <- c("NXPH1", "PAM", "LHX6", "NR2F2", "ADARB2",  "PROX1")
FeaturePlot(combined, features = sel.feature,  pt.size = 0.01, ncol = 3)
VlnPlot(combined.small, features = sel.feature, pt.size = 0, ncol = 3, idents = sel.idents)

################################################################
# Normalize, scale, find variable genes and dimension reduction
################################################################
combined <- NormalizeData(combined.small)
combined <- FindVariableFeatures(combined, selection.method = 'vst', nfeatures = 2000)
hvg <- VariableFeatures(combined)
var_regex = '^HLA-|^IG[HJKL]|^RNA|^MT|^RP' # remove HLA, immunoglobulin, RNA, MT, and RP genes based on HUGO gene names
hvg = grep(var_regex, hvg, invert=T, value=T)

combined %<>%
  ScaleData(vars.to.regress = vars.reg) %>%
  RunPCA(features = hvg) %>%
  FindNeighbors(dims = 1:ndims) %>%
  FindClusters(resolution = res) %>%
  RunUMAP(dims = 1:ndims) %>%
  RunTSNE(dims = 1:ndims)

###############################################
# Regressing out sources of unwanted variation
###############################################
# Cell-Cycle Scoring and Regression
###############################################
# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# calculate cell cycle score for each cell
combined <- CellCycleScoring(combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# Visualize the distribution of cell cycle markers across
RidgePlot(combined, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)

### For each gene, Seurat models the relationship between gene expression and the S and G2M cell cycle scores. The scaled residuals of this model represent a ‘corrected’ expression matrix
combined <- ScaleData(combined, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(marrow))

###########################################################
# Exploration technical batch effects and correct for them
###########################################################
### explore in UMAP or PCA, find the unwanted variation brought by technical BatchVariable1, BatchVariable2
combined <- SCTransform(combined,vars.to.regress = c("BatchVariable1", "BatchVariable2"),verbose = FALSE)

#------------------------------------------------------
#------------------------------------------------------
# DOWNSTREAM ANALYSIS BASED ON QUALITY CONTROLLED DATA
#------------------------------------------------------
#------------------------------------------------------
#####################
# Pseudobulk RNA seq 
#####################
plotBulkHeatmap <- function(x){
  combined.tmp <- ScaleData(combined, features = x)
  x <- intersect(x, GetAssayData(combined.tmp, slot = 'scale.data'))
  mat <- AverageExpression(combined.tmp, features = x, slot = 'scale.data')
  mat1 <- mat$RNA
  re <- pheatmap::pheatmap(mat1, angle_col = 45,  border = NA)
  return(re)
}

### heatmap, clusters vs averaging genne expressions
gn <- c('SOX9', 'HOXP', 'MKI67', 'EOMES', 'NEUROD2', 'SATB2', 'NR4A2', 'GAD1', 'GAD2')
plotBulkHeatmap(gn)

###############################
# Gene Differential Expression 
###############################
### different from the algorithms used for Pseudobulk RNA data, for all the clusters defined
mk <- FindMarkers(combined, min.pct = 0.3, min.diff.pct = 0.2, only.pos = T)

####################################
# Functional Annotation (KEGG & GO)
####################################
### Transfer gene symbol into entrez id (matching the database entries we are using)
geneid.ls <- deg.ls %>% map(~{
  gene.df <- select(org.Mmu.eg.db,
                    keys = .x,
                    columns = c("ENTREZID", "SYMBOL"),
                    keytype = "SYMBOL")
  
  gene <- gene.df$ENTREZID
  gene <- gene[which(!is.na(gene))]
  gene <- unique(gene)
  
  return(gene)
})

gene.ls <- geneid.ls[c(1, 2, 8)]

# KEGG and GO pathway enrichment tests
compKEGG <- compareCluster(geneCluster   = gene.ls,
                           fun           = "enrichKEGG",
                           pvalueCutoff  = 0.05,
                           pAdjustMethod = "BH", 
                           organism = "mcc")
dotplot(compKEGG, showCategory = 10, title = "KEGG Pathway Enrichment Analysis")

compGO <- compareCluster(geneCluster   = gene.ls,
                         fun           = "enrichGO",
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH", 
                         OrgDb = org.Mmu.eg.db, 
                         ont = 'BP')
dotplot(compGO, showCategory = 10, title = "GO Enrichment Analysis")

########################
# Cell Type Annotation
########################
### singleR and seurat, both packages have multiple choices for cell type annotation:
### Clustering + cell type marker genes; Machine Learning like Random Forrest; Reference-Based Annotation

# Using existing references "celldex"
# Assigns labels to cells based on the reference samples with the highest Spearman rank correlations, using only the marker genes between pairs of labels to focus on the relevant differences between cell types
ref <- BlueprintEncodeData()
pred <- SingleR(test=srat, ref=ref, labels=ref$label.main)
plotScoreHeatmap(pred) # heatmap of the per-cell and label scores

# SCENIC package: ROC-AUC to compare the accuracy of cell type identification bsaed on different approaches
srat2 <- ClassifyCells(srat, classifier = rf_classifier)
# Access AUC scores
auc_scores <- srat2$AUC

#################################
# Trajectory Inference (Monocle2)
#################################
cds <- as.CellDataSet(combined) ### change to format of cell data set

### Estimate size factor
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

### Ordering cells
deg <- readRDS('data/Demo_CombinedSeurat_SCT_Preprocess_FilterLQCells_DEGPerCluster_Minpct03Mindiffpct02.rds')
deg <- deg[which(deg$cluster %in% unique(combined$cluster)), ]
sel.gene <- unique(deg$gene)
cds <- monocle::setOrderingFilter(cds, sel.gene)

# dimension reduciton
cds <- monocle::reduceDimension(cds, method = 'DDRTree')

# ordering cells
cds <- monocle::orderCells(cds)

# ordering cells by assigning root nodes
GM_state <- function(cds){
  if (length(unique(cds$State)) > 1){
    T0_counts <- table(cds$State, cds$cluster)[,"RG"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}
cds <- monocle::orderCells(cds, root_state =  GM_state(cds))

### Visualization
# color could by "cluster","State", "Pseudotime"
monocle::plot_cell_trajectory(cds, color_by = "cluster") + facet_wrap(~State)

# per gene, relative expression vs pseudo time 
my_genes <- c("HOPX", "MKI67", "EOMES", "NEUROD2", "SATB2")
cds_subset <- cds[my_genes,]
plot_genes_in_pseudotime(cds_subset, color_by = "cluster")

#####################################
# Gene regulatory networks (SCENIC)
#####################################
# co-expression modules & transcriptioin target database,  such as cisTarget or TRRUST
exprMat <- combined@assays$RNA@data
scenic_results <- run_SCENIC(exprMat, tf_database, cisTarget)
plot_SCENIC(scenic_results)

###########################################
# Motif Deviations and Footprinting (ArchR)
###########################################
addArchRThreads(threads = 16) # parallel threads
proj <- addBgdPeaks(proj) # add bg peak
proj <- addDeviationsMatrix(
  ArchRProj = proj,
  peakAnnotation = "Motif",
  force = TRUE,
  binarize = TRUE
) # add deviations scores

plotVarDev <- getVarDeviations(proj, name = "MotifMatrix", plot = TRUE)
motifs <- c("PAX6", "EOMES", "NEUROD2", "DLX2")
markerMotifs <- getFeatures(proj, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
markerMotifs 

### Footprint
seFoot <- getFootprints(
  ArchRProj = proj, 
  positions = motifPositions[markerMotifs], 
  groupBy = "predictedGroup_Un"
)

plotFootprints(
  seFoot = seFoot,
  ArchRProj = proj, 
  normMethod = "Subtract",
  plotName = "Footprints-Subtract-Bias",
  addDOC = FALSE,
  smoothWindow = 5
)
