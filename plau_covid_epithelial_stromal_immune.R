### scRNA-seq analysis
###
### lung parenchima samples from control individuals,
### individuals with lung fibrosis and individulas diseased from severe Covid19
###
### The datasets were published in 2021
### in DOI: 10.1126/scitranslmed.abe4282 (COVID-19 samples) and
### DOI: 10.1126/sciadv.aba1972 (fibrosis samples)
### Processed matrices were downloaded from GEO: GSE158127
### 
### Human cell cycle markers were downloaded from
### https://hbctraining.github.io/scRNA-seq/lessons/06_SC_SCT_and_integration.html
###
### Anna KARPUKHINA, CNRS UMR9018, Institut Gustave Roussy
### 2022

library(Seurat)
library(SeuratDisk)
library(Matrix)
library(ggplot2)
library(tidyverse)
library(SingleCellExperiment)
library(scales)
library(cowplot)
library(RCurl)
library(DESeq2)
library(pheatmap)
library(dplyr)
library(scater)
library(purrr)
library(Matrix.utils)
library(edgeR)
library(magrittr)
library(reshape2)
library(S4Vectors)
library(tibble)
library(apeglm)
library(ggsignif)
library(ggpubr)

################################################################################
## Epithelial ##
data <- Convert("data/GSE158127_01epithelial.h5ad", ".h5seurat")
epithelial <- LoadH5Seurat(data)
epithelial@meta.data$condition <- NA
epithelial@meta.data$condition[which(str_detect(epithelial@meta.data$Diagnosis,"^Control"))] <- "control"
epithelial@meta.data$condition[which(str_detect(epithelial@meta.data$Diagnosis,"^COVID"))] <- "COVID-19"
epithelial@meta.data$condition[which(str_detect(epithelial@meta.data$Diagnosis,"^IPF"))] <- "IPF/Other PF"
epithelial@meta.data$condition[which(str_detect(epithelial@meta.data$Diagnosis,"^Other"))] <- "IPF/Other PF"

################################################################################
### Figure 1 ###
DimPlot(epithelial,reduction = "umap", group.by = "Cluster", split.by = "condition")
DimPlot(epithelial, reduction = "umap", group.by = "Cell.Type")

################################################################################
### Figure 2 ###
FeaturePlot(epithelial,reduction = "umap", split.by = "condition",
            features = c("PLAU", "PLAUR"))
################################################################################
### Figure 3 ###
VlnPlot(epithelial, split.by = "condition", features = c("PLAU", "PLAUR"))
VlnPlot(epithelial, features = c("PLAU", "PLAUR"), group.by = "Cluster", split.by = "condition")
VlnPlot(epithelial, features = c("PLAU", "PLAUR"), split.by = "Cluster", group.by = "condition")

test <- FetchData(epithelial, vars = c("PLAU", "condition", "Cluster"))
compare_means(PLAU ~ condition, data = test[test$Cluster == "KRT17+ KRT5-", ], method = "t.test", alternative="greater", p.adjust.method = "fdr")
compare_means(PLAU ~ condition, data = test[test$Cluster == "SCGB3A2+", ], method = "t.test", alternative="greater", p.adjust.method = "fdr")
compare_means(PLAU ~ condition, data = test[test$Cluster == "AT2", ], method = "t.test", alternative="greater", p.adjust.method = "fdr")
compare_means(PLAU ~ condition, data = test[test$Cluster == "Basal", ], method = "t.test", alternative="less", p.adjust.method = "fdr")
compare_means(PLAU ~ condition, data = test[test$Cluster == "Ciliated", ], method = "t.test", alternative="greater", p.adjust.method = "fdr")
compare_means(PLAU ~ condition, data = test[test$Cluster == "Club", ], method = "t.test", alternative="greater", p.adjust.method = "fdr")
compare_means(PLAU ~ condition, data = test[test$Cluster == "Differentiating ciliated", ], method = "t.test", alternative="greater", p.adjust.method = "fdr")

summary(test[test$Cluster == "KRT17+ KRT5-", ][test[test$Cluster == "KRT17+ KRT5-", ]$condition == "control",]$PLAU)
summary(test[test$Cluster == "KRT17+ KRT5-", ][test[test$Cluster == "KRT17+ KRT5-", ]$condition == "COVID-19",]$PLAU)
summary(test[test$Cluster == "KRT17+ KRT5-", ][test[test$Cluster == "KRT17+ KRT5-", ]$condition == "IPF/Other PF",]$PLAU)

summary(test[test$Cluster == "SCGB3A2+", ][test[test$Cluster == "SCGB3A2+", ]$condition == "control",]$PLAU)
summary(test[test$Cluster == "SCGB3A2+", ][test[test$Cluster == "SCGB3A2+", ]$condition == "COVID-19",]$PLAU)
summary(test[test$Cluster == "SCGB3A2+", ][test[test$Cluster == "SCGB3A2+", ]$condition == "IPF/Other PF",]$PLAU)

summary(test[test$Cluster == "AT2", ][test[test$Cluster == "AT2", ]$condition == "control",]$PLAU)
summary(test[test$Cluster == "AT2", ][test[test$Cluster == "AT2", ]$condition == "COVID-19",]$PLAU)
summary(test[test$Cluster == "AT2", ][test[test$Cluster == "AT2", ]$condition == "IPF/Other PF",]$PLAU)

summary(test[test$Cluster == "Basal", ][test[test$Cluster == "Basal", ]$condition == "control",]$PLAU)
summary(test[test$Cluster == "Basal", ][test[test$Cluster == "Basal", ]$condition == "COVID-19",]$PLAU)
summary(test[test$Cluster == "Basal", ][test[test$Cluster == "Basal", ]$condition == "IPF/Other PF",]$PLAU)

summary(test[test$Cluster == "Ciliated", ][test[test$Cluster == "Ciliated", ]$condition == "control",]$PLAU)
summary(test[test$Cluster == "Ciliated", ][test[test$Cluster == "Ciliated", ]$condition == "COVID-19",]$PLAU)
summary(test[test$Cluster == "Ciliated", ][test[test$Cluster == "Ciliated", ]$condition == "IPF/Other PF",]$PLAU)

summary(test[test$Cluster == "Club", ][test[test$Cluster == "Club", ]$condition == "control",]$PLAU)
summary(test[test$Cluster == "Club", ][test[test$Cluster == "Club", ]$condition == "COVID-19",]$PLAU)
summary(test[test$Cluster == "Club", ][test[test$Cluster == "Club", ]$condition == "IPF/Other PF",]$PLAU)

summary(test[test$Cluster == "Differentiating ciliated", ][test[test$Cluster == "Differentiating ciliated", ]$condition == "control",]$PLAU)
summary(test[test$Cluster == "Differentiating ciliated", ][test[test$Cluster == "Differentiating ciliated", ]$condition == "COVID-19",]$PLAU)
summary(test[test$Cluster == "Differentiating ciliated", ][test[test$Cluster == "Differentiating ciliated", ]$condition == "IPF/Other PF",]$PLAU)

test <- FetchData(epithelial, vars = c("PLAUR", "condition", "Cluster"))
compare_means(PLAUR ~ condition, data = test[test$Cluster == "KRT17+ KRT5-", ], method = "t.test", alternative="greater", p.adjust.method = "fdr")
compare_means(PLAUR ~ condition, data = test[test$Cluster == "SCGB3A2+", ], method = "t.test", alternative="greater", p.adjust.method = "fdr")
compare_means(PLAUR ~ condition, data = test[test$Cluster == "AT2", ], method = "t.test", alternative="greater", p.adjust.method = "fdr")
compare_means(PLAUR ~ condition, data = test[test$Cluster == "Basal", ], method = "t.test", alternative="less", p.adjust.method = "fdr")
compare_means(PLAUR ~ condition, data = test[test$Cluster == "Ciliated", ], method = "t.test", alternative="greater", p.adjust.method = "fdr")
compare_means(PLAUR ~ condition, data = test[test$Cluster == "Club", ], method = "t.test", alternative="greater", p.adjust.method = "fdr")
compare_means(PLAUR ~ condition, data = test[test$Cluster == "Differentiating ciliated", ], method = "t.test", alternative="greater", p.adjust.method = "fdr")

summary(test[test$Cluster == "KRT17+ KRT5-", ][test[test$Cluster == "KRT17+ KRT5-", ]$condition == "control",]$PLAUR)
summary(test[test$Cluster == "KRT17+ KRT5-", ][test[test$Cluster == "KRT17+ KRT5-", ]$condition == "COVID-19",]$PLAUR)
summary(test[test$Cluster == "KRT17+ KRT5-", ][test[test$Cluster == "KRT17+ KRT5-", ]$condition == "IPF/Other PF",]$PLAUR)

summary(test[test$Cluster == "SCGB3A2+", ][test[test$Cluster == "SCGB3A2+", ]$condition == "control",]$PLAUR)
summary(test[test$Cluster == "SCGB3A2+", ][test[test$Cluster == "SCGB3A2+", ]$condition == "COVID-19",]$PLAUR)
summary(test[test$Cluster == "SCGB3A2+", ][test[test$Cluster == "SCGB3A2+", ]$condition == "IPF/Other PF",]$PLAUR)

summary(test[test$Cluster == "AT2", ][test[test$Cluster == "AT2", ]$condition == "control",]$PLAUR)
summary(test[test$Cluster == "AT2", ][test[test$Cluster == "AT2", ]$condition == "COVID-19",]$PLAUR)
summary(test[test$Cluster == "AT2", ][test[test$Cluster == "AT2", ]$condition == "IPF/Other PF",]$PLAUR)

summary(test[test$Cluster == "Basal", ][test[test$Cluster == "Basal", ]$condition == "control",]$PLAUR)
summary(test[test$Cluster == "Basal", ][test[test$Cluster == "Basal", ]$condition == "COVID-19",]$PLAUR)
summary(test[test$Cluster == "Basal", ][test[test$Cluster == "Basal", ]$condition == "IPF/Other PF",]$PLAUR)

summary(test[test$Cluster == "Ciliated", ][test[test$Cluster == "Ciliated", ]$condition == "control",]$PLAUR)
summary(test[test$Cluster == "Ciliated", ][test[test$Cluster == "Ciliated", ]$condition == "COVID-19",]$PLAUR)
summary(test[test$Cluster == "Ciliated", ][test[test$Cluster == "Ciliated", ]$condition == "IPF/Other PF",]$PLAUR)

summary(test[test$Cluster == "Club", ][test[test$Cluster == "Club", ]$condition == "control",]$PLAUR)
summary(test[test$Cluster == "Club", ][test[test$Cluster == "Club", ]$condition == "COVID-19",]$PLAUR)
summary(test[test$Cluster == "Club", ][test[test$Cluster == "Club", ]$condition == "IPF/Other PF",]$PLAUR)

summary(test[test$Cluster == "Differentiating ciliated", ][test[test$Cluster == "Differentiating ciliated", ]$condition == "control",]$PLAUR)
summary(test[test$Cluster == "Differentiating ciliated", ][test[test$Cluster == "Differentiating ciliated", ]$condition == "COVID-19",]$PLAUR)
summary(test[test$Cluster == "Differentiating ciliated", ][test[test$Cluster == "Differentiating ciliated", ]$condition == "IPF/Other PF",]$PLAUR)

################################################################################
## Stromal ##

data <- Convert("data/GSE158127_03stromal.h5ad", ".h5seurat")
stromal <- LoadH5Seurat(data)
stromal@meta.data$condition <- NA
stromal@meta.data$condition[which(str_detect(stromal@meta.data$Diagnosis,"^Control"))] <- "control"
stromal@meta.data$condition[which(str_detect(stromal@meta.data$Diagnosis,"^COVID"))] <- "COVID-19"
stromal@meta.data$condition[which(str_detect(stromal@meta.data$Diagnosis,"^IPF"))] <- "IPF/Other PF"
stromal@meta.data$condition[which(str_detect(stromal@meta.data$Diagnosis,"^Other"))] <- "IPF/Other PF"

################################################################################
### Figure 1 ###
DimPlot(stromal,reduction = "umap", group.by = "Cluster", split.by = "condition")
DimPlot(stromal, reduction = "umap", group.by = "Cell.Type")

################################################################################
### Figure 2 ###
FeaturePlot(stromal,reduction = "umap", split.by = "condition",
            features = c("PLAU", "PLAUR"))
################################################################################
### Figure 3 ###
VlnPlot(stromal, split.by = "condition", features = c("PLAU", "PLAUR"))
VlnPlot(stromal, features = c("PLAU", "PLAUR"), group.by = "Cluster", split.by = "condition")

test <- FetchData(stromal, vars = c("PLAU", "condition", "Cluster"))
compare_means(PLAU ~ condition, data = test[test$Cluster == "Fibroblasts-1", ], method = "t.test", alternative="greater", p.adjust.method = "fdr")
compare_means(PLAU ~ condition, data = test[test$Cluster == "Fibroblasts-2", ], method = "t.test", alternative="greater", p.adjust.method = "fdr")
compare_means(PLAU ~ condition, data = test[test$Cluster == "Fibroblasts-3", ], method = "t.test", alternative="greater", p.adjust.method = "fdr")
compare_means(PLAU ~ condition, data = test[test$Cluster == "Fibroblasts-4", ], method = "t.test", alternative="greater", p.adjust.method = "fdr")
compare_means(PLAU ~ condition, data = test[test$Cluster == "Mesothelial", ], method = "t.test", alternative="two.sided", p.adjust.method = "fdr")
compare_means(PLAU ~ condition, data = test[test$Cluster == "Smooth Muscle Cells", ], method = "t.test", alternative="greater", p.adjust.method = "fdr")

summary(test[test$Cluster == "Fibroblasts-1", ][test[test$Cluster == "Fibroblasts-1", ]$condition == "control",]$PLAU)
summary(test[test$Cluster == "Fibroblasts-1", ][test[test$Cluster == "Fibroblasts-1", ]$condition == "COVID-19",]$PLAU)
summary(test[test$Cluster == "Fibroblasts-1", ][test[test$Cluster == "Fibroblasts-1", ]$condition == "IPF/Other PF",]$PLAU)

summary(test[test$Cluster == "Fibroblasts-2", ][test[test$Cluster == "Fibroblasts-2", ]$condition == "control",]$PLAU)
summary(test[test$Cluster == "Fibroblasts-2", ][test[test$Cluster == "Fibroblasts-2", ]$condition == "COVID-19",]$PLAU)
summary(test[test$Cluster == "Fibroblasts-2", ][test[test$Cluster == "Fibroblasts-2", ]$condition == "IPF/Other PF",]$PLAU)

summary(test[test$Cluster == "Fibroblasts-3", ][test[test$Cluster == "Fibroblasts-3", ]$condition == "control",]$PLAU)
summary(test[test$Cluster == "Fibroblasts-3", ][test[test$Cluster == "Fibroblasts-3", ]$condition == "COVID-19",]$PLAU)
summary(test[test$Cluster == "Fibroblasts-3", ][test[test$Cluster == "Fibroblasts-3", ]$condition == "IPF/Other PF",]$PLAU)

summary(test[test$Cluster == "Fibroblasts-4", ][test[test$Cluster == "Fibroblasts-4", ]$condition == "control",]$PLAU)
summary(test[test$Cluster == "Fibroblasts-4", ][test[test$Cluster == "Fibroblasts-4", ]$condition == "COVID-19",]$PLAU)
summary(test[test$Cluster == "Fibroblasts-4", ][test[test$Cluster == "Fibroblasts-4", ]$condition == "IPF/Other PF",]$PLAU)

summary(test[test$Cluster == "Mesothelial", ][test[test$Cluster == "Mesothelial", ]$condition == "control",]$PLAU)
summary(test[test$Cluster == "Mesothelial", ][test[test$Cluster == "Mesothelial", ]$condition == "COVID-19",]$PLAU)
summary(test[test$Cluster == "Mesothelial", ][test[test$Cluster == "Mesothelial", ]$condition == "IPF/Other PF",]$PLAU)

summary(test[test$Cluster == "Smooth Muscle Cells", ][test[test$Cluster == "Smooth Muscle Cells", ]$condition == "control",]$PLAU)
summary(test[test$Cluster == "Smooth Muscle Cells", ][test[test$Cluster == "Smooth Muscle Cells", ]$condition == "COVID-19",]$PLAU)
summary(test[test$Cluster == "Smooth Muscle Cells", ][test[test$Cluster == "Smooth Muscle Cells", ]$condition == "IPF/Other PF",]$PLAU)

test <- FetchData(stromal, vars = c("PLAUR", "condition", "Cluster"))
compare_means(PLAUR ~ condition, data = test[test$Cluster == "Fibroblasts-1", ], method = "t.test", alternative="greater", p.adjust.method = "fdr")
compare_means(PLAUR ~ condition, data = test[test$Cluster == "Fibroblasts-2", ], method = "t.test", alternative="greater", p.adjust.method = "fdr")
compare_means(PLAUR ~ condition, data = test[test$Cluster == "Fibroblasts-3", ], method = "t.test", alternative="less", p.adjust.method = "fdr")
compare_means(PLAUR ~ condition, data = test[test$Cluster == "Fibroblasts-4", ], method = "t.test", alternative="greater", p.adjust.method = "fdr")
compare_means(PLAUR ~ condition, data = test[test$Cluster == "Mesothelial", ], method = "t.test", alternative="greater", p.adjust.method = "fdr")
compare_means(PLAUR ~ condition, data = test[test$Cluster == "Smooth Muscle Cells", ], method = "t.test", alternative="greater", p.adjust.method = "fdr")

summary(test[test$Cluster == "Fibroblasts-1", ][test[test$Cluster == "Fibroblasts-1", ]$condition == "control",]$PLAUR)
summary(test[test$Cluster == "Fibroblasts-1", ][test[test$Cluster == "Fibroblasts-1", ]$condition == "COVID-19",]$PLAUR)
summary(test[test$Cluster == "Fibroblasts-1", ][test[test$Cluster == "Fibroblasts-1", ]$condition == "IPF/Other PF",]$PLAUR)

summary(test[test$Cluster == "Fibroblasts-2", ][test[test$Cluster == "Fibroblasts-2", ]$condition == "control",]$PLAUR)
summary(test[test$Cluster == "Fibroblasts-2", ][test[test$Cluster == "Fibroblasts-2", ]$condition == "COVID-19",]$PLAUR)
summary(test[test$Cluster == "Fibroblasts-2", ][test[test$Cluster == "Fibroblasts-2", ]$condition == "IPF/Other PF",]$PLAUR)

summary(test[test$Cluster == "Fibroblasts-3", ][test[test$Cluster == "Fibroblasts-3", ]$condition == "control",]$PLAUR)
summary(test[test$Cluster == "Fibroblasts-3", ][test[test$Cluster == "Fibroblasts-3", ]$condition == "COVID-19",]$PLAUR)
summary(test[test$Cluster == "Fibroblasts-3", ][test[test$Cluster == "Fibroblasts-3", ]$condition == "IPF/Other PF",]$PLAUR)

summary(test[test$Cluster == "Fibroblasts-4", ][test[test$Cluster == "Fibroblasts-4", ]$condition == "control",]$PLAUR)
summary(test[test$Cluster == "Fibroblasts-4", ][test[test$Cluster == "Fibroblasts-4", ]$condition == "COVID-19",]$PLAUR)
summary(test[test$Cluster == "Fibroblasts-4", ][test[test$Cluster == "Fibroblasts-4", ]$condition == "IPF/Other PF",]$PLAUR)

summary(test[test$Cluster == "Mesothelial", ][test[test$Cluster == "Mesothelial", ]$condition == "control",]$PLAUR)
summary(test[test$Cluster == "Mesothelial", ][test[test$Cluster == "Mesothelial", ]$condition == "COVID-19",]$PLAUR)
summary(test[test$Cluster == "Mesothelial", ][test[test$Cluster == "Mesothelial", ]$condition == "IPF/Other PF",]$PLAUR)

summary(test[test$Cluster == "Smooth Muscle Cells", ][test[test$Cluster == "Smooth Muscle Cells", ]$condition == "control",]$PLAUR)
summary(test[test$Cluster == "Smooth Muscle Cells", ][test[test$Cluster == "Smooth Muscle Cells", ]$condition == "COVID-19",]$PLAUR)
summary(test[test$Cluster == "Smooth Muscle Cells", ][test[test$Cluster == "Smooth Muscle Cells", ]$condition == "IPF/Other PF",]$PLAUR)

################################################################################
## Immune ##

data <- Convert("data/GSE158127_02immune.h5ad", ".h5seurat")
immune <- LoadH5Seurat(data)
immune@meta.data$condition <- NA
immune@meta.data$condition[which(str_detect(immune@meta.data$Diagnosis,"^Control"))] <- "control"
immune@meta.data$condition[which(str_detect(immune@meta.data$Diagnosis,"^COVID"))] <- "COVID-19"
immune@meta.data$condition[which(str_detect(immune@meta.data$Diagnosis,"^IPF"))] <- "IPF/Other PF"
immune@meta.data$condition[which(str_detect(immune@meta.data$Diagnosis,"^Other"))] <- "IPF/Other PF"

################################################################################
### Figure 1 ###
DimPlot(immune,reduction = "umap", group.by = "Cluster", split.by = "condition", raster=FALSE)
DimPlot(immune, reduction = "umap", group.by = "Cell.Type")

################################################################################
### Figure 2 ###
FeaturePlot(immune,reduction = "umap", split.by = "condition",
            features = c("PLAU", "PLAUR"))
################################################################################
### Figure 3 ###
VlnPlot(immune, split.by = "condition", features = c("PLAU", "PLAUR"))
VlnPlot(immune, features = c("PLAU", "PLAUR"), group.by = "Cluster", split.by = "condition", raster=FALSE)

test <- FetchData(immune, vars = c("PLAU", "condition", "Cluster"))
compare_means(PLAU ~ condition, data = test[test$Cluster == "AM1", ], method = "t.test", alternative="greater", p.adjust.method = "fdr")
compare_means(PLAU ~ condition, data = test[test$Cluster == "AM2", ], method = "t.test", alternative="greater", p.adjust.method = "fdr")
compare_means(PLAU ~ condition, data = test[test$Cluster == "MoM1", ], method = "t.test", alternative="greater", p.adjust.method = "fdr")
compare_means(PLAU ~ condition, data = test[test$Cluster == "MoM3", ], method = "t.test", alternative="greater", p.adjust.method = "fdr")
compare_means(PLAU ~ condition, data = test[test$Cluster == "MoM2", ], method = "t.test", alternative="greater", p.adjust.method = "fdr")
compare_means(PLAU ~ condition, data = test[test$Cluster == "Monocytes", ], method = "t.test", alternative="greater", p.adjust.method = "fdr")

summary(test[test$Cluster == "AM1", ][test[test$Cluster == "AM1", ]$condition == "control",]$PLAU)
summary(test[test$Cluster == "AM1", ][test[test$Cluster == "AM1", ]$condition =="COVID-19",]$PLAU)
summary(test[test$Cluster == "AM1", ][test[test$Cluster == "AM1", ]$condition == "IPF/Other PF",]$PLAU)

summary(test[test$Cluster == "AM2", ][test[test$Cluster == "AM2", ]$condition == "control",]$PLAU)
summary(test[test$Cluster == "AM2", ][test[test$Cluster == "AM2", ]$condition =="COVID-19",]$PLAU)
summary(test[test$Cluster == "AM2", ][test[test$Cluster == "AM2", ]$condition == "IPF/Other PF",]$PLAU)

summary(test[test$Cluster == "MoM1", ][test[test$Cluster == "MoM1", ]$condition == "control",]$PLAU)
summary(test[test$Cluster == "MoM1", ][test[test$Cluster == "MoM1", ]$condition =="COVID-19",]$PLAU)
summary(test[test$Cluster == "MoM1", ][test[test$Cluster == "MoM1", ]$condition == "IPF/Other PF",]$PLAU)

summary(test[test$Cluster == "MoM2", ][test[test$Cluster == "MoM2", ]$condition == "control",]$PLAU)
summary(test[test$Cluster == "MoM2", ][test[test$Cluster == "MoM2", ]$condition =="COVID-19",]$PLAU)
summary(test[test$Cluster == "MoM2", ][test[test$Cluster == "MoM2", ]$condition == "IPF/Other PF",]$PLAU)

summary(test[test$Cluster == "MoM3", ][test[test$Cluster == "MoM3", ]$condition == "control",]$PLAU)
summary(test[test$Cluster == "MoM3", ][test[test$Cluster == "MoM3", ]$condition =="COVID-19",]$PLAU)
summary(test[test$Cluster == "MoM3", ][test[test$Cluster == "MoM3", ]$condition == "IPF/Other PF",]$PLAU)

summary(test[test$Cluster == "Monocytes", ][test[test$Cluster == "Monocytes", ]$condition == "control",]$PLAU)
summary(test[test$Cluster == "Monocytes", ][test[test$Cluster == "Monocytes", ]$condition =="COVID-19",]$PLAU)
summary(test[test$Cluster == "Monocytes", ][test[test$Cluster == "Monocytes", ]$condition == "IPF/Other PF",]$PLAU)

test <- FetchData(immune, vars = c("PLAUR", "condition", "Cluster"))
compare_means(PLAUR ~ condition, data = test[test$Cluster == "AM1", ], method = "t.test", alternative="less", p.adjust.method = "fdr")
compare_means(PLAUR ~ condition, data = test[test$Cluster == "AM2", ], method = "t.test", alternative="less", p.adjust.method = "fdr")
compare_means(PLAUR ~ condition, data = test[test$Cluster == "MoM1", ], method = "t.test", alternative="less", p.adjust.method = "fdr")
compare_means(PLAUR ~ condition, data = test[test$Cluster == "MoM2", ], method = "t.test", alternative="less", p.adjust.method = "fdr")
compare_means(PLAUR ~ condition, data = test[test$Cluster == "MoM3", ], method = "t.test", alternative="greater", p.adjust.method = "fdr")
compare_means(PLAUR ~ condition, data = test[test$Cluster == "Monocytes", ], method = "t.test", alternative="less", p.adjust.method = "fdr")

summary(test[test$Cluster == "AM1", ][test[test$Cluster == "AM1", ]$condition == "control",]$PLAUR)
summary(test[test$Cluster == "AM1", ][test[test$Cluster == "AM1", ]$condition =="COVID-19",]$PLAUR)
summary(test[test$Cluster == "AM1", ][test[test$Cluster == "AM1", ]$condition == "IPF/Other PF",]$PLAUR)

summary(test[test$Cluster == "AM2", ][test[test$Cluster == "AM2", ]$condition == "control",]$PLAUR)
summary(test[test$Cluster == "AM2", ][test[test$Cluster == "AM2", ]$condition =="COVID-19",]$PLAUR)
summary(test[test$Cluster == "AM2", ][test[test$Cluster == "AM2", ]$condition == "IPF/Other PF",]$PLAUR)

summary(test[test$Cluster == "MoM1", ][test[test$Cluster == "MoM1", ]$condition == "control",]$PLAUR)
summary(test[test$Cluster == "MoM1", ][test[test$Cluster == "MoM1", ]$condition =="COVID-19",]$PLAUR)
summary(test[test$Cluster == "MoM1", ][test[test$Cluster == "MoM1", ]$condition == "IPF/Other PF",]$PLAUR)

summary(test[test$Cluster == "MoM2", ][test[test$Cluster == "MoM2", ]$condition == "control",]$PLAUR)
summary(test[test$Cluster == "MoM2", ][test[test$Cluster == "MoM2", ]$condition =="COVID-19",]$PLAUR)
summary(test[test$Cluster == "MoM2", ][test[test$Cluster == "MoM2", ]$condition == "IPF/Other PF",]$PLAUR)

summary(test[test$Cluster == "MoM3", ][test[test$Cluster == "MoM3", ]$condition == "control",]$PLAUR)
summary(test[test$Cluster == "MoM3", ][test[test$Cluster == "MoM3", ]$condition =="COVID-19",]$PLAUR)
summary(test[test$Cluster == "MoM3", ][test[test$Cluster == "MoM3", ]$condition == "IPF/Other PF",]$PLAUR)

summary(test[test$Cluster == "Monocytes", ][test[test$Cluster == "Monocytes", ]$condition == "control",]$PLAUR)
summary(test[test$Cluster == "Monocytes", ][test[test$Cluster == "Monocytes", ]$condition =="COVID-19",]$PLAUR)
summary(test[test$Cluster == "Monocytes", ][test[test$Cluster == "Monocytes", ]$condition == "IPF/Other PF",]$PLAU)

################################################################################
### Quality metrics ### -- same pipeline for all three datasets
seurat_epithelial$nCount_RNA <- seurat_epithelial$No.of.UMIs
seurat_epithelial$nFeature_RNA <- seurat_epithelial$No.of.genes
seurat_epithelial[["percent.mt"]] <- PercentageFeatureSet(seurat_epithelial, pattern="^MT-")
seurat_epithelial$log10GenesPerUMI <- log10(seurat_epithelial$nFeature_RNA)/log10(seurat_epithelial$nCount_RNA)
table(PercentageFeatureSet(object=seurat_epithelial, pattern="^MT-"))

## cell counts per sample
seurat_epithelial@meta.data %>%
  ggplot(aes(x=Sample.Name, fill=Study)) + 
  geom_bar() + theme_classic() + ggtitle("N Cells")

seurat_epithelial@meta.data %>%
  ggplot(aes(x=Sample.Status, fill=Sample.Status)) + 
  geom_bar() + theme_classic() + ggtitle("N Cells")

seurat_epithelial@meta.data %>%
  ggplot(aes(x=Diagnosis, fill=Diagnosis)) + 
  geom_bar() + theme_classic() + ggtitle("N Cells")

#VUILD, TILD - pulmonary fibrosis
#VUHD, THD - control

# UMIs per cell
seurat_epithelial@meta.data %>%
  ggplot(aes(color=Sample.Name, x=nCount_RNA, fill= Sample.Name)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)

seurat_epithelial@meta.data %>%
  ggplot(aes(color=Diagnosis, x=nCount_RNA, fill= Diagnosis)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)

# Genes per cell
seurat_epithelial@meta.data %>%
  ggplot(aes(color=Sample.Name, x=nFeature_RNA, fill= Sample.Name)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)

seurat_epithelial@meta.data %>%
  ggplot(aes(color=Diagnosis, x=nFeature_RNA, fill= Diagnosis)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)

# UMIs vs Genes detected
seurat_epithelial@meta.data %>%
  ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=X..of.mito.genes)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~Diagnosis)

# mitochondrial genes per cell
seurat_epithelial@meta.data %>%
  ggplot(aes(color=Sample.Name, x=X..of.mito.genes, fill=Sample.Name)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 20)

seurat_epithelial@meta.data %>%
  ggplot(aes(color=Diagnosis, x=X..of.mito.genes, fill=Diagnosis)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 20)

# log10genes per UMI
seurat_epithelial@meta.data %>%
  ggplot(aes(color=Sample.Name, x=log10GenesPerUMI, fill=Sample.Name)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.8)

################################################################################
### Filtering ###

seurat_epithelial <- subset(x = seurat_epithelial, subset= (nCount_RNA >= 500) & 
                            (nFeature_RNA >= 250) & 
                            (log10GenesPerUMI > 0.80) & 
                            (X..of.mito.genes < 30))
# UMIs vs Genes detected
seurat_epithelial@meta.data %>%
  ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=X..of.mito.genes)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~Diagnosis)

# filtering genes expressed in less than 40 cells  
counts <- GetAssayData(object = seurat_epithelial, slot = "counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 40
counts <- counts[keep_genes, ]
seurat_epithelial <- CreateSeuratObject(counts, meta.data = seurat_epithelial@meta.data)

seurat_epithelial@meta.data %>%
  ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=X..of.mito.genes)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~Diagnosis)

remove(counts)
remove(nonzero)
remove(keep_genes)
#saveRDS(seurat_epithelial, "results/filtered_seurat.rds")
seurat_epithelial_orig <- readRDS("results/filtered_seurat.rds")

################################################################################
### Checking for cell cycle phases ###

phase <- NormalizeData(seurat_epithelial)
load("data/cycle.rda")
phase <- CellCycleScoring(phase, g2m.features = g2m_genes, s.features = s_genes)
phase <- FindVariableFeatures(phase, selection.method = "vst")
phase <- ScaleData(phase)
phase <- RunPCA(phase)
# PCA colored by cell cycle phase
DimPlot(phase,
        reduction = "pca",
        group.by= "Phase",
        split.by = "Phase")
# no need to regress to cell cycle phase
remove(phase)

options(future.globals.maxSize = 40000 * 1024^2)


DimPlot(seurat_epithelial,reduction = "umap", group.by = "Cluster", split.by = "condition")
DimPlot(seurat_epithelial, reduction = "umap", group.by = "Cell.Type")

