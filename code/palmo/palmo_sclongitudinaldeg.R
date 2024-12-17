# Palmo tutorial 2 data + sclongitudinalDEG

renv::load()
.libPaths()

library(RSpectra)
library(RANN) # Needed by PALMO
library(RColorBrewer)
library(PALMO)
library(Seurat)
library(tidyverse)

out_dir = file.path("results", "palmo_tutorial_2")

# Load data

pbmc <- readRDS("data/AIFI-scRNA-PBMC-FinalData.RDS")
metaData <- pbmc@meta.data
pbmc@meta.data$Sample <- pbmc@meta.data$orig.ident
pbmc@meta.data$celltype <- gsub(" ", "-", pbmc@meta.data$celltype)
pbmc@meta.data$celltype <- gsub("_", "-", pbmc@meta.data$celltype)
#Load annotation datas
load("data/AIFI-Metadata.Rda")

avgGroup <- "celltype"
#Celltypes observed in dataset
cell_type <- sort(unique(pbmc@meta.data$celltype))
#Celltypes selected for analysis consisting atleast >5% of cells in each celltype.
celltype_oi <- c("CD4_Naive","CD4_TEM","CD4_TCM","CD4_CTL","CD8_Naive",
                 "CD8_TEM","CD8_TCM","Treg","MAIT","gdT",
                 "NK", "NK_CD56bright",
                 "B_naive", "B_memory", "B_intermediate",
                 "CD14_Mono","CD16_Mono",
                 "cDC2","pDC")

# Create PALMO object

#Create PALMO object
palmo_obj <- createPALMOobject(anndata=ann, data=pbmc)

#Assign Sample, PTID and Time parameters
palmo_obj <- annotateMetadata(data_object=palmo_obj,
                              sample_column= "Sample",
                              donor_column= "PTID",
                              time_column= "Time")

# Custom code
palmo_obj@filePATH = out_dir

#The expression dataframe columns merged with input annotation dataframe. Only overlapping samples kept. Missing annotations with Sample, Donor/participant, or Time columns are removed from downstream analysis.
palmo_obj <- mergePALMOdata(data_object=palmo_obj, datatype="singlecell")

celltype_oi <- c("NK", "B-naive")
palmo_obj <- sclongitudinalDEG(data_object=palmo_obj, scassay="RNA",
                               group_column="celltype",
                               group_oi = celltype_oi)

saveRDS(palmo_obj, file.path(out_dir, "palmo_obj_sclongitudinalDEG.rds"))
