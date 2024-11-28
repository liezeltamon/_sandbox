# Palmo tutorial 2

renv::load()

library(PALMO)
library(Seurat)
library(tidyverse)

setwd("test")

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
# palmo_obj@filePATH <- "/ceph/project/sims-lab/ltamon/git/_sandbox/results/palmo_tutorial_2"
# palmo_obj

#The expression dataframe columns merged with input annotation dataframe. Only overlapping samples kept. Missing annotations with Sample, Donor/participant, or Time columns are removed from downstream analysis.
palmo_obj <- mergePALMOdata(data_object=palmo_obj, datatype="singlecell")

celltype_oi <- c("NK", "B-naive")
palmo_obj <- sclongitudinalDEG(data_object=palmo_obj, scassay="RNA",
                               group_column="celltype",
                               group_oi = celltype_oi)









palmo_obj <- avgExpCalc(data_object=palmo_obj,
                        assay="RNA", group_column="celltype")

head(palmo_obj@curated[["anndata"]]) #merged annotation data
head(palmo_obj@curated[["data"]]) #scRNA average expression data

##

data_object = palmo_obj
group_column = "celltype"
assay = "RNA"
  
anndata <- data_object@curated$anndata
dataObj <- data_object@curated$SeuratObj

## Add sample group to metadata Define sample group and Calculate
## average expression
dataObj@meta.data$group <- dataObj@meta.data[, group_column]
dataObj@meta.data$Sample_group <- paste(dataObj@meta.data$Sample,
                                        dataObj@meta.data$group,
                                        sep = ":")
dataObj@meta.data$Sample_group <- gsub(" ", "_", dataObj@meta.data$Sample_group)
metaData <- dataObj@meta.data

## Check assay
DefaultAssay(dataObj) <- assay
## Average expression on log-scaled data
dataObj@assays[[assay]]@counts <- dataObj@assays[[assay]]@data

scrna_avgmat <- AverageExpression(object = dataObj, assays = assay,
                                  slot = "counts",
                                  group.by = "Sample_group",
                                  verbose = TRUE)
cn <- data.frame(colnames(scrna_avgmat[[assay]]))
#cn <- gsub("data\\[, 1]", "", as.character(row.names(cn)))
scrna_avgmat <- data.frame(scrna_avgmat[[assay]], check.names = FALSE,
                           stringsAsFactors = FALSE)

#colnames(scrna_avgmat) <- cn
colnames(scrna_avgmat) <- cn$colnames.scrna_avgmat..assay...

message(date(), ": scRNA Average expression finished")

## Keep genes with avgExpression > zero
rowDF <- rowSums(scrna_avgmat)
rowDF <- rowDF[rowDF > 0]
mat <- scrna_avgmat[names(rowDF), ]
message(date(), ": Keeping genes with avg expression >0")

## Create annotation
cn <- data.frame(Sample_group = colnames(mat))
temp <- data.frame(do.call(rbind, strsplit(cn$Sample_group,
                                           split = "_")),
#temp <- data.frame(do.call(rbind, strsplit(cn$Sample_group,
#                                           split = ":")),
                   stringsAsFactors = FALSE)

# cn <- data.frame(cn, Sample = temp$X1, group = temp$X2,
#                  stringsAsFactors = FALSE)
cn <- data.frame(cn, Sample = temp$X2, group = temp$X1,
                 stringsAsFactors = FALSE)
row.names(cn) <- cn$Sample_group
cn <- merge(cn, anndata, by = "Sample", all = TRUE)
cn <- cn[!is.na(cn$Sample_group), ]
row.names(cn) <- cn$Sample_group
ann <- cn
ann[, ncol(ann) + 1] <- ann$group
colnames(ann) <- c(colnames(ann)[1:ncol(ann) - 1], group_column)
ann$Sample_group_i <- paste(ann$group, ann$PTID, sep = ":")
rm(cn)

palmo_obj@curated$anndata <- ann
palmo_obj@curated$data <- mat
palmo_obj@rownames <- row.names(mat)
palmo_obj@colnames <- colnames(mat)

#####
#CV profile
palmo_obj <- cvCalcSCProfile(data_object=palmo_obj,
                             housekeeping_genes=c("GAPDH", "ACTB"),
                             meanThreshold = 0.1,
                             fileName="scrna")

#Features contributing towards donor variations
#Variance decomposition
featureSet <- c("Time","celltype")
palmo_obj <- lmeVariance(data_object=palmo_obj,
                         featureSet=featureSet,
                         meanThreshold=0.1, cl=4,
                         fileName="scrna")
var_decomp <- palmo_obj@result$variance_decomposition
head(var_decomp[,featureSet])
#Variance contributing features
plots <- variancefeaturePlot(vardata=var_decomp, featureSet=featureSet,
                             Residual=F,
                             cols=c("purple", "darkgreen"))
