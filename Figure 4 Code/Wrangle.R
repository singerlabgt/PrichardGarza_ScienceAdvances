# PRELIMINARIES ################################

# Steps to complete before running Wrangle.R:

# Upload fastq files to 10X Genomics Cloud Analysis: https://cloud.10xgenomics.com/cloud-analysis
# Send fastq for analysis in Cell Ranger to create count matrix file (Feature / cell matrix, filtered, TAR file)
# Cell Ranger Count v7.0.1, mm10 2020-A, include introns, skip BAM
# Download the TAR file using the Cloud Command Line Interface (CLI)
# Extract the files from the .tar.gz file using untar() in R or any other method (see below)
# Put the barcodes.tsv, genes.tsv (or features.tsv), and matrix.mtx witin the same directory/folder

# Outputs of Wrangle.R:

# Creates and saves filtered Seurat Objects of each sample from 10x Genomics count matrices
# (You will run the first section then navigate to Azimuth web app for refernce mapping)
# Adds Azimuth reference mapping results to each Seurat Object and saves a new copy
# Saves GSVA gene sets as R data file for later access during analysis (skip if not doing GSVA)
# Creates and saves pseudo-bulk aggregated count matrices for each sample then each cell type
# Saves results of DESeq2 analysis on cell type pseudo-bulk data for 40Hz vs 20Hz and No Flicker

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

if (!require("pacman")) install.packages("pacman")
pacman::p_load(Seurat,tidyverse,rio,DESeq2,readxl,GSVA)



# CREATE SEURAT OBJECTS ##################

## Untar the downloaded count matrix and put all three extracted files in same directory
sapply(1:12,function(x){untar(paste0("10x Genomics Cell Ranger Output/Feature Matrices Filtered TAR/S",x,"/filtered_feature_bc_matrix.tar.gz"),
                              exdir = paste0("10X Genomics Cell Ranger Output/Count Matrices/S",x))}) %>% invisible()


## Create and save Seurat objects for each sample to input to Azimuth ####

dir.create("Seurat Objects",showWarnings = FALSE)
lapply(1:12,function(i){
  cond = case_when(i %in% 1:4 ~ "Control",i %in% 5:8 ~ "20Hz",i %in% 9:12 ~ "40Hz")
  count.matrix = Read10X(data.dir = paste0("10x Genomics Cell Ranger Output/Count Matrices/S",i))
  data.so = CreateSeuratObject(counts = count.matrix, project = "CART") %>%
    AddMetaData(metadata = paste0("Sample",i), col.name = "Sample") %>%
    AddMetaData(metadata = cond, col.name = "Stim")
  
  # Filter out cells with too few (possibly droplets) or too many counts (possibly multiple cells) and cells with too much mito RNA
  data.so[["percent.mt"]] <- PercentageFeatureSet(data.so, pattern = "^mt-") # QC: mitochondrial RNA
  data.so <- subset(data.so, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 20)
  
  data.so = NormalizeData(data.so)
  data.so = FindVariableFeatures(data.so, selection.method = "vst", nfeatures = 2000)
  saveRDS(data.so,file=paste0("Seurat Objects/S",i,".rds"))
}) %>% invisible()

# At this point, navigate to Azimuth https://azimuth.hubmapconsortium.org/
# Use app for mapping to mouse motor cortex reference: https://app.azimuth.hubmapconsortium.org/app/mouse-motorcortex
# Upload the seurat objects one sample at a time and download the resulting cell type predictions and UMAP projections
# In the next step, we load each seurat object and add in the Azimuth labels as metadata and attach the UMAP

## Reference mapping with Azimuth ############

# Add Azimuth data to each sample and save new Seurat Object 
for(i in ind){
  sample.so = readRDS(file=paste0("Seurat Objects/S",i,".rds"))
  az_pred = import(paste0("Azimuth/S",i,"/azimuth_pred.tsv"))
  az_umap_out = import(paste0("Azimuth/S",i,"/azimuth_umap.Rds"))
  az_out = list(umap=az_umap_out,pred.df=az_pred)
  saveRDS(az_out,paste0("Azimuth/S",i,"/az_out.rds"))
  
  # Some cell type labels were combined due to very low numbers and total absence in some samples
  az_pred$predicted.subclass[az_pred$predicted.subclass=="L6 CT"] = "L6 CT/L6b"
  az_pred$predicted.subclass[az_pred$predicted.subclass=="L6b"] = "L6 CT/L6b"
  az_pred$predicted.subclass[az_pred$predicted.subclass=="Sst Chodl"] = "Sst"
  az_pred$predicted.subclass[az_pred$predicted.subclass=="L6 IT Car3"] = "L6 IT"
  az_pred$predicted.subclass[az_pred$predicted.subclass=="Sncg"] = "Vip/Sncg"
  az_pred$predicted.subclass[az_pred$predicted.subclass=="Vip"] = "Vip/Sncg"
  
  sample.so <- sample.so %>%
    AddMetaData(metadata = az_pred$predicted.subclass, col.name = "Subclass") %>%
    AddMetaData(metadata = az_pred$predicted.class, col.name = "Class") %>%
    AddMetaData(metadata = az_pred$predicted.cluster, col.name = "Cluster") %>%
    AddAzimuthResults(paste0("Azimuth/S",i,"/az_out.rds"))
  
  saveRDS(sample.so,file=paste0("Seurat Objects/S",i,"_az.rds"))
  rm(sample.so)
}

cell_types=c("Micro-PVM","L2/3 IT","L6 IT","Astro","Sst","VLMC","Oligo","OPC","L5 IT","Meis2",
             "L5 ET","L5/6 NP","Pvalb","Endo","Lamp5","Peri","Vip/Sncg","L6 CT/L6b")
saveRDS(cell_types,paste0("R Data/cell_types.rds"))  # save this for later


# CREATE PSEUDOBULK DATA ###############

## Create Pseudobulk CSV files for each Sample, Aggregate Counts #####

dir.create("Pseudobulk Data/Aggregated Counts by Sample",recursive = TRUE,showWarnings=FALSE)
pb_samples = lapply(1:12,function(i){
  sample.so = readRDS(paste0("Seurat Objects/S",i,"_az.rds"))
  Idents(sample.so) = "Subclass"
  pb = AggregateExpression(sample.so,slot="counts")[[1]]  # sums up gene counts in all cells by cell type 
  export(pb,paste0("Pseudobulk Data/Aggregated Counts by Sample/pb_S",i,".csv"),row.names=TRUE)  # save count matrix as a csv for later access
  return(as.data.frame(pb))
  })
names(pb_samples) = str_c("S",1:12)

# Using list of sample pseudo-bulk (pb) matrices, create a matrix for each separate cell type
cell_types = readRDS("R Data/cell_types.rds")
dir.create("Pseudobulk Data/Aggregated Counts by Cell Type",recursive = TRUE,showWarnings=FALSE)
pb_ct = lapply(cell_types,function(ct){
  pb = do.call(cbind,lapply(pb_samples,function(x){select(x,ct)}) )
  colnames(pb)=names(pb_samples) # col names are sample labels
  rownames(pb)=pb_list[[1]]$V1  # row names are gene names
  saveRDS(pb,paste0("Pseudobulk Data/Aggregated Counts by Cell Type/pb_",str_replace(ct,pattern = "/",replacement = "-"),".rds"))
  rio::export(pb,paste0("Pseudobulk Data/Aggregated Counts by Cell Type/pb_",str_replace(ct,pattern = "/",replacement = "-"),".csv"),row.names=TRUE)
  return(as.data.frame(pb))
})
names(pb_ct) = cell_types
rm(pb_samples)

## Normalize Pseudobulk data with DESeq2, run and save DESeq2 results by cell type  ###########
# save results in separate folders for 40Hz vs 20Hz and 40Hz vs No Flicker

cell_types=readRDS(paste0("R Data/cell_types.rds"))
metadata=read_xlsx("22131-03 metadata.xlsx") %>%
  mutate(Frequency = factor(Frequency,levels = c("40Hz","20Hz","None")))
dir.create("Pseudobulk Data/DESeq2 Results 40 vs 20",recursive = TRUE,showWarnings=FALSE)
dir.create("Pseudobulk Data/DESeq2 Results 40 vs No Flicker",recursive = TRUE,showWarnings=FALSE)

lapply(1:length(pb_ct),function(i){
  pb = pb_ct[[i]]
  ct = names(pb_ct)[i]
  # DESeq2 Normalization  
  dds <- DESeqDataSetFromMatrix(countData = pb,
                                colData = metadata,
                                design = ~ Frequency)
  dds <- DESeq(dds)
  res20v40 <- results(dds, name="Frequency_20Hz_vs_40Hz")
  resNonev40 <- results(dds, name="Frequency_None_vs_40Hz")
  
  export(as.data.frame(res20v40),paste0("Pseudobulk Data/DESeq2 Results 40 vs 20/",str_replace(ct,pattern = "/",replacement = "-"),"_all.csv"),row.names=TRUE)
  export(as.data.frame(resNonev40),paste0("Pseudobulk Data/DESeq2 Results 40 vs No Flicker/",str_replace(ct,pattern = "/",replacement = "-"),"_all.csv"),row.names=TRUE)
}) %>% invisible()
