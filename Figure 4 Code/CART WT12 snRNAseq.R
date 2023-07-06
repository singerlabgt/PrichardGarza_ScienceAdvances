# PRELIMINARIES ################################

# Run Wrangle.R first to create Seurat Objects and pseudobulk data

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

if (!require("pacman")) install.packages("pacman")
if (!require("BiocManager")) install.packages("BiocManager")
pacman::p_load(tidyverse,rio,gplots,cowplot,RColorBrewer,Seurat,ggpubr,heatmap3,WoodLabFunctions,EnhancedVolcano)

theme_set(theme_cowplot())
stim_palette= c("#bebebe","#1071b9","#ec242a")  # gray for No Flicker, blue for 20Hz, red for 40Hz


# FIGURE PANEL A: UMAP of All Samples from Azimuth Projections ##################

ind=1:12  # index for 12 samples, using in for loop
umap_all = data.frame()  # create blank data frame which will be filled with UMAP projection data for each sample

# load UMAP projection from each sample to umap_all
# UMAP projections are created in Azimuth by mapping our query data onto the reference set (primary mouse motor cortex)
for(i in ind){
  rm(SO)  # remove last seurat object from environment to limit RAM use
  SO = readRDS(file=paste0("Seurat Objects/S",i,"_az.rds"))  # load the seurat object for sample number i
  # save the umap projection coordinates alongside cell type labels, stimulation condition, and sample number
  umap=SO@reductions$umap.proj@cell.embeddings %>% as.data.frame() %>%
    cbind(Stim=SO@meta.data$Stim,Class=SO@meta.data$Class,Subclass=SO@meta.data$Subclass,Sample=SO@meta.data$Sample)
  umap_all = rbind(umap_all,umap)  # add this sample's data to the total data frame
}

# create a ggplot object with the cell type labels at each cluster centroid, to be added to each other ggplot 
# DimPlot makes labeling centroids super easy but limits other capabilities, so I switch to ggplot moving forward
p_labels = DimPlot(SO, reduction = "umap.proj", group.by="Subclass",label = TRUE, repel = TRUE)
rm(SO)  # remove last seurat object from environment to limit RAM use

# Create a UMAP colored by Subclass (cell type)
p1=ggplot(data=umap_all,aes_string(x="UMAP_1",y="UMAP_2",color="Subclass")) +
  geom_point(size=0.2)
p1$layers[[2]] = p_labels$layers[[2]]
# plot(p1)

# Create a UMAP colored by Stimulation
p2=ggplot(data=umap_all %>% arrange(sample(1:nrow(umap_all),nrow(umap_all))),aes_string(x="UMAP_1",y="UMAP_2",color="Stim")) +
  geom_point(size=0.001,alpha=0.1) + 
  scale_color_manual(values=c("#1071b9","#ec242a","#bebebe")) +
  theme(legend.position="none") +
  xlab("UMAP1") +
  ylab("UMAP2")
p2$layers[[2]] = p_labels$layers[[2]]
# plot(p2)

# Create a UMAP colored by sample number
p3=ggplot(data=umap_all,aes_string(x="UMAP_1",y="UMAP_2",color="Sample")) +
  geom_point(size=0.2)
p3$layers[[2]] = p_labels$layers[[2]]
# plot(p3)

# Create a UMAP colored by Class (GABA,gluta,non-neuronal)
p4=ggplot(data=umap_all,aes_string(x="UMAP_1",y="UMAP_2",color="Class")) +
  geom_point(size=0.2)
p4$layers[[2]] = p_labels$layers[[2]]
# plot(p4)


dir.create(paste0("Figures/Panel A - UMAP"),recursive = TRUE, showWarnings = FALSE)
png(paste0("Figures/Panel A - UMAP/UMAP by Subclass.png"),res=1000,units="in",height=4.5,width=7); p1; dev.off()
png(paste0("Figures/Panel A - UMAP/UMAP by Stim.png"),res=1000,units="in",height=4.5,width=5); p2; dev.off()
png(paste0("Figures/Panel A - UMAP/UMAP by Sample.png"),res=1000,units="in",height=4.5,width=5.5); p3; dev.off()
png(paste0("Figures/Panel A - UMAP/UMAP by Class.png"),res=1000,units="in",height=4.5,width=5.5); p4; dev.off()

pdf(paste0("Figures/Panel A - UMAP/UMAP by Subclass.pdf"),height=4.5,width=7); p1; dev.off()
pdf(paste0("Figures/Panel A - UMAP/UMAP by Stim.pdf"),height=4.5,width=5); p2; dev.off()
pdf(paste0("Figures/Panel A - UMAP/UMAP by Sample.pdf"),height=4.5,width=5.5); p3; dev.off()
pdf(paste0("Figures/Panel A - UMAP/UMAP by Class.pdf"),height=4.5,width=5.5); p4; dev.off()



# FIGURE PANEL B: Number of Pseudo-bulk DE Genes by Cell Type #############

cell_types=readRDS(paste0("R Data/cell_types.rds"))
NumDE = data.frame()
DEs = list()
pb_ct_dir = "Pseudobulk Data/DESeq2 Results "

for(i in 1:length(cell_types)){
  ct = cell_types[i]
  safe_ct = str_replace(ct,pattern = "/",replacement = "-") # replace incompatible characters
  DE_20v40 = import(paste0(pb_ct_dir,"40 vs 20/",safe_ct,"_all.csv")) %>%
    filter(padj < 0.05) %>%
    mutate(comparison = "20Hz vs 40Hz") %>%
    mutate(celltype = ct)
  DE_Nonev40 = import(paste0(pb_ct_dir,"40 vs No Flicker/",safe_ct,"_all.csv")) %>%
    filter(padj < 0.05) %>%
    mutate(comparison = "No Flicker vs 40Hz") %>%
    mutate(celltype = ct)
  add = data.frame(ct=ct,num20v40=nrow(DE_20v40),numNonev40=nrow(DE_Nonev40))
  NumDE = rbind(NumDE,add)
  
  DEs$new = rbind(DE_20v40,DE_Nonev40)
  names(DEs)[which(names(DEs)=="new")] = str_replace(ct,pattern = "/",replacement = "-")
}

# save a csv of all DE genes (padj < 0.05) across all cell types and both comparisons
DE_all = do.call(rbind,DEs)
export(DE_all,"Pseudobulk Data/All DE Genes by Cell Type and Comparison.csv")

NumDE = NumDE %>%
  arrange(desc(`num20v40`),desc(`numNonev40`)) %>%
  mutate(ct = factor(ct,levels=ct)) %>%
  gather(key="comparison",value="n",2:3) %>%
  mutate(comp_group = case_when(
    comparison == "num20v40" ~ "40Hz vs 20Hz",
    comparison == "numNonev40" ~ "40Hz vs Light",
  ))

stim_palette= c("#bebebe","#1071b9","#ec242a")  # gray for No Flicker, blue for 20Hz, red for 40Hz

dir.create(paste0("Figures/Panel B - Bar Plot Number of DE Genes by Cell Type"),recursive = TRUE, showWarnings = FALSE)
pdf("Figures/Panel B - Bar Plot Number of DE Genes by Cell Type/DE Genes Bar Plot 40Hz vs 20Hz and Light.pdf",width=11,height=5)
ggplot(NumDE,aes(x=ct,y=n,fill=comp_group)) +
  geom_col(position="dodge",width=0.8) +
  geom_text(aes(label = n), vjust = -0.5,position=position_dodge(width=0.8),size=2.7) +
  scale_fill_manual(values=c(stim_palette[2],stim_palette[1]),name="Comparison:") + 
  theme(axis.text.x = element_text(angle=90,hjust=.95,vjust=0.5)) +
  xlab("Cell Type") +
  ylab("Differentially Expressed Genes")
dev.off()



# FIGURE PANEL C: Volcano plots for cell types of interest  #############

ind = 1:12
cell_types=readRDS(paste0("R Data/cell_types.rds"))

# Define colors to use in volcano plots by stimulation group
volcano_palette = c(`None`="#bebebe",        # grey
                    `No Flicker`="#bebebe",  # grey
                    `20Hz`="#1071b9",        # blue
                    `40Hz`="#ec242a",        # red
                    `ns`="black")

# comp_group = "20Hz"
comp_group = "None"
cell_types_subset = c("Micro-PVM","Astro","L2/3 IT")   # Focusing analysis on three cell types of interest
dir.create(paste0("Figures/Panel C - Pseudobulk Volcanos"),recursive = TRUE, showWarnings = FALSE)


for(ct in cell_types_subset){
  if(comp_group == "None"){comp_label = "No Flicker"; folder="No Flicker/"}
  if(comp_group == "20Hz"){comp_label = comp_group; folder="20/"}
  ct_comp=import(paste0("Pseudobulk Data/DESeq2 Results 40 vs ",folder,str_replace(ct,pattern = "/",replacement = "-"),"_all.csv")) %>%
  drop_na() %>%
  dplyr::rename(gene=V1) %>%
  mutate(logp = -log(padj,10)) %>%
  mutate(elevated = case_when(
    log2FoldChange < 0 & padj < 0.05 ~ "40Hz",
    log2FoldChange > 0 & padj < 0.05 ~ comp_label,
    TRUE ~ "ns"
  )) %>%
  mutate(labels = case_when(
    elevated != "ns" ~ gene,
    TRUE ~ ""
  )) %>%
  mutate(flip_FC = -log2FoldChange) # for flipping axis to show 40Hz on right side

  p_title = ct
  if(ct == "Astro"){p_title = "Astrocytes"}
  if(ct == "Micro-PVM"){p_title = "Microglia / PVM"}
  if(ct == "L2/3 IT"){p_title = "L2/3 IT Neurons"}

p1=ggscatter(data=ct_comp,x="flip_FC",y="logp",title = p_title,
          ylab="p-value (-log10)",xlab="Fold Change (log2)",color="elevated",
          repel = T,label="labels") +
  theme(legend.position = "none",
        plot.title = element_text(hjust=0.5,size=20,face="bold"))
p2=ggplot() +
  geom_hline(yintercept = -log(0.05,10),linetype="dashed",color="gray") +
  annotate("text",y={-0.03*range(ct_comp$logp)[2]}-log(0.05,10),x=0.95*min(ct_comp$flip_FC,na.rm = TRUE),label="p=0.05",color="gray")
p1$layers = c(p2$layers, p1$layers) 
p1$layers[[4]]$aes_params$size=3
p1$layers[[4]]$aes_params$fontface="italic"

cairo_pdf(paste0("Figures/Panel C - Pseudobulk Volcanos/",str_replace(ct,pattern = "/",replacement = "-")," Volcano ",comp_group," vs 40Hz.pdf"),width=5,height=5)
print(ggpar(p1,palette = volcano_palette))
dev.off()

png(paste0("Figures/Panel C - Pseudobulk Volcanos/",str_replace(ct,pattern = "/",replacement = "-")," Volcano ",comp_group," vs 40Hz.png"),width=7,height=7,units="in",res=1000)
print(ggpar(p1,palette = volcano_palette))
dev.off()
}



# FIGURE PANEL D/E: Gene Ontology on Pseudobulk by Cell Type: Prep Gene Lists #############
# create .txt files of gene query and background sets for easy upload to GO PANTHER web app

ind=1:12
cell_types=readRDS(paste0("R Data/cell_types.rds"))
cell_types_subset = c("Micro-PVM","Astro","L2/3 IT")
dir.create(paste0("Pseudobulk Data/Gene Ontology 40 vs 20"),recursive = TRUE, showWarnings = FALSE)
dir.create(paste0("Pseudobulk Data/Gene Ontology 40 vs No Flicker"),recursive = TRUE, showWarnings = FALSE)

for(i in 1:length(cell_types_subset)){
  ct = cell_types_subset[i]
  safe_ct = str_replace(ct,pattern = "/",replacement = "-")
  
  All_20v40 = import(paste0("Pseudobulk Data/DESeq2 Results 40 vs 20/",safe_ct,"_all.csv")) %>%
    drop_na()
  # Gene Background
  export(data.frame(gene=All_20v40$V1),paste0("Pseudobulk Data/Gene Ontology 40 vs 20/",safe_ct," Background.txt"),col.names=FALSE)
  # DE Genes
  DE_20v40 = filter(All_20v40,padj<0.05 & log2FoldChange < 0)
  export(data.frame(gene=DE_20v40$V1),paste0("Pseudobulk Data/Gene Ontology 40 vs 20/",safe_ct," Significantly Elevated in 40Hz.txt"),col.names=FALSE)
  
  
  
  DE_Nonev40 = import(paste0("Pseudobulk Data/DESeq2 Results 40 vs No Flicker/",safe_ct,"_all.csv")) %>%
    drop_na()
  # Gene Background
  export(data.frame(gene=All_20v40$V1),paste0("Pseudobulk Data/Gene Ontology 40 vs No Flicker/",safe_ct," Background.txt"),col.names=FALSE)
  # DE Genes
  DE_20v40 = filter(All_20v40,padj<0.05 & log2FoldChange < 0)
  export(data.frame(gene=DE_20v40$V1),paste0("Pseudobulk Data/Gene Ontology 40 vs No Flicker/",safe_ct," Significantly Elevated in 40Hz.txt"),col.names=FALSE)
}

# conduct GO analysis online via PANTHER 17.0: http://pantherdb.org/webservices/go/overrep.jsp
# save table of results to Pseudobulk Data/Gene Ontology Results 40 vs 20


# FIGURE PANEL D/E: Gene Ontology on Pseudobulk by Cell Type: Bar Graphs for Enriched Processes ############

GO_results_microglia = read_tsv("Pseudobulk Data/Gene Ontology 40 vs 20/Micro-PVM GO Results.txt",skip = 11)
colnames(GO_results_microglia) = c("GO_biological_process",
                                   "set_n_ref",
                                   "set_n_query_actual",
                                   "set_n_query_expected",
                                   "over_under",
                                   "fold_enrichment",
                                   "fishers_p",
                                   "padj_fdr")
GO_results_microglia_filtered = GO_results_microglia %>%
  filter(padj_fdr < 0.25) %>%
  mutate(fold_enrichment = as.numeric(fold_enrichment)) %>%
  arrange(fold_enrichment) %>%
  mutate(GO_biological_process = factor(GO_biological_process,levels=GO_biological_process)) %>%
  mutate(logp = -log(padj_fdr,10))

GO_results_L23 = read_tsv("Pseudobulk Data/Gene Ontology 40 vs 20/L2-3 IT GO Results.txt",skip = 11)
colnames(GO_results_L23) = colnames(GO_results_microglia)
GO_results_L23_filtered = GO_results_L23 %>%
  filter(padj_fdr < 0.25) %>%
  mutate(fold_enrichment = as.numeric(fold_enrichment)) %>%
  arrange(fold_enrichment) %>%
  mutate(GO_biological_process = factor(GO_biological_process,levels=GO_biological_process)) %>%
  mutate(logp = -log(padj_fdr,10))


dir.create(paste0("Figures/Panel D - Bar Plot Upregulated GO Sets in Microglia-PVM"),recursive = TRUE, showWarnings = FALSE)
pdf("Figures/Panel D - Bar Plot Upregulated GO Sets in Microglia-PVM/Microglia Upregulated Processes.pdf",height=3.5,width=12)
ggplot(GO_results_microglia_filtered,aes(x=GO_biological_process,y=fold_enrichment,
                                         label=str_c(signif(fishers_p,1)," / ",signif(padj_fdr,1)),fill=logp)) +
  geom_col(color="black") +
  coord_flip() +
  geom_text(hjust=-0.1) +
  expand_limits(y=23) +
  xlab("") +
  ylab("Fold Enrichment") +
  ggtitle("Upregulated Biological Processes in Microglia / PVM") +
  theme(plot.title = element_text(face="bold",hjust=0.5,size=20)) +
  scale_fill_gradient(low="white",high="purple",name="FDR p (-log10)") +
  annotate("text",x=3,y=0.9*max(GO_results_microglia_filtered$fold_enrichment),label="Fisher's p / FDR-adjusted p")
dev.off()


dir.create(paste0("Figures/Panel E - Bar Plot Upregulated GO Sets in L2-3 IT"),recursive = TRUE, showWarnings = FALSE)
pdf("Figures/Panel E - Bar Plot Upregulated GO Sets in L2-3 IT/L2_3 IT Upregulated Processes.pdf",height=3.5,width=12)
ggplot(GO_results_L23_filtered,aes(x=GO_biological_process,y=fold_enrichment,
                                         label=str_c(signif(fishers_p,1)," / ",signif(padj_fdr,1)),fill=logp)) +
  geom_col(color="black") +
  coord_flip() +
  geom_text(hjust=-0.1) +
  expand_limits(y=4.2) +
  xlab("") +
  ylab("Fold Enrichment") +
  ggtitle("Upregulated Biological Processes in L2/3 IT Neurons") +
  theme(plot.title = element_text(face="bold",hjust=0.5,size=20)) +
  scale_fill_gradient(low="white",high="purple",name="FDR p (-log10)") +
  annotate("text",x=3,y=1*max(GO_results_L23_filtered$fold_enrichment),label="Fisher's p / FDR-adjusted p")
dev.off()

