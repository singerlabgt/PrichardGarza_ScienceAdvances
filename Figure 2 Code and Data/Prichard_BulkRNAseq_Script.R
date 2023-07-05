# ---- Preliminaries ---- 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

options(java.parameters = "-Xmx8g")
library("fastcluster")
library("gplots")
library("dendextend")
library(reshape)
library(matrixStats)
library("readxl")
library("heatmap3")
library("compare")
library("gridGraphics")
library("ggplot2")
library("Rmisc") 
library("ggpubr")
library("GSVA")
library("RColorBrewer")
library("Hmisc")
library("ggfortify")#for PCA
library('EnhancedVolcano')
library('dplyr')
library("gplots")
source("summarySE.R")
library("DESeq2")
library("EnhancedVolcano")
library("plotrix")

#Uppercase String pasting function
r_ucfirst <- function (str) {
  paste(toupper(substring(str, 1, 1)), tolower(substring(str, 2)), sep = "")
}



# ---- Set heatmap3 color bar parameters ---- 
breakBarColors=c(-200,seq(-1.5, 1.5, 0.01),200) #Outside numbers clip outliers. This is for zscoring.
barColors = colorpanel(length(breakBarColors)-1, "snow4", "white", "mediumorchid1")

breakBarColorsCor=c(-200,seq(-1, 1.5, 0.01),200) #Outside numbers clip outliers. This is for zscoring.
barColorsCor = colorpanel(length(breakBarColorsCor)-1, "black", "white", "orange")


# ---- Load gene expression data ---- 

dataIn=read.csv("normalizedCounts_Prichard.csv")
countData=data.matrix(dataIn[,2:(ncol(dataIn))])
rownames(countData)=dataIn$X
metaData=c(rep("20",7), rep("40",7)) #Label  stimulation frequency


#Define condition color bars
colorBar=matrix(,nrow=ncol(countData), ncol=1)
colorBar[which(metaData=="20")]="slateblue1"
colorBar[which(metaData=="40")]="red"


#pdfPath=paste("./FiguresOut/nocluster.pdf",sep='')
#pdf(pdfPath, width=6,height=8,pointsize = 14)
heatmap3(countData, ColSideWidth=1, 
         col=barColors, breaks=breakBarColors, legendfun=function()showLegend(legend=c(NA),col=c(NA),cex=1.5), 
         ColSideColors = c(colorBar),
         Rowv=NA, Colv=NA,  scale="row",
         labRow = NA, labCol = NA) 
#dev.off()

#z-score the data
countDataZ=t(apply(countData,1,scale))
rownames(countDataZ)=rownames(countData)
colnames(countDataZ)=colnames(countData)


# ---- DESeq2 ----
# Round normalized counts
y=round(countData); rownames(y)=toupper(rownames(countData))


#40Hz vs 20Hz
pheno=metaData
coldata=data.frame(sample=colnames(y), condition=pheno, row.names="sample")
coldata$condition = factor(coldata$condition) 
dds = DESeqDataSetFromMatrix(countData = y, colData = coldata, 
                             design = ~ condition)
dds$condition <- relevel(dds$condition, ref = "20")
dds = DESeq(dds)
resultsNames(dds)
res=results(dds, name = "condition_40_vs_20")

write.csv(as.data.frame(res), file="p_deseq2_40_vs_20_Prichard.csv")


# ---- Gene Set Variation Analysis ----

# Read in the C2 Gene Sets
geneSetsMat=read.csv("c2.cp.v7.5.1.symbols.csv", header=FALSE)
geneSetsMat=geneSetsMat[,-2]
geneSetsMat_T=t(geneSetsMat)
colnames(geneSetsMat_T)=geneSetsMat_T[1,]
geneSetsMat_T=geneSetsMat_T[-1,]
geneSets3=as.list(as.data.frame(geneSetsMat_T))
geneSets=lapply(geneSets3, function(x) x[!is.na(x)])


geneSetEnrich=gsva(y, geneSets, min.sz=7, mx.diff=TRUE, kcdf="Poisson")
meanDiff=rowMeans(geneSetEnrich[,which(pheno=="40")])-rowMeans(geneSetEnrich[,which(pheno=="20")])


# ---- Permutation p-value ---- 

#Define number of permutations
  R=1000

#Compute mean differences between groups during permutation   
meanOut1=matrix(NA,nrow=nrow(geneSetEnrich), ncol = R)
for (i in 1:R)
  {
    yShuffle=y
    rownames(yShuffle)=rownames(yShuffle)[sample(1:dim(y)[1])] #random shuffling gene row titles in matrix
    gsvaBoot=gsva(yShuffle,geneSets, min.sz=7, mx.diff=TRUE, kcdf="Poisson") 
    
    
    meanOut1[,i]=rowMeans(gsvaBoot[,metaData=="40"])-rowMeans(gsvaBoot[,metaData=="20"]) 
    print(i)

  }
#Compute mean differences between groups for true gene assignments   
meanTrue1=rowMeans(geneSetEnrich[,metaData=="40"])-rowMeans(geneSetEnrich[,metaData=="20"]) 
      
#Compute p values  
pBoot1=matrix(NA,nrow=nrow(geneSetEnrich), ncol = 1)
for (j in 1:nrow(meanOut1))
{
  pBoot1[j] = mean(abs(meanOut1[j,]) > abs(meanTrue1[j]))
}

p1FDR=p.adjust(pBoot1, "fdr")

dataOut=data.frame("Pathway"=rownames(geneSetEnrich), meanTrue1, pBoot1, p1FDR)
write.csv(dataOut, file="pPermute_pathways_Prichard.csv", row.names=FALSE)


# ---- Make heatmaps and bar plots for selected pathways ----

gsvaList=c("REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM",
           "REACTOME_COMPLEMENT_CASCADE",
           "BIOCARTA_STRESS_PATHWAY",
           "REACTOME_MAP2K_AND_MAPK_ACTIVATION", 
           "REACTOME_ACETYLCHOLINE_NEUROTRANSMITTER_RELEASE_CYCLE",
           "REACTOME_NEUROTRANSMITTER_RELEASE_CYCLE",
           "REACTOME_NEUROTRANSMITTER_RECEPTORS_AND_POSTSYNAPTIC_SIGNAL_TRANSMISSION",
           "WP_PURINERGIC_SIGNALING")


indGSVASelect=match(gsvaList, rownames(geneSetEnrich))
gsvaSelect=geneSetEnrich[indGSVASelect,]


for (indGeneSetPlot in 1:nrow(gsvaSelect))
  
{
  
  dataBar=data.frame(unlist(gsvaSelect[indGeneSetPlot,]), unlist(colorBar))
  
  colnames(dataBar)=c("ES","colors")
  
  #Summarize stats for each group color
  statsOut = summarySE(dataBar, measurevar="ES", groupvars="colors")
  colnames(statsOut)[3]="ES"
  
  #meanBlue4vsOR4[indGeneSetPlot]=statsOut$ES[statsOut$color=="orangered4"]-statsOut$ES[statsOut$color=="skyblue4"]
  
  #Re-order the colors to match the heatmaps
  dataBar$colors = factor(dataBar$colors, levels = c("slateblue1","red"))
  statsOut$colors = factor(statsOut$colors, levels = c("slateblue1","red"))
  
  group.colors <- c(slateblue1 = "slateblue1", red= "red")
  

  figBar=ggplot(dataBar, aes(x=colors, y=ES, fill=colors)) + 
    geom_bar(data=statsOut,position=position_dodge(), stat="identity", colour="black") +
    geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 2)+
    geom_errorbar(data=statsOut,aes(ymin=ES-se, ymax=ES+se),
                  width=.4, size=1.5,                    # Width of the error bars
                  position=position_dodge(.9))+
    geom_crossbar(data=statsOut, aes(ymin = ES, ymax = ES),
                  size=0.5,col="black", width = .5)+
    
    theme(panel.background = element_rect(fill = 'white', colour = 'black'), 
          text = element_text(size=4),
          axis.text.x = element_text(face = "plain", color = "black", 
                                     size = 18, angle=0),
          axis.text.y = element_text(face = "plain", color = "black", 
                                     size = 18),
          axis.title=element_text(size=14,face="plain"),
          legend.position = "none")+
    xlab("")+
    ylab("Enrichment Score")+
    ggtitle(rownames(gsvaSelect)[indGeneSetPlot])+
    scale_fill_manual(values=c("#1A71B9", "#EC2426"))+
    ylim(-0.6,0.6)+
    scale_x_discrete(labels=c("20Hz", "40Hz"))+
    geom_hline(yintercept=0)
  
  pngPath=paste("./AnnotatedFigs/geneSets_Cluster_", indGeneSetPlot,".png",sep='')
  png(pngPath, width=3.5,height=4,units="in",res=600)
  print(figBar)
  dev.off()
  #graphics.off()
  
  pdfPath=paste("./AnnotatedFigs/geneSets_Cluster_", indGeneSetPlot,".pdf",sep='')
  pdf(pdfPath, width=3.5,height=4,pointsize = 12)
  print(figBar)
  dev.off()
  #graphics.off()
  
  
  #plot heatmaps for each gene set 
    #Get index of gene set associated with enrichment
    indGetGS=which(match(names(geneSets),rownames(gsvaSelect)[indGeneSetPlot])==1)
    
    
    indGS_inData=match(as.vector(unlist(geneSets[indGetGS])),rownames(y))
    indGS_inData=indGS_inData[!is.na(indGS_inData)]
  
  
    if (length(indGS_inData)>1){
      pngPath=paste("./AnnotatedFigs/geneSetHeatMap_Cluster_",indGeneSetPlot,".pdf",sep='')
      pdf(pngPath,width=6,height=6)
      #png(pngPath, height=700, width=650, pointsize=12)
      print({
        
        blah=t(apply(y[indGS_inData,], 1, scale)) 
        rownames(blah)=rownames(countData[indGS_inData,]) #use lower case gene symbols for mouse
        blah=blah[!is.na(blah[,1]),]
        
    
        purpleBlueDiff=rowMeans(blah[,colorBar=="red"])-rowMeans(blah[,colorBar=="slateblue1"])
        indBlah=sort(purpleBlueDiff, index.return=TRUE)$ix
       
        colorBar2=colorBar
        colorBar2[colorBar=="slateblue1"]="#1A71B9"
          colorBar2[colorBar=="red"]="#EC2426"
          
        heatmap3(blah[indBlah,], ColSideColors =colorBar2, ColSideWidth=1, ColSideLabs=NA,
                 col=barColors, breaks=breakBarColors,legendfun=function()showLegend(legend=c(NA),col=c(NA),cex=1.5), 
                 Rowv=NA, Colv=NA,  scale="row",
                 labRow = row.names(blah)[indBlah], labCol = NA, cexRow=0.15,
                 main=rownames(gsvaSelect)[indGeneSetPlot]) 
      })
      dev.off
    
    
    geneTable=data.frame("gene"=rownames(blah),blah)
    colnames(geneTable)=c("gene", rep("20",7), rep("40",7))
    
    if (nrow(blah)>20){
      geneTableTranspose=t(blah[indBlah[c(c(1:10),c((nrow(blah)-9):nrow(blah)))],])
    }else{
      geneTableTranspose=t(blah[indBlah,])
    }  
    
    newCol=c(rep("20",7), rep("40",7))
    geneTableTranspose=data.frame("gene"=newCol,geneTableTranspose)
    #colnames(geneTableTranspose)[1]="gene"
    geneTableTransposeMelted=melt(geneTableTranspose, id.vars = "gene")
    
    colnames(geneTableTransposeMelted)=c("variable", "gene", "value")

    statsOut=geneTableTransposeMelted %>% group_by(variable, gene) %>% summarise(se=std.error(value), value=mean(value))

    figBarGenes=ggplot(geneTableTransposeMelted, aes(x=gene, y=value, fill=variable)) + 
       geom_col(data=statsOut,position=position_dodge(0.5), width=0.5, colour="black") +
      #geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 2)+
       geom_errorbar(data=statsOut,aes(ymin=value-se, ymax=value+se),
                   width=.4, size=1,                    # Width of the error bars
                    position=position_dodge(0.45))+
      # geom_crossbar(data=statsOut, aes(ymin = ES, ymax = ES),
      #             size=0.5,col="black", width = .5)+

      theme(panel.background = element_rect(fill = 'white', colour = 'black'), 
            text = element_text(size=4),
            axis.text.x = element_text(face = "plain", color = "black", 
                                       size = 18, angle=90, hjust=1, vjust=0.5),
            axis.text.y = element_text(face = "plain", color = "black", 
                                       size = 18),
            axis.title=element_text(size=14,face="plain"),
            legend.position = "none")+
      xlab("")+
      ylab("z-score")+
      ggtitle(rownames(gsvaSelect)[indGeneSetPlot])+
      scale_fill_manual(values=c("#1A71B9", "#EC2426"))+
      #ylim(-0.6,0.6)+
      #scale_x_discrete(labels=c("20Hz", "40Hz"))+
      geom_hline(yintercept=0)
    
    print({figBarGenes})
    
    pngPath=paste("./AnnotatedFigs/geneSetHeatMap_Bar_",indGeneSetPlot,".png",sep='')
    png(pngPath,width=8,height=3.25,units="in",res=600)
    print({figBarGenes})
    dev.off()
    
    pdfPath=paste("./AnnotatedFigs/geneSetHeatMap_Bar_",indGeneSetPlot,".pdf",sep='')
    pdf(pdfPath,width=8,height=3.25)
    print({figBarGenes})
    dev.off()
    
    inddeseq=match(toupper(rownames(blah[indBlah,])),rownames(res))
    
    
    deseqOutput=res[inddeseq,]
    
    rownames(deseqOutput)=r_ucfirst(rownames(deseqOutput)) #reset the gene symbols to mouse with first letter upper case
    
    write.csv(deseqOutput, paste("./AnnotatedFigs/geneSetDESEQ2_",indGeneSetPlot,"_",gsvaList[indGeneSetPlot],".csv",sep='') )
    
    }
  graphics.off()
}





