#!/usr/bin/env Rscript

source("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")
suppressPackageStartupMessages(library("scran"))
suppressPackageStartupMessages(library("scater"))

loginfo(" begin ...")

#lname <- load("/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/gsnap_out/quantification/scater/liver.P5.SCESet.RData")
lname <- load("/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/gsnap_out/quantification/scater/liver.P5.SCESet.withCountTPM.RData")
#lname <- load("/WPS1/zhenglt/work/proj_xy/integrated/cross.patient/exp/lung.P6.SCESet.afterClustering.RData")

sce.noMAIT <- sce[,sce$majorCluster!="CD8_Cluster3-SLC4A10"]
#sce.noMAIT <- sce[c("6347","7127","4082","308","920","915","3856"),sce$majorCluster!="CD8_Cluster3-SLC4A10"]

sce.plot <- sce.noMAIT["3164",]
#sce.plot <- sce.noMAIT[,grepl("CD4",sce.noMAIT$leafCluster)]
featureNames(sce.plot)<-featureData(sce.plot)$geneSymbol
apply(tpm(sce.plot),1,summary)
apply(tpm(sce.plot),1,function(x){ sum(x>3)/length(x) })

pdf("/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/cross.patient/sce.NR4A1.noMAIT.pdf",width=10,height=8)
p1 <- plotExpression(sce.plot, "NR4A1", x = "sampleType", show_median=T, exprs_values = "tpm",log2_values = T) + 
    labs(list(y = "log2(TPM+1)"))
p2 <- plotExpression(sce.plot, "NR4A1", x = "majorCluster", show_median=T, exprs_values = "tpm",log2_values = T) +
    theme(axis.text.x  = element_text(angle=45, vjust=0.5)) +
    labs(list(y = "log2(TPM+1)"))
multiplot(p1, p2,cols = 1)
dev.off()



