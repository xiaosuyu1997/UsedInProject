#!/usr/bin/env Rscript

args <- commandArgs(T)
if(length(args)<4)
{
    cat("usage: getCEFFromRData.R <out.dir> <sample.id> <in.rdata.file> <het.gene.file/N>\n")
    q()
}

out.dir <- args[1]
sample.id <- args[2]
in.rdata.file <- args[3]
arg4 <- args[4]

source("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")
suppressPackageStartupMessages(library("R.utils"))
lenv <- loadToEnv(in.rdata.file)
vstMat.blind <- lenv[["vstMat.blind"]]
myDesign <- data.frame(patient=lenv[["myDesign"]]$patient,sample=rownames(lenv[["myDesign"]]),libType=lenv[["myDesign"]]$libType,sampleType=lenv[["myDesign"]]$sampleType,stringsAsFactors = F)

if(file.exists(arg4))
{
    het.gene.file <- arg4
    het.gene.vec <- as.character(read.table(het.gene.file,header = T)$x)
    N <- length(het.gene.vec)
}else
{
    N <- as.numeric(arg4)
    rowVar <- apply(vstMat.blind,1,var)
    het.gene.vec <- rownames(vstMat.blind)[order(rowVar,decreasing = TRUE)[seq_len(min(N, length(rowVar)))]]
}
#out.dir <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/BackSPIN/OUT/P1118"
#sample.id <- "P1118"
#in.rdata.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/clustering/DESeq2.OUT.phase01-29.multiCore.postTCellMarkerFilter/P1118/DESeq2.RData"
#het.gene.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/P1118.P2/OUT.scLVM/P1118/P1118.var.geneID.txt"


dir.create(out.dir,recursive = T,showWarnings = F)
out.rowann <- data.frame(geneID=het.gene.vec,gname=entrezToXXX(het.gene.vec),hh="")
out.colann <- t(data.frame(cellType=myDesign$sampleType,sample=myDesign$sample))
out.mat <- vstMat.blind[het.gene.vec,]
##write.table(out.df,sprintf("%s/%s.vst.txt",out.dir,sample.id),sep = "\t",row.names = F,quote = F)
##write.table(myDesign,sprintf("%s/%s.design.txt",out.dir,sample.id),sep = "\t",row.names = F,quote = F)

out.con <- file(sprintf("%s/%s.cef",out.dir,sample.id),"w")
cat(sprintf("CEF\t1\t%d\t%d\t%d\t%d\t0\n",ncol(out.rowann)-1,nrow(out.colann),nrow(out.mat),ncol(out.mat)),file=out.con)
cat(sprintf("Genome\thg19\n"),file=out.con)
cat(sprintf("\t\t"),file=out.con)
write.table(out.colann["cellType",,drop=F],file=out.con,quote = F,sep = "\t",row.names = T,col.names = F)
cat(sprintf("\t\t"),file=out.con)
write.table(out.colann["sample",,drop=F],file=out.con,quote = F,sep = "\t",row.names = T,col.names = F)
write.table(t(names(out.rowann)[-ncol(out.rowann)]),file=out.con,quote = F,sep = "\t",row.names = F,col.names = F)
write.table(cbind(out.rowann,signif(out.mat,8)),file=out.con,quote = F,sep = "\t",row.names = F,col.names = F)
close(out.con)
