#!/usr/bin/env Rscript


args <- commandArgs(T)
if(length(args)<3)
{
    cat("plot.clustering.R <rdata.file> <out.dir> <sample.id>\n")
    q()
}

suppressPackageStartupMessages(library("R.utils"))

rdata.file <- args[1]
out.dir <- args[2]
sample.id <- args[3]

lenv <- loadToEnv(rdata.file)
vsd.blind <- lenv[["vsd.blind"]]
vstMat.blind <- lenv[["vstMat.blind"]]
myDesign <- lenv[["myDesign"]]

source("/lustre1/zeminz_pkuhpc/zhenglt/02.pipeline/cancer/lib/myFunc.R")

print(paste0("args[1]: ",args[1]))
print(paste0("args[2]: ",args[2]))
print(paste0("args[3]: ",args[3]))
qcAndViz(vsd.blind,vstMat.blind,myDesign,out.dir,extra=paste0(".",sample.id))


