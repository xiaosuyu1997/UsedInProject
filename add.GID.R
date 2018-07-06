#!/usr/bin/env Rscript

source("/Share/BP/zhenglt/02.pipeline/cancer/lib/myFunc.R")
#source("/lustre1/zeminz_pkuhpc/02.pipeline/cancer/lib/myFunc.R")

#args <- commandArgs(T)
#if(length(args)<2)
#{
#	cat("usage: add.GSymbol.R <in.file> <out.file>\n")
#	q()
#}
#in.file <- args[1]
#out.file <- args[2]
in.file <- file("stdin")
out.file <- ""

in.table <- read.table(in.file,header=T,check.names = F,stringsAsFactors = F,row.names = 1,sep="\t")
##gname <- entrezToXXX(rownames(in.table))
geneID <- XXXToEntrez(rownames(in.table))
out.table <- data.frame(geneID=geneID,symbol=rownames(in.table))
rownames(out.table) <- rownames(in.table)
out.table <- cbind(out.table,in.table)
write.table(out.table,out.file,sep="\t",quote = F,row.names = F,col.names = T)

