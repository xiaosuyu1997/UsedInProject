#!/usr/bin/env Rscript

args <- commandArgs(T)
if(length(args)<2)
{
    cat("cal.cor.R <in.rdata> <output.prefix>\n")
    q()
}
rd.file <- args[1]
out.prefix <- args[2]

lname <- load(rd.file)
lname

suppressPackageStartupMessages(library("WGCNA"))
suppressPackageStartupMessages(library("DESeq2"))

cal.sample.cor <- function(mat,out.prefix)
{
	expr.cor <- cor(mat)
	diag(expr.cor) <- NA
	expr.cor.summary <- t(apply(expr.cor,1,function(x){ c(mean(x,na.rm = T),median(x,na.rm = T),min(x,na.rm = T),max(x,na.rm = T)) } ))
	colnames(expr.cor.summary) <- c("mean","mediaan","min","max")
	expr.cor.summary <- data.frame(sampleID=rownames(expr.cor.summary),expr.cor.summary)
	print(expr.cor.summary)
	write.table(expr.cor.summary,paste0(out.prefix,".sample.cor.txt"),row.names = F,sep="\t",quote = F)
	pdf(paste0(out.prefix,".sample.cor.pdf"),width = 10,height = 8)
	plot(hclust(as.dist(1-expr.cor),"average"),xlab="",sub = "")
	dev.off()
}

vsd <- varianceStabilizingTransformation(dds, blind=T)
vstMat <- assay(vsd)
cal.sample.cor(vstMat,out.prefix)

