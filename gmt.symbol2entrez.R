#!/usr/bin/env Rscript

##source("/Share/BP/zhenglt/02.pipeline/cancer/lib/myFunc.R")
source("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")

#args <- commandArgs(T)
#if(length(args)<2)
#{
#	cat("usage: add.GSymbol.R <in.file> <out.file>\n")
#    q()
#}
#in.file <- args[1]
#out.file <- args[2]
in.file <- file("stdin")
out.file <- ""

library("org.Hs.eg.db")
#no_col <- max(count.fields(in.file, sep = "\t"))
#in.table <- read.table(in.file,sep="\t",fill=TRUE,col.names=1:no_col,check.names = F,stringsAsFactors = F,header = F)

id.mapping.list <- unlist(as.list(org.Hs.egSYMBOL2EG))

gmt.dat <- readGMT(in.file)

gset.gid <- sapply(gmt.dat$gSet,function(x){ id.mapping.list[x] })

out.dat <- sapply(seq_along(gset.gid),function(i){ 
                      x <- gset.gid[[i]]
                      paste0(c(names(gset.gid)[i],gmt.dat$gLink[i],x),collapse="\t") })

writeLines(out.dat)


