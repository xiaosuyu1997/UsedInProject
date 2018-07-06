#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-i", "--infile", type="character", required=TRUE, help="infile")
parser$add_argument("-o", "--outprefix", type="character", required=TRUE, help="outprefix")
parser$add_argument("-s", "--sample", type="character", default="SAMPLE", help="sample id  [default %(default)s]")
parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="whether verbose mode [default %(default)s]")
args <- parser$parse_args()
print(args)

in.file <- args$infile
out.prefix <- args$outprefix
sample.id <- args$sample

#in.file <- "/WPS1/zhenglt/work/proj_xy/integrated/cross.patient/exp/lung.P9.scran.RData"
#out.prefix <- "/WPS1/zhenglt/work/proj_xy/integrated/cross.patient/exp/lung.P9.scran.pool.fltGeneT.fltSampleT"
#sample.id <- "lung.P9"

#source("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")
#library("dplyr")
library("scater")
#library("BiocParallel")
library("scran")
#library("plotrix")
#library("statmod")
#library("limSolve")

lname <- load(in.file)

### counts
out.df <- rowData(sce.norm)[,c("geneID","geneSymbol")]
out.df <- cbind(as.data.frame(out.df),counts(sce.norm))
conn <- gzfile(sprintf("%s.counts.tab.gz",out.prefix),open = "w")
write.table(out.df,file = conn,row.names = F,col.names = T,sep = "\t",quote = F)
close(conn)

### tpm
out.df <- rowData(sce.norm)[,c("geneID","geneSymbol")]
out.df <- cbind(as.data.frame(out.df),assay(sce.norm,"tpm"))
conn <- gzfile(sprintf("%s.tpm.tab.gz",out.prefix),open = "w")
write.table(out.df,file = conn,row.names = F,col.names = T,sep = "\t",quote = F)
close(conn)

### norm_exprs (log(count/size_factor+1))
out.df <- rowData(sce.norm)[,c("geneID","geneSymbol")]
out.df <- cbind(as.data.frame(out.df),assay(sce.norm,"norm_exprs"))
conn <- gzfile(sprintf("%s.norm_exprs.tab.gz",out.prefix),open = "w")
write.table(out.df,file = conn,row.names = F,col.names = T,sep = "\t",quote = F)
close(conn)

### exprs (centered by patient, based on norm_exprs)
out.df <- rowData(sce.norm)[,c("geneID","geneSymbol")]
out.df <- cbind(as.data.frame(out.df),assay(sce.norm,"centered_norm_exprs"))
conn <- gzfile(sprintf("%s.exprs.tab.gz",out.prefix),open = "w")
write.table(out.df,file = conn,row.names = F,col.names = T,sep = "\t",quote = F)
close(conn)




