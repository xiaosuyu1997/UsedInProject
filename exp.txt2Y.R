#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-i", "--inFile", type="character", required=TRUE, help="infile")
parser$add_argument("-o", "--outputPrefix", type="character", required=TRUE, help="output prefix")
args <- parser$parse_args()
print(args)

in.file <- args$inFile
out.prefix <- args$outputPrefix

#in.file <- "/WPS1/zhenglt/work/proj_xy/integrated/cross.patient/exp/lung.P6.TPM.txt.gz"
#out.prefix <- "/WPS1/zhenglt/work/proj_xy/integrated/cross.patient/exp/lung.P6.TPM.RData"

Y <- read.table(in.file,header = T,sep = "\t",check.names = F,stringsAsFactors = F)
rownames(Y) <- Y[,1]
Y <- as.matrix(Y[,c(-1,-2)])
Y <- log2(Y+1)
save(Y,file = sprintf("%s.Y.RData",out.prefix))
