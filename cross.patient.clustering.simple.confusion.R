#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-a", "--aFile", type="character", required=TRUE, help="a file; two columns(sample, densClust)")
parser$add_argument("-b", "--bFile", type="character", required=TRUE, help="b file; contains sample, majorCluster")
parser$add_argument("-o", "--outprefix", type="character", required=TRUE, help="outprefix")
parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="whether verbose mode [default %(default)s]")
args <- parser$parse_args()

print(args)

a.file <- args$aFile
b.file <- args$bFile
out.prefix <- args$outprefix

#a.file <- "OUT.tsne.density.var1500.CD8.v0/L1C1/colonC.CD8.simple.k2.clust.txt"
#b.file <- "/WPSnew/zhenglt/work/proj_zl/integrate0924/cross.patient.b20170523/clusteringResult/colonC.clustering.final.forPlot.allNoMAIT.nGene.txt"

a.table <- read.table(a.file,header = T,check.names = F,sep="\t",stringsAsFactors = F)
b.table <- read.table(b.file,header = T,check.names = F,sep="\t",stringsAsFactors = F)
##b.table <- subset(b.table,stype=="CD8" & invariantTCR=="diverse")
colnames(a.table) <- c("sample","densClust")
library("dplyr")
c.table <- left_join(a.table,b.table)
out <- unclass(table(c.table[,c("densClust","majorCluster")]))
out.df <- cbind(data.frame(densClust=sprintf("C%s",rownames(out)),stringsAsFactors = F),
                as.data.frame(out,stringsAsFactors=F))
write.table(out.df,sprintf("%s.cons.txt",out.prefix),sep="\t",quote = F,row.names = F)
