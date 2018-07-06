#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-o", "--outputPrefix", type="character", required=TRUE, help="output prefix")
parser$add_argument("-r", "--resListFile", type="character", required=TRUE, help="res.list.file")
parser$add_argument("-p", "--pattern", type="character", default="^KEGG_", help="pattern [default %(default)s]")
parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="verbose mode [default %(default)s]")
##parser$add_argument("-s", "--sample", type="character", default="SAMPLE", help="sample id  [default %(default)s]")
###parser$add_argument("-d", "--pdfWidth", type="integer", default=10, help="pdf width [default %(default)s]")
###parser$add_argument("-u", "--pdfHeight", type="integer", default=8, help="pdf heihgt [default %(default)s]")
args <- parser$parse_args()
print(args)

out.prefix <- args$outputPrefix
###plist.file <- args$plistFile
res.list.file <- args$resListFile
###clusterColorFile <- args$clusterColorFile
pattern.str <- args$pattern

###source("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")
###suppressPackageStartupMessages(library("ComplexHeatmap"))
###suppressPackageStartupMessages(library("circlize"))
###suppressPackageStartupMessages(library("gridBase"))
###suppressPackageStartupMessages(library("dendextend"))
###suppressPackageStartupMessages(library("RColorBrewer"))
###suppressPackageStartupMessages(library("gplots"))

#out.prefix <- "/WPS1/zhenglt/work/proj_xy/integrated/cross.patient.rename/pathway_list/CD8.KEGG"
###plist.file <- "/WPS1/zhenglt/work/proj_xy/integrated/cross.patient.rename/pathway_list/CD8.KEGG.plist"
#res.list.file <- "/WPS1/zhenglt/work/proj_xy/integrated/cross.patient.rename/pathway_list/CD8.cp.res.list"
###clusterColorFile <- "/WPS1/zhenglt/work/proj_xy/integrated/cross.patient.rename/curated.clusterColor.txt"
#pattern.str <- "KEGG_"

#plist <- read.table(plist.file,header = F,stringsAsFactors = F,check.names = F)$V1
res.list <- read.table(res.list.file,header = F,stringsAsFactors = F,check.names = F)
for(i in seq_len(nrow(res.list))){
    .i.table <- read.table(res.list[i,2],header = T,stringsAsFactors = F,check.names = F,sep="\t")
    rownames(.i.table) <- .i.table[,1]
    .i.table <- .i.table[grepl(pattern.str,rownames(.i.table),perl = T),]
    .i.table[,"p.adj"] <- p.adjust(.i.table[,"p.value"],method = "BH")
    .i.table <- .i.table[order(.i.table[,"p.adj"],.i.table[,"p.value"]),]
    write.table(.i.table,file = sprintf("%s.%s.txt",out.prefix,res.list[i,1]),sep = "\t",row.names = F,col.names = T,quote = F)
}
out.list <- t(sapply(seq_len(nrow(res.list)),function(i){
    c(res.list[i,1],sprintf("%s.%s.txt",out.prefix,res.list[i,1]))
}))
write.table(out.list,file = sprintf("%s.list",out.prefix),sep = "\t",row.names = F,col.names = F,quote = F)









