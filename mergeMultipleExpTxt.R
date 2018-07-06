#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-i", "--inputFile", type="character", required=TRUE, help="input file. format: <id>{tab}<file>")
parser$add_argument("-o", "--outputFile", type="character", required=TRUE, help="output file")
parser$add_argument("-s", "--sample", type="character", default="SAMPLE", help="sample id  [default %(default)s]")
parser$add_argument("-a", "--annCol", type="character", default="1,2", help="annotation column(s) [default %(default)s]")
args <- parser$parse_args()

print(args)

inputFile <- args$inputFile
outputFile <- args$outputFile
sample.id <- args$sample

annCol <- as.numeric(unlist(strsplit(args$annCol,",")))

fTable <- read.table(inputFile,header = F,stringsAsFactors = F,check.names = F)
flist <- fTable$V2
out.table <- NULL
ann <- NULL
for(i in seq_along(flist)){
    in.table <- read.table(flist[i],header = T,sep = "\t",stringsAsFactors = F,check.names = F)
    rownames(in.table) <- in.table[,1]
    if(i==1) { 
        ann <- in.table[,annCol,drop=F] 
        out.table <- in.table[,-annCol,drop=F]
    }else{
        out.table <- cbind(out.table,in.table[,-annCol,drop=F])
    }
    cat(sprintf("%s: (%s, %s)\n",fTable$V1[i],nrow(in.table),ncol(in.table[,-annCol,drop=F])))
}
if(grepl("gz$",outputFile,perl=T)){
    conn <- gzfile(outputFile,"w")
}else{
    conn <- file(outputFile,"w")
}
write.table(cbind(ann,out.table),conn,row.names = F,sep = "\t",quote = F)
close(conn)
cat("----------\n")
cat(sprintf("%s: (%s, %s)\n",outputFile,nrow(out.table),ncol(out.table)))
