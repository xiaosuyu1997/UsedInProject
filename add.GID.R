#!/usr/bin/env Rscript

source("/WPSnew/zhenglt/02.pipeline/cancer/lib/myFunc.R")
#source("/lustre1/zeminz_pkuhpc/02.pipeline/cancer/lib/myFunc.R")


suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-a", "--geneType", type="character", default="SYMBOL", help="one of SYMBOL, ENSG [default %(default)s]")
parser$add_argument("-b", "--species", type="character", default="human", help="one of human, mouse [default %(default)s]")
parser$add_argument("-c", "--column", type="integer", default=1, help="column index for gene [default %(default)s]")
parser$add_argument("-i", "--inFile", type="character", required=TRUE, help="input file; - for stdin")
parser$add_argument("-o", "--outFile", type="character", required=TRUE, help="output file; - for stdout")
parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="whether verbose mode [default %(default)s]")
args <- parser$parse_args()

#print(args)
args.geneType <- args$geneType
args.species <- args$species
args.column <- args$column
args.inFile <- args$inFile
args.outFile <- args$outFile

#args <- commandArgs(T)
#if(length(args)<2)
#{
#	cat("usage: add.GSymbol.R <in.file> <out.file>\n")
#	q()
#}
#in.file <- args[1]
#out.file <- args[2]
if(args.inFile=="-"){
    in.file <- file("stdin")
}else{
    in.file <- file(args.inFile)
}
if(args.outFile=="-"){
    out.file <- ""
}else{
    out.file <- args.outFile
}

in.table <- read.delim(in.file,header=T,check.names = F,stringsAsFactors = F,row.names = args.column,sep="\t")
##gname <- entrezToXXX(rownames(in.table))
geneID <- XXXToEntrez(rownames(in.table),type=args.geneType,species=args.species)
out.table <- data.frame(geneID=geneID,symbol=rownames(in.table))
rownames(out.table) <- rownames(in.table)
out.table <- cbind(out.table,in.table)
write.table(out.table,out.file,sep="\t",quote = F,row.names = F,col.names = T)

if(args.inFile!="-"){
    close(in.file)
}
if(args.outFile!="-"){
    close(out.file)
}

