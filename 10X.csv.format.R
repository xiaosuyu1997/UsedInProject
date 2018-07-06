#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-i", "--inputFile", type="character", required=TRUE, help="input file")
parser$add_argument("-o", "--outPrefix", type="character", required=TRUE, help="output prefix")
parser$add_argument("-s", "--sample", type="character", default="SAMPLE", help="sample id  [default %(default)s]")
parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="whether verbose mode [default %(default)s]")
args <- parser$parse_args()

print(args)
inputFile <- args$inputFile
outPrefix <- args$outPrefix
sample.id <- args$sample

in.file <- inputFile
out.prefix <- outPrefix

#in.file <- "PR10_ZZM_10X-PBMC.csv.gz"
#out.prefix <- "PR10_ZZM_10X-PBMC.csv.formated"

in.table <- read.csv(in.file,header = T,check.names=F,stringsAsFactors=F)
rownames(in.table) <- in.table[,1]
in.table <- as.matrix(in.table[,-1])
colnames(in.table) <- sprintf("%s.%s",sample.id,colnames(in.table))

ID.mapping.table <- read.table("/DBS/DB_temp/zhangLab/ensemble/release88/Homo_sapiens.GRCh38.cdna.all.ID.mapping",header = F,check.names = F,stringsAsFactors = F)
colnames(ID.mapping.table) <- c("tid","gid","gname","chr","beg","end","strand","gene_biotype","transcript_biotype")
ID.mapping.table$gid <- gsub("\\..+","",ID.mapping.table$gid)

in.table[1:4,1:8]

Y <- in.table
g.GNAME <- ID.mapping.table$gname
names(g.GNAME) <- ID.mapping.table$gid
save(Y,g.GNAME,file = sprintf("%s.RData",out.prefix))

out.table <- data.frame(geneID=rownames(in.table),geneSymbol=g.GNAME[rownames(in.table)],stringsAsFactors = F)
out.table <- cbind(out.table,in.table)

conn <- gzfile(sprintf("%s.txt.gz",out.prefix),"w")
write.table(out.table,conn,sep = "\t",row.names = F,quote = F)
close(conn)

