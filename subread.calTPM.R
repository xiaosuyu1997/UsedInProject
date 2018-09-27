#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-i", "--infile", type="character", required=TRUE, help="infile")
parser$add_argument("-a", "--idMappingFile", type="character",
                    default="/WPSnew/zhenglt/work/proj_fh/data/refdata-cellranger-mm10-1.2.0/genes/ID.mapping",
                    help="id mapping file [default %(default)s]")
parser$add_argument("-s", "--sample", type="character", default="SAMPLE", help="sample [default %(default)s]")
parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="whether verbose mode [default %(default)s]")
args <- parser$parse_args()
print(args)

in.file <- args$infile
id.mapping.file <- args$idMappingFile
sample.id <- args$sample
#in.file <- "/WPSnew/zhenglt/work/proj_fh/byBatch/batch01-04/OUT/B01.C24/STAR/B01.C24.subread.exp.optO"
#id.mapping.file <- "/WPSnew/zhenglt/work/proj_fh/data/refdata-cellranger-mm10-1.2.0/genes/ID.mapping"
#sample.id <- "B01.C24"

library("data.table")

in.table <- fread(in.file)
colnames(in.table)[1] <- "geneID"
colnames(in.table)[ncol(in.table)] <- sample.id

id.mapping.table <- fread(id.mapping.file,header=F)
colnames(id.mapping.table) <- c("geneID","geneSymbol")

o.table <- id.mapping.table[in.table[,c("geneID","Length",sample.id),with=F],,on="geneID"]

#### count
write.table(as.data.frame(o.table[,-c("Length")]),sprintf("%s.RC.txt",in.file),
            row.names = F,sep = "\t",quote = F)

#### TPM
nrc <- o.table[[sample.id]]/o.table[["Length"]]
o.table[[sample.id]] <- nrc*1e6/sum(nrc)
write.table(as.data.frame(o.table[,-c("Length")]),sprintf("%s.TPM.txt",in.file),
            row.names = F,sep = "\t",quote = F)







