#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-i", "--inFile", type="character", required=TRUE, help="infile")
parser$add_argument("-o", "--outprefix", type="character", required=TRUE, help="outprefix")
parser$add_argument("-s", "--seed", type="integer", help="infile")
parser$add_argument("-m", "--minFreq", default=3, type="integer", help="minFreq")
args <- parser$parse_args()
print(args)

in.file <- args$inFile
myseed <- args$seed
out.prefix <- args$outprefix
min.freq <- args$minFreq

#in.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/iterative/recurrent.diffGene.txt"
#myseed <- "1463363866"
#out.prefix <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/iterative/recurrent.diffGene"

in.table <- read.table(in.file,header = F,check.names = F,stringsAsFactors = F)
colnames(in.table) <- c("geneSymbol","recurrence")
print(head(in.table))

library("wordcloud")
library("RColorBrewer")

if(is.null(myseed)){ myseed <- as.integer(Sys.time()) }
set.seed(myseed)
print(myseed)
pdf(sprintf("%s.wordcloud.pdf",out.prefix),width = 10,height = 8)
wordcloud(in.table$geneSymbol,in.table$recurrence,min.freq = min.freq,random.color = F, 
                    colors=colorRampPalette(rev(brewer.pal(n = 7, name =  "RdBu")))(11) )
dev.off()
