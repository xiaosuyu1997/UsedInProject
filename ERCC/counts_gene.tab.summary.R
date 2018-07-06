#!/usr/bin/env Rscript

args <- commandArgs(T)
if(length(args)<2)
{
    cat(sprintf("counts_gene.tab.summary.R <exp.file> <out.prefix>\n"))
    q()
}
#exp.file <- "P1118.phase23.count.withERCC.tab.gz"
#out.prefix <- "nCount.txt"
exp.file <- args[1]
out.prefix <- args[2]

exp.table <- read.table(exp.file,header = T,row.names = 1,sep = "\t",check.names = F,stringsAsFactors = F)
exp.table <- exp.table[,-1,drop=F]
#print(head(exp.table[,1:8]))

exp.ERCC.table <- exp.table[grepl(pattern = "^ERCC-",x = rownames(exp.table),perl = T),,drop=F]
exp.Gene.table <- exp.table[!grepl(pattern = "^ERCC-",x = rownames(exp.table),perl = T),,drop=F]
#print(head(exp.ERCC.table[,1:8]))
#print(head(exp.Gene.table[,1:8]))
nCountERCC <- apply(exp.ERCC.table,2,sum)
nCountGene <- apply(exp.Gene.table,2,sum)
out <- data.frame(sample=names(nCountERCC),nCountERCC=nCountERCC,nCountGene=nCountGene,nCountTotal=nCountERCC+nCountGene)
write.table(out,out.prefix,sep = "\t",quote = F,row.names = F)
