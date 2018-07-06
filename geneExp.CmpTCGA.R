#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-a", "--TCGAFile", type="character", required=TRUE, help="TCGAFile")
parser$add_argument("-b", "--myFile", type="character", required=TRUE, help="myFile")
parser$add_argument("-o", "--outprefix", type="character", required=TRUE, help="outprefix")
parser$add_argument("-s", "--sample", type="character", default="SAMPLE", help="sample id  [default %(default)s]")
parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="whether verbose mode [default %(default)s]")
args <- parser$parse_args()
print(args)

exp.TCGA.file <- args$TCGAFile
exp.my.file <- args$myFile
out.prefix <- args$outprefix
sample.id <- args$sample

#source("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")
#library("dplyr")


#exp.TCGA.file <- "/WPSnew/zhenglt/work/TCR_chunhong/TCGA/quantification/mdata/TPM.LUAD.txt.gz"
#exp.my.file <- "/WPS1/zhenglt/work/proj_xy/integrated/bulkRNAseq/stat/lungC.P5.bulkRNA.someGene.txt"
#out.prefix <- "/WPS1/zhenglt/work/proj_xy/integrated/bulkRNAseq/stat/lungC.P5.bulkRNA.cmp.TCGA/lungC.P5"
#sample.id <- "lungC.P5"

### expression data of TCGA
exp.TCGA.table <- read.table(exp.TCGA.file,header = T,check.names = F,stringsAsFactors = F)
rownames(exp.TCGA.table) <- exp.TCGA.table[,1]
g.GNAME.TCGA <- exp.TCGA.table[,2]
names(g.GNAME.TCGA) <- rownames(exp.TCGA.table)
exp.TCGA.table <- as.matrix(exp.TCGA.table[,c(-1,-2)])

### expression data of mine
exp.my.table <- read.table(exp.my.file,header = T,check.names = F,stringsAsFactors = F)
rownames(exp.my.table) <- exp.my.table[,1]
g.GNAME.my <- exp.my.table[,2]
names(g.GNAME.my) <- rownames(exp.my.table)
exp.my.table <- exp.my.table[,c(-1,-2)]

EPSILON <- 0.1
pdf(sprintf("%s.tpm.dist.pdf",out.prefix),width = 6,height = 6)
par(mar=c(5,5,4,2),cex.lab=1.5,cex.main=1.5)
for(gid in rownames(exp.my.table)){
    print(gid)
    v <- log2(exp.TCGA.table[gid,]+EPSILON)
    v.density <- density(v)
    v.density.pp <- pretty(v.density$x)
    hist(v,breaks=20,col="gray",freq=F,xlim=c(v.density.pp[1],v.density.pp[length(v.density.pp)]),border = NA,main=sprintf("%s (%s)",g.GNAME.TCGA[gid],sample.id),xlab="log2(TPM+0.1)")
    lines(v.density,col="orange")
    abline(v = log2(exp.my.table[gid,]+EPSILON),col="red",lwd=1,lty=2)
}
dev.off()

exp.my.table.ecdf <- t(sapply(rownames(exp.my.table),function(x){
    ecdf(exp.TCGA.table[x,])(exp.my.table[x,])
}))
colnames(exp.my.table.ecdf) <- colnames(exp.my.table)
out.df <- data.frame(geneID=rownames(exp.my.table.ecdf))
out.df$geneSymbol <- g.GNAME.my[out.df$geneID]
out.df <- cbind(out.df,exp.my.table.ecdf)
write.table(out.df,sprintf("%s.tpm.ecdf.txt",out.prefix),quote = F,row.names = F,sep = "\t")

