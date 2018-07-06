#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-i", "--inFile", type="character", required=TRUE, help="infile")
parser$add_argument("-a", "--sampleDescFile", type="character", required=TRUE, help="sample.desc.file")
parser$add_argument("-b", "--geneFile", type="character", required=TRUE, help="gene list file")
parser$add_argument("-o", "--outputPrefix", type="character", required=TRUE, help="output prefix")
parser$add_argument("-s", "--sample", type="character", default="SAMPLE", help="sample id  [default %(default)s]")
parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="verbose mode [default %(default)s]")
parser$add_argument("-w", "--disableLog", action="store_true", default=FALSE, help="disable log transform [default %(default)s]")
args <- parser$parse_args()
print(args)

in.file <- args$inFile
sample.desc.file <- args$sampleDescFile
gene.file <- args$geneFile
out.prefix <- args$outputPrefix
sample.id <- args$sample
mode.verbose <- args$verbose
disable.log <- args$disableLog

##in.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/gsnap_out/quantification/sizeFactor/P0322/P0322.all.countData.sfNormalized.txt.gz"
##sample.desc.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/byTCellType/clusteringResult/P0322/clustering.P0322.iterative.leaf.txt"
##gene.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/byTCellType//geneList/TH.cytotoxic.LTB"
##out.prefix <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/byTCellType/selectedSamples/P0322.TH.LTB.cytotoxic/P0322.violin.crossLeaves"
##sample.id <- "P0322"
##mode.verbose <- FALSE
##disable.log <- FALSE

suppressPackageStartupMessages(library("vioplot"))
suppressPackageStartupMessages(library("plotrix"))
source("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")
loginfo("begin ...")
in.table <- read.table(in.file,header = T,sep = "\t",check.names = F,stringsAsFactors = F,row.names = 1)
g.GNAME <- in.table[,1]
names(g.GNAME) <- rownames(in.table)
in.table <- as.matrix(in.table[,-1])
if(!disable.log){ in.table <- as.matrix(log2(in.table+1)) }

print(dim(in.table))
print(in.table[1:4,1:8])

#### sample data
sample.desc <- read.table(sample.desc.file,header = T,sep = "\t",check.names = F,stringsAsFactors = F)
rownames(sample.desc) <- sample.desc$sample
print(dim(sample.desc))
print(head(sample.desc))
s.f <- intersect(sample.desc$sample,colnames(in.table))
sample.desc <- sample.desc[s.f,,drop=F]
in.table <- in.table[,s.f,drop=F]
#### gene data
gene.desc <- read.table(gene.file,header = T,sep = "\t",check.names = F,stringsAsFactors = F,colClasses = c("character","character"))
g.f <- intersect(gene.desc$geneID,rownames(in.table))

cls <- unique(subset(sample.desc,leafCluster!="uncharacterized",select="leafCluster",drop=T))
col.cls <- auto.colSet(length(cls),name = "Dark2")
for(j in seq_along(g.f)){
    loginfo(sprintf("plot gene: %s",g.GNAME[g.f[j]]))
    pdf(sprintf("%s.%s.pdf",out.prefix,g.GNAME[g.f[j]]),width = 16,height = 8)
    par(mar=c(8.5,7,4,2),cex.lab=1.8,cex.main=1.8)
    ylim <- pretty(range(in.table[g.f[j],]),n=5)
    plot(0,0,type="n",xlim=c(0.5,length(cls)+0.5),ylim=ylim[c(1,length(ylim))],
         xaxt = 'n', xlab ="",ylab="Expression",main=sprintf("%s",g.GNAME[g.f[j]]))
    for(cl in seq_along(cls)){
        tryCatch(vioplot(na.omit(in.table[g.f[j], subset(sample.desc,leafCluster==cls[cl],select="sample",drop=T)]), 
                         at = cl, add = T, col = col.cls[cl]), 
                 error = function(e) { loginfo(sprintf("vioplot() catch exception: gene %s",g.GNAME[g.f[j]])); e })
    }
    staxlab(1,at = seq_along(cls),labels=cls,srt=30, cex=1.1,adj=1,top.line=2)
    ###mtext(sprintf("FDR: %4.4f",ttest.out$ttest.out.sig[g.list.1[j],"p.adj"]),side = 1,line = -2.1,adj = 0.05,cex=1.2)
    ###mtext(sprintf("FC: %4.4f",ttest.out$ttest.out.sig[g.list.1[j],"fc"]),side = 1,line = -1.1,adj = 0.05,cex=1.2)
    dev.off()
}


