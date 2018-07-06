#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-i", "--inFile", type="character", required=TRUE, help="infile")
parser$add_argument("-a", "--sampleDescFile", type="character", required=TRUE, help="sample.desc.file")
parser$add_argument("-b", "--geneDescFile", type="character", required=TRUE, help="gene.desc.file")
parser$add_argument("-o", "--outputPrefix", type="character", required=TRUE, help="output prefix")
parser$add_argument("-s", "--sample", type="character", default="SAMPLE", help="sample id  [default %(default)s]")
parser$add_argument("-y", "--cellTypeColorFile", type="character", help="cellTypeColorFile")
parser$add_argument("-l", "--clonotypeFile", type="character", help="clonotype file")
parser$add_argument("-u", "--u", action="store_true", default=FALSE, help="unique gene [default %(default)s]")
parser$add_argument("-e", "--scale", action="store_true", default=FALSE, help="do scale [default %(default)s]")
args <- parser$parse_args()
print(args)

in.file <- args$inFile
sample.desc.file <- args$sampleDescFile
gene.desc.file <- args$geneDescFile
out.prefix <- args$outputPrefix
sample.id <- args$sample
cellTypeColorFile <- args$cellTypeColorFile
clonotype.file <- args$clonotypeFile
opt.u <- args$u
do.scale <- args$scale

##opt.u <- FALSE
##do.scale <- F
##in.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/OUT.scLVM/P0508/sfIgnoreERCC/P0508.het.countGeneData.sfNormalized"
##sample.desc.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/iterative/P0508.sample.txt"
##gene.desc.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/iterative/P0508.markerGene.slim"
##out.prefix <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/iterative/test.P0508.markerGene.tsne"
##sample.id <- "P0508"
##cellTypeColorFile <- "/WPS1/zhenglt/work/TCR_chunhong/data/CellType.color"
##clonotype.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/tracer/OUT/P0508/filtered_TCR_summary/P0508.summary.cell.reassigneClonotype.methodChunhong.r.txt"

source("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")

suppressPackageStartupMessages(require("ComplexHeatmap"))
suppressPackageStartupMessages(require("circlize"))
suppressPackageStartupMessages(require("gridBase"))
suppressPackageStartupMessages(require("dendextend"))
suppressPackageStartupMessages(require("RColorBrewer"))
suppressPackageStartupMessages(require("gplots"))

### clonotype data
clonotype.data <- read.clonotype(in.file = clonotype.file,ctype.col = "C_strict")
#### exp data
in.table <- read.table(in.file,header = T,sep = "\t",check.names = F,stringsAsFactors = F,row.names = 1)

in.table <- in.table[,-1]
in.table <- as.matrix(log2(in.table+1))

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
gene.desc <- read.table(gene.desc.file,header = F,sep = "\t",check.names = F,stringsAsFactors = F)
colnames(gene.desc) <- c("patient","geneName","markerClass")
gene.desc <- subset(gene.desc,markerClass!="uncharacterized")
print(dim(gene.desc))
print(head(gene.desc))

#### cell type color
sampleTypeColor <- read.SampleTypeColor(cellTypeColorFile)
sampleTypeColor <- sampleTypeColor[intersect(names(sampleTypeColor),unique(sample.desc$sampleType))]
leafClusterColor <- structure(colorRampPalette(brewer.pal(12,"Paired"))(length(unique(sample.desc$leafCluster))),
                              names=unique(sample.desc$leafCluster))
leafClusterColor["uncharacterized"] <- "gray"



tsne.out <- runTSNEAnalysis(in.table,out.prefix,legend=names(leafClusterColor),col.points=leafClusterColor[sample.desc$leafCluster],col.legend=leafClusterColor,pch=16,pch.legend=16,inPDF=TRUE,eps=2.0,dims=2,k=NULL,do.dbscan=F,myseed=NULL,width.pdf = 14,margin.r = 24,legend.inset = -0.55)

tsne.dbscan.out <- runTSNEAnalysis(in.table,sprintf("%s.dbScan",out.prefix),legend=names(sampleTypeColor),col.points=sampleTypeColor[sample.desc$sampleType],col.legend=sampleTypeColor,pch=16,pch.legend=16,inPDF=TRUE,eps=1.45,dims=2,k=NULL,do.dbscan=T,myseed=tsne.out[["30"]]$myseed,width.pdf = 14,margin.r = 24,legend.inset = -0.55)


