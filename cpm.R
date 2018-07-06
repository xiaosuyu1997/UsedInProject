#!/usr/bin/env Rscript


suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-d", "--designFile", type="character", required=TRUE, help="design matrix file")
parser$add_argument("-i", "--inputFile", type="character", required=TRUE, help="input file")
parser$add_argument("-f", "--geneFile", type="character", help="gene list file")
parser$add_argument("-o", "--outDir", type="character", required=TRUE, help="output directory")
parser$add_argument("-s", "--sample", type="character", default="SAMPLE", help="sample id  [default %(default)s]")
parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="whether verbose mode [default %(default)s]")

args <- parser$parse_args()

print(args)

designFile <- args$designFile
inputFile <- args$inputFile
geneFile <- args$geneFile
out.dir <- args$outDir
sample.id <- args$sample

#designFile <- "/WPSnew/zhenglt/work/proj_zl/integrate0924/cross.patient.b20170926/exp/colonC.P10.scran.fltGeneT.fltSampleT.design.txt"
#inputFile <- "/WPSnew/zhenglt/work/proj_zl/integrate0924/cross.patient.b20170926/exp/colonC.P10.countData.txt.gz"
#geneFile <- "/WPSnew/zhenglt/work/proj_zl/integrate0924/cross.patient.b20170926/exp/colonC.P10.scran.pool.fltGeneT.fltSampleT.MethodAVE.RC1.glist"
#out.dir <- "/WPSnew/zhenglt/work/proj_zl/integrate0924/cross.patient.b20170926/exp/cpm"
#sample.id <- "colonC.P10"

##### frequently used code
dir.create(out.dir,recursive = T,showWarnings = F)
out.prefix <- sprintf("%s/%s",out.dir,sample.id)
if(file.exists("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")){
	source("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")
}else{
	source("/lustre1/zeminz_pkuhpc/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")
}
library("data.table")
#suppressPackageStartupMessages(library("stringr"))
#suppressPackageStartupMessages(library("dplyr"))
#suppressPackageStartupMessages(library("limma"))
#suppressPackageStartupMessages(library("sva"))
#suppressPackageStartupMessages(library("R.utils"))

d.table <- fread(designFile)
in.table <- fread(sprintf("gzip -cd %s ",inputFile))
gene.table <- fread(geneFile)
gene.table[,geneID:=as.character(geneID)]
in.table <- cbind(in.table[,c(1,2)],in.table[,d.table$sample,with=F])
in.table[,geneID:=as.character(geneID)]
in.table <- in.table[gene.table[,.(geneID)],,on="geneID"]
in.table[1:4,1:8]
g.GNAME <- in.table[,symbol]
names(g.GNAME) <- in.table[,geneID]
f.na <- is.na(g.GNAME)
g.GNAME[f.na] <- names(g.GNAME)[f.na]
cdat.mtx <- as.matrix(in.table[,c(-1,-2)])
rownames(cdat.mtx) <- in.table[,geneID]
cdat.mtx[1:4,1:8]
##### 

#### cpm
cpm.mtx <- apply(cdat.mtx,2,function(x){ log2(1+x*1e6/sum(x)) })
cpm.mtx[1:4,1:8]
#cpm.mtx.2 <- apply(cdat.mtx,2,function(x){ log2(1+x*1e5/sum(x)) })
#cpm.mtx.2[1:4,1:8]
Y <- cpm.mtx
save(Y,g.GNAME,file=sprintf("%s.cpm.RData",out.prefix))

#### centered by patient
cpm.centered <- my.centerData(Y, as.data.frame(d.table), do.scale=F)
cpm.centered[1:4,1:8]
Y <- cpm.centered
save(Y,g.GNAME,file=sprintf("%s.cpm.centered.RData",out.prefix))







