#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))
parser <- ArgumentParser()
parser$add_argument("-i", "--inFile", type="character", required=TRUE, help="infile")
parser$add_argument("-d", "--designFile", type="character", help="sample desc file")
parser$add_argument("-o", "--outprefix", type="character",required=TRUE, help="output prefix")
parser$add_argument("-s", "--sample", type="character",default="SAMPLE", help="sample id [default %(default)s]")
##parser$add_argument("-a", "--expThreshold", type="double",default="1", help="expression threshold [default %(default)s]")
parser$add_argument("-n", "--nCores", type="integer",default="8", help="ncores [default %(default)s]")
args <- parser$parse_args()
print(args)

library("mclust")
library("scater")
library("scran")
library("R.utils")
library("data.table")
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("doParallel"))

if(file.exists("/WPSnew/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R"))
{
	source("/WPSnew/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")
}else if(file.exists("/WPSnew/zhenglt/work/proj_zl/integrate0924/cross.patient.b20180409/cancer/lib/scRNAToolKit.R")){
    source("/WPSnew/zhenglt/work/proj_zl/integrate0924/cross.patient.b20180409/cancer/lib/scRNAToolKit.R")
}else if(file.exists("/lustre1/zeminz_pkuhpc/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")){
	source("/lustre1/zeminz_pkuhpc/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")
}

exp.file <- args$inFile
d.file <- args$designFile
out.prefix <- args$outprefix
sample.id <- args$sample
##exp.threshold <- args$expThreshold
args.ncores <- args$nCores

#args.ncores <- 8
#exp.file <- "/lustre1/zeminz_pkuhpc/zhenglt/work/proj_zl/integrate0924/cross.patient.b20170523/exp/colonC.P8.scran.RData"
#out.prefix <- "/lustre1/zeminz_pkuhpc/zhenglt/work/proj_zl/integrate0924/cross.patient.b20170523/exp/colonC.P8.exp.bin"
#d.file <- "/WPSnew/zhenglt/work/proj_zl/integrate0924/cross.patient.b20170926/exp/colonC.P10.scran.fltGeneT.fltSampleT.design.txt"

myDesign <- NULL
if(!is.null(d.file) && file.exists(d.file)){
    myDesign <- fread(d.file,key="sample")
}

if(grepl("RData$",exp.file,perl = T)){
    lenv <- loadToEnv(exp.file)
    if("sce.norm" %in% names(lenv)){
        exp.table <- tpm(lenv[["sce.norm"]])
        ###g.GNAME <- fData(sce.norm)$display.name
        g.GNAME <- rowData(lenv[["sce.norm"]])$display.name
        ###names(g.GNAME) <- fData(sce.norm)$geneID
        names(g.GNAME) <- rowData(lenv[["sce.norm"]])$geneID
    }else if("Y" %in% names(lenv)){
        #### Y is data.table, first two column geneID and geneSymbol
        exp.table <- lenv$Y
        g.GNAME <- exp.table[,geneSymbol]
        names(g.GNAME) <- exp.table[,geneID]
        f.na <- is.na(g.GNAME)
        g.GNAME[f.na] <- names(g.GNAME)[f.na]
        exp.table <- as.matrix(exp.table[,c(-1,-2)])
        rownames(exp.table) <- names(g.GNAME)
        rm(lenv)
    }else{
        stop(sprintf("no expression data found!!!"))
    }
}else{
##}else if(grepl("(.txt|.txt.gz)$",exp.file,perl = T)){
    exp.table <- read.table(exp.file,header = T,check.names = F,stringsAsFactors = F,row.names = 1)
    g.GNAME <- exp.table[,1]
    names(g.GNAME) <- rownames(exp.table)
    f.na <- is.na(g.GNAME)
    g.GNAME[f.na] <- names(g.GNAME)[f.na]
    exp.table <- as.matrix(exp.table[,-1,drop=F])
}

if(!is.null(myDesign)){
    f.sample <- intersect(colnames(exp.table),myDesign[,sample])
    exp.table <- exp.table[,f.sample]
}

#oo <- binarizedExp(log2(exp.table["50943",]+0.01),ofile="./IL17F.bin.png",G=2:3,e.TH=log2(5+0.01),e.name="FOXP3",verbose=F)
#oo2 <- binarizedExp(log2(exp.table["3605",]+0.01),ofile="./IL17A.bin.png",G=2:3,e.TH=log2(30+0.01),e.name="IL17A",verbose=T)

registerDoParallel(cores = args.ncores)
exp.bin <- ldply(rownames(exp.table),function(v){ 
			 ##binarizedExp(log2(exp.table[v,]+0.01),G=2:3,e.name=v,verbose = F)
			 binarizedExp(log2(exp.table[v,]),G=2:3,e.name=v,verbose = F,zero.as.low=F)
			},.progress = "none",.parallel=T)
rownames(exp.bin) <- rownames(exp.table)
exp.bin <- as.matrix(exp.bin)
save(exp.bin,g.GNAME,file=sprintf("%s.RData",out.prefix))

out.stat <- apply(exp.bin,2,function(x){ sum(x>0)  })
out.stat.df <- data.frame(sample=names(out.stat),nGene=out.stat,stringsAsFactors = F)
write.table(out.stat.df,sprintf("%s.stat.txt",out.prefix),row.names = F,sep = "\t",quote = F)


