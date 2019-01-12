#!/usr/bin/env Rscript


suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-d", "--designFile", type="character", required=TRUE, help="design matrix file")
parser$add_argument("-i", "--inputFile", type="character", required=TRUE, help="input file")
parser$add_argument("-f", "--geneFile", type="character", help="gene list file")
parser$add_argument("-o", "--outDir", type="character", required=TRUE, help="output directory")
parser$add_argument("-s", "--sample", type="character", default="SAMPLE", help="sample id  [default %(default)s]")
parser$add_argument("-b", "--batch", type="character", help="column in the design file used as batch  [default %(default)s]")
parser$add_argument("-e", "--geneIDType", type="character",default="entrezID",
                    help="geneID type (entrezID, ensemblID) [default %(default)s]")
parser$add_argument("-j", "--measure", type="character", default="norm_exprs", help="measure (tpm, exprs, norm_exprs)[default %(default)s]")
parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="whether verbose mode [default %(default)s]")
parser$add_argument("-l", "--log", action="store_true", default=FALSE, help="whether do log transform [default %(default)s]")
parser$add_argument("-r", "--center", action="store_true", default=FALSE,
                    help="whether center genes by individual [default %(default)s]")
parser$add_argument("-m", "--notFilter", action="store_true", default=FALSE, help="don't filtering genes [default %(default)s]")
parser$add_argument("-a", "--method", type="character", default="limma", help="method (limma, combat, zscore)  [default %(default)s]")

args <- parser$parse_args()

print(args)

designFile <- args$designFile
inputFile <- args$inputFile
geneFile <- args$geneFile
out.dir <- args$outDir
sample.id <- args$sample
args.geneIDType <- args$geneIDType
args.measure <- args$measure
args.log <- args$log
args.center <- args$center
args.notFilter <- args$notFilter
args.method <- args$method
args.batch <- args$batch

print("args.batch:")
#print(args.batch)
print(if(is.null(args.batch)) "NULL" else args.batch)

#designFile <- "/WPS1/zhenglt/work/proj_xy/integrated/clustering.otherMethod/d/lungC.P9.CD8.designUsed.txt"
#inputFile <- "/WPS1/zhenglt/work/proj_xy/integrated/cross.patient/exp/lung.P9.scran.RData"
#geneFile <- NULL
#out.dir <- "/WPS1/zhenglt/work/proj_xy/integrated/clustering.otherMethod/test"
#sample.id <- "colonC.rmIndEff"
#args.geneIDType <- "entrezID"
#args.measure <- "norm_exprs"
#args.log <- F
#args.center <- F
#args.notFilter <- F

##### frequently used code
dir.create(out.dir,recursive = T,showWarnings = F)
out.prefix <- sprintf("%s/%s",out.dir,sample.id)
if(file.exists("/WPSnew/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")){
	source("/WPSnew/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")
}else{
	source("/lustre1/zeminz_pkuhpc/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")
}
suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("dplyr"))
#suppressPackageStartupMessages(library("densityClust"))
suppressPackageStartupMessages(library("limma"))
#suppressPackageStartupMessages(library("sva"))
suppressPackageStartupMessages(library("R.utils"))

g.inputD <- processInput(designFile,cellTypeColorFile=NULL,inputFile=inputFile,args.notFilter,geneFile,args.center,args.log,
                         args.measure=args.measure)
myDesign <- g.inputD$myDesign
sampleTypeColor <- g.inputD$sampleTypeColor
patient.col.list <- g.inputD$patient.col.list
Y <- g.inputD$Y
g.GNAME <- g.inputD$g.GNAME
#### rename rownames
Y.newRowname <- Y
newName <- g.GNAME[rownames(Y)]
names(newName) <- rownames(Y)
f.gname.na <- is.na(newName)
newName[f.gname.na] <- names(newName)[f.gname.na]
rownames(Y.newRowname) <- unname(newName)
#if(args.geneIDType=="ensemblID"){ Y <- Y.newRowname }

##### 

#print(str(Y))
loginfo("Before correction:")
print(range(Y))
print(Y[1:4,1:8])
#print(head(myDesign))
#all(rownames(myDesign)==colnames(Y))

loginfo("After correction:")
if(args.method=="limma"){
    loginfo("use limma")
    Y.cor.limma <- removeBatchEffect(Y,batch=myDesign$patient,batch2=if(is.null(args.batch)) NULL else myDesign[,args.batch])
    print(range(Y.cor.limma))
    print(Y.cor.limma[1:4,1:8])
    Y <- Y.cor.limma
}else if(args.method=="combat"){
    loginfo("use combat")
    pdf(sprintf("%s.diag.pdf",out.prefix),width = 8,height = 6)
    Y.cor.combat = ComBat(dat=Y, batch=myDesign$patient, mod=NULL, par.prior=TRUE, prior.plots=T)
    dev.off()
    print(range(Y.cor.combat))
    print(Y.cor.combat[1:4,1:8])
    Y <- Y.cor.combat
}else if(args.method=="zscore"){
    Y.zscore <- my.centerData(Y, myDesign, do.scale=T)
    print(range(Y.zscore))
    print(Y.zscore[1:4,1:8])
    Y <- Y.zscore
}

save(Y,g.GNAME,file=sprintf("%s.rmIndEff.%s.RData",out.prefix,args.method))

