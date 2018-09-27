#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-d", "--designFile", type="character", required=TRUE, help="design matrix file")
parser$add_argument("-i", "--inputFile", type="character", required=TRUE, help="input file")
parser$add_argument("-x", "--binFile", type="character", help="binarized exp file")
parser$add_argument("-t", "--tpmFile", type="character", help="tpmFile")
parser$add_argument("-f", "--geneFile", type="character", help="gene list file")
parser$add_argument("-c", "--cellTypeColorFile", default="/WPS1/zhenglt/work/TCR_chunhong/data/CellType.color", 
                    type="character", help="cellTypeColorFile [default %(default)s]")
parser$add_argument("-o", "--outDir", type="character", required=TRUE, help="output directory")
parser$add_argument("-s", "--sample", type="character", default="SAMPLE", help="sample id  [default %(default)s]")
parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="whether verbose mode [default %(default)s]")
parser$add_argument("-l", "--log", action="store_true", default=FALSE, help="whether do log transform [default %(default)s]")
parser$add_argument("-r", "--center", action="store_true", default=FALSE, help="whether center genes by individual [default %(default)s]")
##parser$add_argument("-m", "--notFilter", action="store_true", default=FALSE, help="don't filtering genes [default %(default)s]")
parser$add_argument("-m", "--minCell", type="integer", default=5, help="keep genes which have >INT cells with non-zeror expression value [default %(default)s]")
parser$add_argument("-n", "--nKeep", type="integer", help="use top nKeep genes [default %(default)s]")
#### method specific
parser$add_argument("-q", "--F.FDR", type="double", default=0.05, help="threshold of F.FDR [default %(default)s]")
parser$add_argument("-k", "--ncores", type="integer", default=8, help="num of cors to use [default %(default)s]")
parser$add_argument("-e", "--geneIDType", type="character",default="entrezID", help="geneID type (entrezID, ensemblID) [default %(default)s]")
parser$add_argument("-b", "--normExprs", action="store_true", default=FALSE, help="use norm_exprs assay [default %(default)s]")
args <- parser$parse_args()

print(args)
designFile <- args$designFile
inputFile <- args$inputFile
geneFile <- args$geneFile
cellTypeColorFile <- args$cellTypeColorFile
out.dir <- args$outDir
sample.id <- args$sample
args.log <- args$log
args.center <- args$center
nKeep <- args$nKeep
#args.notFilter <- args$notFilter
args.minCell <- args$minCell
args.geneIDType <- args$geneIDType
args.ncores <- args$ncores
args.norm.exprs <- args$normExprs
binFile <- args$binFile
tpmFile <- args$tpmFile
args.F.FDR <- args$F.FDR
args.verbose <- args$verbose

##### frequently used code
dir.create(out.dir,recursive = T,showWarnings = F)
if(file.exists("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")){
    source("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")
}else if(file.exists("/WPSnew/zhenglt/work/proj_zl/integrate0924/cross.patient.b20180409/cancer/lib/scRNAToolKit.R")){
    source("/WPSnew/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")
}
#suppressPackageStartupMessages(require("ComplexHeatmap"))
#suppressPackageStartupMessages(require("circlize"))
#suppressPackageStartupMessages(require("gridBase"))
#suppressPackageStartupMessages(require("dendextend"))
suppressPackageStartupMessages(library("RColorBrewer"))
suppressPackageStartupMessages(library("gplots"))
suppressPackageStartupMessages(library("plotrix"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("Rmisc"))
suppressPackageStartupMessages(library("R.utils"))


###suppressPackageStartupMessages(library("stringr"))
out.prefix <- sprintf("%s/%s.de.aov",out.dir,sample.id)
g.inputD <- processInput(designFile,cellTypeColorFile,inputFile,geneFile = NULL,
                         args.center=args.center,args.log=args.log,args.norm.exprs = args.norm.exprs,args.notFilter=T)
myDesign <- g.inputD$myDesign
sampleTypeColor <- g.inputD$sampleTypeColor
patient.col.list <- g.inputD$patient.col.list
Y <- g.inputD$Y
#if(!args.notFilter)
{
    f <- apply(Y,1,function(x){ nE <- sum(x>0); return( nE > args.minCell & nE/length(x) > 0.01 )  })
    Y <- Y[f,]
}

g.GNAME <- g.inputD$g.GNAME
i.mcls <- which(colnames(myDesign)=="majorCluster")
##### 

#### rename rownames
Y.newRowname <- Y
newName <- g.GNAME[rownames(Y)]
names(newName) <- rownames(Y)
f.gname.na <- is.na(newName)
newName[f.gname.na] <- names(newName)[f.gname.na]
rownames(Y.newRowname) <- unname(newName)

exp.bin <- NULL
if(!is.null(binFile) && file.exists(binFile)){
    lenv <- loadToEnv(binFile)
    exp.bin <- lenv$exp.bin[rownames(Y),colnames(Y),drop=F]
}

exp.tpm <- NULL
if(!is.null(tpmFile) && file.exists(tpmFile)){
    lenv <- loadToEnv(tpmFile)
    exp.tpm <- log2(lenv$Y[rownames(Y),colnames(Y),drop=F]+1)
}

#if(args.geneIDType=="ensemblID"){ Y <- Y.newRowname }

### get differential expressed genes

clusterColor <- structure(colorRampPalette(brewer.pal(12,"Paired"))(length(unique(myDesign$majorCluster))),
                              names=sort(as.character(unique(myDesign$majorCluster))))

print("str(Y):")
print(str(Y))

clusterMarker.res <- my.clusterMarkerGene(Y,
                                          clust=myDesign$majorCluster,
                                          ann.col=sampleTypeColor[myDesign[,"sampleType"]],
                                          clust.col=clusterColor,
                                          out.prefix=sprintf("%s", out.prefix),
                                          n.cores=args.ncores,
                                          sampleType = myDesign[,"sampleType"],
                                          original.labels=F,
                                          gid.mapping=g.GNAME,
                                          sampleTypeColSet = sampleTypeColor,
                                          exp.bin = exp.bin,
                                          exp.tpm=exp.tpm,
                                          F.FDR.THRESHOLD=args.F.FDR,verbose=args.verbose,minCell=args.minCell)

if(args$verbose){
    save(clusterMarker.res,file=sprintf("%s.RData",out.prefix))
}

