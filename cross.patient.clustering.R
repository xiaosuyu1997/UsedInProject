#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-d", "--designFile", type="character", required=TRUE, help="design matrix file")
parser$add_argument("-i", "--inputFile", type="character", required=TRUE, help="input file")
parser$add_argument("-f", "--geneFile", type="character", help="gene list file")
parser$add_argument("-c", "--cellTypeColorFile", default="/WPS1/zhenglt/work/TCR_chunhong/data/CellType.color", 
                    type="character", help="cellTypeColorFile [default %(default)s]")
parser$add_argument("-e", "--clonotypeFile", type="character", help="clonotype file")
parser$add_argument("-o", "--outDir", type="character", required=TRUE, help="output directory")
parser$add_argument("-s", "--sample", type="character", default="SAMPLE", help="sample id  [default %(default)s]")
parser$add_argument("-a", "--geneIDType", type="character", default="entrez", help="entrez, ensemblID  [default %(default)s]")
###parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="whether verbose mode [default %(default)s]")
parser$add_argument("-l", "--log", action="store_true", default=FALSE, help="whether do log transform [default %(default)s]")
parser$add_argument("-r", "--center", action="store_true", default=FALSE, help="whether center genes by individual [default %(default)s]")
parser$add_argument("-y", "--onlyY", action="store_true", default=FALSE, help="only center the data [default %(default)s]")
parser$add_argument("-k", "--notSC3", action="store_true", default=FALSE, help="don't run SC3 pipeline [default %(default)s]")
parser$add_argument("-m", "--notFilter", action="store_true", default=FALSE, help="don't filtering genes [default %(default)s]")
parser$add_argument("-n", "--nKeep", type="integer", default=1500, help="use top nKeep genes [default %(default)s]")
args <- parser$parse_args()

print(args)

designFile <- args$designFile
inputFile <- args$inputFile
geneFile <- args$geneFile
cellTypeColorFile <- args$cellTypeColorFile
clonotype.file <- args$clonotypeFile
out.dir <- args$outDir
sample.id <- args$sample
args.log <- args$log
args.center <- args$center
nKeep <- args$nKeep
args.notSC3 <- args$notSC3
args.notFilter <- args$notFilter
args.geneIDType  <- args$geneIDType

#### TEST DATA 
###designFile <- "/WPSnew/zhenglt/work/proj_zl/integrate0924/cross.patient/OUT.CD8/colonC.CD8.designUsed.txt"
###inputFile <- "/WPSnew/zhenglt/work/proj_zl/integrate0924/cross.patient/exp/colonC.P5.scran.RData"
###geneFile <- NULL
###cellTypeColorFile <- "/WPS1/zhenglt/work/TCR_chunhong/data/CellType.color"
###clonotype.file <- NULL
###out.dir <- "/WPSnew/zhenglt/work/proj_zl/integrate0924/cross.patient/test"
###sample.id <- "colon"
###args.log <- F
###args.center <- F
###nKeep <- 1500
###args.notSC3 <- T
###args.notFilter <- F

dir.create(out.dir,recursive = T,showWarnings = F)
source("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")

g.inputD <- processInput(designFile,cellTypeColorFile,inputFile,args.notFilter,geneFile,args.center,args.log)
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
if(args.geneIDType=="ensemblID"){ Y <- Y.newRowname }


loginfo("... all samples.")

doit <- function(dat.plot,extra="",myseed=NULL){
    pca.res <- runPCAAnalysis(dat.plot,sprintf("%s/%s%s.het.PCA",out.dir,sample.id,extra),
                   myDesign[,"sampleType"],sampleTypeColor[names(sampleTypeColor) %in% unique(as.character(myDesign$sampleType))],
                   ntop=NULL,main=sample.id,gid.mapping=g.GNAME)

    tsne.res <- runTSNEAnalysis(dat.plot,sprintf("%s/%s%s.het.tSNE",out.dir,sample.id,extra),
                    col.points = sampleTypeColor[as.character(myDesign$sampleType)],
                    legend=c(names(sampleTypeColor)[names(sampleTypeColor) %in% unique(as.character(myDesign$sampleType))],
                             unique(myDesign$libType)),
                    col.legend=c(sampleTypeColor[names(sampleTypeColor) %in% unique(as.character(myDesign$sampleType))],
                                 rep("black",length(unique(myDesign$libType)))),
                    pch=(as.numeric(as.factor(myDesign$libType))-1) %% 26,
                    pch.legend=c(rep(20,sum(names(sampleTypeColor) %in% unique(as.character(myDesign$sampleType)))),
                                 (seq_along(unique(myDesign$libType))-1) %% 26),
                    myseed=myseed
                    )

    ##for(nn in unique(c(50,100,150,200,250,300,350,400,450,500,1000,nKeep,2000,100000)))
    for(nn in unique(c(50,100)))
    {
    runHierarchicalClusteringAnalysis(dat.plot,
                    sprintf("%s/%s%s.het.hclustering.n%s",out.dir,sample.id,extra,ifelse(nn>99999,"All",nn)),
                    myDesign[,"sampleType"],
                    sampleTypeColor[names(sampleTypeColor) %in% unique(as.character(myDesign$sampleType))],
                    clonotype.col=NULL,
                    patient.col.list = patient.col.list,
                    ntop=nn,
                    complexHeatmap.use=TRUE,do.cuttree=T,
                    verbose=TRUE,main="Variable Genes",gid.mapping=g.GNAME)
    }
    if(!args.notSC3){
        sc3.res <- runSC3Analysis(dat.plot,sprintf("%s/%s%s.SC3",out.dir,sample.id,extra),
                       myDesign[sname,"sampleType"],
                       sampleTypeColor[names(sampleTypeColor) %in% unique(as.character(myDesign[sname,"sampleType"]))],
                       do.log.scale=FALSE,n.cores=10,gid.mapping=g.GNAME,data.forDE=Y)
    }
    ###save(sc3.res,file=sprintf("%s/%s.SC3.RData",out.dir,sample.id))
}

Y.sd <- apply(Y, 1, sd)
Y.sd.sort <- sort(Y.sd,decreasing=TRUE)
Y.sd.df <- data.frame(geneID=names(Y.sd.sort),geneSymbol=g.GNAME[names(Y.sd.sort)])
Y.sd.df <- cbind(Y.sd.df,Y.sd.sort)
write.table(Y.sd.df,sprintf("%s/%s.Y.sd.sort.txt",out.dir,sample.id),row.names = F,sep = "\t",quote = F)

g.f <- head(names(Y.sd.sort),n=nKeep)
doit(Y[g.f,],extra=".centered")
###doit(Y[names(Y.sd.sort)[1:nKeep],],extra=".centered.correct",myseed = 1471060113)
###doit(Y[names(Y.sd.sort)[1:nKeep],],extra=".centered.problem",myseed = 1471060113)

#save.image(sprintf("%s/%s.all.RData",out.dir,sample.id))

