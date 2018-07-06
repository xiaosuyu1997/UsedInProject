#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-o", "--outDir", type="character", required=TRUE, help="output directory")
parser$add_argument("-s", "--sample", type="character", default="SAMPLE", help="sample id  [default %(default)s]")
parser$add_argument("-c", "--cellTypeColorFile", type="character", required=TRUE, help="cellTypeColorFile")
parser$add_argument("-g", "--variableGeneFile", type="character", required=TRUE, help="variable genes' file")
parser$add_argument("-r", "--rdataFile", type="character", required=TRUE, help="rdata file")
parser$add_argument("-d", "--designFile", type="character", required=TRUE, help="design matrix file")
parser$add_argument("-e", "--outlierFile", type="character", help="outlier file")
##parser$add_argument("-b", "--clonotypeFile", type="character", required=FALSE, help="clonotype file")
parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="whether verbose mode [default %(default)s]")
args <- parser$parse_args()

print(args)

out.dir <- args$outDir
sample.id <- args$sample
cellType.color.file <- args$cellTypeColorFile
variableGene.file <- args$variableGeneFile
rdata.file <- args$rdataFile
design.file <- args$designFile
verboseMode <- args$verbose
outlierFile <- args$outlierFile
#### ------ TEST DATA ------
#out.dir <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/consensusClusterPlus/test"
#sample.id <- "P0205"
#cellType.color.file <- "/WPS1/zhenglt/work/TCR_chunhong/data/CellType.color"
#variableGene.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/P1118.P2/OUT.scLVM/P0205/sfEndo/P0205.fitTechNoise.ERCC.counts.variableGenes.txt"
#rdata.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/clustering/DESeq2.OUT.phase01-29.multiCore.postTCellMarkerFilter/P0205/DESeq2.RData"
#design.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/P1118.P2/sample.design/P0205.design.ERCCFlt.txt"
#verboseMode <- TRUE
#outlierFile <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/consensusClusterPlus/test/P0205.outlier.txt"
#### ------ TEST DATA ------


source("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")

suppressPackageStartupMessages(library("ConsensusClusterPlus"))
suppressPackageStartupMessages(library("R.utils"))
lenv <- loadToEnv(rdata.file)
vsd.blind <- lenv[["vsd.blind"]]
vstMat.blind <- lenv[["vstMat.blind"]]

loginfo("process begin...")
dir.create(out.dir,showWarnings = F,recursive = T)
myDesign<-read.table(design.file,header=T,row.names="sample",check.names=F,colClasses=c("factor","character","factor","factor"))
sampleTypeColor <- read.SampleTypeColor(cellType.color.file)
variableGene.table <- read.table(variableGene.file,header = T,sep = "\t",check.names = F,stringsAsFactors = F)
myGFilter <- as.character(variableGene.table[variableGene.table$is_het,"geneID"])
mySFilter <- NULL
if(!is.null(outlierFile) && file.exists(outlierFile)){
    try({
        mySFilter <- read.table(outlierFile,header = F,check.names = F,stringsAsFactors = F)[,1]
        mySFilter <- rownames(myDesign)[!rownames(myDesign) %in% mySFilter]
    })
}

runConsensusClusterPlus <- function(X,out.dir,ntop=NULL,sfilter=NULL,plotType="pdf"){
    dir.create(out.dir,showWarnings = F,recursive = T)
    if(!is.null(sfilter)) { X <- X[,sfilter,drop=F] }
    if(!is.null(ntop)) { 
        rowVar <- apply(X,1,var)
        select <- order(rowVar,decreasing = TRUE)[seq_len(min(ntop, length(rowVar)))]
        X <- X[select,]
    }
    #rowCV <- apply(X,1,function(v){ sd(v)/mean(v) })
    #f.rowCV <- rowCV > 0.8
    #X <- X[f.rowCV,]
    rowVar <- apply(X,1,var)
    f.rowVar <- rowVar > 2 
    X <- X[f.rowVar,]

    results <- ConsensusClusterPlus(X,title=out.dir,seed=12345,plot=plotType,
                                    maxK=15,reps=100,pItem=0.9,pFeature=1.0,
                                    clusterAlg="hc",distance="pearson",verbose=TRUE)
    icl = calcICL(results,title=out.dir,plot=plotType)
    ret <- list(res=results,icl=icl)
    return(ret)
}

res <- runConsensusClusterPlus(vstMat.blind[myGFilter,],sprintf("%s/withoutSFilter",out.dir))

### outlier in the first round run
#write.table(names(which((res[["res"]][[20]][["consensusClass"]]!=1))),"P0205.outlier.txt",row.names=F,col.names=F,quote=F)
#res.r2 <- runConsensusClusterPlus(vstMat.blind[myGFilter,][1:500,],out.dir,sfilter=mySFilter)
res.r2 <- runConsensusClusterPlus(vstMat.blind[myGFilter,],sprintf("%s/withSFilter",out.dir),sfilter=mySFilter)
#res.r2 <- runConsensusClusterPlus(vstMat.blind[myGFilter,grepl("^TTS",colnames(vstMat.blind))],sprintf("%s/TTR",out.dir),ntop=500)

dir.create(sprintf("%s/hclustering",out.dir),showWarnings = F,recursive = T)
runHierarchicalClusteringAnalysis(vstMat.blind[myGFilter,],sprintf("%s/hclustering/%s.NULL",out.dir,sample.id),myDesign$sampleType,sampleTypeColor,k=1,clonotype.col=NULL,ntop=NULL,complexHeatmap.use=TRUE,verbose=FALSE)
runHierarchicalClusteringAnalysis(vstMat.blind[myGFilter,],sprintf("%s/hclustering/%s.3000",out.dir,sample.id),myDesign$sampleType,sampleTypeColor,k=1,clonotype.col=NULL,ntop=3000,complexHeatmap.use=TRUE,verbose=FALSE)
runHierarchicalClusteringAnalysis(vstMat.blind[myGFilter,],sprintf("%s/hclustering/%s.2500",out.dir,sample.id),myDesign$sampleType,sampleTypeColor,k=1,clonotype.col=NULL,ntop=2500,complexHeatmap.use=TRUE,verbose=FALSE)
runHierarchicalClusteringAnalysis(vstMat.blind[myGFilter,],sprintf("%s/hclustering/%s.2000",out.dir,sample.id),myDesign$sampleType,sampleTypeColor,k=1,clonotype.col=NULL,ntop=2000,complexHeatmap.use=TRUE,verbose=FALSE)
runHierarchicalClusteringAnalysis(vstMat.blind[myGFilter,],sprintf("%s/hclustering/%s.1000",out.dir,sample.id),myDesign$sampleType,sampleTypeColor,k=1,clonotype.col=NULL,ntop=1000,complexHeatmap.use=TRUE,verbose=FALSE)
runHierarchicalClusteringAnalysis(vstMat.blind[myGFilter,],sprintf("%s/hclustering/%s.500",out.dir,sample.id),myDesign$sampleType,sampleTypeColor,k=1,clonotype.col=NULL,ntop=500,complexHeatmap.use=TRUE,verbose=FALSE)

X <- vstMat.blind[myGFilter,]
rowVar <- apply(X,1,var)
select <- order(rowVar,decreasing = TRUE)
X <- X[select,]
runHierarchicalClusteringAnalysis(X[1:25,],sprintf("%s/hclustering/%s.test.1.25",out.dir,sample.id),myDesign$sampleType,sampleTypeColor,k=1,clonotype.col=NULL,ntop=NULL,complexHeatmap.use=TRUE,verbose=FALSE)
runHierarchicalClusteringAnalysis(X[26:50,],sprintf("%s/hclustering/%s.test.26.50",out.dir,sample.id),myDesign$sampleType,sampleTypeColor,k=1,clonotype.col=NULL,ntop=NULL,complexHeatmap.use=TRUE,verbose=FALSE)
runHierarchicalClusteringAnalysis(X[1:50,],sprintf("%s/hclustering/%s.test.1.50",out.dir,sample.id),myDesign$sampleType,sampleTypeColor,k=1,clonotype.col=NULL,ntop=NULL,complexHeatmap.use=TRUE,verbose=FALSE)
runHierarchicalClusteringAnalysis(X[51:100,],sprintf("%s/hclustering/%s.test.51.100",out.dir,sample.id),myDesign$sampleType,sampleTypeColor,k=1,clonotype.col=NULL,ntop=NULL,complexHeatmap.use=TRUE,verbose=FALSE)
runHierarchicalClusteringAnalysis(X[101:200,],sprintf("%s/hclustering/%s.test.101.200",out.dir,sample.id),myDesign$sampleType,sampleTypeColor,k=1,clonotype.col=NULL,ntop=NULL,complexHeatmap.use=TRUE,verbose=FALSE)
runHierarchicalClusteringAnalysis(X[201:300,],sprintf("%s/hclustering/%s.test.201.300",out.dir,sample.id),myDesign$sampleType,sampleTypeColor,k=1,clonotype.col=NULL,ntop=NULL,complexHeatmap.use=TRUE,verbose=FALSE)
runHierarchicalClusteringAnalysis(X[301:400,],sprintf("%s/hclustering/%s.test.301.400",out.dir,sample.id),myDesign$sampleType,sampleTypeColor,k=1,clonotype.col=NULL,ntop=NULL,complexHeatmap.use=TRUE,verbose=FALSE)
runHierarchicalClusteringAnalysis(X[401:500,],sprintf("%s/hclustering/%s.test.401.500",out.dir,sample.id),myDesign$sampleType,sampleTypeColor,k=1,clonotype.col=NULL,ntop=NULL,complexHeatmap.use=TRUE,verbose=FALSE)
runHierarchicalClusteringAnalysis(X[501:1000,],sprintf("%s/hclustering/%s.test.501.1000",out.dir,sample.id),myDesign$sampleType,sampleTypeColor,k=1,clonotype.col=NULL,ntop=NULL,complexHeatmap.use=TRUE,verbose=FALSE)
runHierarchicalClusteringAnalysis(X[1001:1500,],sprintf("%s/hclustering/%s.test.1001.1500",out.dir,sample.id),myDesign$sampleType,sampleTypeColor,k=1,clonotype.col=NULL,ntop=NULL,complexHeatmap.use=TRUE,verbose=FALSE)
runHierarchicalClusteringAnalysis(X[1501:2000,],sprintf("%s/hclustering/%s.test.1501.2000",out.dir,sample.id),myDesign$sampleType,sampleTypeColor,k=1,clonotype.col=NULL,ntop=NULL,complexHeatmap.use=TRUE,verbose=FALSE)
runHierarchicalClusteringAnalysis(X[2001:2500,],sprintf("%s/hclustering/%s.test.2001.2500",out.dir,sample.id),myDesign$sampleType,sampleTypeColor,k=1,clonotype.col=NULL,ntop=NULL,complexHeatmap.use=TRUE,verbose=FALSE)
runHierarchicalClusteringAnalysis(X[2501:3000,],sprintf("%s/hclustering/%s.test.2501.3000",out.dir,sample.id),myDesign$sampleType,sampleTypeColor,k=1,clonotype.col=NULL,ntop=NULL,complexHeatmap.use=TRUE,verbose=FALSE)

save.image(file=sprintf("%s/%s.ConsensusClusterPlus.RData",out.dir,sample.id))

loginfo("process end...")
