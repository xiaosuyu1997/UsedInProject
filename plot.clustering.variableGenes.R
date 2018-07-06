#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-i", "--inputFile", type="character", required=TRUE, help="input file in these two format: RData contain vstMat.blind obj; TPM in txt format")
parser$add_argument("-d", "--designFile", type="character", required=TRUE, help="design matrix file;required if --inputFile is in txt format")
parser$add_argument("-c", "--cellTypeColorFile", type="character", required=TRUE, help="cell type color file")
parser$add_argument("-o", "--outDir", type="character", required=TRUE, help="output directory")
parser$add_argument("-s", "--sampleID", type="character",default="SAMPLE", required=TRUE, help="sample id [default %(default)s]")
parser$add_argument("-l", "--clonotypeFile", type="character", help="clonotype file")
parser$add_argument("-e", "--variableGeneFile", type="character", help="variable gene file")

args <- parser$parse_args()
print(args)
rdata.file <- args$inputFile
cellType.color.file  <- args$cellTypeColorFile
out.dir <- args$outDir
sample.id <- args$sampleID
design.file <- args$designFile
clonotype.file <- args$clonotypeFile
variable.gene.file <- args$variableGeneFile

suppressPackageStartupMessages(library("R.utils"))
myDesign<-read.table(design.file,header=T,row.names="sample",check.names=F,colClasses=c("factor","character","factor","factor"))

###-------------  TEST -------------    
#rdata.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/clustering/DESeq2.OUT.phase01-28.multiCore.postTCellMarkerFilter/P1202/DESeq2.RData"
#out.dir <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/clustering/test/P1202"
#sample.id  <- "P1202"
#myDesign <- read.table("/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/clustering/marker/P1202.filter.by.marker.design.txt",header=T,row.names="sample",check.names=F,colClasses=c("factor","character","factor","factor"))
#clonotype.file <- NULL
#cellType.color.file  <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/clustering/CellType.color"
#rdata.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/clustering/DESeq2.OUT.phase01-21.multiCore.postTCellMarkerFilter/P0205/DESeq2.RData"
#out.dir <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/clustering/test/P0205"
#sample.id <- "P0205"
#myDesign <- read.table("/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/clustering/marker/P0205.filter.by.marker.design.txt",header=T,row.names="sample",check.names=F,colClasses=c("factor","character","factor","factor"))
#clonotype.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/tracer/OUT/P0205/filtered_TCR_summary/P0205.summary.cell.reassigneClonotype.txt"
###-------------  TEST -------------    

lenv <- loadToEnv(rdata.file)
vsd.blind <- lenv[["vsd.blind"]]
vstMat.blind <- lenv[["vstMat.blind"]]

#source("/Share/BP/zhenglt/02.pipeline/cancer/lib/myFunc.R")
source("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")
#source("/lustre1/zeminz_pkuhpc/zhenglt/02.pipeline/cancer/lib/myFunc.R")

loginfo("process begin...")

dir.create(out.dir,showWarnings = F,recursive = T)
sampleTypeColor <- read.SampleTypeColor(cellType.color.file)
#####qcAndVizMat(vstMat.blind,myDesign,out.dir,colSet=sampleTypeColor,intgroup=c("sampleType"),ntop=500,extra=paste0(".",sample.id),sfilter=NULL, gfilter=NULL,complexHeatmap.use=TRUE)
loginfo("... all samples.")
#qcAndVizMat(vstMat.blind,myDesign,out.dir,colSet=sampleTypeColor,intgroup=c("sampleType"),ntop=2000,extra=paste0(".",sample.id),sfilter=NULL, gfilter=NULL)
clonotype.share.data <- read.clonotype(in.file = clonotype.file,ctype.col = "C_share")
clonotype.strict.data <- read.clonotype(in.file = clonotype.file,ctype.col = "C_strict")
vst.data.used <- vstMat.blind
if(!is.null(variable.gene.file))
{
    variable.gene.list <- as.character(read.table(variable.gene.file,header = T,stringsAsFactors = F,check.names = F)[,"x"])
    vst.data.used <- vst.data.used[variable.gene.list,]
}

sname <- intersect(rownames(myDesign),colnames(vst.data.used))
sc3.res <- runSC3Analysis(vst.data.used[,sname],sprintf("%s/%s.SC3",out.dir,sample.id),
               myDesign[sname,"sampleType"],sampleTypeColor[names(sampleTypeColor) %in% unique(as.character(myDesign[sname,"sampleType"]))],
               do.log.scale=FALSE,n.cores=4)
save(sc3.res,file=sprintf("%s/%s.SC3.RData",out.dir,sample.id))
q()

qcAndVizMat(vst.data.used,myDesign,out.dir,
            colSet=sampleTypeColor,intgroup=c("sampleType"),ntop=500,extra=paste0(".",sample.id),
            sfilter=NULL, gfilter=NULL,
            complexHeatmap.use=TRUE,
            clonotype.col = list(share=clonotype.share.data,strict=clonotype.strict.data),
            runNMF=FALSE)
#qcAndVizMat(vst.data.used,myDesign,out.dir,colSet=sampleTypeColor,intgroup=c("sampleType"),ntop=500,extra=paste0(".",sample.id,".C_share"),sfilter=NULL, gfilter=NULL,complexHeatmap.use=TRUE,clonotype.col = clonotype.share.data,runNMF=TRUE)
#qcAndVizMat(vst.data.used,myDesign,out.dir,colSet=sampleTypeColor,intgroup=c("sampleType"),ntop=500,extra=paste0(".",sample.id,".C_strict"),sfilter=NULL, gfilter=NULL,complexHeatmap.use=TRUE,clonotype.col = clonotype.strict.data,runNMF=TRUE)
loginfo("process end.")
#invisible(sapply(levels(myDesign$sampleType),function(x){
#                 loginfo(paste0("... ",x))
#                 s.out.dir <- paste0(out.dir,"/",x)
#                 dir.create(s.out.dir,showWarnings = F,recursive = T)
#                 qcAndVizMat(vstMat.blind,myDesign,s.out.dir,colSet=sampleTypeColor,intgroup=c("sampleType"),ntop=500,extra=sprintf(".%s.%s.%s",sample.id,x,"C_share"),sfilter=rownames(subset(myDesign,sampleType==x)), gfilter=NULL,complexHeatmap.use=TRUE,clonotype.col = clonotype.share.data)
#                 qcAndVizMat(vstMat.blind,myDesign,s.out.dir,colSet=sampleTypeColor,intgroup=c("sampleType"),ntop=500,extra=sprintf(".%s.%s.%s",sample.id,x,"C_strict"),sfilter=rownames(subset(myDesign,sampleType==x)), gfilter=NULL,complexHeatmap.use=TRUE,clonotype.col = clonotype.strict.data)
#                 #qcAndVizMat(vstMat.blind,myDesign,s.out.dir,colSet=sampleTypeColor,intgroup=c("sampleType"),ntop=500,extra=sprintf(".%s.%s",sample.id,x),sfilter=rownames(subset(myDesign,sampleType==x)), gfilter=NULL)
#}))

loginfo("process end.")

