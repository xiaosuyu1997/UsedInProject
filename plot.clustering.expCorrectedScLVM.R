#!/usr/bin/env Rscript


args <- commandArgs(T)
if(length(args)<7)
{
    cat("plot.clustering.R <out.dir> <sample.id> <YFile> <cellNameFile> <geneIDFile> <designFile> <cellType.color.file> [isConvergedFile]\n")
    q()
}

out.dir <- args[1]
sample.id <- args[2]
YFile <- args[3]
cellNameFile <- args[4]
geneIDFile <- args[5]
designFile <- args[6]
cellType.color.file  <- args[7]
isConvergedFile <- args[8]

### for test only
#out.dir <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/test/post"
#sample.id <- "P0205"
#YFile <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/OUT.scLVM/P0205/sfEndo/P0205.vars.Ycorr.txt"
#cellNameFile <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/OUT.scLVM/P0205/sfEndo/P0205.var.cellNames.txt"
#geneIDFile <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/OUT.scLVM/P0205/sfEndo/P0205.var.geneID.txt"
#designFile <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/OUT.scLVM/P0205/sfEndo/P0205.designUsed.txt"
#cellType.color.file <- "/WPS1/zhenglt/work/TCR_chunhong/data/CellType.color"
###
dir.create(out.dir,showWarnings = F,recursive = T)

source("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")
#source("/lustre1/zeminz_pkuhpc/zhenglt/02.pipeline/cancer/lib/myFunc.R")

Y <- read.table(YFile,header = F)
geneID <- read.table(geneIDFile,header = T)$x
cellName <- read.table(cellNameFile,header = T)$x
rownames(Y) <- geneID
colnames(Y) <- cellName

if(!is.null(isConvergedFile) && file.exists(isConvergedFile)){
    isConverged<-as.logical(read.table(isConvergedFile,header=F)$V1)
    Y <- Y[isConverged,]
}

Y.df <- data.frame(geneID=rownames(Y),geneName=entrezToXXX(rownames(Y)))
Y.df <- cbind(Y.df,Y)
write.table(Y.df,file = sprintf("%s/Ycorr.df.txt",out.dir),quote = F,sep = "\t",row.names = F,col.names = T)
myDesign<-read.table(designFile,header=T,row.names="sample",check.names=F,colClasses=c("factor","character","factor","factor"))

sampleTypeColor <- read.SampleTypeColor(cellType.color.file)
### check Y
pdf(sprintf("%s/check.Ycorr.%s.pdf",out.dir,sample.id),width = 8,height = 8)
par(mar=c(5,5,4,8),cex.lab=1.5,cex.main=1.5,cex.main=1.5)
x.mean <- apply(Y,1,mean)
y.var <- apply(Y,1,var)
y.sd <- apply(Y,1,sd)
y.cv2 <- apply(Y,1,function(xx){ var(xx)/(mean(xx)^2)} )
plot(x.mean,y.cv2,xlab="Mean",log="xy",pch=20,cex=0.3,ylab=expression("CV"^2),main="")
plot(x.mean,y.sd,xlab="Mean",pch=20,cex=0.3,ylab="SD",main="")
plot(x.mean,y.var,xlab="Mean",pch=20,cex=0.3,ylab="Var",main="")
dev.off()

loginfo("... all samples.")

sname <- intersect(rownames(myDesign),colnames(Y))

runPCAAnalysis(Y[,sname],sprintf("%s/%s.het.scLVMCorrected.PCA",out.dir,sample.id),
               myDesign[sname,"sampleType"],sampleTypeColor[names(sampleTypeColor) %in% unique(as.character(myDesign[sname,"sampleType"]))],
               ntop=NULL,clonotype.col=NULL,main=sample.id)

runTSNEAnalysis(Y[,sname],sprintf("%s/%s.het.scLVMCorrected.tSNE",out.dir,sample.id),
                names(sampleTypeColor)[names(sampleTypeColor) %in% unique(as.character(myDesign[sname,"sampleType"]))],
                sampleTypeColor[as.character(myDesign[sname,"sampleType"])],sampleTypeColor)
for(nn in c(50,100,200,300,500,1000,2000,9999999))
{
    runHierarchicalClusteringAnalysis(Y[,sname],
                    sprintf("%s/%s.het.scLVMCorrected.hclustering.n%s",out.dir,sample.id,ifelse(nn<99999,nn,"All")),
                    myDesign[sname,"sampleType"],sampleTypeColor[names(sampleTypeColor) %in% unique(as.character(myDesign[sname,"sampleType"]))],
                    clonotype.col=NULL,
                    ntop=nn,
                    complexHeatmap.use=TRUE,do.cuttree=T,
                    verbose=FALSE,main="Variable Genes")
}

sc3.res <- runSC3Analysis(Y[,sname],sprintf("%s/%s.SC3",out.dir,sample.id),
               myDesign[sname,"sampleType"],sampleTypeColor[names(sampleTypeColor) %in% unique(as.character(myDesign[sname,"sampleType"]))],
               do.log.scale=FALSE,n.cores=8)
save(sc3.res,file=sprintf("%s/%s.SC3.RData",out.dir,sample.id))


#runNMFAnalysis(Y[,sname],sprintf("%s/%s.het.scLVMCorrected.NMF",out.dir,sample.id),
#               myDesign[sname,"sampleType",drop=F],
#               list(sampleType=sampleTypeColor))
######qcAndVizMat(Y,myDesign,out.dir,colSet=sampleTypeColor,intgroup=c("sampleType"),ntop=100000,extra=paste0(".",sample.id),sfilter=NULL, gfilter=NULL)

loginfo("... ...")
save(Y,myDesign,file = sprintf("%s/%s.Y.RData",out.dir,sample.id))

