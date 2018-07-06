#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-i", "--inFile", type="character", required=TRUE, help="infile")
parser$add_argument("-a", "--sampleDescFile", type="character", required=TRUE, help="sample.desc.file")
parser$add_argument("-b", "--geneFile", type="character", required=TRUE, help="gene list file")
parser$add_argument("-o", "--outputPrefix", type="character", required=TRUE, help="output prefix")
parser$add_argument("-s", "--sample", type="character", default="SAMPLE", help="sample id  [default %(default)s]")
parser$add_argument("-y", "--cellTypeColorFile", type="character", default="/WPS1/zhenglt/work/TCR_chunhong/data/CellType.color", 
                    help="cellTypeColorFile [default %(default)s]")
parser$add_argument("-w", "--clusterColorFile", type="character", 
                    help="clusterColorFile [default %(default)s]")
parser$add_argument("-f", "--clonotypeFile", type="character", help="clonotype file")
parser$add_argument("-t", "--rowSort", action="store_true", default=FALSE, help="row sort [default %(default)s]")
parser$add_argument("-l", "--log", action="store_true", default=FALSE, help="whether do log transform [default %(default)s]")
parser$add_argument("-r", "--center", action="store_true", default=FALSE, 
                    help="whether center genes by individual [default %(default)s]")
parser$add_argument("-m", "--scoreMatrix", type="character", help="score matrix used for annotating samples")
parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="verbose mode [default %(default)s]")
parser$add_argument("-e", "--expSort", action="store_true", default=FALSE, help="sort by exp in each gene [default %(default)s]")
####parser$add_argument("-z", "--binFile", type="character", help="binarized expression file")
parser$add_argument("-d", "--pdfWidth", type="integer", default=25, help="pdf width [default %(default)s]")
args <- parser$parse_args()
print(args)


in.file <- args$inFile
sample.desc.file <- args$sampleDescFile
gene.file <- args$geneFile
out.prefix <- args$outputPrefix
sample.id <- args$sample
cellTypeColorFile <- args$cellTypeColorFile
clusterColorFile <- args$clusterColorFile
clonotype.file <- args$clonotypeFile
args.rowSort <- args$rowSort
args.log <- args$log
args.center <- args$center
scoreMatrix.file <- args$scoreMatrix
mode.verbose <- args$verbose
#args.rowSort <- FALSE
pdf.width <- args$pdfWidth
exp.sort <- args$expSort
#bin.exp.file <- args$binFile

###in.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/cross.patient/clusteringResult/liver.tSNE.colorByCluster.CD8.liver.CD8.Y.RData"
#in.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/cross.patient/tpm/tpmDat.liver.txt.gz"
#sample.desc.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/cross.patient/clusteringResult/liver.leaf.noUncharacterized.slim.addMajor.CD8.sort.txt"
#gene.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/cross.patient/SLC4A10/SLC4A10.speGene.list"
#out.prefix <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/cross.patient/SLC4A10/test.figure1C"
#sample.id <- "CD8"
#cellTypeColorFile <- "/WPS1/zhenglt/work/TCR_chunhong/data/CellType.color"
#clonotype.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/cross.patient/clonotype/clonotype.liver.txt"
#args.log <- FALSE
#args.center <- FALSE
#mode.verbose <- TRUE
#scoreMatrix.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/cross.patient/gsva/gsva.bindea.liver.txt.gz"
#args.rowSort <- TRUE
#clusterColorFile <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/cross.patient/clusteringResult/liver.tSNE.colorByCluster.useAOV.n3000.liver.majorClusterColor.txt"

source("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")

suppressPackageStartupMessages(require("ComplexHeatmap"))
suppressPackageStartupMessages(require("circlize"))
suppressPackageStartupMessages(require("gridBase"))
suppressPackageStartupMessages(require("dendextend"))
suppressPackageStartupMessages(require("RColorBrewer"))
suppressPackageStartupMessages(require("gplots"))
suppressPackageStartupMessages(library("vioplot"))
suppressPackageStartupMessages(library("plotrix"))

loginfo("begin ...")


#### sample data
sample.desc <- read.table(sample.desc.file,header = T,sep = "\t",check.names = F,stringsAsFactors = F)
rownames(sample.desc) <- sample.desc$sample
print(dim(sample.desc))
print(head(sample.desc))
print(str(sample.desc))
##q()

### clonotype data
clonotype.data <- read.clonotype(in.file = clonotype.file,ctype.col = "C_strict")
#### exp data
if(grepl("RData$",in.file,perl=T)){
    suppressPackageStartupMessages(library("R.utils"))
    lenv <- loadToEnv(in.file)
    Y <- lenv[["Y"]]
}else{
    in.table <- read.table(in.file,header = T,sep = "\t",stringsAsFactors = F,check.names = F)
    rownames(in.table) <- in.table[,1]
    Y <- in.table[,c(-1,-2)]
}
sname <- intersect(rownames(sample.desc),colnames(Y))
sample.desc <- sample.desc[sname,,drop=F]
Y <- Y[,sname,drop=F]

f <- apply(Y,1,function(x){ nE <- sum(x>0); return( nE > 5 & nE/length(x) > 0.01 )  })
Y <- Y[f,]
if(args.log) { Y <- log2(Y+1) }
if(args.center){
    Y.new <- c()
    for(pp in unique(sample.desc$patient)){
        Y.block <- t(scale(t(Y[,subset(sample.desc,patient==pp,select="sample",drop=T)]),center = T,scale = F))
        Y.new <- cbind(Y.new,Y.block)
        ##print(apply(Y.block[1:4,],1,mean))
        ##print(apply(Y.block[1:4,],1,sd))
    }
    Y <- Y.new
    Y <- Y[,sname,drop=F]
}
print(dim(Y))
print(Y[1:4,1:6])

g.GNAME <- entrezToXXX(rownames(Y))
names(g.GNAME) <- rownames(Y)

if(mode.verbose){
    save(Y,g.GNAME,sample.desc,file=sprintf("%s.Y.RData",out.prefix))
}
### binarized exp data
##if(!is.null(bin.exp.file) && file.exists(bin.exp.file)){
##    bin.exp.table <- read.table(bin.exp.file,header = T,sep = "\t",check.names = F,stringsAsFactors = F,row.names = 1)
##    bin.exp.table <- bin.exp.table[,-1]
##}

#### gene data
gene.desc <- read.table(gene.file,header = T,sep = "\t",check.names = F,stringsAsFactors = F,
                        colClasses = c("character","character","character"))

g.f <- intersect(gene.desc$geneID,rownames(Y))
Y.inGeneList <- Y[g.f,,drop=F]
### cell type color
sampleTypeColor <- read.SampleTypeColor(cellTypeColorFile)
sampleTypeColor <- sampleTypeColor[ names(sampleTypeColor) %in% unique(sample.desc[colnames(Y),"sampleType"]) ]

### major cluster color
nMajor <- length(unique(sample.desc$majorCluster))
if(!is.null(clusterColorFile) && file.exists(clusterColorFile)){
    majorClusterColor <- read.SampleTypeColor(clusterColorFile)
    majorClusterColor <- majorClusterColor[ names(majorClusterColor) %in% unique(sample.desc[colnames(Y),"majorCluster"]) ]
}else{
    majorClusterColor <- structure(colorRampPalette(brewer.pal(nMajor,"Paired"))(nMajor),
                              names=unique(sample.desc$majorCluster))
}
print(majorClusterColor)
### reorder samples by leaf cluster than hclust
#dat.plot <- c()
### row order
#if(args.rowSort){
#    loginfo("sort row by hclust")
#    dat.plot <- Y.inGeneList
#}else{
#    loginfo("sort row by Group then hclust")
#    for(g.cls in unique(gene.desc$Group))
#    {
#        g.desc <- subset(gene.desc,Group==g.cls)
#        dat.t <- Y.inGeneList[g.desc$geneID,,drop=F]
#        if(nrow(dat.t)>2){
#            tryCatch({
#                hc.t <- hclust(dist(dat.t),"complete")
#                dat.t <- dat.t[hc.t$order,,drop=F]
#            },error=function(e) e)
#        }
#        dat.plot <- rbind(dat.plot,dat.t)
#    }
#}
### column order
dat.plot <- Y.inGeneList
t.dat.plot <- dat.plot
dat.plot <- c()
for(s.cls in (unique(sample.desc$majorCluster)))
{
    s.desc <- subset(sample.desc,majorCluster==s.cls)
    ###dat.t <- Y.inGeneList[,rownames(s.desc),drop=F]
    dat.t <- t.dat.plot[,rownames(s.desc),drop=F]
    if(ncol(dat.t)>2){
        tryCatch({
            hc.t <- hclust(dist(t(dat.t)), "complete")
            dat.t <- dat.t[,hc.t$order,drop=F]
            if(!is.null(clonotype.data)){
                ### clonotype
                clonotype.d.t <- clonotype.data$ctypeII[colnames(dat.t)]
                clonotype.d.t[is.na(clonotype.d.t)] <- "N0"
                names(clonotype.d.t) <- colnames(dat.t)
                ### sampleType
                sampleType.d.t <- factor(sample.desc[colnames(dat.t),"sampleType"],
                                         levels=c("PTC","NTC","TTC","PTH","NTH","TTH","PTR","NTR","TTR"))
                ### 
                ###dat.t <- dat.t[,names(sort(clonotype.d.t))]
                dat.t <- dat.t[,order(clonotype.d.t,sampleType.d.t)]
            }
        },error = function(e) e)
    }
    if(!is.null(dat.plot)) { dat.plot <- cbind(dat.plot,dat.t) 
    }else{ dat.plot <- dat.t }
}

#####
nMax=20
gid.list <- rownames(dat.plot)
gname.list <- g.GNAME[gid.list]
nUsed <- min(nMax,length(gid.list))
nIdx <- ncol(dat.plot)
wTrack <- 1
nTrack <- 1
pdf(sprintf("%s.barplot.pdf",out.prefix),width=pdf.width,height=8)
layout(mat = matrix(seq_len(nUsed+2),ncol = 1,byrow = T),heights = c(1,rep(1,nUsed+1),1))
par(mar=c(0,12,0,2),cex.lab=1.5)

plot(seq_len(nIdx),rep(nTrack*wTrack+1.5,nIdx),ylim=c(0,nTrack*wTrack+1.5),xlab="",ylab="",xaxt="n",yaxt="n",axes=0,type="n")
#rect(seq_len(nIdx)-0.5,rep(0.5,length(nIdx)),seq_len(nIdx)+0.5,rep(1*wTrack+0.5,length(nIdx)),
#     col=sampleTypeColor[sample.desc[colnames(dat.plot),"sampleType"]],border=NA)
#ctCol <- clonotype.data$colII[clonotype.data$ctypeII[colnames(dat.plot)]]
#ctCol[is.na(ctCol)] <- "gray"
#rect(seq_len(nIdx)-0.5,rep(1*wTrack+0.5,length(nIdx)),seq_len(nIdx)+0.5,rep(2*wTrack+0.5,length(nIdx)),
#     col=ctCol,border=NA)
#rect(seq_len(nIdx)-0.5,rep(2*wTrack+0.5,length(nIdx)),seq_len(nIdx)+0.5,rep(3*wTrack+0.5,length(nIdx)),
#     col=majorClusterColor[sample.desc[colnames(dat.plot),"majorCluster"]],border=NA)

rect(seq_len(nIdx)-0.5,rep(0*wTrack+0.5,length(nIdx)),seq_len(nIdx)+0.5,rep(1*wTrack+0.5,length(nIdx)),
     col=majorClusterColor[sample.desc[colnames(dat.plot),"majorCluster"]],border=NA)

par(mar=c(0,12,0,2),cex.lab=1.5)
for(i in seq_len(nUsed)){
    plot(seq_along(dat.plot[gid.list[i],]),dat.plot[gid.list[i],],type="h",
         col=majorClusterColor[sample.desc[colnames(dat.plot),"majorCluster"]],
         xlab="",main="",xaxt="n",yaxt="n",las=2, ylab="",ylim=c(0,14),lwd=0.1)
    staxlab(2,at = 7,labels = gname.list[i],srt=0,adj=1,cex=2)
}
dev.off()

loginfo("end.")
