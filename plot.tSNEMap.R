#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-i", "--inFile", type="character", required=TRUE, help="infile")
parser$add_argument("-a", "--sampleDescFile", type="character", required=TRUE, help="sample.desc.file")
parser$add_argument("-b", "--geneFile", type="character", required=TRUE, help="gene list file (\"ALL\" for all genes pass the internal filter)")
parser$add_argument("-o", "--outputPrefix", type="character", required=TRUE, help="output prefix")
parser$add_argument("-s", "--sample", type="character", default="SAMPLE", help="sample id  [default %(default)s]")
parser$add_argument("-y", "--cellTypeColorFile", type="character", default="/WPS1/zhenglt/work/TCR_chunhong/data/CellType.color", 
                    help="cellTypeColorFile [default %(default)s]")
parser$add_argument("-w", "--clusterColorFile", type="character", 
                    help="clusterColorFile [default %(default)s]")
#parser$add_argument("-f", "--clonotypeFile", type="character", help="clonotype file")
parser$add_argument("-l", "--log", action="store_true", default=FALSE, help="whether do log transform [default %(default)s]")
parser$add_argument("-r", "--center", action="store_true", default=FALSE, 
                    help="whether center genes by individual [default %(default)s]")
parser$add_argument("-n", "--nMax", type="integer", help="max gene to plot")
#parser$add_argument("-p", "--perGene", action="store_true", default=FALSE, help="one figure per gene [default %(default)s]")
parser$add_argument("-t", "--tsneFile", type="character", required=TRUE, help="tsne file")
parser$add_argument("-c", "--commonScale", action="store_true", default=FALSE, help="use common legend scale for geneOnMap [default %(default)s]")
parser$add_argument("-d", "--pdfWidth", type="integer", default=8, help="pdf width [default %(default)s]")
parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="verbose mode [default %(default)s]")
parser$add_argument("-u", "--normExprs", action="store_true", default=FALSE, help="use norm_exprs assay [default %(default)s]")
args <- parser$parse_args()
print(args)

in.file <- args$inFile
sample.desc.file <- args$sampleDescFile
gene.file <- args$geneFile
out.prefix <- args$outputPrefix
sample.id <- args$sample
cellTypeColorFile <- args$cellTypeColorFile
clusterColorFile <- args$clusterColorFile
##clonotype.file <- args$clonotypeFile
args.log <- args$log
args.center <- args$center
nMax <- args$nMax
##args.perGene <- args$perGene
args.commonScale <- args$commonScale
pdf.width <- args$pdfWidth
mode.verbose <- args$verbose
tsne.file <- args$tsneFile
args.norm.exprs <- args$normExprs

#in.file <- "/WPS1/zhenglt/work/proj_xy/integrated/cross.patient/exp/lung.P9.scran.RData"
#sample.desc.file <- "/WPS1/zhenglt/work/proj_xy/integrated/cross.patient/clusteringResult/lungC.all.leaf.addMajor.final.forPlot.txt"
#gene.file <- "ALL"
#out.prefix <- "/WPS1/zhenglt/work/proj_xy/integrated/cross.patient/violin/lungC"
#sample.id <- "lungC T Cells"
#cellTypeColorFile <- "/WPS1/zhenglt/work/TCR_chunhong/data/CellType.color"
#clusterColorFile <- "/WPS1/zhenglt/work/proj_xy/integrated/cross.patient/curated.clusterColor.txt"
#clonotype.file <- "/WPS1/zhenglt/work/proj_xy/integrated/tcr/TCRinfo.zyy.formatted.txt"
#args.rowSort <- F
#args.log <- F
#args.center <- F
#args.perGene <- T
#tsne.RData <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/cross.patient/new.clusteringResult/liver.tSNE.colorByCluster.useAOV.n1000.obj.RData"

source("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")

suppressPackageStartupMessages(require("ComplexHeatmap"))
suppressPackageStartupMessages(require("circlize"))
suppressPackageStartupMessages(require("gridBase"))
suppressPackageStartupMessages(require("dendextend"))
suppressPackageStartupMessages(require("RColorBrewer"))
suppressPackageStartupMessages(require("gplots"))
suppressPackageStartupMessages(library("vioplot"))
suppressPackageStartupMessages(library("plotrix"))
suppressPackageStartupMessages(library("vioplot"))
suppressPackageStartupMessages(library("R.utils"))
suppressPackageStartupMessages(require("fields"))
loginfo("begin ...")

g.inputD <- processInput(designFile = sample.desc.file,cellTypeColorFile = cellTypeColorFile ,inputFile = in.file,
                         args.notFilter=T,geneFile = NULL, args.center,args.log,args.norm.exprs = args.norm.exprs)
myDesign <- g.inputD$myDesign
myDesign$sampleType <- paste(myDesign$stype,sapply(strsplit(myDesign$sampleType,""),function(x){x[1]}),sep=".")
sample.desc <- myDesign
###sampleTypeColor <- g.inputD$sampleTypeColor
patient.col.list <- g.inputD$patient.col.list
Y <- g.inputD$Y
g.GNAME <- g.inputD$g.GNAME
print("info of Y:")
print(dim(Y))
print(Y[1:4,1:8])

tenv <- loadToEnv(tsne.file)
if("tsne.out" %in% names(tenv)){
    tsne.out <- tenv[["tsne.out"]]
}else if("tsne.res" %in% names(tenv)){
    tsne.out <- tenv[["tsne.res"]]
}else{
    cat("No tsne.out or tsne.res found !!!!!\n")
    q()
}

#### sample data
#sample.desc <- read.table(sample.desc.file,header = T,sep = "\t",check.names = F,stringsAsFactors = F)
#rownames(sample.desc) <- sample.desc$sample
#print(dim(sample.desc))
#print(head(sample.desc))
#print(str(sample.desc))
#
##### tsne map file
#tenv <- loadToEnv(tsne.file)
#tsne.out <- tenv[["tsne.out"]]
#
#### clonotype data
####clonotype.data <- read.clonotype(in.file = clonotype.file,ctype.col = "C_strict")
##### exp data
#if(grepl("\\.scran\\.RData$",in.file,perl=T)){
#    lenv <- loadToEnv(in.file)
#    Y <- exprs(lenv[["sce.norm"]])
#    args.log <- F
#    args.center <- F
#}else if(grepl("RData$",in.file,perl=T)){
#    lenv <- loadToEnv(in.file)
#    Y <- lenv[["Y"]]
#}else{
#    in.table <- read.table(in.file,header = T,sep = "\t",stringsAsFactors = F,check.names = F)
#    rownames(in.table) <- in.table[,1]
#    Y <- in.table[,c(-1,-2)]
#    sname <- intersect(rownames(sample.desc),colnames(Y))
#    sample.desc <- sample.desc[sname,,drop=F]
#    Y <- Y[,sname,drop=F]
#}
#
#f <- apply(Y,1,function(x){ nE <- sum(x>0); return( nE > 5 & nE/length(x) > 0.01 )  })
#Y <- Y[f,]
#if(args.log) { Y <- log2(Y+1) }
#if(args.center){
#    Y.new <- c()
#    for(pp in unique(sample.desc$patient)){
#        Y.block <- t(scale(t(Y[,subset(sample.desc,patient==pp,select="sample",drop=T)]),center = T,scale = F))
#        Y.new <- cbind(Y.new,Y.block)
#        ##print(apply(Y.block[1:4,],1,mean))
#        ##print(apply(Y.block[1:4,],1,sd))
#    }
#    Y <- Y.new
#    Y <- Y[,sname,drop=F]
#}
#print(dim(Y))
#print(Y[1:4,1:6])
#
#g.GNAME <- entrezToXXX(rownames(Y))
#names(g.GNAME) <- rownames(Y)
#g.GNAME.f <- which(is.na(g.GNAME))
#g.GNAME[g.GNAME.f] <- names(g.GNAME)[g.GNAME.f]
#
#if(mode.verbose){
#    ##save(Y,g.GNAME,sample.desc,file=sprintf("%s.Y.RData",out.prefix))
#}

#### gene data
if(gene.file=="ALL")
{
    ####g.f <- head(rownames(Y),n=20)
    g.f <- rownames(Y)
}else{
    gene.desc <- read.delim(gene.file,header = T,sep = "\t",check.names = F,stringsAsFactors = F)
                            ###colClasses = c("character","character","character"))
    gene.desc$geneID <- as.character(gene.desc$geneID)
    g.f <- intersect(gene.desc$geneID,rownames(Y))
}
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

### tSNE map
### save the coordinate
shared.slist <- intersect(colnames(Y),rownames(tsne.out[["30"]]$Rtsne.res$Y))
Y <- Y[,shared.slist,drop=F]
out.tsne.coor.df <- as.data.frame(tsne.out[["30"]]$Rtsne.res$Y[shared.slist,,drop=F])
colnames(out.tsne.coor.df) <- c("tSNE1","tSNE2")
out.tsne.coor.df <- cbind(data.frame(cellID=shared.slist,stringsAsFactors = F),out.tsne.coor.df)
write.table(out.tsne.coor.df,sprintf("%s.coord.txt",out.prefix),row.names = F,sep = "\t",quote = F)

gOnMapDir <- sprintf("%s.GeneOnSNEMap",out.prefix)
dir.create(gOnMapDir,showWarnings = F,recursive = T)
for(i in seq_along(g.f)){
    gid <- g.f[i]
    gname <- g.GNAME[gid]
    firstC <- unlist(strsplit(gname,""))[1]
    dir.create(sprintf("%s/%s",gOnMapDir,firstC),showWarnings = F,recursive = T)
    loginfo(sprintf("try to plot gene: %s",gname))
    if(gid %in% rownames(Y)){
        pdf(file=sprintf("%s/%s/geneOnSNE.%s.pdf",gOnMapDir,firstC,gname),width=9,height=8)
        #layout(matrix(c(1,2),ncol=2),widths =c(0.7,0.3))
        par(mar=c(5,5,4,8),cex.lab=1.5,cex.main=3.0,cex.axis=1.5)
        plot(tsne.out[["30"]]$Rtsne.res$Y[,1],tsne.out[["30"]]$Rtsne.res$Y[,2], t='n',pch=16,col="lightgray", main=sprintf("%s",gname),xlab="Dim1",ylab="Dim2")
        s.f <- order(Y[gid,],decreasing=F)
        Y.level <- pretty(Y[gid,],n=10)
        ppalette <- brewer.pal(9,"YlOrRd")
        #ppalette <- c("yellow","red","black")
        if(args.commonScale){
            t.lo <- -2
            t.hi <- 12
            t.vec <- Y[gid,s.f]
            t.vec[t.vec < t.lo] <- t.lo
            t.vec[t.vec > t.hi] <- t.hi
            gid.color <- rgb( colorRamp( ppalette )( (t.vec-t.lo)/(t.hi-t.lo) )/255 )
        }else{
            gid.color <- as.character(cut(Y[gid,s.f],
                     breaks=quantile(c(Y.level[1],Y.level[length(Y.level)]), seq(0,1,0.01)),
                     labels=addalpha(colorRampPalette(ppalette)(100),alpha = 1)))
        }
        points(tsne.out[["30"]]$Rtsne.res$Y[s.f,],col=gid.color,pch=20,cex=1.0)
        if(args.commonScale){
            image.plot(zlim=c(t.lo,t.hi),legend.only=TRUE,col=colorRampPalette(ppalette)(100),legend.width=2.5,legend.mar=5.0)
        }else{
            image.plot(zlim=c(Y.level[1],Y.level[length(Y.level)]),legend.only=TRUE,col=colorRampPalette(ppalette)(100),legend.width=2.5,legend.mar=5.0)
        }
        dev.off()
    }
}


#####

loginfo("end.")
