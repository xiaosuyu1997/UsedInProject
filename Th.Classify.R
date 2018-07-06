#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(require("RColorBrewer"))
suppressPackageStartupMessages(require("plotrix"))

parser <- ArgumentParser()
parser$add_argument("-i", "--inFile", type="character", required=TRUE, help="infile")
parser$add_argument("-a", "--sampleDescFile", type="character", required=TRUE, help="sample.desc.file")
parser$add_argument("-b", "--geneFile", type="character", required=TRUE, help="gene list file")
parser$add_argument("-o", "--outputPrefix", type="character", required=TRUE, help="output prefix")
parser$add_argument("-k", "--groupBy", type="character",default="leafCluster", help="group by")
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
##parser$add_argument("-m", "--scoreMatrix", type="character", help="score matrix used for annotating samples")
parser$add_argument("-n", "--nMax", type="integer", help="max gene to plot")
parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="verbose mode [default %(default)s]")
parser$add_argument("-e", "--expSort", action="store_true", default=FALSE, help="sort by exp in each gene [default %(default)s]")
####parser$add_argument("-z", "--binFile", type="character", help="binarized expression file")
parser$add_argument("-d", "--pdfWidth", type="integer", default=25, help="pdf width [default %(default)s]")
args <- parser$parse_args()
print(args)

args.groupBy <- args$groupBy
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
nMax <- args$nMax

#gene.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/cross.patient/TF.ThType.list"
#in.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/cross.patient/tpm/tpmDat.liver.Y.RData"
#sample.desc.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/cross.patient/new.clusteringResult/liver.leaf.addMajor.CD4.final.txt"
#out.prefix <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/cross.patient/ThType/liver.P5.ThType"
#sample.id <- "CD4"
#cellTypeColorFile <- "/WPS1/zhenglt/work/TCR_chunhong/data/CellType.color"
#clonotype.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/cross.patient/clonotype/clonotype.liver.txt"
#args.log <- FALSE
#args.center <- FALSE
#mode.verbose <- TRUE
#clusterColorFile <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/cross.patient/new.clusteringResult/predefined.majorClusterColor.txt"

source("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")

#suppressPackageStartupMessages(require("ComplexHeatmap"))
#suppressPackageStartupMessages(require("circlize"))
#suppressPackageStartupMessages(require("gridBase"))
#suppressPackageStartupMessages(require("dendextend"))
#suppressPackageStartupMessages(require("RColorBrewer"))
#suppressPackageStartupMessages(require("gplots"))
#suppressPackageStartupMessages(library("vioplot"))
#suppressPackageStartupMessages(library("plotrix"))
#suppressPackageStartupMessages(library("vioplot"))

loginfo("begin ...")


#### sample data
sample.desc <- read.table(sample.desc.file,header = T,sep = "\t",check.names = F,stringsAsFactors = F)
rownames(sample.desc) <- sample.desc$sample
print(dim(sample.desc))
print(head(sample.desc))
print(str(sample.desc))
##q()

### clonotype data
#clonotype.data <- read.clonotype(in.file = clonotype.file,ctype.col = "C_strict")
#### exp data
suppressPackageStartupMessages(library("R.utils"))
if(grepl("\\.scran\\.RData$",in.file,perl=T)){
    lenv <- loadToEnv(in.file)
    Y <- exprs(lenv[["sce.norm"]])
    args.log <- F
    args.center <- F
}else if(grepl("RData$",in.file,perl=T)){
    lenv <- loadToEnv(in.file)
    Y <- lenv[["Y"]]
}else{
    in.table <- read.table(in.file,header = T,sep = "\t",stringsAsFactors = F,check.names = F)
    rownames(in.table) <- in.table[,1]
    Y <- in.table[,c(-1,-2)]
}
### rename "sampleType"
sample.desc$sampleType <- paste(sample.desc$stype,sapply(strsplit(sample.desc$sampleType,""),function(x){x[1]}),sep=".")
sample.desc <- subset(sample.desc,stype=="CD4")
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
gene.desc <- read.delim(gene.file,header = T,sep = "\t",check.names = F,stringsAsFactors = F,
                        colClasses = c("character","character","character"))

expT <- 4
ThClass <- c("Th1","Th2","Th17","Treg","TfH")
names(ThClass) <- c("TBX21","GATA3","RORC","FOXP3","BCL6")
g.f <- intersect(gene.desc$geneID,rownames(Y))
Y.inGeneList <- t(Y[g.f,,drop=F])
colnames(Y.inGeneList) <- g.GNAME[colnames(Y.inGeneList)]
Y.inGeneList.Z <- scale(Y.inGeneList)
varName <- colnames(Y.inGeneList)
ttt<-t(apply(Y.inGeneList,1,function(x,expT=4){
                            if(sum(x>expT)==0) { return(c("NA","NA")) }
                            note <- "mult"
                            if(sum(x>expT)==1){  note <- "single" }
                            idX <- which.max(x)
                            return(c(ThClass[varName[idX]],note))
                        }))
colnames(ttt) <- c("ThClass","Note")
out.df <- data.frame(geneID=rownames(Y.inGeneList),stringsAsFactors = F,check.names = F)
out.df <- cbind(out.df,Y.inGeneList)
out.df <- cbind(out.df,ttt)
s.f <- intersect(rownames(out.df),rownames(sample.desc))
out.df <- cbind(sample.desc[s.f,],out.df[s.f,])
table(out.df$ThClass,out.df$Note)

ff <- out.df$Note=="single"
for(cc in c("majorCluster","sampleType")){
    dat.plot <- table(out.df$ThClass[ff],out.df[ff,cc])
    dat.plot <- dat.plot[rownames(dat.plot)!="NA",]
    dat.plot.a <- apply(dat.plot,2,function(x){x/sum(x)})
    ccColor <- brewer.pal(5,"Set2")
    pdf(sprintf("%s.ThTypeSubtype.%s.pdf",out.prefix,cc),width=10,height=8)
    par(mar=c(12,6,4,8),cex.lab=2,cex.axis=2,cex.axis=2)
    xx <- barplot(dat.plot.a,beside = F,xaxt="n",col=ccColor)
    staxlab(1,at = xx,labels=colnames(dat.plot.a),srt=if(length(dat.plot) >= 5) 45 else 0, cex=1.3)
    legend("topright",legend=rownames(dat.plot.a),inset=c(-0.205,0),xpd=T,cex=1.5,fill=ccColor)
    dev.off()
}

write.table(out.df,file = sprintf("%s.txt",out.prefix),sep = "\t",row.names = F,quote = F)
head(out.df[out.df$Note=="mult",c(-1,-3,-4,-5)])
#sum(out.df$GATA3>expT)
#sum(out.df$GATA3>expT & out.df$TBX21>expT)
#sum(out.df$GATA3>expT & out.df$RORC>expT)
#sum(out.df$GATA3>expT & out.df$FOXP3>expT)
#sum(out.df$GATA3>expT & out.df$BCL6>expT)
#sum(out.df$FOXP3>expT)
#sum(out.df$TBX21>expT)
#sum(out.df$RORC>expT)
#sum(out.df$BCL6>expT)

q()
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
##nMax=10
gid.list <- rownames(dat.plot)
gname.list <- g.GNAME[gid.list]
if(is.null(nMax)) 
{ 
    nUsed <- length(gid.list) 
}else{ 
    nUsed <- min(nMax,length(gid.list)) 
}
nIdx <- ncol(dat.plot)
wTrack <- 1
nTrack <- 1
uClus <- unique(sample.desc$majorCluster)

pdf(sprintf("%s.vioplot.pdf",out.prefix),width=pdf.width,height=nUsed*1+2)
layout(mat = matrix(seq_len(nUsed+2),ncol = 1,byrow = T),heights = c(1,rep(1,nUsed+1),8))
par(mar=c(5,16,0,4),cex.lab=1.5)
plot(seq_along(uClus),rep(nTrack*wTrack+1.5,length(uClus)),
     xlim=c(0.5,length(uClus)+0.5), ylim=c(0,nTrack*wTrack+1.5),xlab="",ylab="",xaxt="n",yaxt="n",axes=0,type="n")
rect(seq_along(uClus)-0.5,rep(0*wTrack+0.5,length(uClus)),seq_along(uClus)+0.5,rep(1*wTrack+0.5,length(uClus)),
     col=majorClusterColor[uClus],border=NA)
#text(seq_along(uClus)-0.5,rep(1*wTrack+0.5+0.2,length(uClus)),labels=uClus,adj=0,cex=1.0)
staxlab(1,at=seq_along(uClus),labels=uClus,srt=30)

par(mar=c(0,16,0,4),cex.lab=1.5)
for(i in seq_len(nUsed)){
    plot(0,0,type="n",xlim=c(0.5,length(uClus)+0.5), ylim=c(0,14), 
         xaxt = 'n',yaxt="n", xlab ="",ylab="",main="",axes=0)
    par(bty="n")
    for(jj in seq_along(uClus))
    {
        s.desc <- subset(sample.desc,majorCluster==uClus[jj])
        d.vec <- unlist(dat.plot[gid.list[i],rownames(s.desc),drop=T])
        ###d.col <- rgb(colorRamp(brewer.pal(9,"YlOrRd"))(mean(d.vec)/14)/255)
        t.lo <- 0
        t.hi <- 14
        ###d.vec[d.vec < t.lo] <- t.lo
        ###d.vec[d.vec > t.hi] <- t.hi
        cat(sprintf("gene: %s, mean: %4.4f, median: %4.4f\n",gname.list[i],mean(d.vec),median(d.vec)))
        t.mean <- mean(d.vec)
        if(t.mean > t.hi) { 
            t.mean <- t.hi
        }else if(t.mean < t.lo) { 
            t.mean <- t.lo 
        }
        d.col <- rgb(colorRamp(c("yellow","red","black"))(t.mean/t.hi)/255)
        tryCatch(vioplot(d.vec,at=jj,add = T,drawRect=F,col=d.col,border=NA),
                 error = function(e) { loginfo(sprintf("vioplot() catch exception: gene %s",gname.list[i])); e })
    }
    staxlab(2,at = 7,labels = gname.list[i],srt=0,adj=1,cex=3)
}
par(mar=c(4,16,1,4),cex.lab=3,cex.axis=2.5,tck=-0.5)
##image(matrix(1:100),col=colorRampPalette(brewer.pal(9,"YlOrRd"))(100), ylab="",yaxt="n",xaxt="n")
image(matrix(1:100),col=colorRampPalette(c("yellow","red","black"))(100), ylab="",yaxt="n",xaxt="n")
axis(1,seq(0,1,1/7),labels = seq(0,14,14/7),padj=1.2)
dev.off()

loginfo("end.")
