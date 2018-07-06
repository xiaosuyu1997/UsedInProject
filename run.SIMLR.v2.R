#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-d", "--designFile", type="character", required=TRUE, help="design matrix file")
parser$add_argument("-i", "--inputFile", type="character", required=TRUE, help="input file")
parser$add_argument("-c", "--cellTypeColorFile", default="/WPS1/zhenglt/work/TCR_chunhong/data/CellType.color", 
                    type="character", help="cellTypeColorFile [default %(default)s]")
parser$add_argument("-o", "--outDir", type="character", required=TRUE, help="output directory")
parser$add_argument("-s", "--sample", type="character", default="SAMPLE", help="sample id  [default %(default)s]")
parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="whether verbose mode [default %(default)s]")
parser$add_argument("-l", "--log", action="store_true", default=FALSE, help="whether do log transform [default %(default)s]")
parser$add_argument("-r", "--center", action="store_true", default=FALSE,
                    help="whether center genes by individual [default %(default)s]")
parser$add_argument("-m", "--notFilter", action="store_true", default=FALSE, help="don't filtering genes [default %(default)s]")
parser$add_argument("-n", "--nKeep", type="integer", default=1500, help="use top nKeep genes [default %(default)s]")
###parser$add_argument("-b", "--binFile", type="character", help="binarized exp file [default %(default)s]")
###parser$add_argument("-u", "--iterMax", type="integer", default=4, help="number of iter [default %(default)s]")
parser$add_argument("-k", "--k", type="integer",default=2, help="specified optimal k [default %(default)s]")
parser$add_argument("-a", "--large", action="store_true", default=FALSE, help="large mode? [default %(default)s]")
###parser$add_argument("-x", "--densityPar", type="character", help="specified densityPar file [default %(default)s]")
parser$add_argument("-e", "--geneIDType", type="character",default="entrezID",
                    help="geneID type (entrezID, ensemblID) [default %(default)s]")
parser$add_argument("-j", "--measure", type="character", help="measure (tpm, exprs, norm_exprs)[default %(default)s]")
#### method specific
parser$add_argument("-t", "--myseed", type="character", help="seed [default %(default)s]")
args <- parser$parse_args()

print(args)

designFile <- args$designFile
inputFile <- args$inputFile
out.dir <- args$outDir
sample.id <- args$sample
dir.create(out.dir,showWarnings = F,recursive = T)
out.prefix <- sprintf("%s/%s",out.dir,sample.id)
cellTypeColorFile <- args$cellTypeColorFile
args.notFilter <- args$notFilter
args.center <- args$center
args.log <- args$log
args.measure <- args$measure
args.myseed <- args$myseed
args.k <- args$k
args.large <- args$large
args.nKeep <- args$nKeep

suppressPackageStartupMessages(library("SIMLR"))
suppressPackageStartupMessages(library("igraph"))
suppressPackageStartupMessages(library("RColorBrewer"))
suppressPackageStartupMessages(library("ks"))
suppressPackageStartupMessages(library("fields"))
source("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")

#inputFile <- "/WPS1/zhenglt/work/proj_xy/integrated/cross.patient/exp/lung.P9.scran.RData"
#designFile <- "/WPS1/zhenglt/work/proj_xy/integrated/clustering.otherMethod/d/lungC.P9.CD8.designUsed.txt"
#out.prefix <- "/WPS1/zhenglt/work/proj_xy/integrated/clustering.otherMethod/OUT.SIMLR/test.CD8/test.CD8"
#cellTypeColorFile <- "/WPS1/zhenglt/work/TCR_chunhong/data/CellType.color"
#args.notFilter <- F
#args.center <- F
#args.log <- F
#args.measure <- NULL
#args.myseed <- "123456"
#args.k <- 6

dir.create(dirname(out.prefix),showWarnings = F,recursive = T)

g.inputD <- processInput(designFile,cellTypeColorFile,inputFile,args.notFilter,geneFile=NULL,args.center,args.log,
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

#### select variable genes
Y.sd <- apply(Y, 1, sd)
Y.sd.sort <- sort(Y.sd,decreasing=TRUE)
Y.sd.df <- data.frame(geneID=names(Y.sd.sort),geneSymbol=g.GNAME[names(Y.sd.sort)])
Y.sd.df <- cbind(Y.sd.df,Y.sd.sort)
write.table(Y.sd.df,sprintf("%s.Y.sd.sort.txt",out.prefix),row.names = F,sep = "\t",quote = F)
g.f <- head(names(Y.sd.sort),n=args.nKeep)

## ----SIMLR_Large_Scale_run, warning=FALSE--------------------------------
set.seed(args.myseed)
loginfo("SIMLR begin")
if(args.large){
    SIMLR.res = SIMLR_Large_Scale(X = Y[g.f,], c = args.k)
}else{
    SIMLR.res = SIMLR(X = Y[g.f,], c = args.k)
}
loginfo("SIMLR end")

out.cls.df <- data.frame(sample=colnames(Y),
                         cluster=SIMLR.res$y$cluster,
                         stringsAsFactors = F)
write.table(out.cls.df,sprintf("%s.cluster.txt",out.prefix),row.names = F,quote = F,sep = "\t")

clusterColor <- structure(colorRampPalette(brewer.pal(12,"Paired"))(length(unique(SIMLR.res$y$cluster))),
                          names=sprintf("%s",sort(as.character(unique(SIMLR.res$y$cluster)))))

.clust <- SIMLR.res$y$cluster
.gene.table <- my.clusterMarkerGene(Y, clust=.clust, ann.col=sampleTypeColor[myDesign$sampleType],
                                            clust.col=clusterColor,
                                            out.prefix=sprintf("%s.k%d.clust", out.prefix,args.k),
                                            n.cores=NULL,
                                            original.labels=F,
                                            sampleType = myDesign$sampleType,
                                            sampleTypeColSet = sampleTypeColor,gid.mapping = g.GNAME)

##### visualization by tSNE
## a reference tSNE
tsne.res.cls <- runTSNEAnalysis(Y[head(g.f,n=1500),,drop=F],sprintf("%s.k%d.viz",out.prefix,args.k),
                                col.points = clusterColor[sprintf("%s",.clust)],
                                legend=sprintf("C%s",names(clusterColor)),
                                col.legend=clusterColor,preSNE = NULL,
                                myseed=args.myseed,do.dbscan = F,do.scale = T)

## tSNE using DE genes
if(!is.null(.gene.table) && nrow(.gene.table$aov.res$aov.out.sig)>30){
    tsne.res.cls.de <- runTSNEAnalysis(Y[rownames(.gene.table$aov.res$aov.out.sig),],
                                       sprintf("%s.k%d.de.viz",out.prefix,args.k),
                                       col.points = clusterColor[as.character(.clust)],
                                       legend=sprintf("C%s",names(clusterColor)),
                                       col.legend=clusterColor,preSNE = NULL,
                                       myseed=args.myseed,do.dbscan = F,do.scale = T)
}

#### plot on SIMLR space
pdf(file=sprintf("%s.k%d.SIMLR.pdf",out.prefix,args.k),width=10,height=6)
##par(mar=c(5,5,4,margin.r),cex.lab=1.5,cex.main=1.5)
layout(mat = matrix(c(1,2),ncol=2),widths = c(0.6,0.4))
par(mar=c(5,5,4,2))
    ### color by predefined sample type
plot(SIMLR.res$ydata, t='n', main="",xlab="SIMLR component 1",ylab="SIMLR component 2")
points(SIMLR.res$ydata,col=clusterColor[as.character(.clust)],pch=20)
par(mar=c(5,0,4,2))
plot.new()
legend("left",legend=sprintf("C%s",names(clusterColor)),
       fill = NULL,xpd = NA,cex=1.2,pch=16,border =NA,col = clusterColor)
### density map
layout(mat = matrix(c(1,2),ncol=2),widths = c(0.6,0.4))
par(mar=c(6,6,6,6))
.density <- kde(SIMLR.res$ydata)
.zz <- c(10,20,30,40,50,60,70,80,90)
plot(.density,display="filled.contour2", cont=.zz,xlab="Dim1", ylab="Dim2")
image.plot(zlim=c(0,.zz[length(.zz)]),legend.only=TRUE, col = c("transparent", rev(heat.colors(length(.zz)))),
           axis.args=list( at=.zz, labels=sprintf("%s%%",100-.zz)), legend.width=2.5,legend.mar=4.0)
dev.off()

save.image(file=sprintf("%s.all.RData",out.prefix))


#### ----nmi <- performance <- large <- scale-----------------------------------------
##nmi_2 = compare(ZeiselAmit$true_labs[,1], example_large_scale$y$cluster, method="nmi")
##print(nmi_2)
##
#### ----image <- large <- scale, fig.show='hide', fig.width=5, fig.height=5,results='hide'----
##plot(example_large_scale$ydata, 
##     col = c(topo.colors(ZeiselAmit$n_clust))[ZeiselAmit$true_labs[,1]], 
##     xlab = "SIMLR component 1", 
##     ylab = "SIMLR component 2", 
##     pch = 20, 
##     main="SIMILR 2D visualization for ZeiselAmit")
##
#### ----SIMLR_Feature_Ranking_run, results='hide'---------------------------
###set.seed(11111)
###ranks = SIMLR_Feature_Ranking(A=BuettnerFlorian$results$S,X=BuettnerFlorian$in_X)
##
##
