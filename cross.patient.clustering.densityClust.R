#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-d", "--designFile", type="character", required=TRUE, help="design matrix file")
parser$add_argument("-i", "--inputFile", type="character", required=TRUE, help="input file")
parser$add_argument("-f", "--geneFile", type="character", help="gene list file")
parser$add_argument("-c", "--cellTypeColorFile", default="/WPS1/zhenglt/work/TCR_chunhong/data/CellType.color", 
                    type="character", help="cellTypeColorFile [default %(default)s]")
parser$add_argument("-o", "--outDir", type="character", required=TRUE, help="output directory")
parser$add_argument("-s", "--sample", type="character", default="SAMPLE", help="sample id  [default %(default)s]")
parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="whether verbose mode [default %(default)s]")
parser$add_argument("-l", "--log", action="store_true", default=FALSE, help="whether do log transform [default %(default)s]")
parser$add_argument("-r", "--center", action="store_true", default=FALSE, help="whether center genes by individual [default %(default)s]")
parser$add_argument("-m", "--notFilter", action="store_true", default=FALSE, help="don't filtering genes [default %(default)s]")
parser$add_argument("-n", "--nKeep", type="integer", default=1500, help="use top nKeep genes [default %(default)s]")
#### method specific
parser$add_argument("-y", "--transform", type="character", default="none", help="do transformation(none, pca, eigenmap) [default %(default)s]")
parser$add_argument("-a", "--rho", type="double", help="rho for densityClust [default %(default)s]")
parser$add_argument("-b", "--delta", type="double", help="delta for densityClust [default %(default)s]")
parser$add_argument("-e", "--objFile", type="character", help="obj RData file contain tsne result")
#parser$add_argument("-k", "--kbatch", type="character", default="2,3,4,5,6,7,8,9,10", help="kbatch [default %(default)s]")
#parser$add_argument("-w", "--method", type="character", default="kmeans", help="clustering method (kmeans, hclust) [default %(default)s]")
parser$add_argument("-t", "--myseed", type="character", help="seed [default %(default)s]")
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
#args.notSC3 <- args$notSC3
args.notFilter <- args$notFilter

args.transform <- args$transform
args.kbatch <- args$kbatch
args.rho <- args$rho
args.delta <- args$delta
#if(!is.null(args.kbatch)) { args.kbatch <- as.numeric(unlist(strsplit(args.kbatch,","))) }
#args.method <- args$method
args.myseed <- args$myseed
args.objFile <- args$objFile

###args.objFile <- "/WPSnew/zhenglt/work/proj_zl/integrate0924/cross.patient.b20170523/clusteringResult/colonC.CD8.tSNE.colorByCluster.n20000.importantGeneByRF.AUC65/colonC.CD8.useAOV.n20000.obj.RData"
###inputFile <- "/WPSnew/zhenglt/work/proj_zl/integrate0924/cross.patient.b20170523/exp/colonC.P8.scran.RData"
###designFile="/WPSnew/zhenglt/work/proj_zl/integrate0924/cross.patient.b20170523/clusteringResult/colonC.clustering.final.forPlot.CD8.txt"
###cellTypeColorFile="/WPS1/zhenglt/work/TCR_chunhong/data/CellType.color"
###args.notFilter=T
###geneFile=NULL
###args.center=F
###args.log=F

##### frequently used code
dir.create(out.dir,recursive = T,showWarnings = F)
out.prefix <- sprintf("%s/%s.densityClust",out.dir,sample.id)
source("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")
suppressPackageStartupMessages(library("densityClust"))
suppressPackageStartupMessages(library("R.utils"))
suppressPackageStartupMessages(library("RColorBrewer"))
suppressPackageStartupMessages(library("ks"))
suppressPackageStartupMessages(library("fields"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("plotrix"))

g.inputD <- processInput(designFile,cellTypeColorFile,inputFile,args.notFilter,geneFile,args.center,args.log)
myDesign <- g.inputD$myDesign
sampleTypeColor <- g.inputD$sampleTypeColor
patient.col.list <- g.inputD$patient.col.list
Y <- g.inputD$Y
g.GNAME <- g.inputD$g.GNAME
##### 

if(!is.null(args.objFile) && file.exists(args.objFile)){
    lenv <- loadToEnv(args.objFile)
    ###lname <- load(args.objFile)
    tsne.Y <- lenv$tsne.out$`30`$Rtsne.res$Y
    dist.obj <- dist(tsne.Y)
    clust.res <- densityClust::densityClust(dist.obj, gaussian = T)
    if(is.null(args.rho)){
        args.rho <- quantile(clust.res$rho, probs = 0.95)
    }
    if(is.null(args.delta)){
        args.delta <- quantile(clust.res$delta, probs = 0.95)
    }

    set.seed(1000)
    pdf(sprintf("%s.decision.pdf",out.prefix),width = 5,height = 5)
    clust.res <- densityClust::findClusters(clust.res, rho = args.rho, delta = args.delta,plot=T)
    dev.off()

    clusterColor <- structure(colorRampPalette(brewer.pal(12,"Paired"))(length(unique(clust.res$clusters))),
                              names=sort(as.character(unique(clust.res$clusters))))
    clusterColor["0"] <- "gray"

    ### color by clusters
    plot.tsne.points(tsne.Y[names(clust.res$rho),,drop=F],
                 sprintf("%s.points.colByCluster.peakT.pdf",out.prefix),
                 tsne.col.points=clusterColor[as.character(clust.res$clusters)],
                 col.tsne.legend=clusterColor,
                 tsne.legend=names(clusterColor),
                 pch=20,nclusters=length(clusterColor),peak = names(clust.res$peaks))
    plot.tsne.points(tsne.Y[names(clust.res$rho),,drop=F],
                 sprintf("%s.points.colByCluster.peakF.pdf",out.prefix),
                 tsne.col.points=clusterColor[as.character(clust.res$clusters)],
                 col.tsne.legend=clusterColor,
                 tsne.legend=names(clusterColor),
                 pch=20,nclusters=length(clusterColor),peak = NULL)

    ### density
    .density <- kde(tsne.Y)

    pdf(sprintf("%s.density.pdf",out.prefix),width = 5,height = 5)
    par(mar=c(5,4,5,6))
    .zz <- c(10,20,30,40,50,60,70,80,90)
    plot(.density,display="filled.contour2", cont=.zz,xlab="Dim1", ylab="Dim2")
    image.plot(zlim=c(0,.zz[length(.zz)]),legend.only=TRUE, col = c("transparent", rev(heat.colors(length(.zz)))),
               axis.args=list( at=.zz, labels=sprintf("%s%%",100-.zz)), legend.width=2.0,legend.mar=4.5)
    plot(.density,display="filled.contour2", cont=.zz,xlab="Dim1", ylab="Dim2")
    points(tsne.Y[names(clust.res$peaks),,drop=F],pch=3,cex=2,col="black")
    image.plot(zlim=c(0,.zz[length(.zz)]),legend.only=TRUE, col = c("transparent", rev(heat.colors(length(.zz)))),
               axis.args=list( at=.zz, labels=sprintf("%s%%",100-.zz)), legend.width=2.0,legend.mar=4.5)
    dev.off()

    ### marker genes
    clusterMarker.res <- my.clusterMarkerGene(Y[,names(clust.res$rho),drop=F],
                                            clust=clust.res$clusters,
                                            ann.col=sampleTypeColor[myDesign[names(clust.res$rho),"sampleType"]],
                                            clust.col=clusterColor,
                                            out.prefix=sprintf("%s.marker", out.prefix),
                                            n.cores=8,
                                            sampleType = myDesign[names(clust.res$rho),"sampleType"],
                                            sampleTypeColSet = sampleTypeColor)
    
    out.df <- left_join(myDesign,data.frame(sample=names(clust.res$rho),clustersByDensity=sprintf("C%02d",clust.res$clusters),stringsAsFactors = F))
    write.table(out.df,sprintf("%s.cellInfo.txt",out.prefix),row.names = F,sep = "\t",quote = F)

    for(locType in c("P","N","T","")){
        out.df.block <- out.df[grepl(pattern = sprintf("^%s",locType),out.df$sampleType,perl = T),]
        dat.plot <- out.df.block[,c("clustersByDensity","patient")] %>% table %>% apply(.,MARGIN = 2,FUN = function(x){ x/sum(x) })
        .patient.order <- c("P0701","P1207","P1212","P1228","P0909","P0413","P1012","P0215")
        .patient.order <- intersect(.patient.order,unique(out.df.block$patient))
        .patient.order <- c(.patient.order,setdiff(unique(out.df.block$patient),.patient.order))
        dat.plot <- dat.plot[,.patient.order]
        .clusterColor <- clusterColor
        names(.clusterColor) <- sprintf("C%02d",as.numeric(names(.clusterColor)))
        pdf(sprintf("%s.distByPatients.loc%s.pdf",out.prefix,if(locType=="") "All" else locType),width=8,height=8)
        par(mar=c(5,6,4,18),cex.lab=2,cex.main=2,cex.axis=1.5)
        xx <- barplot(dat.plot,beside = F,xlab="",ylab="Clusters",main="",col=.clusterColor[rownames(dat.plot)],xaxt="n")
        staxlab(1,at = xx,labels=colnames(dat.plot),srt=if(length(dat.plot) >= 5) 45 else 0, cex=1.3)
        title(main=sprintf("%s",if(locType=="") "All" else locType),line=1)
        legend("topright",legend=rownames(dat.plot),fill=.clusterColor[rownames(dat.plot)],
               inset = c(-1.05,0),xpd = T,horiz = F,cex=1.2)
        dev.off()
    }



}else{
    #### select variable genes
    Y.sd <- apply(Y, 1, sd)
    Y.sd.sort <- sort(Y.sd,decreasing=TRUE)
    Y.sd.df <- data.frame(geneID=names(Y.sd.sort),geneSymbol=g.GNAME[names(Y.sd.sort)])
    Y.sd.df <- cbind(Y.sd.df,Y.sd.sort)
    write.table(Y.sd.df,sprintf("%s/%s.Y.sd.sort.txt",out.dir,sample.id),row.names = F,sep = "\t",quote = F)

    #### rename rownames
    Y.newRowname <- Y
    newName <- g.GNAME[rownames(Y)]
    names(newName) <- rownames(Y)
    f.gname.na <- is.na(newName)
    newName[f.gname.na] <- names(newName)[f.gname.na]
    rownames(Y.newRowname) <- unname(newName)

    g.f <- head(names(Y.sd.sort),n=nKeep)
    #doit.kmeans(Y,g.f,extra="",myseed=NULL,data.forDE=Y.newRowname)

    ### pca or other transformation
    if(args.transform=="pca"){
        g.out.prefix <- sprintf("%s/%s.%s.pca",out.dir,sample.id,args.method)
        g.do.pca <- T
    }else if(args.transform=="none"){
        g.out.prefix <- sprintf("%s/%s.%s.none",out.dir,sample.id,args.method)
        g.do.pca <- F
    }

    clustering.out <- runSimpleClusteringAnalysis(Y[g.f,],g.out.prefix,myDesign$sampleType,colSet,B=100,
                      legend=c(names(sampleTypeColor)[names(sampleTypeColor) %in% unique(as.character(myDesign$sampleType))],
                               unique(myDesign$libType)),
                      col.legend=c(sampleTypeColor[names(sampleTypeColor) %in% unique(as.character(myDesign$sampleType))],
                                   rep("black",length(unique(myDesign$libType)))),
                      col.points = sampleTypeColor[as.character(myDesign$sampleType)],
                      k.batch = args.kbatch,do.scale=F,do.pca = g.do.pca,data.forDE = Y.newRowname,method=args.method,myseed = args.myseed)
                      ##k.batch = args.kbatch,do.scale=T,do.pca = g.do.pca,data.forDE = NULL,method=args.method)

    #### To do: refine

}

save.image(file = sprintf("%s.all.RData",out.prefix))

