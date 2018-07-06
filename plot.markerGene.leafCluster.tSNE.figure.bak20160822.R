#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-i", "--inFile", type="character", required=TRUE, help="infile")
parser$add_argument("-a", "--sampleDescFile", type="character", required=TRUE, help="sample.desc.file")
parser$add_argument("-b", "--geneDescFile", type="character", help="gene.desc.file")
parser$add_argument("-o", "--outputPrefix", type="character", required=TRUE, help="output prefix")
parser$add_argument("-s", "--sample", type="character", default="SAMPLE", help="sample id  [default %(default)s]")
parser$add_argument("-e", "--seed", type="integer", help="random number seed ")
parser$add_argument("-n", "--nKeep", type="integer", default=3000, help="number of genes used [default %(default)s]")
parser$add_argument("-y", "--cellTypeColorFile", type="character",
                    default="/WPS1/zhenglt/work/TCR_chunhong/data/CellType.color", 
                    help="cellTypeColorFile [default %(default)s]")
parser$add_argument("-w", "--newickFile", type="character",help="newickFile")
parser$add_argument("-q", "--save", action="store_true", default=FALSE, help="save object in .RData file [default %(default)s]")
parser$add_argument("-l", "--log", action="store_true", default=FALSE, help="whether do log transform [default %(default)s]")
parser$add_argument("-r", "--center", action="store_true", default=FALSE, help="whether center genes by individual [default %(default)s]")
#parser$add_argument("-l", "--clonotypeFile", type="character", help="clonotype file")
#parser$add_argument("-u", "--u", action="store_true", default=FALSE, help="unique gene [default %(default)s]")
#parser$add_argument("-e", "--scale", action="store_true", default=FALSE, help="do scale [default %(default)s]")
args <- parser$parse_args()
print(args)

in.file <- args$inFile
sample.desc.file <- args$sampleDescFile
gene.desc.file <- args$geneDescFile
out.prefix <- args$outputPrefix
sample.id <- args$sample
cellTypeColorFile <- args$cellTypeColorFile
rn.seed <- args$seed
newick.file <- args$newickFile
args.log <- args$log
args.center <- args$center
nKeep <- args$nKeep
args.save <- args$save
#clonotype.file <- args$clonotypeFile

##newick.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/byTCellType/clusteringResult/P0322/P0322.Newick.txt"
##in.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/OUT.scLVM/P0322/sfIgnoreERCC/P0322.het.countGeneData.sfNormalized"
##sample.desc.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/byTCellType/clusteringResult/P0322/clustering.P0322.iterative.leaf.txt"
##gene.desc.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/byTCellType/clusteringResult/P0322/clustering.P0322.SC3.marker.list"
##out.prefix <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/byTCellType/clusteringResult/P0322/P0322.tSNE.colorByCluster"
##sample.id <- "P0322"
##cellTypeColorFile <- "/WPS1/zhenglt/work/TCR_chunhong/data/CellType.color"
##clonotype.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/tracer/OUT/P0322/filtered_TCR_summary/P0322.summary.cell.reassigneClonotype.methodChunhong.r.txt"

source("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")

suppressPackageStartupMessages(require("ComplexHeatmap"))
suppressPackageStartupMessages(require("circlize"))
suppressPackageStartupMessages(require("gridBase"))
suppressPackageStartupMessages(require("dendextend"))
suppressPackageStartupMessages(require("RColorBrewer"))
suppressPackageStartupMessages(require("gplots"))
suppressPackageStartupMessages(require("plotrix"))
suppressPackageStartupMessages(require("fields"))
suppressPackageStartupMessages(require("magrittr"))
suppressPackageStartupMessages(require("DECIPHER"))

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
if(args.save){
    save(Y,file=sprintf("%s.%s.Y.RData",out.prefix,sample.id))
    loginfo("Y data saved.")
}

gname <- entrezToXXX(rownames(Y))

#### gene data
if(!is.null(gene.desc.file) && file.exists(gene.desc.file)){
    gene.desc <- read.table(gene.desc.file,header = T,sep = "\t",check.names = F,stringsAsFactors = F)
    gene.desc$geneID <-XXXToEntrez(gene.desc$geneName)
    #colnames(gene.desc) <- c("patient","geneName","markerClass")
    #gene.desc <- subset(gene.desc,markerClass!="uncharacterized")
    print(dim(gene.desc))
    print(head(gene.desc))
}
#### cell type color
sampleTypeColor <- read.SampleTypeColor(cellTypeColorFile)
sampleTypeColor <- sampleTypeColor[intersect(names(sampleTypeColor),unique(sample.desc$sampleType))]
leafClusterColor <- structure(colorRampPalette(brewer.pal(12,"Paired"))(length(unique(sample.desc$leaf))),
                              names=unique(sample.desc$leafCluster))
leafClusterColor["uncharacterized"] <- "lightgray"

Y.sd <- apply(Y, 1, sd)
Y.sd.sort <- sort(Y.sd,decreasing=TRUE)
Y.sd.df <- data.frame(geneID=names(Y.sd.sort),geneSymbol=entrezToXXX(names(Y.sd.sort)))
Y.sd.df <- cbind(Y.sd.df,Y.sd.sort)
write.table(Y.sd.df,sprintf("%s.%s.Y.sd.sort.txt",out.prefix,sample.id),row.names = F,sep = "\t",quote = F)

###nKeep <- 3000

loginfo("tSNE begin.")
###tsne.out <- runTSNEAnalysis(Y,out.prefix,
tsne.out <- runTSNEAnalysis(Y[head(names(Y.sd.sort),n=nKeep),],out.prefix,
                            legend=names(leafClusterColor),
                            col.points=leafClusterColor[sample.desc$leafCluster],
                            col.legend=leafClusterColor,
                            pch=16,pch.legend=16,
                            inPDF=TRUE,eps=2.0,dims=2,k=NULL,
                            do.dbscan=F,myseed=rn.seed,width.pdf = 14,margin.r = 24,legend.inset = -0.55)

tsne.dbscan.out <- runTSNEAnalysis(Y[head(names(Y.sd.sort),n=nKeep),],sprintf("%s.dbScan",out.prefix),
                                   legend=names(sampleTypeColor),
                                   col.points=sampleTypeColor[sample.desc$sampleType],
                                   col.legend=sampleTypeColor,
                                   pch=16,pch.legend=16,
                                   inPDF=TRUE,eps=1.00,dims=2,k=NULL,
                                   do.dbscan=T,myseed=rn.seed,preSNE=tsne.out[["30"]]$Rtsne.res,width.pdf = 14,
                                   margin.r = 24,legend.inset = -0.55)

loginfo("tSNE done.")
if(args.save){ save(tsne.out,tsne.dbscan.out,sample.desc,file=sprintf("%s.obj.RData",out.prefix)) }
### integrated figure
do.plot.integratedFigure <- function(newick.file=NULL){
    if(is.null(newick.file)) {
        pdf(file=sprintf("%s.integrated.pdf",out.prefix),width=16,height=8)
    }else{
        pdf(file=sprintf("%s.integrated.newick.pdf",out.prefix),width=16,height=8)
    }
    layout(matrix(c(1,2,3,4,5,5,5,5),ncol=4,byrow = T),widths = c(0.55,0.15,0.30,0.25),heights = c(0.85,0.15))
    par(mar=c(5,5,4,0),cex.lab=1.5,cex.main=1.5)
    plot(tsne.out[["30"]]$Rtsne.res$Y[,1],tsne.out[["30"]]$Rtsne.res$Y[,2], t='p',pch=16,col="lightgray", main="BarnesHutSNE",xlab="Dim1",ylab="Dim2",cex.lab=2.0,cex.axis=1.5,cex.main=2.0)
    f.points <- !grepl("uncharacterized",sample.desc$leafCluster,perl=T)
    points(tsne.out[["30"]]$Rtsne.res$Y[f.points,],col=leafClusterColor[sample.desc$leafCluster[f.points]],pch=16)
    #legend("bottom",horiz = T,legend=names(leafClusterColor),fill = NULL,inset = -0.01,xpd = NA,cex=1.5,pch=16,border =NA,col = leafClusterColor)
    opar <- par(mar=c(5,11,4,8),cex.lab=1.5,cex.main=1.5)
    if(!is.null(newick.file) && file.exists(newick.file)){
        hc.merged <- ReadDendrogram(newick.file,convertBlanks = F)
    }else{
        ### hc on leafClusters
        leaf.cluster <- unique(sample.desc$leafCluster)
        leaf.cent <- NULL
        for(i in seq_along(leaf.cluster)){
            leaf.cent <- cbind(leaf.cent,rowMeans(Y[,sample.desc[colnames(Y),"leafCluster"]==leaf.cluster[i],drop=F]))
        }
        colnames(leaf.cent) <- leaf.cluster
        f <- grepl("^CD8",leaf.cluster,perl = T)
        if(sum(f)>1){
            hc.CD8.leafCluster <- as.dendrogram(hclust(dist(t(leaf.cent[,f,drop=F]))^2, method = "complete", 
                                 members = table(sample.desc$leafCluster)[leaf.cluster[f]]))
        }else{
            hc.CD8.leafCluster <- NULL
        }
        f <- grepl("^CD4",leaf.cluster,perl = T)
        if(sum(f)>1){
            hc.CD4.leafCluster <- as.dendrogram(hclust(dist(t(leaf.cent[,f,drop=F]))^2, method = "complete", 
                                 members = table(sample.desc$leafCluster)[leaf.cluster[f]]))
        }else{
            hc.CD4.leafCluster <- NULL
        }
        #f <- grepl("^CD8",leaf.cluster,perl = T)
        #hc.CD8.leafCluster <- as.dendrogram(hclust(as.dist(1-cor(leaf.cent[,f,drop=F],method="spearman")), method = "complete", 
        #                         members = table(sample.desc$leafCluster)[leaf.cluster[f]]))
        #f <- grepl("^CD4",leaf.cluster,perl = T)
        #hc.CD4.leafCluster <- as.dendrogram(hclust(as.dist(1-cor(leaf.cent[,f,drop=F],method="spearman")), method = "complete", 
        #                         members = table(sample.desc$leafCluster)[leaf.cluster[f]]))
        if(!is.null(hc.CD8.leafCluster) && is.null(hc.CD4.leafCluster)){
            hc.merged <- hc.CD8.leafCluster
        }else if(is.null(hc.CD8.leafCluster) && !is.null(hc.CD4.leafCluster)){
            hc.merged <- hc.CD4.leafCluster
        }else if(!is.null(hc.CD8.leafCluster) && !is.null(hc.CD4.leafCluster)){
            hc.merged <- merge(hc.CD8.leafCluster,hc.CD4.leafCluster)
        }else{
            loginfo("both hc.CD8.leafCluster and hc.CD4.leafCluster are NULL! error will occur!")
        }
    }
    #hc.merged <- set(hc.merged,"branches_k_col",leafClusterColor[leaf.cluster[!grepl("uncharacterized",leaf.cluster)][order.dendrogram(hc.merged)]])
    hc.merged <- set(hc.merged,"branches_k_col",leafClusterColor[labels(hc.merged)])
    hc.merged <- set(hc.merged,"branches_lwd",3)
    par(mar=c(5,2,4,0),cex.lab=1.5,cex.main=1.5)
    plot(hc.merged,horiz=T,yaxt="n",leaflab="none")
    ##plot(hc.merged,horiz=T,yaxt="n",leaflab="none",edgePar = list(lwd = 1.5:2))

    ### barplot describing cell origin
    dat.plot <- table(sample.desc[,c("sampleType","leafCluster")])
    dat.plot <- dat.plot[,labels(hc.merged)]
    ########dat.plot <- apply(dat.plot,2,function(x){x/sum(x)})
    ###par(mar=c(5,11,4,8),cex.lab=1.5,cex.main=1.5)
    par(mar=c(5,1,4,4),cex.lab=1.5,cex.main=1.5)
    xx <- barplot(dat.plot,horiz = T,beside=F,col=sampleTypeColor[rownames(dat.plot)],yaxt="n",
                  sub = "Cell Origin",cex.axis = 1.5,cex.sub = 2.0)
    ####staxlab(2,at = xx,labels=colnames(dat.plot),srt=0, cex=1.0,adj=1,top.line=0.5)
    #text(rep(par("usr")[1]-10,length(xx)),xx,labels=colnames(dat.plot), cex=0.95,adj=1,xpd=NA)
    points(x=rep(-20,length(xx)),y=xx,cex=3,col=leafClusterColor[colnames(dat.plot)],pch=16,xpd=NA)
    legend("bottomright",legend=rownames(dat.plot),fill=sampleTypeColor[rownames(dat.plot)],
           border=sampleTypeColor[rownames(dat.plot)],cex=1.5,inset=c(-0.08,0.05),xpd=NA,horiz = F)
    ###legend("topright",legend=rownames(dat.plot),fill=sampleTypeColor[rownames(dat.plot)],
    ###       border=sampleTypeColor[rownames(dat.plot)],cex=1.5,inset=c(-0.30,0),xpd=NA)
    
    ### barplot describing patient
    dat.plot <- table(sample.desc[,c("patient","leafCluster")])
    dat.plot <- dat.plot[,labels(hc.merged)]
    par(mar=c(5,0,4,4),cex.lab=1.5,cex.main=1.5)
    patientColor <- auto.colSet(n=nrow(dat.plot),name="Accent")
    xx <- barplot(dat.plot,horiz = T,beside=F,col=patientColor,yaxt="n",
                  sub = "Patient",cex.axis = 1.5,cex.sub = 2.0)
    legend("bottomright",legend=rownames(dat.plot),fill=patientColor,
           border=patientColor,cex=1.5,inset=c(-0.13,0.05),xpd=NA)

    ### legend
    par(mar=c(1,1,2,1),cex.lab=1.5,cex.main=1.5)
    plot.new()
    ##idx.leaf <- ceiling(length(colnames(dat.plot))/2)
    leafCluster.name <- colnames(dat.plot)
    #leafCluster.name <- paste0(leafCluster.name,strrep("x",max(nchar(leafCluster.name))-nchar(leafCluster.name)))
    print(leafCluster.name)
    ii <- 0
    mm <- 3
    nn <- 6
    for(yy in head(1-seq(0,1,1/mm),n=mm)){ 
        for(xx in head(seq(0,1,1/nn),n=nn)) {
            ii <- ii+1
            text(x = xx+0.018,y = yy,labels = leafCluster.name[ii],cex=1.15,adj = c(0,1),offset = 4)
            points(xx+0.012,yy-0.050,cex=3,pch=16,xpd=NA,col=leafClusterColor[leafCluster.name[ii]])
        } 
    }
    dev.off()
}
do.plot.integratedFigure()
#do.plot.integratedFigure(newick.file)

if(!is.null(gene.desc.file) && file.exists(gene.desc.file)){
    gOnMapDir <- sprintf("%s.GeneOnSNEMap",out.prefix)
    dir.create(gOnMapDir,showWarnings = F,recursive = T)
    for(i in seq_len(nrow(gene.desc))){
        gid <- gene.desc[i,"geneID"]
        gname <- gene.desc[i,"geneName"]
        loginfo(sprintf("try to plot gene: %s",gname))
        if(gid %in% rownames(Y)){
            pdf(file=sprintf("%s/geneOnSNE.%s.pdf",gOnMapDir,gname),width=9,height=8)
            #layout(matrix(c(1,2),ncol=2),widths =c(0.7,0.3))
            par(mar=c(5,5,4,6),cex.lab=1.5,cex.main=1.5)
            plot(tsne.out[["30"]]$Rtsne.res$Y[,1],tsne.out[["30"]]$Rtsne.res$Y[,2], t='p',pch=16,col="lightgray", main=sprintf("%s",gname),xlab="Dim1",ylab="Dim2")
            gid.color <- as.character(cut(Y[gid,s.f],
                         breaks=quantile(0:15, seq(0,1,0.01)),
                         labels=colorRampPalette(brewer.pal(9,"YlOrRd"))(100)))
            points(tsne.out[["30"]]$Rtsne.res$Y,col=gid.color,pch=16)
            image.plot(zlim=c(0,15),legend.only=TRUE,col=colorRampPalette(brewer.pal(9,"YlOrRd"))(100))
            dev.off()
        }
    }
}
