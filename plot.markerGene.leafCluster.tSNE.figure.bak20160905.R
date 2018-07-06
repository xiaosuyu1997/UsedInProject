#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-i", "--inFile", type="character", required=TRUE, help="infile")
parser$add_argument("-a", "--sampleDescFile", type="character", required=TRUE, help="sample.desc.file")
parser$add_argument("-b", "--geneDescFile", type="character", help="gene.desc.file")
parser$add_argument("-o", "--outputPrefix", type="character", required=TRUE, help="output prefix")
parser$add_argument("-x", "--excludeGeneFile", type="character", help="gene (to be excluded) list file")
parser$add_argument("-s", "--sample", type="character", default="SAMPLE", help="sample id  [default %(default)s]")
parser$add_argument("-e", "--seed", type="integer", help="random number seed ")
parser$add_argument("-n", "--nKeep", type="integer", default=3000, help="number of genes used [default %(default)s]")
parser$add_argument("-y", "--cellTypeColorFile", type="character",
                    default="/WPS1/zhenglt/work/TCR_chunhong/data/CellType.color", 
                    help="cellTypeColorFile [default %(default)s]")
parser$add_argument("-w", "--newickFile", type="character",help="newickFile")
parser$add_argument("-q", "--save", action="store_true", default=FALSE, help="save object in .RData file [default %(default)s]")
parser$add_argument("-l", "--log", action="store_true", default=FALSE, help="whether do log transform [default %(default)s]")
parser$add_argument("-r", "--center", action="store_true", default=FALSE, 
                    help="whether center genes by individual [default %(default)s]")
parser$add_argument("-f", "--clonotypeFile", type="character", help="clonotype file")
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
excludeGeneFile <- args$excludeGeneFile
clonotype.file <- args$clonotypeFile

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
##suppressPackageStartupMessages(require("MASS"))
suppressPackageStartupMessages(require("ks"))
suppressPackageStartupMessages(require("reshape2"))
suppressPackageStartupMessages(require("ggplot2"))

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
majorClusterColor <- structure(colorRampPalette(brewer.pal(12,"Paired"))(length(unique(sample.desc$majorCluster))),
                              names=unique(sample.desc$majorCluster))
majorClusterColor["uncharacterized"] <- "lightgray"
##leafClusterColor <- structure(colorRampPalette(brewer.pal(12,"Paired"))(length(unique(sample.desc$leaf))),
##                              names=unique(sample.desc$leafCluster))
##leafClusterColor["uncharacterized"] <- "lightgray"
majorClusterColor.df <- data.frame(sampleType=names(majorClusterColor),color=majorClusterColor)
write.table(majorClusterColor.df,sprintf("%s.%s.majorClusterColor.txt",out.prefix,sample.id),row.names = F,sep = "\t",quote = F)


Y.sd <- apply(Y, 1, sd)
Y.sd.sort <- sort(Y.sd,decreasing=TRUE)
Y.sd.df <- data.frame(geneID=names(Y.sd.sort),geneSymbol=entrezToXXX(names(Y.sd.sort)))
Y.sd.df <- cbind(Y.sd.df,Y.sd.sort)
write.table(Y.sd.df,sprintf("%s.%s.Y.sd.sort.txt",out.prefix,sample.id),row.names = F,sep = "\t",quote = F)

###nKeep <- 3000
g.f <- head(names(Y.sd.sort),n=nKeep)
if(!is.null(excludeGeneFile)){
    loginfo(sprintf("exclude genes in file: %s",excludeGeneFile))
    excludeGeneV <- read.table(excludeGeneFile,header = F,stringsAsFactors = F,check.names = F)$V1
    excludeGeneV.gid <- XXXToEntrez(excludeGeneV)
    excludeGeneV.f <- is.na(excludeGeneV.gid)
    excludeGeneV.gid[excludeGeneV.f] <- excludeGeneV[excludeGeneV.f]
    g.f <- setdiff(g.f,excludeGeneV.gid)
}

loginfo(sprintf("Total %d genes will be used.",length(g.f)))
loginfo("tSNE begin.")
###tsne.out <- runTSNEAnalysis(Y,out.prefix,
tsne.out <- runTSNEAnalysis(Y[g.f,],out.prefix,
                            legend=names(majorClusterColor),
                            col.points=majorClusterColor[sample.desc$majorCluster],
                            col.legend=majorClusterColor,
                            pch=16,pch.legend=16,
                            inPDF=TRUE,eps=2.0,dims=2,k=NULL,
                            do.dbscan=F,myseed=rn.seed,width.pdf = 12,margin.r = 20,legend.inset = -0.50)

tsne.dbscan.out <- runTSNEAnalysis(Y[g.f,],sprintf("%s.dbScan",out.prefix),
                                   legend=names(sampleTypeColor),
                                   col.points=sampleTypeColor[sample.desc$sampleType],
                                   col.legend=sampleTypeColor,
                                   pch=16,pch.legend=16,
                                   inPDF=TRUE,eps=1.00,dims=2,k=NULL,
                                   do.dbscan=T,myseed=rn.seed,preSNE=tsne.out[["30"]]$Rtsne.res,width.pdf = 10,
                                   margin.r = 8,legend.inset = -0.21)

loginfo("tSNE done.")
if(args.save){ save(tsne.out,tsne.dbscan.out,sample.desc,file=sprintf("%s.obj.RData",out.prefix)) }
### integrated figure
do.plot.integratedFigure <- function(newick.file=NULL,perplexity=30){
    ###if(is.null(newick.file)) {
    ###    pdf(file=sprintf("%s.integrated.tsne.pdf",out.prefix),width=8,height=8)
    ###}else{
    ###    pdf(file=sprintf("%s.integrated.tsne.newick.pdf",out.prefix),width=8,height=8)
    ###}
    pdf(file=sprintf("%s.integrated.perplexity%d.tsne.contour.pdf",out.prefix,perplexity),width=8,height=8)
    ###layout(matrix(c(1,2,3,4,5,5,5,5),ncol=4,byrow = T),widths = c(0.55,0.15,0.30,0.25),heights = c(0.85,0.15))
    par(mar=c(5,5,4,2),cex.lab=1.5,cex.main=1.5)
    ###plot(tsne.out[["30"]]$Rtsne.res$Y[,1],tsne.out[[as.character(perplexity)]]$Rtsne.res$Y[,2], 
    plot(tsne.out[[as.character(perplexity)]]$Rtsne.res$Y, 
         t='p',pch=16,col="lightgray", main="BarnesHutSNE",xlab="Dim1",ylab="Dim2",cex.lab=2.0,cex.axis=1.5,cex.main=2.0)
    for(m.cls in unique(sample.desc$majorCluster)){
        f.points <- sample.desc$majorCluster==m.cls
        ###tsne.density <- kde2d(tsne.out[["30"]]$Rtsne.res$Y[f.points,1],tsne.out[["30"]]$Rtsne.res$Y[f.points,2],n=50)
        ###contour(tsne.density,xlab="",ylab="",add=TRUE,nlevels = 3,drawlabels=F)
        tsne.density <- kde(tsne.out[[as.character(perplexity)]]$Rtsne.res$Y[f.points,])
        plot(tsne.density, cont=c(80),lwd=4,col=majorClusterColor[m.cls],add=T,drawlabels=F)

    }
    dev.off()

    pdf(file=sprintf("%s.integrated.perplexity%d.tsne.colPoints.pdf",out.prefix,perplexity),width=8,height=8)
    par(mar=c(5,5,4,2),cex.lab=1.5,cex.main=1.5)
    ###plot(tsne.out[["30"]]$Rtsne.res$Y[,1],tsne.out[[as.character(perplexity)]]$Rtsne.res$Y[,2], 
    plot(tsne.out[[as.character(perplexity)]]$Rtsne.res$Y, 
         t='p',pch=16,col="lightgray", main="BarnesHutSNE",xlab="Dim1",ylab="Dim2",cex.lab=2.0,cex.axis=1.5,cex.main=2.0)
    f.points <- !grepl("uncharacterized",sample.desc$majorCluster,perl=T)
    points(tsne.out[[as.character(perplexity)]]$Rtsne.res$Y[f.points,],
           col=majorClusterColor[sample.desc$majorCluster[f.points]],pch=16)
    dev.off()
    ####legend("bottom",horiz = T,legend=names(leafClusterColor),fill = NULL,inset = -0.01,xpd = NA,cex=1.5,pch=16,border =NA,col = leafClusterColor)
    #opar <- par(mar=c(5,11,4,8),cex.lab=1.5,cex.main=1.5)
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
    ###
    mapping.leafToMajor <- unique(sample.desc[,c("leafCluster","majorCluster")])
    mapping.leafToMajor.v <- mapping.leafToMajor$majorCluster
    names(mapping.leafToMajor.v) <- mapping.leafToMajor$leafCluster

    #hc.merged <- set(hc.merged,"branches_k_col",leafClusterColor[leaf.cluster[!grepl("uncharacterized",leaf.cluster)][order.dendrogram(hc.merged)]])
    hc.merged <- set(hc.merged,"branches_k_col",majorClusterColor[mapping.leafToMajor.v[labels(hc.merged)]])
    hc.merged <- set(hc.merged,"branches_lwd",3)
   
    pdf(file=sprintf("%s.integrated.perplexity%d.tree.pdf",out.prefix,perplexity),width=10,height=8)
    layout(matrix(c(1,2,3),ncol=3,byrow = T),widths = c(0.15,0.30,0.25),heights = c(1))
    par(mar=c(6,2,4,0),cex.lab=1.5,cex.main=1.5)
    plot(hc.merged,horiz=T,yaxt="n",leaflab="none")
    ##plot(hc.merged,horiz=T,yaxt="n",leaflab="none",edgePar = list(lwd = 1.5:2))

    ### barplot describing cell origin
    dat.plot <- table(sample.desc[,c("sampleType","leafCluster")])
    dat.plot <- dat.plot[,labels(hc.merged)]
    ########dat.plot <- apply(dat.plot,2,function(x){x/sum(x)})
    ###par(mar=c(5,11,4,8),cex.lab=1.5,cex.main=1.5)
    par(mar=c(6,1,4,4),cex.lab=1.5,cex.main=1.5)
    xx <- barplot(dat.plot,horiz = T,beside=F,col=sampleTypeColor[rownames(dat.plot)],yaxt="n",
                  sub = "Cell Origin",cex.axis = 1.5,cex.sub = 2.0)
    ####staxlab(2,at = xx,labels=colnames(dat.plot),srt=0, cex=1.0,adj=1,top.line=0.5)
    #text(rep(par("usr")[1]-10,length(xx)),xx,labels=colnames(dat.plot), cex=0.95,adj=1,xpd=NA)
    points(x=rep(-20,length(xx)),y=xx,cex=3,col=majorClusterColor[mapping.leafToMajor.v[colnames(dat.plot)]],pch=16,xpd=NA)
    legend("bottomright",legend=rownames(dat.plot),fill=sampleTypeColor[rownames(dat.plot)],
           border=sampleTypeColor[rownames(dat.plot)],cex=1.5,inset=c(-0.08,0.05),xpd=NA,horiz = F)
    ###legend("topright",legend=rownames(dat.plot),fill=sampleTypeColor[rownames(dat.plot)],
    ###       border=sampleTypeColor[rownames(dat.plot)],cex=1.5,inset=c(-0.30,0),xpd=NA)
   
    ### scater: leafCluster ~ sampleType, by patient
    for(pp in unique(sample.desc$patient)){
        f.CD8 <- grepl("^CD8",sample.desc$leafCluster,perl = T)
        f.CD4 <- grepl("^CD4",sample.desc$leafCluster,perl = T)
        pdf(file=sprintf("%s.integrated.perplexity%d.%s.scatter.pdf",out.prefix,perplexity,pp),
            width=8,height=5*(as.numeric(sum(f.CD8)>0) + as.numeric(sum(f.CD4)>0)) )
        par(mar=c(5,5,4,2),cex.lab=1.5,cex.main=1.5)
        pushViewport(viewport(layout=grid.layout( as.numeric(sum(f.CD8)>0) + as.numeric(sum(f.CD4)>0),1)))
        if(sum(f.CD8)>0){
            dat.plot <- table(subset(sample.desc,patient==pp & f.CD8, select=c("sampleType","leafCluster")))
            dat.plot <- melt(dat.plot)
            dat.plot <- subset(dat.plot,value>0)
            print(ggplot(dat.plot, aes(sampleType, leafCluster)) + geom_point(aes(size = value)) + scale_size_area(max_size = 12) +
                  theme_minimal() + ggtitle(sprintf("%s (CD8)",pp)) + xlab("") + ylab(""), 
                  vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
        }
        if(sum(f.CD4)>0){
            dat.plot <- table(subset(sample.desc,patient==pp & f.CD4,select=c("sampleType","leafCluster")))
            dat.plot <- melt(dat.plot)
            dat.plot <- subset(dat.plot,value>0)
            print(ggplot(dat.plot, aes(sampleType, leafCluster)) + geom_point(aes(size = value)) + scale_size_area(max_size = 12) +
                  theme_minimal() + ggtitle(sprintf("%s (CD4)",pp)) + xlab("") + ylab(""), 
                  vp=viewport(layout.pos.row = 1+as.numeric(sum(f.CD8)>0), layout.pos.col = 1))
        }
        dev.off()
    }
    ### bar plot descrie subtype in TTC and TTH(R)
    s.table <- aggregate(sample~patient+sampleType+majorCluster,data=sample.desc,FUN=length)
    s.table$patient <- factor(s.table$patient,
                                    levels=c("P0322", "P1116", "P0508", "P0205", "P0407"))
    s.table$majorCluster <- factor(s.table$majorCluster,
                                    levels=c("CD8_RGS1.exhausted","CD8_RGS1.GZMK","CD8_SLC4A10","CD8_CX3CR1","CD8_naive",
                                             "CD4_GZMA.cytotoxic","CD4_GZMA.main","CD4_GZMA.exhausted","CD4_P.other","CD4_P.PTR","CD4_FOXP3-CTLA4"))
    dat.plot <- s.table[grepl("^TTC",s.table$sampleType,perl=T) & s.table$majorCluster!="CD8_naive",,drop=F]
    if(nrow(dat.plot)>0){
        dat.plot <- acast(dat.plot,majorCluster~patient,value.var="sample")
        dat.plot[is.na(dat.plot)] <- 0
        ##apply(dat.plot,2,sum)
        dat.plot <- apply(dat.plot,2,function(x){ x/sum(x) } )
        print(dat.plot)
        dat.plot.df <- data.frame(subtype=rownames(dat.plot))
        dat.plot.df <- cbind(dat.plot.df,dat.plot)
        write.table(dat.plot.df,sprintf("%s.integrated.perplexity%d.CD8.subtypeDist.txt",out.prefix,perplexity),quote = F,sep = "\t",row.names = F)

        pdf(sprintf("%s.integrated.perplexity%d.CD8.subtypeDist.pdf",out.prefix,perplexity),width=8,height=8)
        par(mar=c(5,6,4,11),cex.lab=2,cex.main=2,cex.axis=1.5)
        xx <- barplot(dat.plot,beside = F,xlab="",ylab="Cell Subtype",main="CD8+ T cell",col=majorClusterColor[rownames(dat.plot)],xaxt="n")
        #staxlab(1,at = xx,labels=names(dat.plot),srt=if(length(dat.plot) >= 5) 45 else 0, cex=1.3,adj=0.5,top.line=2)
        staxlab(1,at = xx,labels=colnames(dat.plot),srt=if(length(dat.plot) >= 5) 45 else 0, cex=1.3)
        legend("topright",legend=rownames(dat.plot),fill=majorClusterColor[rownames(dat.plot)],inset = c(-0.45,0),xpd = T,horiz = F)
        dev.off()
    }
    dat.plot <- s.table[grepl("^(TTH|TTR)",s.table$sampleType,perl=T) & !grepl("^CD4_P",s.table$majorCluster,perl=T),,drop=F]
    if(nrow(dat.plot)>0){
        dat.plot <- acast(dat.plot,majorCluster~patient,value.var="sample",fun.aggregate=sum)
        dat.plot[is.na(dat.plot)] <- 0
        ##apply(dat.plot,2,sum)
        dat.plot <- apply(dat.plot,2,function(x){ x/sum(x) } )
        print(dat.plot)
        dat.plot.df <- data.frame(subtype=rownames(dat.plot))
        dat.plot.df <- cbind(dat.plot.df,dat.plot)
        write.table(dat.plot.df,sprintf("%s.integrated.perplexity%d.CD8.subtypeDist.txt",out.prefix,perplexity),quote = F,sep = "\t",row.names = F)

        pdf(sprintf("%s.integrated.perplexity%d.CD4.subtypeDist.pdf",out.prefix,perplexity),width=8,height=8)
        par(mar=c(5,6,4,11),cex.lab=2,cex.main=2,cex.axis=1.5)
        xx <- barplot(dat.plot,beside = F,xlab="",ylab="Cell Subtype",main="CD4+ T cell",col=majorClusterColor[rownames(dat.plot)],xaxt="n")
        #staxlab(1,at = xx,labels=names(dat.plot),srt=if(length(dat.plot) >= 5) 45 else 0, cex=1.3,adj=0.5,top.line=2)
        staxlab(1,at = xx,labels=colnames(dat.plot),srt=if(length(dat.plot) >= 5) 45 else 0, cex=1.3)
        legend("topright",legend=rownames(dat.plot),fill=majorClusterColor[rownames(dat.plot)],inset = c(-0.45,0),xpd = T,horiz = F)
        dev.off()
    }

    ### barplot describing patient
    dat.plot <- table(sample.desc[,c("patient","leafCluster")])
    dat.plot <- dat.plot[,labels(hc.merged)]
    par(mar=c(6,0,4,4),cex.lab=1.5,cex.main=1.5)
    patientColor <- auto.colSet(n=nrow(dat.plot),name="Accent")
    xx <- barplot(dat.plot,horiz = T,beside=F,col=patientColor,yaxt="n",
                  sub = "Patient",cex.axis = 1.5,cex.sub = 2.0)
    legend("bottomright",legend=rownames(dat.plot),fill=patientColor,
           border=patientColor,cex=1.5,inset=c(-0.13,0.05),xpd=NA)
    dev.off()

    ### legend
    pdf(file=sprintf("%s.integrated.perplexity%d.legend.pdf",out.prefix,perplexity),width=16,height=2)
    par(mar=c(1,1,2,1),cex.lab=1.5,cex.main=1.5)
    plot.new()
    ##idx.leaf <- ceiling(length(colnames(dat.plot))/2)
    #leafCluster.name <- colnames(dat.plot)
    #leafCluster.name <- paste0(leafCluster.name,strrep("x",max(nchar(leafCluster.name))-nchar(leafCluster.name)))
    majorCluster.name <- unique(mapping.leafToMajor.v[colnames(dat.plot)])
    print(majorCluster.name)
    ii <- 0
    mm <- 3
    nn <- 4
    for(yy in head(1-seq(0,1,1/mm),n=mm)){ 
        for(xx in head(seq(0,1,1/nn),n=nn)) {
            ii <- ii+1
            text(x = xx+0.030,y = yy,labels = majorCluster.name[ii],cex=1.6,adj = c(0,1),offset = 4)
            points(xx+0.012,yy-0.050,cex=4,pch=16,xpd=NA,col=majorClusterColor[majorCluster.name[ii]])
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
