#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-i", "--inFile", type="character", required=TRUE, help="infile")
parser$add_argument("-a", "--sampleDescFile", type="character", required=TRUE, help="sample.desc.file")
parser$add_argument("-b", "--geneFile", type="character", required=TRUE, help="gene list file")
parser$add_argument("-o", "--outDir", type="character", required=TRUE, help="output dir")
parser$add_argument("-s", "--sample", type="character", default="SAMPLE", help="sample id  [default %(default)s]")
parser$add_argument("-t", "--expThreshold", type="double", default="3", help="expression threshold  [default %(default)s]")
parser$add_argument("-n", "--ncores", type="integer", default="8", help="num of cores  [default %(default)s]")
parser$add_argument("-k", "--groupBy", type="character",default="leafCluster", help="group by")
#parser$add_argument("-y", "--cellTypeColorFile", type="character", default="/WPS1/zhenglt/work/TCR_chunhong/data/CellType.color", 
#                    help="cellTypeColorFile [default %(default)s]")
parser$add_argument("-w", "--clusterColorFile", type="character", 
                    help="clusterColorFile [default %(default)s]")
parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="verbose mode [default %(default)s]")
args <- parser$parse_args()
print(args)

in.file <- args$inFile
sample.desc.file <- args$sampleDescFile
gene.file <- args$geneFile
out.dir <- args$outDir
sample.id <- args$sample
ncores <- args$ncores
##cellTypeColorFile <- args$cellTypeColorFile
clusterColorFile <- args$clusterColorFile
mode.verbose <- args$verbose
expT <- args$expThreshold
args.groupBy <- args$groupBy

source("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")

########## TEST INPUT ######
#in.file <- "/WPSnew/zhenglt/work/proj_h/byPhase/phase01/stat/LCPD5.phase01.TPM.tab.gz"
#sample.desc.file <- "/WPSnew/zhenglt/work/proj_h/byPhase/phase01/sample.design.phase01.LCPD5.txt"
#gene.file <- "/WPSnew/zhenglt/work/proj_h/byPhase/phase01/qc/marker/markerList.list"
#out.dir <- "/WPSnew/zhenglt/work/proj_h/byPhase/phase01/qc/marker/test"
#sample.id <- "LCPD5"
#ncores <- 8
#args.groupBy <- "sampleType"
#clusterColorFile <- NULL
#in.file <- "/WPS1/zhenglt/work/proj_xy/integrated/cross.patient/exp/lung.P9.scran.pool.fltGeneT.fltSampleT.tpm.tab.gz"
#sample.desc.file  <- "/WPS1/zhenglt/work/proj_xy/integrated/cross.patient/clusteringResult/lungC.all.leaf.addMajor.final.forPlot.txt"
#gene.file <- "/WPS1/zhenglt/work/proj_xy/integrated/cross.patient/cytokine/geneSet.fromGO.list"
#out.dir <- "/WPS1/zhenglt/work/proj_xy/integrated/cross.patient/cytokine/test/P1202"
#sample.id <- "P1202"
#clusterColorFile <- "/WPS1/zhenglt/work/proj_xy/integrated/cross.patient/clusteringResult/lungC.tSNE.colorByCluster.n1000/lungC.useAOV.n1000.lungC.majorClusterColor.txt"
#mode.verbose <- T
#expT <- 3
#args.groupBy <- "majorCluster"

########## TEST      ######

loginfo(paste0("process sample ",sample.id))
loginfo(paste0("mkdir ",out.dir," ..."))
dir.create(out.dir,recursive = T,showWarnings = F)

### expression data
in.table <- read.table(in.file,sep="\t",stringsAsFactors = F,header = T,row.names = 1,check.names = F)
g.GNAME <- in.table[,1]
names(g.GNAME) <- rownames(in.table)
in.table <- as.matrix(in.table[,-1])
### marker list
marker.file.list <- read.table(gene.file,sep="\t",stringsAsFactors = F,header = F,row.names = 1,check.names = F)
marker.list <- lapply(seq_len(nrow(marker.file.list)),function(i){
           .ginfo <- read.table(marker.file.list[i,1],sep = "\t",stringsAsFactors = F,header = T,check.names = F)
           if(!"geneID" %in% colnames(.ginfo)){
            .ginfo <- read.table(marker.file.list[i,1],sep = "\t",stringsAsFactors = F,header = F,check.names = F)
            colnames(.ginfo) <- c("geneID","geneSymbol")
           }
           .ginfo$geneID <- as.character(.ginfo$geneID)
           rownames(.ginfo) <- .ginfo$geneID
           .ginfo })
names(marker.list) <- rownames(marker.file.list)
### sample desc file
myDesign<-read.table(sample.desc.file,header=T,check.names=F,stringsAsFactors = )
rownames(myDesign) <- myDesign$sample
f.sample <- intersect(rownames(myDesign),colnames(in.table))
in.table <- in.table[,f.sample]
myDesign <- myDesign[f.sample,]
### group color
clusterColor <- read.SampleTypeColor(clusterColorFile)

suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(colorspace))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(gridBase))
suppressPackageStartupMessages(library(fields))
suppressPackageStartupMessages(require("plyr"))
suppressPackageStartupMessages(require("doParallel"))
#library(GetoptLong)


myRowBarAnnotation <- function (x, gp = gpar(fill = "#CCCCCC"), ...) 
{
  require(gridBase)
  x = x
  factor = 0.6
  #data_scale = range(x)
  #print(data_scale)
  #data_scale = data_scale + c(-0.05, 0.05) * (data_scale[2] - data_scale[1])
  #print(data_scale)
  data_scale=c(0,1)
  breaks = grid.pretty(data_scale)

  row = function(index) {
    n = length(index)
    if (n != length(x)) {
      stop(paste0("Length of index should be ", length(x)))
    }
    pushViewport(viewport(xscale = data_scale, yscale = c(0.5, n + 0.5)))
    #grid.rect()
    grid.rect(x = data_scale[1], y = seq_along(index), 
              width = x[index] - data_scale[1], 
              height = 1 * factor, just = "left", 
              default.units = "native",gp=gp)
    grid.xaxis(at = breaks, label = breaks, main = FALSE, gp = gpar(fontsize = 15))
    upViewport()
  }
}

plotC<-function(X,designM,outDir,colSet=NULL,n=50,extra="",filter=NULL, gfilter=NULL,X2=NULL,X2.ann=NULL,sort.by.prevalence=FALSE,...)
{
  if(!is.null(filter))
  {
    X <- X[,filter,drop=F]
    designM <- designM[filter,,drop=F]
  }
  ## make the X and designM consistant
  X <- X[,rownames(designM),drop=F]
  if(!is.null(gfilter))
  {
    X <- X[rownames(X) %in% gfilter[!is.na(gfilter)],,drop=F]
  }
  
  if(is.null(colSet))
  {
      annCol <- auto.colSet(n = length(unique(designM[,args.groupBy])))
      names(annCol) <- unique(designM[,args.groupBy])
  }else
  {
      annCol <- colSet[names(colSet) %in% levels(designM[,args.groupBy])]
  }
  print(annCol)
  ## heatmap by most variable genes
  #if(is.null(n))
  #{
  #  n=nrow(X)
  #}
  #rowVar <- apply(X,1,var)
  #select <- order(rowVar,decreasing = T)[1:n]
  #dat.plot <-X[select,] 
  dat.plot <- X
  dat.plot.export <- data.frame(geneID=rownames(dat.plot),geneSymbol=entrezToXXX(rownames(dat.plot)))
  dat.plot.export <- cbind(dat.plot.export,dat.plot)
  rownames(dat.plot) <- entrezToXXX(rownames(dat.plot))
  m <- ncol(dat.plot)
  n <- nrow(dat.plot)  
  write.table(dat.plot.export,paste0(outDir,"/Marker.n",n,extra,".txt",sep=""),sep="\t",row.names = F,quote = F)
  
  ### density plot
  ##dat.plot <- log10(dat.plot + 1e-6)
  dat.plot <- log2(dat.plot + 1)
  cat(paste("generate figure: ",outDir,"/Marker.n",n,extra,".densityplot.pdf\n",sep=""))
  pdf(paste(outDir,"/Marker.n",n,extra,".densityplot.pdf",sep=""),width=8,height=8)
  par(mar=c(5,6,4,2),cex.lab=2)
  sapply(rownames(dat.plot), function(x){
    dat.plot.gene <- dat.plot[x,]
    dd <- list()
    yy <- c()
    dd.all <- density(dat.plot.gene)
    yy <- c(yy, max(dd.all[["y"]]))
    for(i in seq_along(names(annCol)))
    {
      dat.plot.gene.i <- dat.plot[x, rownames(designM)[designM[,args.groupBy]==names(annCol)[i]]]
      if(length(dat.plot.gene.i) > 2)
      {
          dd[[names(annCol)[i]]] <- density(dat.plot.gene.i)
          yy <- c(yy, max(dd[[names(annCol)[i]]][["y"]]))
      }
    }
    ##plot(density(dat.plot.gene),lwd=2,xlab=expression(log[10](RPKM + epsilon)),main=x,cex.lab=2,ylim=c(0,max(yy)))
    ##plot(density(dat.plot.gene),lwd=2,xlab=expression(log[10](Exp + epsilon)),main=x,cex.lab=2,ylim=c(0,max(yy)))
    plot(density(dat.plot.gene),lwd=2,xlab=expression(log[2](Exp + 1)),main=x,cex.lab=2,ylim=c(0,max(yy)))
    lapply(names(dd),function(n){
        lines(dd[[n]],col=annCol[n],lwd=2)
    }) 
    legend("topright",legend=names(annCol),fill=annCol)
  })
  dev.off()
  #return(0)
  
  #dat.plot[dat.plot < 1]=0
  #dat.plot[dat.plot >= 1]=1
  ## log scale
  #eFilter <- (dat.plot < -5.99)
  eFilter <- (dat.plot < 0)
  #dat.plot[eFilter]=0
  #dat.plot[!eFilter]=1

  if(sort.by.prevalence){
    p.order <- order(apply(dat.plot,1,function(x){ sum(x>expT)  }),decreasing = T)
    dat.plot <- dat.plot[p.order,,drop=F]
  }

  ### simple barplot
  cat(paste("generate figure: ",outDir,"/Marker.n",n,extra,".barplot.pdf\n",sep=""))
  pdf(paste(outDir,"/Marker.n",n,extra,".barplot.pdf",sep=""),width=8,height=6)
  par(mar=c(6,4,3,2))
  for(i in seq_along(names(annCol)))
  {
    dat.plot.i <- (dat.plot[, rownames(designM)[designM[,args.groupBy]==names(annCol)[i]], drop = FALSE])
    if(ncol(dat.plot.i)>2)
    {
        ann.bar.dat <- apply(dat.plot.i,1,function(x){ sum(x>expT)/length(x)})
        xx <- barplot(ann.bar.dat,col=annCol[i],ylim=c(0,1),cex.sub=1.2,sub=paste0("NumOfSample=",ncol(dat.plot.i)),ylab=names(annCol)[i],add = F, xaxt="n")
        text(xx,par("usr")[3]-0.02, srt = 60, adj = 1, labels = names(ann.bar.dat), xpd = NA,cex=34/max(nrow(dat.plot.i),34))
    }
  }
  dev.off()
  #return(0)
  ### Big heatmap
  if(!is.null(X2) && is.null(X2.ann))
  {
      X2 <- X2[rownames(designM),,drop = FALSE]
  }
  cat(paste("generate figure: ",outDir,"/Marker.n",n,extra,".pdf\n",sep=""))
  pdf(paste(outDir,"/Marker.n",n,extra,".pdf",sep=""),width=length(annCol)*8.5,height=13)
  par(mar=c(6,4,5,8)+0.1)
  plot.new()
  legend("bottom",legend=names(annCol),fill=annCol,border=annCol,cex=2.5,horiz=T,inset=-0.1,xpd=T)
  image.plot(zlim=c(0,12),legend.only=TRUE,col=colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),legend.width=3.0,legend.mar=7.0,axis.args=list(cex.axis=2.5))
  vps <- baseViewports()
  pushViewport(vps$inner, vps$figure, vps$plot)
  hh <- NULL
  for(i in seq_along(names(annCol)))
  {
    dat.plot.i <- dat.plot[, rownames(designM)[designM[,args.groupBy]==names(annCol)[i]], drop = FALSE]
    if(ncol(dat.plot.i)>2)
    {
        ann.bar.dat <- apply(dat.plot.i,1,function(x){ sum(x>expT)/length(x)})

        top_annotation_height <- unit(1.5, "cm")
        if(!is.null(X2) && is.null(X2.ann))
        {
            annDF <- data.frame(sampleType=rep(names(annCol)[i],ncol(dat.plot.i)))
            annDF <- cbind(annDF,X2[colnames(dat.plot.i),,drop = FALSE])
            annColList <- list(sampleType=annCol)
            invisible(lapply(colnames(X2),function(x){
                                            extra.ann.type <- unique(X2[,x])
                                            annColList[[x]] <<- brewer.pal(length(extra.ann.type),"Dark2")
                                            names(annColList[[x]]) <<- extra.ann.type
                                            }))
            ha.col <- HeatmapAnnotation(df = annDF,
                                    col = annColList, show_legend = (i==length(annCol)))
            top_annotation_height <-  unit(1.5*(ncol(X2)+1), "cm")
        }else if(is.null(X2) && is.null(X2.ann))
        {
            annDF <- data.frame(sampleType=rep(names(annCol)[i],ncol(dat.plot.i)))
            annColList <- list(sampleType=annCol)
            ha.col <- HeatmapAnnotation(df = annDF, col = annColList, show_legend = F)
            ###ha.col <- HeatmapAnnotation(df = annDF, col = annColList, show_legend = (i==length(annCol)))
        }
        #ha.row <- HeatmapAnnotation(barplot = myRowBarAnnotation(ann.bar.dat, gp = gpar(fill = "#008000")),
        #                           which = "row", width = unit(6, "cm"))
        ### !!!  surprising things !!! should reverse the ann.bar.dat in version ComplexHeatmap1.6.0
        ha.row <- HeatmapAnnotation(barplot = myRowBarAnnotation(ann.bar.dat[length(ann.bar.dat):1], gp = gpar(fill = "#008000")),
                                   which = "row", width = unit(6, "cm"))

        ###test.out[[i]]<<-dat.plot.i
        ###test.out.dend[[i]] <<- as.dendrogram(hclust(dist(t(dat.plot.i))))
        # col=c("1"="red","0"="blue"),
        z.lo <- 0
        z.hi <- 12
        dat.plot.i.dend <- as.dendrogram(hclust(dist(t(dat.plot.i))))
        if(attr(dat.plot.i.dend,"height")==0){ dat.plot.i.dend <- F }
        ht <- Heatmap(dat.plot.i,
                      col = colorRamp2(seq(z.lo,z.hi,length=100), 
                                     colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100), space="LAB"),
                      column_names_gp = gpar(fontsize = 12*55/m),row_names_gp = gpar(fontsize = 12*55/max(n,25)),
                      column_dend_height = unit(6, "cm"), row_dend_width = unit(6, "cm"),
                      cluster_columns = dat.plot.i.dend,
                      show_heatmap_legend = F , row_names_max_width = unit(10,"cm"),
                      heatmap_legend_param = list(grid_width = unit(1.5, "cm"), 
                                                  grid_height = unit(1.5, "cm"),
                                                  title_gp = gpar(fontsize = 24, fontface = "bold"),
                                                  label_gp = gpar(fontsize = 20), color_bar = "continuous"),
                      top_annotation_height = top_annotation_height,
                      top_annotation = ha.col,...)
        if(is.null(hh))
        {
          hh <- ht+ha.row
        }else
        {
          hh <- hh+ht+ha.row
        }
    }
  }
  draw(hh, newpage=FALSE)
  dev.off()
}

print(in.table[1:4,1:8])
print(head(myDesign))
registerDoParallel(cores = ncores)
ret <- ldply(seq_len(length(marker.list)),function(i){
        dir.create(sprintf("%s/%s",out.dir,names(marker.list)[i]),showWarnings = F,recursive = T)
        plotC(in.table,myDesign,sprintf("%s/%s",out.dir,names(marker.list)[i]),colSet=clusterColor,n=NULL,gfilter=marker.list[[i]][,"geneID"],extra=paste0(".",sample.id),
              cluster_rows = F,row_names_side = "left",X2 = NULL,X2.ann = NULL,sort.by.prevalence=T)
    },.progress = "none",.parallel=T)

loginfo("end.")

