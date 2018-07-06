#!/usr/bin/env Rscript


args <- commandArgs(T)
if(length(args)<6)
{
    cat("plot.Marker.fromTab.R <cellType.color file> <sample id> <design file> <output.prefix> <in.tab> <marker list> [<in.tab2> <marker list2>]/[extra data]\n")
    q()
}

cellType.color.file  <- args[1]
sampleID <- args[2]
design.file <- args[3]
out.dir <- args[4]
tab.file <- args[5]
marker.file <- args[6]
if(length(args)==8)
{
    tab2.file <- args[7]
    marker2.file <- args[8]
}else if(length(args)==7)
{
    tab2.file <- args[7]
    marker2.file <- NULL
}else
{
    tab2.file <- NULL
    marker2.file <- NULL
}
expT <- 3

source("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")

########## TEST INPUT ######

#cellType.color.file  <- "/WPS1/zhenglt/work/TCR_chunhong/data/CellType.color"
#sampleID <- "P1202"
#design.file <- "/WPS1/zhenglt/work/TCR_chunhong/phase31/sample.design.phase31.P1202.txt"
#out.dir <- "/WPS1/zhenglt/work/TCR_chunhong/phase31/qc/marker/test"
#tab.file <- "/WPS1/zhenglt/work/TCR_chunhong/phase31/stat/P1202.phase31.RPKM.tab.gz" 
#marker.file <- "/WPS1/zhenglt/work/TCR_chunhong/phase31/qc/marker/RPKM.filter0/P1202/P1202.marker.TCellType/Marker.n7.TCellType.txt"
#tab2.file <- NULL
#marker2.file  <- NULL
########## TEST      ######

loginfo(paste0("process sample ",sampleID))
loginfo(paste0("mkdir ",out.dir," ..."))
dir.create(out.dir,recursive = T,showWarnings = F)

in.table <- read.table(tab.file,sep="\t",stringsAsFactors = F,header = T,row.names = 1,check.names = F)
in.table <- as.matrix(in.table[,-1])
marker.list <- read.table(marker.file,sep="\t",stringsAsFactors = F,header = F,row.names = 1,check.names = F)
names(marker.list) <- c("MName")
print(str(in.table))

if(!is.null(tab2.file) && !is.null(marker2.file))
{
  in2.table <- read.table(tab2.file,sep="\t",stringsAsFactors = F,header = T,row.names = 1,check.names = F)
  in2.table <- as.matrix(in2.table[,-1])
  marker2.list <- read.table(marker2.file,sep="\t",stringsAsFactors = F,header = F,row.names = 1,check.names = F)
  names(marker2.list) <- c("MName")
}else if(!is.null(tab2.file) && is.null(marker2.file))
{
  #in2.table <- read.table(tab2.file,sep="\t",stringsAsFactors = F,header = T,row.names = 1,check.names = F)
  #in2.table <- t(in2.table)
  in2.table <- read.table(tab2.file,header=T,check.names = F,stringsAsFactors = F,row.names=1)
  marker2.list <- NULL
}else
{
    in2.table <- NULL
    marker2.list <- NULL
}


###myDesign<-read.table(design.file,header=T,row.names="sample",check.names=F,colClasses=c("factor","character","factor","character"))
myDesign<-read.table(design.file,header=T,row.names="sample",check.names=F)
myDesign$sampleType <- factor(myDesign$sampleType,levels=unique(myDesign$sampleType))
sampleTypeColor <- read.SampleTypeColor(cellType.color.file)

suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(colorspace))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(gridBase))
suppressPackageStartupMessages(library(fields))
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
      sampleType.n.levels <- length(unique(designM$sampleType))
      #sampleType.n.levels <- length(levels(designM$sampleType))
      if(length(levels(designM$sampleType)) <= 9)
      {
        colSet <- brewer.pal(sampleType.n.levels,"Set1")
        #colSet <- brewer.pal(sampleType.n.levels,"Paired") ### n==12
        annCol <<- structure(colSet, 
                            names =levels(designM$sampleType) [levels(designM$sampleType) %in% unique(as.character(designM$sampleType))])
      }else
      {
        colSet <- rainbow(sampleType.n.levels)
        annCol <- structure(colSet, 
                            names = levels(designM$sampleType) [levels(designM$sampleType) %in% unique(as.character(designM$sampleType))])
      }
  }else
  {
      annCol <- colSet[names(colSet) %in% levels(designM$sampleType)]
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
  if(!is.null(X2) && !is.null(X2.ann))
  {
    ## select the genes in X2.ann only; make the X2 and designM consistant
    X2 <- X2[rownames(X2.ann),rownames(designM),drop=F]
    tf <- data.frame(geneID=rownames(X2),geneSymbol=X2.ann$MName)
    tf <- cbind(tf,X2)
    if(all.equal(colnames(dat.plot.export),colnames(tf)) != TRUE)
    {
        stop("dat.plot.export and tf must have the same sample order")
    }
    dat.plot.export <- rbind(dat.plot.export,tf)
    rownames(X2) <- X2.ann$MName
    dat.plot <- rbind(dat.plot,X2)
  }
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
      dat.plot.gene.i <- dat.plot[x, rownames(designM)[designM$sampleType==names(annCol)[i]]]
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
    dat.plot.i <- (dat.plot[, rownames(designM)[designM$sampleType==names(annCol)[i]], drop = FALSE])
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
    #print(names(annCol)[i])
    #print(dim(dat.plot))
    dat.plot.i <- dat.plot[, rownames(designM)[designM$sampleType==names(annCol)[i]], drop = FALSE]
    if(ncol(dat.plot.i)>2)
    {
        ann.bar.dat <- apply(dat.plot.i,1,function(x){ sum(x>expT)/length(x)})

        top_annotation_height <- unit(1.5, "cm")
        if(!is.null(X2) && is.null(X2.ann))
        {
            annDF <- data.frame(sampleType=rep(names(annCol)[i],ncol(dat.plot.i)))
            #print(dim(X2))
            #print(dim(dat.plot.i))
            annDF <- cbind(annDF,X2[colnames(dat.plot.i),,drop = FALSE])
            annColList <- list(sampleType=annCol)
            invisible(lapply(colnames(X2),function(x){
                                            extra.ann.type <- unique(X2[,x])
                                            annColList[[x]] <<- brewer.pal(length(extra.ann.type),"Dark2")
                                            names(annColList[[x]]) <<- extra.ann.type
                                            }))
            #print(annColList)
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

#plotC(in.table,myDesign,out.dir,n=NULL,gfilter=rownames(marker.list),extra=paste0(".",sampleID))
plotC(in.table,myDesign,out.dir,colSet=sampleTypeColor,n=NULL,gfilter=rownames(marker.list),extra=paste0(".",sampleID),
      cluster_rows = F,row_names_side = "left",X2 = in2.table,X2.ann = marker2.list,sort.by.prevalence=T)
