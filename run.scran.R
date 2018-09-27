#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-i", "--infile", type="character", required=TRUE, help="infile")
parser$add_argument("-t", "--tpmfile", type="character", required=TRUE, help="tpmfile")
parser$add_argument("-o", "--outprefix", type="character", required=TRUE, help="outprefix")
parser$add_argument("-e", "--ensembl", action="store_true", default=FALSE, help="whether geneID is ensembl id [default %(default)s]")
parser$add_argument("-s", "--sample", type="character", default="SAMPLE", help="sample id  [default %(default)s]")
parser$add_argument("-m", "--method", type="character", default="pool", help="normalization method  [default %(default)s]")
parser$add_argument("-p", "--quickClustMethod", type="character", default="hclust", help="quick cluster method (igraph or hclust) [default %(default)s]")
parser$add_argument("-d", "--designFile", required=TRUE, type="character", help="sample desc file")
parser$add_argument("-a", "--filterGene", action="store_true", default=FALSE, help="whether filterGene [default %(default)s]")
parser$add_argument("-r", "--filterGeneMethod", type="character", default="AVE", help="method used to filterGene (AVE, RC, TPM) [default %(default)s]")
parser$add_argument("-u", "--nmad", type="integer", default=3, help="median - [INT]*mad as lower threshold for filtering by read counts.[default %(default)s]")
parser$add_argument("-k", "--numOfCells", type="integer", default=10, 
                    help="keep only genes expressed in more than [INT] cells. Used only when filterGeneMethod is RC [default %(default)s]")
parser$add_argument("-n", "--numOfReads", type="integer", default=0,
                    help="genes are considered expressed if reads/TPM larger than [INT]. Used only when filterGeneMethod is RC/TPM [default %(default)s]")
parser$add_argument("-w", "--sizeFactor", type="double", default=0.1, help="size factor threshold [default %(default)s]")
parser$add_argument("-b", "--filterSample", action="store_true", default=FALSE, help="whether filterSample [default %(default)s]")
parser$add_argument("-c", "--outlier", type="character", help="predefined outlier list, [default %(default)s]")
parser$add_argument("-z", "--sizes", type="character",default="20,40,60,80,100", help="sizes used for scran's sf estimation, [default %(default)s]")
parser$add_argument("-q", "--quickCluster", action="store_true", default=FALSE, help="whether do quick clustering, used for scran's sf estimation [default %(default)s]")
parser$add_argument("-j", "--species", type="character",default="human", help="species (human, mouse) [default %(default)s]")
parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="whether verbose mode [default %(default)s]")
args <- parser$parse_args()
print(args)

in.file <- args$infile
tpm.file <- args$tpmfile
out.prefix <- args$outprefix
args.ensembl <- args$ensembl
sample.id <- args$sample
arg.method <- args$method
sample.desc.file <- args$designFile
arg.filterGene <- args$filterGene
arg.filterSample <- args$filterSample
pre.outlier.file <- args$outlier
arg.sizes <- as.integer(unlist(strsplit(x = args$sizes,split = ",")))
arg.do.quickCluster <- args$quickCluster
arg.numOfCells <- args$numOfCells
arg.numOfReads <- args$numOfReads
arg.filterGeneMethod <- args$filterGeneMethod
####args.species <- "human"
args.nmad <- args$nmad
args.species <- args$species
args.quickClustMethod <- args$quickClustMethod
args.sizeFactor <- args$sizeFactor

source("/WPSnew/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")
library("data.table")
library("dplyr")
library("scater")
library("BiocParallel")
library("scran")
library("plotrix")
library("statmod")
library("limSolve")
library("R.utils")
library("tictoc")

dir.create(dirname(out.prefix),showWarnings = F,recursive = T)

#in.file <- "/WPSnew/zhenglt/work/proj_h/evalTech.Yiyezhan/phase01/stat/LCPD3-3.count.tab.gz"
#tpm.file <- "/WPSnew/zhenglt/work/proj_h/evalTech.Yiyezhan/phase01/stat/LCPD3-3.TPM.tab.gz"
#out.prefix <- "/WPSnew/zhenglt/work/proj_h/evalTech.Yiyezhan/phase01/norm/LCPD3-3.scran"
#sample.id <- "LCPD3-3"
#sample.desc.file <- "/WPSnew/zhenglt/work/proj_h/evalTech.Yiyezhan/phase01/norm/LCPD3-3.d.txt"
###sample.desc.file <- "/WPS1/zhenglt/work/proj_xy/integrated/sample.design/P0616A.marker.MT.sample.desc.txt"
#arg.method <- "pool"
#arg.filterGene <- T
#arg.filterSample <- T
#pre.outlier.file <- "Marker:/WPSnew/zhenglt/work/proj_h/evalTech.Yiyezhan/phase01/qc/MT.marker/all.marker.exp.ann.outlier.list"
#arg.do.quickCluster <- T
#arg.sizes <- c(20,40,60,80)
#args.ensembl <- F

##pre.outlier.list <- NULL

# plot size factor distribution
plotSizeFactorDist <- function(sf,out.prefix,sample.id.toHighlight=NULL)
{
    pdf(sprintf("%s.%s",out.prefix,"sizeFactor.pdf"),width = 10,height = 6)
    par(cex.lab=1.5,mar=c(5,5,4,2)+0.1)
    cdata <- sort(sf)
    ccol <- rep("darkblue",length(sf))
    names(ccol) <- names(cdata)
    if(!is.null(sample.id.toHighlight)) {
        ccol[sample.id.toHighlight] <- "red"
    }
    barplot(cdata,col=ccol,border=NA,xlab="cell index",ylab="size factor",xaxt="n")
    box(which = "inner",lwd=4)
    par(new=TRUE, oma=c(1,7,1,1),mar=c(5,5,1,2)+0.1)
    layout(matrix(4:1,2))
    b.midpoint <- barplot(cdata[1:10],col=ccol[1:10],border=NA,xlab="",ylab="",xaxt="n")
    text(b.midpoint,y=-0.01, srt = 45, adj = 1, labels = names(cdata)[1:10], xpd = TRUE,cex=1.0)
    box(which = "figure")
    dev.off()
}


getPredefinedOutlier <- function(file)
{
    if(is.null(file)){ return(NULL) }
    ret.list <- list()
    tt <- unlist(strsplit(x = file,split = ";",perl = T))
    for(i in seq_along(tt)){
        .t <- unlist(strsplit(tt[i],split=":"))
        ret.list[[.t[1]]] <- read.table(.t[2],header = F,check.names = F,stringsAsFactors = F)$V1
    }
    return(ret.list)
}
####pre.outlier.list <- getPredefinedOutlier(NULL)
pre.outlier.list <- getPredefinedOutlier(pre.outlier.file)
#if(!("Marker" %in% names(pre.outlier.list))){ pre.outlier.list[["Marker"]] <- c() }
#if(!("MT" %in% names(pre.outlier.list))){ pre.outlier.list[["MT"]] <- c() }
#print("pre.outlier.list:")
#print(str(pre.outlier.list))

## sample description
myDesign <- read.table(sample.desc.file,header = T,check.names = F,stringsAsFactors = F)
rownames(myDesign) <- myDesign$sample
## count
if(grepl(".RData$",in.file,perl = T)){
    lenv <- loadToEnv(in.file)
    exp.in <- lenv$Y
    g.GNAME <- lenv$g.GNAME
}else{
    #exp.in <- read.table(in.file,header = T,check.names = F,stringsAsFactors = F)
    exp.table <- fread(sprintf("gzip -cd %s",in.file))
    exp.in <- as.matrix(exp.table[,c(-1,-2)])
    rownames(exp.in) <- exp.table[[1]]
    g.GNAME <- exp.table[[2]]
    names(g.GNAME) <- rownames(exp.in)
    ###exp.in <- as.matrix(exp.in[,c(-1,-2)])
}
## TPM
if(grepl(".RData$",tpm.file,perl = T)){
    renv <- loadToEnv(tpm.file)
    exp.tpm <- renv$Y
}else{
    ###exp.tpm <- read.table(tpm.file,header = T,check.names = F,stringsAsFactors = F)
    exp.tpm.table <- fread(sprintf("gzip -cd %s",tpm.file))
    exp.tpm <- as.matrix(exp.tpm.table[,c(-1,-2)])
    rownames(exp.tpm) <- exp.tpm.table[[1]]
    ###exp.tpm <- exp.tpm[,c(-1,-2)]
}
##
f.sample <- intersect(rownames(myDesign),colnames(exp.in))
f.sample <- intersect(f.sample,colnames(exp.tpm))
exp.in <- exp.in[,f.sample]
myDesign <- myDesign[f.sample,]
exp.tpm <- exp.tpm[,f.sample]
#head(exp.in[,1:6])
#head(exp.tpm[,1:6])
#head(myDesign)
print(all.equal(colnames(exp.in),rownames(myDesign)))
print(all.equal(colnames(exp.tpm),rownames(myDesign)))

### normalization
doNormalization <- function(dat,tpmData,filterGene=T,filterSample=F,method="pool",
                            outlier.list=NULL,do.quickCluster=F,sizes=c(20, 40, 60, 80, 100),sample.desc.df=NULL)
{
  require("BiocParallel")
  require("scran")
  require("SingleCellExperiment")
  require("DESeq2")
  require("VennDiagram")
  #pDat <- NULL
  #if(!is.null(sample.desc.df)){ pDat= new("AnnotatedDataFrame", data = sample.desc.df) }
  #sce <- newSCESet(countData=dat, phenoData=pDat)
  #tpm(sce) <- as.matrix(tpmData)
  #featureData(sce)[["geneID"]] <- rownames(sce)
  #featureData(sce)[["geneSymbol"]] <- g.GNAME[rownames(sce)]

  sce <- SingleCellExperiment(assays = list(counts = dat, tpm=as.matrix(tpmData)))
  if(!is.null(sample.desc.df)){ colData(sce) <- DataFrame(sample.desc.df) }
  rowData(sce) <- DataFrame("geneID"=rownames(sce),"geneSymbol"=g.GNAME[rownames(sce)])

  ### gene symbol for plot
  ##.gname <- featureData(sce)$geneSymbol
  .gname <- rowData(sce)$geneSymbol
  dup.idx <- which(duplicated(.gname) | is.na(.gname))
  .gname[dup.idx] <- rownames(sce)[dup.idx]
  rowData(sce)[["display.name"]] <- .gname
  if(args.ensembl){
    rowData(sce)[["ensemblID"]]=rownames(sce)
  }else{
    rowData(sce)[["ensemblID"]]=entrezToXXX(rownames(sce),type = "ENSG")
  }
  print(dim(sce))

  ### QC
  sce <- calculateQCMetrics(sce)
  write.table(colData(sce)[,c("total_counts","total_features")],
              sprintf("%s.%s.fltGene%s.fltSample%s.sce.QCMetrics.txt",out.prefix,arg.method,
                      ifelse(arg.filterGene,"T","F"),ifelse(arg.filterSample,"T","F")),row.names = T,sep = "\t",quote = F)

  print(colnames(colData(sce)))
  print(sce)
  hist(sce$total_counts/1e6, xlab="Library sizes (millions)", main=sample.id, breaks=50, col="grey80", ylab="Number of cells")
  hist(sce$total_features, xlab="Number of expressed genes", main=sample.id, breaks=50, col="grey80", ylab="Number of cells")
  fontsize <- theme(axis.text=element_text(size=8), axis.title=element_text(size=14))
  #save(sce,fontsize,file="tmp.RData")
  #print(plotPCA(sce, pca_data_input="pdata") + 
  #      labs(title=sample.id) + 
  #      theme(plot.title = element_text(size = 18, hjust = 0.5)) + fontsize)
  ### filtering samples
  libsize.drop <- isOutlier(sce$total_counts, nmads=args.nmad, type="lower", log=TRUE)
  feature.drop <- isOutlier(sce$total_features, nmads=args.nmad, type="lower", log=TRUE)
  libsize.med <- median(log10(sce$total_counts), na.rm = TRUE)
  libsize.mad <- mad(log10(sce$total_counts), center = libsize.med, na.rm = TRUE)
  print("library size:")
  print(libsize.med)
  print(libsize.mad)
  ##hist(log10(sce$total_counts),breaks = max(length(sce$total_counts)/100,1) )
  hist(log10(sce$total_counts))
  abline(v=libsize.med-5*libsize.mad,lty=2)
  abline(v=libsize.med-4*libsize.mad,lty=2)
  abline(v=libsize.med-3*libsize.mad,lty=2)

  sample.usage.df <- data.frame(sample=colnames(sce),stringsAsFactors = F,libsize.drop=libsize.drop,feature.drop=feature.drop)
  print("outlier.list:")
  print(str(outlier.list))
  if(!is.null(outlier.list)){
      for(i in seq_along(outlier.list)){
          sample.usage.df[[ names(outlier.list)[[i]] ]] <- sample.usage.df$sample %in% outlier.list[[i]]
      }
  }
  print(str(sample.usage.df))
  print(head(sample.usage.df))

  sample.usage.df$size.factor.drop <- F
  idxLastFeatureForFilter <- ncol(sample.usage.df)
  sample.usage.df[["final.drop"]] <- apply(sample.usage.df[,2:idxLastFeatureForFilter],1,function(x){ sum(x)>0 })
  if(filterSample){ sce <- sce[,!sample.usage.df$final.drop] }

  ##### check expression distribution and other QC
  ave.counts <- rowMeans(counts(sce))
  ###
  hist(log10(ave.counts), breaks = 100, main = sample.id, col = "grey80", xlab = expression(Log[10]~"average count"))
  abline(v=log10(1),col="blue", lwd=2, lty=2)
  ###
  sce.plot <- sce
  rownames(sce.plot) <- rowData(sce)[["display.name"]]
  print(plotQC(sce.plot, type = "highest-expression", n=50) +  ylab(sprintf("Feature (%s)",sample.id)) + fontsize)
  numcells <- nexprs(sce, byrow=TRUE)
  ##alt.keep <- numcells >= 10
  ##print(sum(alt.keep))
  smoothScatter(log10(ave.counts), numcells, main=sample.id, 
                xlab=expression(Log[10]~"average count"), 
                ylab="Number of expressing cells")

  ### filtering genes
  if(arg.filterGeneMethod=="AVE"){
      keep <- ave.counts >= 1
      cat(sprintf("Number of genes with average count > 1: %d\n",sum(keep)))
  }else if(arg.filterGeneMethod=="RC"){
      keep <- apply(counts(sce),1,function(x){ nE <- sum(x>arg.numOfReads); return(nE>arg.numOfCells) })
      cat(sprintf("Number of genes expressed(read count > %d) in more than %d cells: %d\n",arg.numOfReads,arg.numOfCells,sum(keep)))
  }else if(arg.filterGeneMethod=="TPM"){
      keep <- apply(tpm(sce),1,function(x){ nE <- sum(x>arg.numOfReads); return(nE>arg.numOfCells) })
      cat(sprintf("Number of genes expressed(TPM > %d) in more than %d cells: %d\n",arg.numOfReads,arg.numOfCells,sum(keep)))
  }
  
  if(filterGene){ sce <- sce[keep,] }
  #### normalization
  if(method=="pool"){
      ##sce <- computeSumFactors(sce,sizes=c(20,30,40),positive=T)
      if(do.quickCluster){
        #"igraph" or "hclust"
        tic("quick clustering ")
        clusters <- quickCluster(sce,method=args.quickClustMethod)
        toc()
        cat("[quick cluster]: \n")
        print(table(clusters))
        tic("compute sizefactor ")
        sce <- computeSumFactors(sce,cluster=clusters,sizes=sizes,positive=TRUE)
        toc()
      }else{
        tic("compute sizefactor ")
        sce <- computeSumFactors(sce,sizes=sizes,positive=TRUE)
        toc()
      }
      sf <- sizeFactors(sce)
      threshold.sf <- args.sizeFactor
      ##threshold.sf <- 0.1
      ##threshold.sf <- 0.05
      MAX.Iter <- 0
      i.Iter <- 0

      plotSizeFactorDist(sf,out.prefix=sprintf("%s.i.%d",out.prefix,i.Iter),sample.id.toHighlight=NULL)
      write.table(data.frame(sample=colnames(sce),sf=sf),sprintf("%s.i.%d.txt",out.prefix,i.Iter),
                  row.names = F,sep = "\t",quote = F)

      while(filterSample && sum(sf<=threshold.sf)>0 && i.Iter<MAX.Iter){
        i.Iter <- i.Iter + 1
        warning(sprintf("number of size factors <=%4.2f: %d",threshold.sf,sum(sf<=threshold.sf)))
        drop.sample.bySF <- colnames(sce)[which(sf<=threshold.sf)]
        cat(sprintf("remove cells with low size factor:"))
        print(drop.sample.bySF)
        sample.usage.df[drop.sample.bySF,"size.factor.drop"] <- T
        sample.usage.df[["final.drop"]] <- apply(sample.usage.df[,2:idxLastFeatureForFilter],1,function(x){ sum(x)>0 })
        sce <- sce[,sf>threshold.sf]
        ### re-calculate size factors
        if(do.quickCluster){
            clusters <- quickCluster(sce)
            ####clusters <- clusters[sf>threshold.sf]
            sce <- computeSumFactors(sce,cluster=clusters,sizes=sizes)
        }else{
            sce <- computeSumFactors(sce,sizes=sizes)
        }
        sf <- sizeFactors(sce)
        plotSizeFactorDist(sf,out.prefix=sprintf("%s.i.%d",out.prefix,i.Iter),sample.id.toHighlight=NULL)
        write.table(data.frame(sample=colnames(sce),sf=sf),sprintf("%s.i.%d",out.prefix,i.Iter), row.names = F,sep = "\t",quote = F)
      }
      cat(sprintf("cells with size factor lower than %4.2f is dropped\n",threshold.sf))
      sce <- sce[,sf>threshold.sf]
      tic("normalize using size factor ")
      sce <- normalize(sce)
      toc()
  }else if(method=="DESeq2"){
    sf <- estimateSizeFactorsForMatrix(counts(sce))
    sizeFactors(sce) <- sf
    #sce <- normalize(sce)
    norm_exprs_mat <- log2(t( t(counts(sce)) / sf )+1)
    norm_exprs(sce) <- norm_exprs_mat
    exprs(sce) <- norm_exprs_mat
  }
  colData(sce)[,"size_factor"] <- sizeFactors(sce)
  ##
  if("MT" %in% colnames(sample.usage.df) && "Marker" %in% colnames(sample.usage.df)){
      venn.plot <- venn.diagram(list("libsize"=subset(sample.usage.df,libsize.drop==TRUE,select="sample",drop=T),
                                     "feature"=subset(sample.usage.df,feature.drop==TRUE,select="sample",drop=T),
                                     "MT"=subset(sample.usage.df,MT==TRUE,select="sample",drop=T),
                                     "Marker"=subset(sample.usage.df,Marker==TRUE,select="sample",drop=T)),filename = NULL,
                                 cat.cex=1.5,cex=1.5)
  }else if("MT" %in% colnames(sample.usage.df)){
      venn.plot <- venn.diagram(list("libsize"=subset(sample.usage.df,libsize.drop==TRUE,select="sample",drop=T),
                                     "feature"=subset(sample.usage.df,feature.drop==TRUE,select="sample",drop=T),
                                     "MT"=subset(sample.usage.df,MT==TRUE,select="sample",drop=T)),filename = NULL,
                                 cat.cex=1.5,cex=1.5)
  }else{
      venn.plot <- venn.diagram(list("libsize"=subset(sample.usage.df,libsize.drop==TRUE,select="sample",drop=T),
                                     "feature"=subset(sample.usage.df,feature.drop==TRUE,select="sample",drop=T)),filename = NULL,
                                 cat.cex=1.5,cex=1.5)
  }
  opar <- par(xpd=T)
  plot.new()
  grid.draw(venn.plot)
  par(opar)
  ##
  print(summary(sizeFactors(sce)))
  plot(sce$total_counts/1e6,sizeFactors(sce),log="xy", main=sample.id, xlab="Library size (millions)", ylab="Size factor")
  return(list("sce"=sce,"sampleUsage"=sample.usage.df))  
}

assignCellCycel <- function(sce,species="human")
{
    human.pairs <- readRDS(system.file("exdata", sprintf("%s_cycle_markers.rds",species), package="scran"))
    ### need ENSGXXX 
    f <- !is.na(rowData(sce)[,"ensemblID"]) & !duplicated(rowData(sce)[,"ensemblID"])
    sce.tmp <- sce[f,]
    #rownames(sce.tmp) <- fData(sce.tmp)[,"ensemblID"]
    rownames(sce.tmp) <- gsub("\\.\\d+$","",rowData(sce.tmp)[,"ensemblID"],perl=T)
    ##assigned <- cyclone(sce, human.pairs,assay="tpm")
    assigned <- cyclone(sce.tmp, human.pairs,assay.type="counts")
    colData(sce)[,"cyclePhase"] <- "S"
    colData(sce)[assigned$scores$G1 > 0.5,"cyclePhase"] <- "G1"
    colData(sce)[assigned$scores$G2M > 0.5,"cyclePhase"] <- "G2M"
    colData(sce)[assigned$scores$G1 > 0.5 & assigned$scores$G2M > 0.5,"cyclePhase"] <- "unknown"
    
    cat("[cyclone]:\n")
    print(table(colData(sce)[,"cyclePhase"]))
    plot(assigned$score$G1, assigned$score$G2M, col="lightblue", pch=16,xlab="G1",ylab="G2M",main=sprintf("CellCycle (%s)",sample.id))
    abline(h=0.5,lty=2,lwd=1.2)
    abline(v=0.5,lty=2,lwd=1.2)
    return(sce)
}

centerByPatient <- function(sce)
{
    Y.new <- c()
    for(pp in unique(sce$patient)){
        f <- colnames(sce)[sce$patient==pp]
        #Y.block <- t(scale(t(exprs(sce)[,f,drop=F]),center = T,scale = F))
        Y.block <- t(scale(t(assay(sce,"norm_exprs")[,f,drop=F]),center = T,scale = F))
        Y.new <- cbind(Y.new,Y.block)
    }
    Y <<- Y.new[,colnames(sce)]
    #assayData(sce)[["centered_norm_exprs"]] <- Y
    assay(sce,"centered_norm_exprs") <- Y
    ###exprs(sce) <- Y
    return(sce)
}

selectOutlier <- function(sce,out.prefix,sf.fit=0.25)
{
    dat.plot <- data.frame(y=sizeFactors(sce),x=sce$total_counts/1e6)
    lm.out <- lm(y~x,subset(dat.plot,y>sf.fit))
    #print(summary(lm.out))
    lm.coeff <- coefficients(lm.out)
    lm.confint <- confint(lm.out)
    pdf(sprintf("%s.dignosis.pdf",out.prefix),width=6,height = 6)
    par(mar=c(5,5,4,2),cex.lab=1.5,cex.main=1.5)
    plot(dat.plot$x,dat.plot$y, log="xy", main=sample.id, xlab="Library size (millions)", ylab="Size factor")
    plot(dat.plot$x,dat.plot$y, log="xy", main=sample.id, xlab="Library size (millions)", ylab="Size factor")
    rect(xleft = c(0.150,0.150,0.800),ybottom = c(0.005,0.450,2.0),xright = c(0.500,0.500,1.200),ytop = c(0.200,0.750,2.2),lwd=1.5,border=c("red","blue","green"),lty=2)
    sce.list <- list()
    sce.list[["red"]] <- sce[,sce$total_counts>150000 & sce$total_counts<500000 & sizeFactors(sce)>0.005 & sizeFactors(sce)<0.200]
    sce.list[["blue"]] <- sce[,sce$total_counts>150000 & sce$total_counts<500000 & sizeFactors(sce)>0.450 & sizeFactors(sce)<0.750]
    sce.list[["green"]] <- sce[,sce$total_counts>800000 & sce$total_counts<1200000 & sizeFactors(sce)>2.0 & sizeFactors(sce)<2.2]
    test.sce <<- sce.list
    #dens.plot <- list("red"=list(),"blue"=list(),"green"=list())
    dens.plot <- list()
    max.y.vec <- c()
    max.x.vec <- c()
    tpm.median <- c()
    #.breaks <<- NULL
    .breaks <<- 0:20
    for(i in seq_along(sce.list)){
        .i.dens <- c()
        for(j in seq_len(ncol(sce.list[[i]]))){
            .f.gene <- tpm(sce.list[[i]])[,j]>0
            .x <- log2(tpm(sce.list[[i]])[.f.gene,j]+1)
            .x[.x > 20] <- 20
            if(!is.null(.breaks)){
                .hist.out<-hist(.x,plot=F,breaks=.breaks)
            }else{
                .hist.out<-hist(.x,plot=F,breaks=11)
                .breaks <<- .hist.out$breaks
            }
            .i.dens <- cbind(.i.dens,.hist.out$density)
            #dens.plot[[i]][[j]] <- density(log2(tpm(sce.list[[i]])[.f.gene,j]+1))
            #max.y.vec <- append(max.y.vec, max(dens.plot[[i]][[j]]$y))
            #max.x.vec <- append(max.x.vec, max(dens.plot[[i]][[j]]$x))
            #tpm.median <- append(tpm.median,median(tpm(sce.list[[i]])[.f.gene,j]))
        }
        dens.plot[[i]] <- .i.dens
    }
    #print(.breaks)
    ee<<-t(sapply(dens.plot,function(x){
               apply(x,1,mad)
                }))
    dd<<-t(sapply(dens.plot,function(x){
               apply(x,1,median)
                }))
    rownames(dd) <<- names(sce.list)
    colnames(dd) <<- sprintf("[%d,%d]",seq_len(ncol(dd))-1,seq_len(ncol(dd)))
    bb<<-barplot(dd,beside = T,col=names(sce.list),xaxt="n",border=NA,main="Histogram of TPM")
    staxlab(1,at = apply(bb,2,mean),labels=colnames(dd),srt=if(ncol(dd) > 5) 45 else 0, cex=0.8,adj=1.2,top.line=2)
    dev.off()
    return(sce.list)
}

plot.CV2andMean <- function(sce,plot=T)
{
    sf <- sizeFactors(sce)
    if(sum(sf<=0)>0){
        warning(sprintf("number of size factors <=0: %d",sum(sf<=0)))
    }
    sce <- sce[,sf>0]
    sf <- sizeFactors(sce)
    nCountsEndo <- t( t(counts(sce)) / sf )
    meansEndo <- rowMeans(nCountsEndo,na.rm = T)
    ##varsEndo <- rowVars( nCountsEndo )
    varsEndo <- apply(nCountsEndo,1,var)
    cv2Endo <- varsEndo / meansEndo^2
    #Do fitting of technical noise
    mincv2 = 0.3
    quan=0.8
    #normalised counts (with size factor)
    minMeanForFitA <- unname( quantile( meansEndo[ which( cv2Endo > mincv2 ) ], quan ) )
    useForFitA <- meansEndo >= minMeanForFitA
    fitA <- glmgam.fit( cbind( a0 = 1, a1tilde = 1/meansEndo[useForFitA] ), cv2Endo[useForFitA] )
    #4. Transform to log-space and propagate error
    eps=1
    LogNcountsEndo=log10(nCountsEndo+eps)
    dLogNcountsEndo=1/((meansEndo+eps)*log(10))
    var_techEndo=(coefficients(fitA)["a0"] + coefficients(fitA)["a1tilde"]/meansEndo)*meansEndo^2
    LogVar_techEndo=(dLogNcountsEndo*sqrt(var_techEndo))^2 #error propagation 
    
    if(plot==TRUE){
        #plot fit 
        par(cex.lab=2.0,cex.main=2.0,cex.axis=1.8,mar=c(5,6,4,2))
        plot( meansEndo, cv2Endo, log="xy", col=1+2*useForFitA, pch=20, cex=0.3, xlab = 'Means', ylab = expression("CV" ^ 2))
        xg <- 10^seq( -3, 5, length.out=100 )
        lines( xg, coefficients(fitA)["a0"] + coefficients(fitA)["a1tilde"]/xg, col='#0000FFF0' )
        legend('bottomleft',c('Genes used for fit', 'Fit baseline variation'),
               pch=c(20, NA),lty =c(NA,1),col=c('green','blue'),cex=1.8)
        title(expression("CV" ^ 2 * " ~ Mean (using endogeneous genes)"))
    }
    res = list()
    res$fit = fitA
    res$techNoiseLog = LogVar_techEndo
}


pdf(sprintf("%s.%s.fltGene%s.fltSample%s.Norm.pdf",out.prefix,arg.method,ifelse(arg.filterGene,"T","F"),ifelse(arg.filterSample,"T","F")),width=6,height = 6)
par(mar=c(5,5,4,2),cex.lab=1.5,cex.main=1.5)
sce.norm.list <- doNormalization(exp.in,exp.tpm,
                                 filterGene=arg.filterGene, filterSample=arg.filterSample,
                                 method=arg.method, outlier.list=pre.outlier.list,
                                 do.quickCluster = arg.do.quickCluster,sizes=arg.sizes,sample.desc.df = myDesign)
dev.off()
sce.norm <- sce.norm.list$sce

## cell cycle assignment
pdf(sprintf("%s.%s.fltGene%s.fltSample%s.cellCycle.pdf",out.prefix,arg.method,
            ifelse(arg.filterGene,"T","F"),
            ifelse(arg.filterSample,"T","F")),width=6,height = 6)
par(mar=c(5,5,4,2),cex.lab=1.5,cex.main=1.5)
tic("assign cellCycle ")
sce.norm <- assignCellCycel(sce.norm,species=args.species)
toc()
dev.off()

myDesign.extend <- left_join(myDesign,sce.norm.list$sampleUsage)
myDesign.extend <- left_join(myDesign.extend,
                             as.data.frame(colData(sce.norm))[,c("sample","total_counts","total_features","size_factor","cyclePhase")])

write.table(myDesign.extend,sprintf("%s.%s.fltGene%s.fltSample%s.design.Extend.txt",out.prefix,arg.method,
                                    ifelse(arg.filterGene,"T","F"),ifelse(arg.filterSample,"T","F")),
            row.names = F,quote = F,sep = "\t")
##write.table(myDesign[subset(myDesign.extend,final.drop==FALSE,select="sample",drop=T),],
write.table(colData(sce.norm)[,c(colnames(myDesign),"total_counts","total_features","size_factor","cyclePhase")],
            sprintf("%s.%s.fltGene%s.fltSample%s.design.txt",out.prefix,arg.method,
                                    ifelse(arg.filterGene,"T","F"),ifelse(arg.filterSample,"T","F")),
            row.names = F,quote = F,sep = "\t")

## fit cv2 ~ mean
pdf(sprintf("%s.%s.fltGene%s.fltSample%s.CV2Mean.pdf",out.prefix,arg.method,ifelse(arg.filterGene,"T","F"),ifelse(arg.filterSample,"T","F")),width=8,height = 8)
plot.CV2andMean(sce.norm.list[["sce"]],plot=T)
dev.off()

#save.image("/WPSnew/zhenglt/work/proj_fh/byBatch/batch01-04/test.RData")
#load("/WPSnew/zhenglt/work/proj_fh/byBatch/batch01-04/test.RData")

## center by patient
##assay(sce.norm,"norm_exprs") <- exprs(sce.norm)
### version difference
assay(sce.norm,"norm_exprs") <- logcounts(sce.norm)
sce.norm <- centerByPatient(sce.norm)
cat(sprintf("save to %s.RData\n",out.prefix))

save(sce.norm,file = sprintf("%s.RData",out.prefix))

