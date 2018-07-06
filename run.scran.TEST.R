#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-i", "--infile", type="character", required=TRUE, help="infile")
parser$add_argument("-t", "--tpmfile", type="character", required=TRUE, help="tpmfile")
parser$add_argument("-o", "--outprefix", type="character", required=TRUE, help="outprefix")
parser$add_argument("-s", "--sample", type="character", default="SAMPLE", help="sample id  [default %(default)s]")
parser$add_argument("-d", "--designFile", required=TRUE, type="character", help="sample desc file")
args <- parser$parse_args()
print(args)

in.file <- args$infile
tpm.file <- args$tpmfile
out.prefix <- args$outprefix
sample.id <- args$sample
sample.desc.file <- args$designFile

library("dplyr")
library("scater")
library("BiocParallel")
library("scran")
library("plotrix")
library("statmod")
library("limSolve")

#in.file <- "/WPS1/zhenglt/work/proj_xy/integrated/quantification/P0616A.count.tab.gz"
#tpm.file <- "/WPS1/zhenglt/work/proj_xy/integrated/quantification/P0616A.TPM.tab.gz"
#out.prefix <- "./P0616A.scran"
#sample.id <- "P0616A"
#sample.desc.file <- "/WPS1/zhenglt/work/proj_xy/integrated/sample.design/P0616A.above0K.sample.desc.txt"

## sample description
myDesign <- read.table(sample.desc.file,header = T,check.names = F,stringsAsFactors = F)
rownames(myDesign) <- myDesign$sample
## count
exp.in <- read.table(in.file,header = T,check.names = F,stringsAsFactors = F)
rownames(exp.in) <- exp.in[,1]
g.GNAME <- exp.in[,2]
names(g.GNAME) <- rownames(exp.in)
exp.in <- exp.in[,c(-1,-2)]
## TPM
exp.tpm <- read.table(tpm.file,header = T,check.names = F,stringsAsFactors = F)
rownames(exp.tpm) <- exp.tpm[,1]
exp.tpm <- exp.tpm[,c(-1,-2)]
##
f.sample <- intersect(rownames(myDesign),colnames(exp.in))
f.sample <- intersect(f.sample,colnames(exp.tpm))
exp.in <- exp.in[,f.sample]
myDesign <- myDesign[f.sample,]
exp.tpm <- exp.tpm[,f.sample]
head(exp.in[,1:6])
head(exp.tpm[,1:6])
head(myDesign)
all.equal(colnames(exp.in),rownames(myDesign))
all.equal(colnames(exp.tpm),rownames(myDesign))

doNormalizationDESeq2 <- function(dat){
    ##
    .sf <- estimateSizeFactorsForMatrix(dat)
    ##
    total_counts <- apply(dat,2,sum)
    plot(.sf, total_counts/1e6, log="xy", main=sample.id, ylab="Library size (millions)", xlab="Size factor")
    ###
    dat.norm <- t( t(dat) / .sf )
    return(log2(dat.norm+1))
}
### normalization
doNormalization <- function(dat,tpmData,filterGene=T,filterSample=F,method="pool")
{
  require("BiocParallel")
  require("scran")
  require("DESeq2")
  sce <- newSCESet(countData=dat)
  tpm(sce) <- as.matrix(tpmData)
  featureData(sce)[["geneSymbol"]] <- g.GNAME
  ### gene symbol for plot
  .gname <- featureData(sce)$geneSymbol
  dup.idx <- which(duplicated(.gname) | is.na(.gname))
  .gname[dup.idx] <- rownames(sce)[dup.idx]
  featureData(sce)[["display.name"]] <- .gname
  print(dim(sce))
  sce <- calculateQCMetrics(sce)
  print(colnames(pData(sce)))
  hist(sce$total_counts/1e6, xlab="Library sizes (millions)", main=sample.id, breaks=20, col="grey80", ylab="Number of cells")
  hist(sce$total_features, xlab="Number of expressed genes", main=sample.id, breaks=20, col="grey80", ylab="Number of cells")
  libsize.drop <- isOutlier(sce$total_counts, nmads=3, type="lower", log=TRUE)
  feature.drop <- isOutlier(sce$total_features, nmads=3, type="lower", log=TRUE)
  print(sum(libsize.drop | feature.drop))
  if(filterSample){ sce <- sce[,!(libsize.drop | feature.drop)] }
  print(data.frame(ByLibSize=sum(libsize.drop), ByFeature=sum(feature.drop), Remaining=ncol(sce)))
  ###           ByLibSize ByFeature Remaining
  ## Samples         2         2       268
  fontsize <- theme(axis.text=element_text(size=8), axis.title=element_text(size=14))
  print(plotPCA(sce, pca_data_input="pdata") + 
        labs(title=sample.id) + 
        theme(plot.title = element_text(size = 18, hjust = 0.5)) + 
        fontsize)
  ave.counts <- rowMeans(counts(sce))
  keep <- ave.counts >= 1
  print(sum(keep))
  ###
  hist(log10(ave.counts), breaks = 100, main = sample.id, col = "grey80", xlab = expression(Log[10]~"average count"))
  abline(v=log10(1),col="blue", lwd=2, lty=2)
  ###
  sce.plot <- sce
  featureNames(sce.plot) <- featureData(sce)[["display.name"]]
  print(plotQC(sce.plot, type = "highest-expression", n=50) +  ylab(sprintf("Feature (%s)",sample.id)) + fontsize)
  numcells <- nexprs(sce, byrow=TRUE)
  alt.keep <- numcells >= 10
  print(sum(alt.keep))
  smoothScatter(log10(ave.counts), numcells, main=sample.id, 
                xlab=expression(Log[10]~"average count"), 
                ylab="Number of expressing cells")
  if(filterGene){ sce <- sce[keep,] }
  if(method=="pool"){
      sce <- computeSumFactors(sce,sizes=c(20,30,40),positive=T)
      #sce <- computeSumFactors(sce)
      sce <- normalize(sce)
  }else if(method=="DESeq2"){
    sf <- estimateSizeFactorsForMatrix(counts(sce))
    sizeFactors(sce) <- sf
    #sce <- normalize(sce)
    norm_exprs_mat <- log2(t( t(counts(sce)) / sf )+1)
    norm_exprs(sce) <- norm_exprs_mat
    exprs(sce) <- norm_exprs_mat
  }
  summary(sizeFactors(sce))
  plot(sizeFactors(sce), sce$total_counts/1e6, log="xy", main=sample.id, ylab="Library size (millions)", xlab="Size factor")
  return(sce)  
}

arg.filterGene <- T
arg.filterSample <- T
pdf(sprintf("%s.pool.fltGene%s.fltSample%s.pdf",out.prefix,ifelse(arg.filterGene,"T","F"),ifelse(arg.filterSample,"T","F")),width=6,height = 6)
par(mar=c(5,5,4,2),cex.lab=1.5,cex.main=1.5)
sce.pool <- doNormalization(exp.in,exp.tpm,filterGene=arg.filterGene,filterSample=arg.filterSample,method="pool")
dev.off()

pdf(sprintf("%s.DESeq2.fltGene%s.fltSample%s.pdf",out.prefix,ifelse(arg.filterGene,"T","F"),ifelse(arg.filterSample,"T","F")),width=6,height = 6)
par(mar=c(5,5,4,2),cex.lab=1.5,cex.main=1.5)
sce.DESeq2 <- doNormalization(exp.in,exp.tpm,filterGene=arg.filterGene,filterSample=arg.filterSample,method="DESeq2")
dev.off()

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

plot.CV2andMean <- function(sce,plot=T){
    sf <- sizeFactors(sce)
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

## compare DESeq2 and pool
pdf(sprintf("%s.sfCmp.fltSample%s.pdf",out.prefix,ifelse(arg.filterSample,"T","F")),width=6,height = 6)
par(mar=c(5,5,4,2),cex.lab=1.5,cex.main=1.5)
plot(sizeFactors(sce.pool),sizeFactors(sce.DESeq2), log="xy", main=sample.id, xlab="Pool method", ylab="DESeq2")
dev.off()
## fit cv2 ~ mean
pdf(sprintf("%s.FitByDESeq2.fltSample%s.pdf",out.prefix,ifelse(arg.filterSample,"T","F")),width=8,height = 8)
plot.CV2andMean(sce.DESeq2,plot=T)
dev.off()
pdf(sprintf("%s.FitByPool.fltSample%s.pdf",out.prefix,ifelse(arg.filterSample,"T","F")),width=8,height = 8)
plot.CV2andMean(sce.pool,plot=T)
dev.off()

save(sce.DESeq2,sce.pool,file = sprintf("%s.RData",out.prefix))
## check TPM distribution of outliers
#outlier.DESeq2 <- selectOutlier(sce.DESeq2,sprintf("%s.DESeq2.fltGene%s",out.prefix,ifelse(arg.filterGene,"T","F")),sf.fit=0.25)

#selectOutlier(sce.pool,sprintf("%s.pool.fltGene%s",out.prefix,ifelse(arg.filterGene,"T","F")),sf.fit=0.25)

#apply(counts(sce.DESeq2)[,which(sizeFactors(sce.DESeq2)<0.01 & sce.DESeq2$total_counts>150000 )],2,sum)
#apply(counts(sce.DESeq2)[,which(sizeFactors(sce.DESeq2)<0.01 & sce.DESeq2$total_counts>150000 )],2,function(x){ median(x[x>0]) })
#slist.good <- which(sizeFactors(sce.DESeq2)>1 & 
#                                sizeFactors(sce.DESeq2)<1.2 & 
#                                sce.DESeq2$total_counts>250000 & 
#                                sce.DESeq2$total_counts<450000 )
#apply(counts(sce.DESeq2)[,slist.good],2,sum)
#apply(counts(sce.DESeq2)[,slist.good],2,function(x){ median(x[x>0]) })

doHeatmap <- function(sce)
{
  require(gplots)
  require(graphics)
  require(dynamicTreeCut)
  var.fit <- trendVar(sce, trend="loess", use.spikes=FALSE, span=0.2)
  var.out <- decomposeVar(sce, var.fit)
  plot(var.out$mean, var.out$total, pch=16, cex=0.6, xlab="Mean log-expression", ylab="Variance of log-expression")
  o <- order(var.out$mean)
  lines(var.out$mean[o], var.out$tech[o], col="dodgerblue", lwd=2)
  
  var.out <- var.out[order(var.out$bio),]
  hvg.out <- var.out[which(var.out$FDR <= 0.05 & var.out$bio >= 0.5),]
  hvg.out <- hvg.out[order(hvg.out$bio, decreasing=TRUE),]
  nrow(hvg.out)
  ##[1] 198
  print(head(hvg.out))
  plotExpression(sce, rownames(hvg.out)[1:min(10,nrow(hvg.out))]) + fontsize
  
  var.cor <- correlatePairs(sce, subset.row=rownames(hvg.out))
  head(var.cor)
  sig.cor <- (var.cor$p.value < 1e-4 & var.cor$rho > 0.4)
  summary(sig.cor)
  head(var.cor[sig.cor,])
  ###
  chosen <- unique(c(var.cor$gene1[sig.cor], var.cor$gene2[sig.cor]))
  #chosen <- head(rownames(hvg.out),n = 50)
  norm.exprs <- exprs(sce)[chosen,,drop=FALSE]
  rownames(norm.exprs) <- fData(sce)[chosen,"display.name"]
  heat.vals <- norm.exprs - rowMeans(norm.exprs)
  heat.out <- heatmap.2(heat.vals, col=bluered, symbreak=TRUE, trace='none', cexRow=0.6,margins = c(8,6) )
  print(plotPCA(sce, exprs_values="exprs", colour_by="total_features", feature_set=chosen) + fontsize)
  set.seed(100)
  out5 <- plotTSNE(sce, exprs_values="exprs", perplexity=5, colour_by="total_features", feature_set=chosen) + fontsize + ggtitle("Perplexity = 5")
  out10 <- plotTSNE(sce, exprs_values="exprs", perplexity=10, colour_by="total_features", feature_set=chosen) + fontsize + ggtitle("Perplexity = 10")
  print(out5)
  print(out10)
  #multiplot(out5, out10, cols=2)
  #out20 <- plotTSNE(sce, exprs_values="exprs", perplexity=20, colour_by="total_features", feature_set=chosen) + fontsize + ggtitle("Perplexity = 20")
  #multiplot(out5, out10, out20, cols=3)
  
  my.dist <- dist(t(norm.exprs))
  my.tree <- hclust(my.dist, method = "ward.D2")
  my.clusters <- unname(cutreeDynamic(my.tree, distM=as.matrix(my.dist), verbose=0))+1
  clust.col <- rainbow(length(unique(my.clusters)))
  heatmap.2(heat.vals, col = bluered, symbreaks = TRUE, trace = 'none', cexRow = 0.3, 
            ColSideColors = clust.col[my.clusters], Colv = as.dendrogram(my.tree), 
            breaks = seq(-5,5,length.out = 21),margins = c(8,6))
  return(hvg.out)
}

#pdf(sprintf("%s.heatmap.pool.pdf",out.prefix),width=8,height = 8)
#hvg.pool <- doHeatmap(sce.pool)
#dev.off()

#pdf(sprintf("%s.heatmap.DESeq2.pdf",out.prefix),width=8,height = 8)
#hvg.DESeq2 <- doHeatmap(sce.DESeq2)
#dev.off()


