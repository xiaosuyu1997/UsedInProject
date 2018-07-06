#!/usr/bin/env Rscript


suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-i", "--pcaFile", type="character", required=TRUE, help="pcaFile")
parser$add_argument("-o", "--outputPrefix", type="character", required=TRUE, help="output prefix")
parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="verbose mode [default %(default)s]")
args <- parser$parse_args()
print(args)

pca.obj.file <- args$pcaFile
out.prefix <- args$outputPrefix

source("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")
library("R.utils")
library("ggplot2")
library("ggExtra")
library("mvoutlier")
library("extremevalues")
#pca.obj.file <- "OUT.PCA.TSNE.usingMKI67hiGList/CD8/colonC.sampleType.PCA.tmp.RData"
#out.prefix <- "OUT.PCA.TSNE.usingMKI67hiGList/CD8/post/colonC.sampleType.PCA"

dir.create(dirname(out.prefix),showWarnings = F,recursive = T)

sampleTypeColorFile <- "/WPS1/zhenglt/work/TCR_chunhong/data/CellType.color"
sampleTypeColor <- read.SampleTypeColor(sampleTypeColorFile)

lenv <- loadToEnv(pca.obj.file)
names(lenv)
lenv$d.ind.coord %>% .[1:4,1:8]

d.ind.coord.df <- lenv$d.ind.coord
d.ind.coord <- as.matrix(d.ind.coord.df[,c(-1,-2)])
d.ind.coord.df[1:4,1:8]
d.ind.coord[1:4,1:8]
sampleType.uniq <- as.character(d.ind.coord.df$sampleType %>% unique)
sampleTypeColor <- sampleTypeColor[names(sampleTypeColor) %in% sampleType.uniq]

plot.tsne.points(d.ind.coord,out.file=sprintf("%s.coord.PC1-PC2.pdf",out.prefix),
                 tsne.col.points=sampleTypeColor[d.ind.coord.df$sampleType],
                 col.tsne.legend=sampleTypeColor,tsne.legend=names(sampleTypeColor),pch=20, peak=NULL,main="PCA",x.dim=c(1,2))

pdf(sprintf("%s.coord.PC1-PC2.mvoutlier.pdf",out.prefix),width=7,height = 7)
uni.res <- uni.plot(d.ind.coord[,1:2], symb=FALSE, quan=1, alpha=0.00005)
dev.off()

##### outlier detection
### other optoins from package outliers:
### dixon.test()
### chisq.out.test()
K <- getOutliers(d.ind.coord.df[,"Dim.1"],method="I",distribution="normal")
L <- getOutliers(d.ind.coord.df[,"Dim.1"],method="II",distribution="normal")
pdf(sprintf("%s.outlier.pdf",out.prefix),width = 10,height = 6)
opar <- par(mfrow=c(1,2))
outlierPlot(d.ind.coord.df[,"Dim.1"],K,mode="qq")
outlierPlot(d.ind.coord.df[,"Dim.1"],L,mode="residual")
dev.off()
par(opar)
print("K$mu,K$sigma,K$limit")
print(c(K$mu,K$sigma,K$limit))
d.ind.coord.df$outlier <- seq_len(nrow(d.ind.coord.df)) %in% K$iRight
######

### my hijack
#rob <- robustbase::covMcd(d.ind.coord[,1:2], alpha = 1)
#xarw <- arw(d.ind.coord[,1:2], rob$center, rob$cov, alpha = 0.00005)
#f.outlier <- uni.res$md> min(sqrt(xarw$cn))
#print("num of outlier by my hijack:")
#sum(f.outlier)
#print("num of outlier by uni.plot:")
#uni.res$outliers %>% sum
#d.ind.coord.df$outlier <- (f.outlier[rownames(d.ind.coord.df)] & (d.ind.coord.df[,"Dim.1"]>5 | d.ind.coord.df[,"Dim.2"]>0) )
p <- ggplot(d.ind.coord.df, aes_string('Dim.1','Dim.2')) + 
        geom_point(aes(colour=outlier),show.legend=F) + scale_colour_manual(values = c("#D95F02","#1B9E77")) +
        #scale_colour_brewer(palette = "Dark2") + 
        theme_bw(15)
pdf(sprintf("%s.coord.PC1-PC2.density.pdf",out.prefix),width=7,height = 7)
pp <- ggExtra::ggMarginal(p, type = 'density', margins = 'both', size = 6, col = '#FF0000', fill = '#FFA500')
print(pp)
dev.off()

out.df <- d.ind.coord.df[,c("sampleType", "sampleName", "Dim.1", "Dim.2","outlier")]
colnames(out.df) <- gsub("sampleName","sample",colnames(out.df))
write.table(out.df,sprintf("%s.add.outlier.txt",out.prefix),row.names = F,quote = F,sep = "\t")


