#!/usr/bin/env Rscript


suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-i", "--expFile", type="character", required=TRUE, help="expFile")
parser$add_argument("-g", "--geneFile", type="character", required=TRUE, help="geneFile")
parser$add_argument("-d", "--designFile", type="character", required=TRUE, help="designFile")
parser$add_argument("-o", "--outputPrefix", type="character", required=TRUE, help="output prefix")
parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="verbose mode [default %(default)s]")
args <- parser$parse_args()
print(args)

exp.file <- args$expFile
g.file <- args$geneFile
d.file <- args$designFile
out.prefix <- args$outputPrefix

if(file.exists("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R"))
{
    source("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")
}else{
    source("/WPSnew/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")
}

library("R.utils")
library("ggplot2")
library("ggExtra")
library("ggpubr")
library("mvoutlier")
library("reshape2")
library("extremevalues")
library("SingleCellExperiment")

#exp.file <- "/WPSnew/zhenglt/work/proj_zl/integrate0924/cross.patient.b20170926/exp/colonC.P10.scran.RData"
#g.file <- "/WPSnew/zhenglt/work/proj_zl/integrate0924/cross.patient.b20170926/functional/MKI67/prof.sig.zyy.list"
#d.file <- "/WPSnew/zhenglt/work/proj_zl/integrate0924/cross.patient.b20170926/clusteringResult/colonC.clustering.final.all.txt"
#out.prefix <- "/WPSnew/zhenglt/work/proj_zl/integrate0924/cross.patient.b20170926/functional/MKI67/test.byAve"

##sampleTypeColorFile <- "/WPS1/zhenglt/work/TCR_chunhong/data/CellType.color"

dir.create(dirname(out.prefix),showWarnings = F,recursive = T)

lenv <- loadToEnv(exp.file)
names(lenv)
if("Y" %in% names(lenv)){
    Y <- lenv$Y
}else{
    #Y <- exprs(lenv$sce.norm)
    Y <- assay(lenv$sce.norm,"centered_norm_exprs")
}
#Y[1:4,1:8]

g.list <- read.table(g.file,header = T,sep = "\t",stringsAsFactors = F)
g.list$geneID <- as.character(g.list$geneID)
rownames(g.list) <- g.list$geneID
g.list <- g.list[intersect(g.list$geneID,rownames(Y)),]

print("g.list")
print(g.list)

myDesign <- read.table(d.file,header = T,sep = "\t",stringsAsFactors = F)
rownames(myDesign) <- myDesign$sample
f.s <- intersect(rownames(myDesign),colnames(Y))
myDesign <- myDesign[f.s,]
Y <- Y[,f.s]

##sampleTypeColor <- read.SampleTypeColor(sampleTypeColorFile)
myDesign$prol.score <- apply(Y[rownames(g.list),rownames(myDesign)],2,mean)

##### outlier detection (use method I at last)
K <- getOutliers(myDesign$prol.score,method="I",distribution="normal")
L <- getOutliers(myDesign$prol.score,method="II",distribution="normal")
pdf(sprintf("%s.outlier.pdf",out.prefix),width = 10,height = 6)
opar <- par(mfrow=c(1,2))
outlierPlot(myDesign$prol.score,K,mode="qq")
outlierPlot(myDesign$prol.score,L,mode="residual")
dev.off()
par(opar)
print("K$mu,K$sigma,K$limit")
print(c(K$mu,K$sigma,K$limit))
myDesign$prol.cls <- "lo"
myDesign$prol.cls[K$iRight] <- "hi"
######

##### plot #####
#prol.score.x <- seq(-2,6,0.001)
##prol.score.x <- seq(K$yMin,K$yMax,0.001)
prol.score.pretty <- pretty(myDesign$prol.score)
prol.score.x <- seq(prol.score.pretty[1], prol.score.pretty[length(prol.score.pretty)],0.001)
prol.score.forPlot <- data.frame(x=prol.score.x,y=dnorm(prol.score.x,mean = K$mu,sd = K$sigma))
p <- ggplot(myDesign, aes(prol.score)) + geom_density(colour="black") + theme_bw() +
    geom_line(data=prol.score.forPlot,aes(x=x,y=y),colour="red") +
    geom_vline(xintercept = K$limit,linetype=2)
ggsave(filename = sprintf("%s.density.pdf",out.prefix),width = 4,height = 3)

myDesign %>% str
myDesign %>% head
prol.cls.count <- unclass(table(myDesign[,c("majorCluster","prol.cls")]))
prol.cls.freq <- t(apply(prol.cls.count,1,function(x){ x/sum(x) }))
print(prol.cls.freq)

dat.plot <- melt(prol.cls.freq)
colnames(dat.plot) <- c("majorCluster","prol.cls","freq")
dat.plot$prol.cls <- factor(dat.plot$prol.cls,levels=c("lo","hi"))

p <- ggbarplot(data = dat.plot,x="majorCluster",y="freq",
               #fill="MKI67",palette = c("#00AFBB", "#E7B800")) +
               fill="prol.cls",palette = "Paired") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file=sprintf("%s.freqInCluster.pdf",out.prefix),width = 6,height = 4)

write.table(myDesign,file = sprintf("%s.prolCls.txt",out.prefix),row.names = F,quote = F,sep = "\t")

