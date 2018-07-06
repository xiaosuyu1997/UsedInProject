#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-d", "--designFile", type="character", required=TRUE, help="design matrix file")
parser$add_argument("-i", "--inputFile", type="character", required=TRUE, help="input file")
parser$add_argument("-j", "--yExpFile", type="character", help="exp file for y")
parser$add_argument("-f", "--geneFile", type="character", help="gene list file")
parser$add_argument("-c", "--cellTypeColorFile", default="/WPS1/zhenglt/work/TCR_chunhong/data/CellType.color", 
                    type="character", help="cellTypeColorFile [default %(default)s]")
parser$add_argument("-o", "--outDir", type="character", required=TRUE, help="output directory")
parser$add_argument("-s", "--sample", type="character", default="SAMPLE", help="sample id  [default %(default)s]")
parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="whether verbose mode [default %(default)s]")
parser$add_argument("-l", "--log", action="store_true", default=FALSE, help="whether do log transform [default %(default)s]")
parser$add_argument("-r", "--center", action="store_true", default=FALSE, help="whether center genes by individual [default %(default)s]")
parser$add_argument("-m", "--notFilter", action="store_true", default=FALSE, help="don't filtering genes [default %(default)s]")
parser$add_argument("-n", "--nKeep", type="integer", help="use top nKeep genes [default %(default)s]")
#### method specific
parser$add_argument("-q", "--noRF", action="store_true", default=FALSE, help="don't run randomForest [default %(default)s]")
parser$add_argument("-x", "--XLabelFile", type="character",
                    help="x label file, contain \"sample\",\"majorCluster\" columns [default %(default)s]")
parser$add_argument("-y", "--YFile", type="character",
                    help="y file, contain \"sample\" column [default %(default)s]")
parser$add_argument("-k", "--ncores", type="integer", default=8, help="num of cors to use [default %(default)s]")
parser$add_argument("-t", "--myseed", type="character", help="seed [default %(default)s]")
#parser$add_argument("-a", "--fracClustered", type="double",default=0.85, help="minimum fraction of clustered samples [default %(default)s]")
parser$add_argument("-e", "--geneIDType", type="character",default="entrezID", help="geneID type (entrezID, ensemblID) [default %(default)s]")
parser$add_argument("-b", "--normExprs", action="store_true", default=FALSE, help="use norm_exprs assay [default %(default)s]")
args <- parser$parse_args()

print(args)
designFile <- args$designFile
inputFile <- args$inputFile
geneFile <- args$geneFile
cellTypeColorFile <- args$cellTypeColorFile
###clonotype.file <- args$clonotypeFile
out.dir <- args$outDir
sample.id <- args$sample
args.log <- args$log
args.center <- args$center
nKeep <- args$nKeep
args.notFilter <- args$notFilter
args.geneIDType <- args$geneIDType

args.x.label.file <- args$XLabelFile
args.y.file <- args$YFile
args.ncores <- args$ncores
args.myseed <- args$myseed
args.norm.exprs <- args$normExprs
args.y.exp.file <- args$yExpFile
#args.kmax <- args$kmax
#args.iterMax <- args$iterMax

##### frequently used code
dir.create(out.dir,recursive = T,showWarnings = F)

source("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")
#suppressPackageStartupMessages(require("ComplexHeatmap"))
#suppressPackageStartupMessages(require("circlize"))
#suppressPackageStartupMessages(require("gridBase"))
#suppressPackageStartupMessages(require("dendextend"))
suppressPackageStartupMessages(require("RColorBrewer"))
suppressPackageStartupMessages(require("gplots"))
#suppressPackageStartupMessages(library("vioplot"))
suppressPackageStartupMessages(library("plotrix"))
#suppressPackageStartupMessages(library("vioplot"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("Rmisc"))
suppressPackageStartupMessages(library("varSelRF"))
suppressPackageStartupMessages(library("R.utils"))

#out.dir <- "/WPSnew/zhenglt/work/proj_zl/integrate0924/cross.patient.b20170523/randomForest/test"
#sample.id <- "colonC.Th"
#designFile <- "/WPSnew/zhenglt/work/proj_zl/integrate0924/cross.patient.b20170523/OUT.Th/colonC.Th.iter.leaf.txt"
#cellTypeColorFile <- "/WPS1/zhenglt/work/TCR_chunhong/data/CellType.color"
#inputFile <- "/WPSnew/zhenglt/work/proj_zl/integrate0924/cross.patient.b20170523/exp/colonC.P8.scran.fltMethodTPM.3rc10cell.RData"
#args.notFilter <- F
#geneFile <- NULL
#args.center <- F
#args.log <- F
#args.norm.exprs <- F
#args.geneIDType <- "entrezID"
#args.x.label.file <- NULL
#args.y.file <- NULL
##args.x.label.file <- "/WPSnew/zhenglt/work/proj_zl/integrate0924/cross.patient.b20170523/randomForest/test/colonC.Th.x.labels.list"
##args.y.file <- "/WPSnew/zhenglt/work/proj_zl/integrate0924/cross.patient.b20170523/randomForest/test/colonC.Th.y.list"
#args.ncores <- 8
#nKeep <- NULL
#args.myseed <- 123456

###suppressPackageStartupMessages(library("stringr"))
out.prefix <- sprintf("%s/%s.pred",out.dir,sample.id)
g.inputD <- processInput(designFile,cellTypeColorFile,inputFile,args.notFilter,geneFile = NULL,
                         args.center,args.log,args.norm.exprs = args.norm.exprs)
myDesign <- g.inputD$myDesign
sampleTypeColor <- g.inputD$sampleTypeColor
patient.col.list <- g.inputD$patient.col.list
Y <- g.inputD$Y
g.GNAME <- g.inputD$g.GNAME
i.mcls <- which(colnames(myDesign)=="majorCluster")
if(!is.null(args.y.exp.file)){
    lenv <- loadToEnv(args.y.exp.file)
    Y.y <- lenv$Y
}else{
    Y.y <- NULL
}
##### 

#### rename rownames
Y.newRowname <- Y
newName <- g.GNAME[rownames(Y)]
names(newName) <- rownames(Y)
f.gname.na <- is.na(newName)
newName[f.gname.na] <- names(newName)[f.gname.na]
rownames(Y.newRowname) <- unname(newName)
if(args.geneIDType=="ensemblID"){ Y <- Y.newRowname }

####
if(!is.null(args.x.label.file) && file.exists(args.x.label.file)){
    x.label.df <- read.table(args.x.label.file,header = T,sep = "\t",check.names = F,stringsAsFactors = F)
    rownames(x.label.df) <- x.label.df$sample
    f.s <- intersect(colnames(Y),rownames(x.label.df))
    x.label.df <- x.label.df[f.s,,drop=F]

}else{
    x.label.df <- myDesign[,c("sample","majorCluster")]
}
x.label.df <- subset(x.label.df,majorCluster!="uncharacterized")

if(!is.null(args.y.file) && file.exists(args.y.file)){
    y.df <- read.table(args.y.file,header = T,sep = "\t",check.names = F,stringsAsFactors = F)
    rownames(y.df) <- y.df$sample
    if(!is.null(Y.y)){
        f.s <- intersect(colnames(Y.y),rownames(y.df))
    }else{
        f.s <- intersect(colnames(Y),rownames(y.df))
    }
    y.df <- y.df[f.s,,drop=F]
    print("y.df:")
    print(str(y.df))
}

### get differential expressed genes

val.res.sil.df.flt <- NULL
tsne.Y <- NULL
aov.res <- NULL
if(is.null(geneFile) || !file.exists(geneFile)){
    #aov.res <- runMultiGroupSpecificGeneTest(Y[,x.label.df$sample,drop=F],grps=x.label.df$majorCluster,
    #                                     out.prefix,mod=NULL,FDR.THRESHOLD=0.05,FC.THRESHOLD=1,verbose=F,n.cores=args.ncores)

    clusterColor <- structure(colorRampPalette(brewer.pal(12,"Paired"))(length(unique(x.label.df$majorCluster))),
                              names=sort(as.character(unique(x.label.df$majorCluster))))
    clusterColor["uncharacterized"] <- "gray"

    clusterMarker.res <- my.clusterMarkerGene(Y[,rownames(x.label.df)],
                                            clust=x.label.df$majorCluster,
                                            ann.col=sampleTypeColor[myDesign[rownames(x.label.df),"sampleType"]],
                                            clust.col=clusterColor,
                                            out.prefix=sprintf("%s.predMCls", out.prefix),
                                            n.cores=args.ncores,
                                            sampleType = myDesign[rownames(x.label.df),"sampleType"],
                                            sampleTypeColSet = sampleTypeColor)

    aov.res <- clusterMarker.res$aov.res
    if(nrow(clusterMarker.res$gene.table)>10){
        #tsne.res <- runTSNEAnalysis(Y[rownames(clusterMarker.res$gene.table),rownames(myDesign),drop=F],
        tsne.res <- runTSNEAnalysis(Y[rownames(clusterMarker.res$aov.res$aov.out.sig) %>% head(n=500),rownames(myDesign),drop=F],
                                    sprintf("%s.predMCls.tsne",out.prefix),
                                    col.points = sampleTypeColor[myDesign[,"sampleType"]],
                                    legend=names(sampleTypeColor),
                                    col.legend=sampleTypeColor,
                                    myseed=args.myseed,do.dbscan = F,do.scale = F,distance.metric = NULL)

        tsne.Y <- tsne.res$`30`$Rtsne.res$Y
        plot.tsne.points(tsne.Y,
                         sprintf("%s.predMCls.tsne.pdf",out.prefix),
                         tsne.col.points=clusterColor[myDesign$majorCluster],
                         col.tsne.legend=clusterColor,
                         tsne.legend=names(clusterColor),
                         pch=20,nclusters=length(clusterColor))

        ### re-evaluation
        val.res <- my.clusteringValidation(t(tsne.Y)[,rownames(myDesign),drop=F],
                                       clust=myDesign$majorCluster,
                                       out.prefix=sprintf("%s.predMCls.val.00",out.prefix),
                                       n.cores=args.ncores,dist.obj=NULL)

        val.res.sil.df <- data.frame(sample=rownames(myDesign),stringsAsFactors = F,
                                cluster.old=val.res$sil[,1],
                                neighbor=val.res$sil[,2],
                                sil_width=val.res$sil[,3],
                                cluster=val.res$sil[,1],
                                majorCluster.old=myDesign$majorCluster,
                                majorCluster=myDesign$majorCluster)

        val.res.sil.df$cluster[val.res.sil.df$sil_width<0] <- 0
        val.res.sil.df$majorCluster[val.res.sil.df$sil_width<0] <- "uncharacterized"
        rownames(val.res.sil.df) <- val.res.sil.df$sample
        val.res.sil.df.flt <- subset(val.res.sil.df,sil_width>0)
        write.table(right_join(myDesign[,-i.mcls],val.res.sil.df),
                    sprintf("%s.predMCls.sil.txt",out.dir,sample.id),row.names = F,quote = F,sep="\t")

        plot.tsne.points(tsne.Y[rownames(val.res.sil.df),,drop=F],
                         sprintf("%s.predMCls.tsne.gray00.pdf",out.prefix),
                         tsne.col.points=clusterColor[val.res.sil.df$majorCluster],
                         col.tsne.legend=clusterColor,
                         tsne.legend=names(clusterColor),
                         pch=20,nclusters=length(clusterColor))

        plot.tsne.points(tsne.Y[rownames(val.res.sil.df.flt),,drop=F],
                         sprintf("%s.predMCls.tsne.rm00.pdf",out.prefix),
                         tsne.col.points=clusterColor[val.res.sil.df.flt$majorCluster],
                         col.tsne.legend=clusterColor,
                         tsne.legend=names(clusterColor),
                         pch=20,nclusters=length(clusterColor))

        clusterMarker.res.rm00 <- my.clusterMarkerGene(Y[,rownames(val.res.sil.df.flt)],
                                            clust=val.res.sil.df.flt$majorCluster,
                                            ann.col=sampleTypeColor[myDesign[rownames(val.res.sil.df.flt),"sampleType"]],
                                            clust.col=clusterColor,
                                            out.prefix=sprintf("%s.predMCls.rm00", out.prefix),
                                            n.cores=args.ncores,
                                            sampleType = myDesign[rownames(val.res.sil.df.flt),"sampleType"],
                                            sampleTypeColSet = sampleTypeColor)
        aov.res <- clusterMarker.res.rm00$aov.res
    }
    if(is.null(args.x.label.file) || !file.exists(args.x.label.file)){
        x.label.df <- val.res.sil.df.flt[,c("sample","majorCluster")]
    }
    if(is.null(args.y.file) || !file.exists(args.y.file)){
        y.df <- subset(val.res.sil.df,majorCluster=="uncharacterized",select="sample")
    }
    
    gene.table.sortByF <- clusterMarker.res$gene.table[order(clusterMarker.res$gene.table$F,decreasing = T),]
    if(is.null(nKeep)){ 
        #g.f <- rownames(aov.res$aov.out.sig) 
        g.f <- rownames(subset(gene.table.sortByF,AUC>0.65 & score.q.value <0.01))
    }else{ 
        g.f <- head(rownames(aov.res$aov.out.sig),n=nKeep) 
    }
}else{
    g.f <- as.character(read.table(geneFile,header = T,check.names = F,stringsAsFactors = F,sep = "\t")$geneID)
    g.f <- intersect(g.f,rownames(Y))
}
loginfo(sprintf("length of g.f: %d\n",length(g.f)))
##q()
if(args$noRF){
    save.image(file = sprintf("%s/%s.RData",out.dir,sample.id))
}

#### random forest
ntree <- 5000
ntreeIterat <- 2000
#ntree <- 500
#ntreeIterat <- 100

#### using important genes
rfsel <- varSelRF(t(Y[g.f,x.label.df$sample]),
                      as.factor(x.label.df$majorCluster), 
                      ntree=ntree, ntreeIterat=ntreeIterat, vars.drop.frac=0.2,
                      whole.range = FALSE,keep.forest = T)

print(g.GNAME[rfsel$selected.vars])
out.importantGene.df <- data.frame(geneID=rfsel$selected.vars,
                                   geneSymbol=g.GNAME[rfsel$selected.vars],stringsAsFactors = F,
                                   importance=rfsel$initialImportances[rfsel$selected.vars,1])
out.importantGene.df <- out.importantGene.df[order(out.importantGene.df$importance,decreasing = T),]
write.table(out.importantGene.df,sprintf("%s.rfsel.importantGenes.txt",out.prefix),row.names = F,sep = "\t",quote = F)
print(rfsel)
print(rfsel$rf.model)

pdf(sprintf("%s.varSelRF.pdf",out.prefix),width = 8,height = 6)
plot(rfsel)
dev.off()

pdf(sprintf("%s.varSelRF.confusion.pdf",out.prefix),width = 8,height = 8)
dat.plot <- rfsel$rf.model$confusion[,-ncol(rfsel$rf.model$confusion),drop=F]
#heatmap.2(dat.plot,col=bluered(100), Rowv = T, Colv = T, scale="row",
#heatmap.2(dat.plot,col=colorRampPalette(brewer.pal(9,"Blues"))(20), Rowv = T, Colv = T, scale="row",
heatmap.2(dat.plot,col=colorRampPalette(brewer.pal(9,"Blues"))(9), Rowv = T, Colv = T, scale="row",
          density.info="none", dendrogram="none", keysize=1.2,
          trace="none", margin=c(10, 10),cexRow=1,cexCol=1)
dev.off()

.tmp.df <- data.frame(majorCluster=rownames(rfsel$rf.model$confusion),
                      ncells=apply(dat.plot,1,sum),
                      oob.error=rfsel$rf.model$confusion[,ncol(rfsel$rf.model$confusion)],stringsAsFactors = F)

if(!is.null(Y.y)){
    pred.new.rfsel <- predict(rfsel$rf.model, newdata = t(Y.y[rfsel$selected.vars,y.df$sample]),type = "prob")
}else{
    pred.new.rfsel <- predict(rfsel$rf.model, newdata = t(Y[rfsel$selected.vars,y.df$sample]),type = "prob")
}


cls <- colnames(pred.new.rfsel)
pred.new.rfsel.vec <-  apply(pred.new.rfsel,1,function(x){ cls[which.max(x)] })
names(pred.new.rfsel.vec) <- rownames(pred.new.rfsel)

###pred.new.rfsel <- structure(as.character(pred.new.rfsel),names=names(pred.new.rfsel))

if(!is.null(args.y.file) && file.exists(args.y.file)){
    pred.new.rfsel.df <- cbind(y.df, data.frame(majorCluster=pred.new.rfsel.vec,stringsAsFactors = F))
}else{
    pred.new.rfsel.df <- cbind(myDesign[y.df$sample,-i.mcls],
                         data.frame(majorCluster=pred.new.rfsel.vec,stringsAsFactors = F))
}

print(table(pred.new.rfsel.df[,c("patient","majorCluster")]))
write.table(pred.new.rfsel.df,sprintf("%s.rfsel.txt",out.prefix),row.names = F,sep = "\t",quote = F)

if(!is.null(val.res.sil.df.flt)){
    plot.tsne.points(tsne.Y[c(rownames(val.res.sil.df.flt),pred.new.rfsel.df$sample),,drop=F],
                     sprintf("%s.predMCls.tsne.rm00.pred00.pdf",out.prefix),
                     tsne.col.points=clusterColor[c(val.res.sil.df.flt$majorCluster,pred.new.rfsel.df$majorCluster)],
                     col.tsne.legend=clusterColor,
                     tsne.legend=names(clusterColor),
                     pch=20,nclusters=length(clusterColor))
}

#### using all genes from g.f
#rf.model <- randomForest(t(Y[g.f,x.label.df$sample]),
#                         as.factor(x.label.df$majorCluster),
#                         importance=T,proximity=T,ntree=ntree)
#print(rf.model)
#### prediction on trainning set
##.pred.train <- predict(rf.model,t(Y[g.f,x.label.df$sample]))
##print(table(observed=as.factor(x.label.df$majorCluster),
##                predicted=.pred.train))
#### prediction on new set
#
#if(!is.null(Y.y)){
#    pred.new <- predict(rf.model,t(Y.y[g.f,y.df$sample]))
#}else{
#    pred.new <- predict(rf.model,t(Y[g.f,y.df$sample]))
#}
#
#if(!is.null(args.y.file) && file.exists(args.y.file)){
#    pred.new.df <- cbind(y.df, data.frame(majorCluster=pred.new,stringsAsFactors = F))
#}else{
#    pred.new.df <- cbind(myDesign[y.df$sample,-i.mcls],
#                     data.frame(majorCluster=pred.new,stringsAsFactors = F))
#}
#
#print(table(pred.new.df[,c("patient","majorCluster")]))
#write.table(pred.new.df,sprintf("%s.txt",out.prefix),row.names = F,sep = "\t",quote = F)

save.image(file = sprintf("%s/%s.RData",out.dir,sample.id))


