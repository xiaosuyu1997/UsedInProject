#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("ROCR"))

parser <- ArgumentParser()
parser$add_argument("-i", "--inFile", type="character", required=TRUE, help="infile")
parser$add_argument("-a", "--sampleDescFile", type="character", required=TRUE, help="sample.desc.file")
parser$add_argument("-b", "--geneFile", type="character", help="gene list file (\"ALL\" for all genes pass the internal filter)")
parser$add_argument("-o", "--outputPrefix", type="character", required=TRUE, help="output prefix")
parser$add_argument("-s", "--sample", type="character", default="SAMPLE", help="sample id  [default %(default)s]")
parser$add_argument("-m", "--majorCluster", type="character", default="MAIT", help="majorCluster  [default %(default)s]")
parser$add_argument("-e", "--assay", type="character", default="norm_exprs", help="assay [default %(default)s]")
parser$add_argument("-p", "--patient", action="store_true", default=FALSE, help="whether add patient in the model [default %(default)s]")
parser$add_argument("-n", "--controlCluster", type="character", default="Th", help="controlCluster  [default %(default)s]")
parser$add_argument("-c", "--logFC", type="double", default="1", help="logFC threshold(limma)  [default %(default)s]")
parser$add_argument("-d", "--fdr", type="double", default="0.05", help="fdr threshold(limma)  [default %(default)s]")
parser$add_argument("-l", "--log", action="store_true", default=FALSE, help="whether do log transform [default %(default)s]")
parser$add_argument("-r", "--center", action="store_true", default=FALSE, 
                    help="whether center genes by individual [default %(default)s]")
parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="verbose mode [default %(default)s]")
args <- parser$parse_args()
print(args)

in.file <- args$inFile
sample.desc.file <- args$sampleDescFile
gene.file <- args$geneFile
out.prefix <- args$outputPrefix
sample.id <- args$sample
args.majorCluster <- args$majorCluster
###cellTypeColorFile <- args$cellTypeColorFile
###clusterColorFile <- args$clusterColorFile
args.log <- args$log
args.center <- args$center
args.logFC <- args$logFC
args.fdr <- args$fdr
args.controlCluster <- args$controlCluster
##mode.verbose <- args$verbose
if(args.controlCluster=="\"\"") { 
    args.controlCluster <- "" 
}else{
    args.controlCluster <- gsub(pattern = "\"",replacement ="", args.controlCluster,perl = T)
}
args.patient <- args$patient
args.assay <- args$assay

#in.file <- "/WPSnew/zhenglt/work/proj_zl/integrate0924/cross.patient/exp/colonC.P5.scran.RData"
#sample.desc.file <- "/WPSnew/zhenglt/work/proj_zl/integrate0924/cross.patient/clusteringResult/colonC.CD8.clustering.forPlot.final.txt"
#gene.file <- "/WPSnew/zhenglt/work/proj_zl/integrate0924/cross.patient/marker/colonC.MAIT.m.tab"
#out.prefix <- "/WPSnew/zhenglt/work/proj_zl/integrate0924/cross.patient/test/colonC.marker"
#sample.id <- "colonC"
#args.majorCluster <- "MAIT"
###cellTypeColorFile <- "/WPS1/zhenglt/work/TCR_chunhong/data/CellType.color"
###clusterColorFile <- "/WPSnew/zhenglt/work/proj_zl/integrate0924/cross.patient/colonC.majorClusterColor.txt"
###args.log <- F
###args.center <- F

dir.create(dirname(out.prefix),showWarnings = F,recursive = T)
source("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")

suppressPackageStartupMessages(require("ComplexHeatmap"))
suppressPackageStartupMessages(require("circlize"))
suppressPackageStartupMessages(require("gridBase"))
suppressPackageStartupMessages(require("dendextend"))
suppressPackageStartupMessages(require("RColorBrewer"))
suppressPackageStartupMessages(require("gplots"))
suppressPackageStartupMessages(library("vioplot"))
suppressPackageStartupMessages(library("plotrix"))
suppressPackageStartupMessages(library("vioplot"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("Rmisc"))
suppressPackageStartupMessages(library("varSelRF"))
suppressPackageStartupMessages(library("scran"))
suppressPackageStartupMessages(library("SingleCellExperiment"))

loginfo("begin ...")


#### sample data
sample.desc <- read.table(sample.desc.file,header = T,sep = "\t",check.names = F,stringsAsFactors = F)
rownames(sample.desc) <- sample.desc$sample
print(dim(sample.desc))
print(head(sample.desc))
print(str(sample.desc))
##q()

### clonotype data
#clonotype.data <- read.clonotype(in.file = clonotype.file,ctype.col = "C_strict")
clonotype.data <- NULL
#### exp data

suppressPackageStartupMessages(library("R.utils"))
##if(grepl("\\.scran\\.RData$",in.file,perl=T)){
if(grepl("\\.RData$",in.file,perl=T) && grepl("\\.scran",in.file,perl=T)){
    lenv <- loadToEnv(in.file)
    #Y <- exprs(lenv[["sce.norm"]])
    ###Y <- assay(lenv[["sce.norm"]],"centered_norm_exprs")
    Y <- assay(lenv[["sce.norm"]],args.assay)
    args.log <- F
    args.center <- F
    ##g.GNAME <- fData(lenv[["sce.norm"]])[,"geneSymbol"]
    g.GNAME <- rowData(lenv[["sce.norm"]])[,"geneSymbol"]
    names(g.GNAME) <- rownames(Y)
    print(str(Y))
}else if(grepl("RData$",in.file,perl=T)){
    suppressPackageStartupMessages(library("R.utils"))
    lenv <- loadToEnv(in.file)
    Y <- lenv[["Y"]]
    g.GNAME <- entrezToXXX(rownames(Y))
    names(g.GNAME) <- rownames(Y)
}else{
    in.table <- read.table(in.file,header = T,sep = "\t",stringsAsFactors = F,check.names = F)
    rownames(in.table) <- in.table[,1]
    Y <- in.table[,c(-1,-2)]
    g.GNAME <- entrezToXXX(rownames(Y))
    names(g.GNAME) <- rownames(Y)
}
g.GNAME.f <- which(is.na(g.GNAME))
g.GNAME[g.GNAME.f] <- names(g.GNAME)[g.GNAME.f]

Y <- Y[names(g.GNAME)[which(!duplicated(g.GNAME))],]

sname <- intersect(rownames(sample.desc),colnames(Y))
sample.desc <- sample.desc[sname,,drop=F]
Y <- Y[,sname,drop=F]

#f <- apply(Y,1,function(x){ nE <- sum(x>0); return( nE > 5 & nE/length(x) > 0.01 )  })
#Y <- Y[f,]
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

### major cluster color
#nMajor <- length(unique(sample.desc$majorCluster))
#if(!is.null(clusterColorFile) && file.exists(clusterColorFile)){
#    majorClusterColor <- read.SampleTypeColor(clusterColorFile)
#    majorClusterColor <- majorClusterColor[ names(majorClusterColor) %in% gsub("\\..$","",unique(sample.desc[colnames(Y),"majorCluster"]),perl=T)]
#}else{
#    majorClusterColor <- structure(colorRampPalette(brewer.pal(nMajor,"Paired"))(nMajor),
#                              names=unique(sample.desc$majorCluster))
#}
#print("majorClusterColor:")
#print(majorClusterColor)

#' @importFrom ROCR prediction performance
#' @importFrom stats aggregate wilcox.test
## code from SC3
getAUC <- function(gene, labels) 
{
    score <- rank(gene)
    # Get average score for each cluster
    ms <- aggregate(score ~ labels, FUN = mean)
    # Get cluster with highest average score
    posgroup <- ms[ms$score == max(ms$score), ]$labels
    # Return negatives if there is a tie for cluster with highest average score
    # (by definition this is not cluster specific)
    if(length(posgroup) > 1) {
        return (c(-1,-1,1))
    }
    # Create 1/0 vector of truths for predictions, cluster with highest
    # average score vs everything else
    truth <- as.numeric(labels == posgroup)
    #Make predictions & get auc using RCOR package.
    pred <- prediction(score,truth)
    val <- unlist(performance(pred,"auc")@y.values)
    pval <- suppressWarnings(wilcox.test(score[truth == 1],
                                         score[truth == 0])$p.value)
    return(c(val,posgroup,pval))
}

run.limma.normalizedExp <- function(inptX,Grp,sample.descD,T.fdr=0.05,T.logFC)
{
	suppressPackageStartupMessages(library("limma"))
	suppressPackageStartupMessages(library("BiocParallel"))
    if(!all(colnames(inptX)==rownames(sample.descD))){
        cat("input data inconsistant!\n")
    }
    ##Grp <- sapply(sample.descD$majorCluster,function(x){ x==target.grp })
    sample.descD$Grp <- Grp
    if(args.patient){
	    design <- model.matrix(~patient+Grp,data=sample.descD)
    }else{
	    design <- model.matrix(~Grp,data=sample.descD)
    }
	#design <- model.matrix(~Grp)
	colnames(design)[length(colnames(design))]<-"II"

	register(MulticoreParam(8))
	fit <- lmFit(inptX, design)
	fit <- eBayes(fit)
	all.table  <- topTable(fit, coef = "II", n = Inf, sort = "p", p = 1)
	#out.resSig <- topTable(fit, coef = "II", n = Inf, sort = "p", p = fdr)

    .Grp.mean <- t(apply(inptX,1,function(x){
                             .mexp <- aggregate(x ~ Grp, FUN = mean)
                             structure(.mexp[,2],names=sprintf("mean.Grp%s",.mexp[,1])) }))

	all.table <- cbind(data.frame(geneID=rownames(all.table),stringsAsFactors = F),
	      data.frame(geneSymbol=entrezToXXX(rownames(all.table)),stringsAsFactors = F),
	      all.table,.Grp.mean[rownames(all.table),])
	out.resSig <- subset(all.table,adj.P.Val<T.fdr & abs(logFC)>T.logFC)

    geneAUCs <- as.data.frame(t(sapply(rownames(out.resSig),function(x){ getAUC(inptX[x,],Grp) })))
    colnames(geneAUCs) <- c("AUC","cluster","score.p.value")
    out.resSig <- cbind(out.resSig,geneAUCs)
    #all(rownames(geneAUCs)==rownames(out.resSig))
	#write.table(all.table,paste0(output.prefix,".all"),row.names = F,quote = F,sep = "\t")
	#write.table(out.resSig,paste0(output.prefix,".sig"),row.names = F,quote = F,sep = "\t")
	list(all=all.table,sig=out.resSig)
	#head(all.table)
	#head(out.resSig)
}

makeGroup <- function(x,n){
    .g <- rep(-1,length(x))
    names(.g) <- n
    .g[grepl(pattern = sprintf("%s",args.controlCluster),x,perl = T)] <- 0
    .g[grepl(pattern = sprintf("%s",args.majorCluster),x,perl = T)] <- 1
    #.g[x==args.majorCluster] <- 1
    .g[.g!=-1]
}
group <- makeGroup(sample.desc$majorCluster,rownames(sample.desc))
print("group:")
print(str(group))
print(table(group))
print(table(sample.desc$majorCluster))

loginfo("begin varSelRF.")
if(!is.null(gene.file) && file.exists(gene.file)){
    gene.table <- read.table(gene.file,header = T,sep = "\t",check.names = F,stringsAsFactors = F)
    rownames(gene.table) <- gene.table[,1]
    .f.g <- group
    if(length(group) > 2*sum(group==1)){
        .f.g <- group[c(sample(names(group)[which(group==0)],size = sum(group==1),replace = F),names(group)[group==1])]
    }
    rfsel <- varSelRF(t(Y[rownames(gene.table),names(.f.g)]),
                      as.factor(group[names(.f.g)]), 
                      ntree=10000, ntreeIterat=2000, vars.drop.frac=0.2)
    rf.model <- randomForest(t(Y[rownames(gene.table),names(.f.g)]),
                             as.factor(group[names(.f.g)]),
                             importance=T,proximity=T,ntree=10000)
    print(g.GNAME[rfsel$selected.vars])
    print(rfsel)
    print(rf.model)
    .pred.train <- predict(rf.model,t(Y[rownames(gene.table),names(.f.g)]))
    print(table(observed=as.factor(group[names(.f.g)]),
                predicted=.pred.train))
    ###pred <- prediction(.pred.train,as.factor(group[names(.f.g)]))
    ###val <- unlist(performance(pred,"auc")@y.values)
}else{
    out.limma <- run.limma.normalizedExp(inptX=Y[,names(group)],Grp = group,sample.descD=sample.desc[names(group),],T.fdr = args.fdr,T.logFC = args.logFC)
    subset(out.limma$sig,logFC>args.logFC) %>% 
        write.table(file = sprintf("%s.limma.fdr%.2f.logFC%.1f.sig.txt",out.prefix,args.fdr,args.logFC),
                    quote = F,sep = "\t",row.names = F)
    if(args.logFC!=2){
        subset(out.limma$sig,logFC>2) %>% 
            write.table(file = sprintf("%s.limma.fdr%.2f.logFC%.1f.sig.txt",out.prefix,args.fdr,2),
                        quote = F,sep = "\t",row.names = F)
    }
    if(args$verbose){
        out.limma$all %>% write.table(file = sprintf("%s.limma.fdr%.2f.logFC%.1f.all.txt",out.prefix,args.fdr,args.logFC),quote = F,sep = "\t",row.names = F)
    }
#    rfsel <- varSelRF(t(Y[rownames(out.limma$sig),]),as.factor(group), ntree=10000, ntreeIterat=2000, vars.drop.frac=0.2)
}

loginfo("end varSelRF.")

#### tmp
#out.df <- data.frame(geneID=rownames(Y),geneSymbol=g.GNAME[rownames(Y)],stringsAsFactors = F)
#write.table(out.df,"/WPSnew/zhenglt/work/proj_zl/integrate0924/cross.patient/exp/colonC.P5.expressedGene.list",row.names = F,quote = F,sep = "\t")


