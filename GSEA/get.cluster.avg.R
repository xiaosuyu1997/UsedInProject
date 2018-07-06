#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-i", "--inFile", type="character", required=TRUE, help="infile")
parser$add_argument("-o", "--outPrefix", type="character", required=TRUE, help="outprefix")
parser$add_argument("-d", "--dFile", type="character", required=TRUE, help="design file")
parser$add_argument("-n", "--ncores", type="integer", default=4, help="number of cpu cores used [default %(default)s]")
args <- parser$parse_args()
print(args)

inputFile <- args$inFile
out.preifx <- args$outPrefix
designFile <- args$dFile
n.cores <- args$ncores

#inputFile <- "/WPSnew/zhenglt/work/proj_zl/integrate0924/cross.patient.b20170523/exp/colonC.P8.scran.RData"
#out.preifx <- "/WPSnew/zhenglt/work/proj_zl/integrate0924/cross.patient.b20170523/functional/ssGSEA/colonC.P8.cluster.avg"
#designFile <- "/WPSnew/zhenglt/work/proj_zl/integrate0924/cross.patient.b20170523/clusteringResult/colonC.clustering.final.forPlot.all.txt"
#n.cores <- 8
cellTypeColorFile <- NULL

dir.create(dirname(out.preifx),recursive = T,showWarnings = F)

if(file.exists("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.")){
    source("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")
}else{
    source("/WPSnew/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")
}

suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("doParallel"))

g.inputD <- processInput(designFile,cellTypeColorFile,inputFile,args.notFilter=F,geneFile = NULL,
                         args.center=F,args.log=F,args.norm.exprs = F)
myDesign <- g.inputD$myDesign
sampleTypeColor <- g.inputD$sampleTypeColor
patient.col.list <- g.inputD$patient.col.list
Y <- g.inputD$Y
g.GNAME <- g.inputD$g.GNAME

f <- apply(Y,1,function(x){ nE <- sum(x>0); return( nE > 5 & nE/length(x) > 0.01 & sum(is.na(x))==0 )  })
Y <- Y[f,]
print("dim of Y")
print(dim(Y))
print("all(colnames(Y)==rownames(myDesign)):")
print(all(colnames(Y)==rownames(myDesign)))

registerDoParallel(cores = n.cores)

avg.exp <- ldply(rownames(Y),function(v){
                     tryCatch({
                     .res <- aggregate(Y[v,],by=list(myDesign$majorCluster),FUN=mean)
                     .res.2 <- aggregate(Y[v,],by=list(myDesign$majorCluster),FUN=sd)
                     structure(c(.res[,2],.res.2[,2]),names=c(sprintf("avg.%s",.res[,1]),sprintf("sd.%s",.res[,1])))
                     } ,error=function(x){ print(v); print(e) } )
			},.progress = "none",.parallel=T)

rownames(avg.exp) <- as.character(rownames(Y))

avg.exp.df <- data.frame(geneID=rownames(avg.exp),
                         geneSymbol=g.GNAME[rownames(avg.exp)],
                         stringsAsFactors = F)
avg.exp.df <- cbind(avg.exp.df,avg.exp)
write.table(avg.exp.df,sprintf("%s.txt",out.preifx),row.names = F,sep = "\t",quote = F)
save(avg.exp,g.GNAME,file=sprintf("%s.RData",out.preifx))

