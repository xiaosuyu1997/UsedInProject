#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-d", "--designFile", type="character", required=TRUE, help="design matrix file")
parser$add_argument("-i", "--inputFile", type="character", required=TRUE, help="input file")
parser$add_argument("-c", "--cellTypeColorFile", default="/WPS1/zhenglt/work/TCR_chunhong/data/CellType.color", 
                    type="character", help="cellTypeColorFile [default %(default)s]")
parser$add_argument("-e", "--clonotypeFile", type="character", help="clonotype file")
parser$add_argument("-o", "--outDir", type="character", required=TRUE, help="output directory")
parser$add_argument("-s", "--sample", type="character", default="SAMPLE", help="sample id  [default %(default)s]")
###parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="whether verbose mode [default %(default)s]")
parser$add_argument("-l", "--log", action="store_true", default=FALSE, help="whether do log transform [default %(default)s]")
parser$add_argument("-r", "--center", action="store_true", default=FALSE, help="whether center genes by individual [default %(default)s]")
parser$add_argument("-n", "--nKeep", type="integer", default=1500, help="use top nKeep genes [default %(default)s]")
parser$add_argument("-k", "--nClusters", type="character", default=4, help="number of clusters [default %(default)s]")
args <- parser$parse_args()

print(args)

designFile <- args$designFile
inputFile <- args$inputFile
cellTypeColorFile <- args$cellTypeColorFile
clonotype.file <- args$clonotypeFile
out.dir <- args$outDir
sample.id <- args$sample
args.log <- args$log
args.center <- args$center
nKeep <- args$nKeep
nClusters <- args$nClusters

#### TEST DATA 
###nKeep <- 1500
###designFile <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/cross.patient/test/test.design.001.txt"
###inputFile <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/cross.patient/OUT.centered/liver.TC/liver.all.RData"
###cellTypeColorFile <- "/WPS1/zhenglt/work/TCR_chunhong/data/CellType.color"
###clonotype.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/cross.patient/clonotype/clonotype.liver.txt"
###out.dir <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/cross.patient/test/OUT.SIMLR"
###sample.id <- "liver"
###args.log <- FALSE
###args.center <- FALSE
###nClusters <- "2,3,4,5,6,7,8"

dir.create(out.dir,recursive = T,showWarnings = F)
source("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")

myDesign<-read.table(designFile,header=T,check.names=F,stringsAsFactors = F)
rownames(myDesign) <- myDesign$sample
sampleTypeColor <- read.SampleTypeColor(cellTypeColorFile)

clonotype.strict.data <- read.clonotype(in.file = clonotype.file,ctype.col = "C_strict")
patient.col.list <- patientColorListFromMyDesign(myDesign)

if(grepl("RData$",inputFile,perl=T)){
    suppressPackageStartupMessages(library("R.utils"))
    loginfo(sprintf("loading RData: %s ...",inputFile))
    lenv <- loadToEnv(inputFile)
    Y <- lenv[["Y"]]
    ###obj.scdn <- lenv[["obj.scdn"]]
    ###Y <- obj.scdn@normalized.endo
}else{
    loginfo(sprintf("reading file: %s ...",inputFile))
    in.table <- read.table(inputFile,header = T,sep = "\t",stringsAsFactors = F,check.names = F)
    rownames(in.table) <- in.table[,1]
    Y <- in.table[,c(-1,-2)]
}

sname <- intersect(rownames(myDesign),colnames(Y))
myDesign <- myDesign[sname,,drop=F]
Y <- Y[,sname,drop=F]

f <- apply(Y,1,function(x){ nE <- sum(x>0); return( nE > 5 & nE/length(x) > 0.01 )  })
Y <- Y[f,]
if(args.log) { Y <- log2(Y+1) }
if(args.center){
    Y.new <- c()
    for(pp in unique(myDesign$patient)){
        Y.block <- t(scale(t(Y[,subset(myDesign,patient==pp,select="sample",drop=T)]),center = T,scale = F))
        Y.new <- cbind(Y.new,Y.block)
        print(apply(Y.block[1:4,],1,mean))
        print(apply(Y.block[1:4,],1,sd))
    }
    Y <- Y.new
}

Y.sd <- apply(Y, 1, sd)
Y.sd.sort <- sort(Y.sd,decreasing=TRUE)
Y.sd.df <- data.frame(geneID=names(Y.sd.sort),geneSymbol=entrezToXXX(names(Y.sd.sort)))
Y.sd.df <- cbind(Y.sd.df,Y.sd.sort)
write.table(Y.sd.df,sprintf("%s/%s.Y.sd.sort.txt",out.dir,sample.id),row.names = F,sep = "\t",quote = F)

loginfo("dimention of data to be used in SIMLR:")
print(dim(Y[names(Y.sd.sort)[1:nKeep],]))

loginfo("... run SIMLR ...")
res.simlr <- list()
for(kk in as.numeric(unlist(strsplit(nClusters,",")))){
    loginfo(sprintf("k=%s",kk))
    res.simlr[[as.character(kk)]] <- run.SIMLR(Y[names(Y.sd.sort)[1:nKeep],],n.clusters=kk,myseed=NULL,my.title=sample.id,out.prefix=sprintf("%s/%s.k%s",out.dir,sample.id,kk),
                           col.points = sampleTypeColor[as.character(myDesign$sampleType)],
                           col.legend=sampleTypeColor[names(sampleTypeColor) %in% unique(as.character(myDesign$sampleType))],
                           pch=(as.numeric(as.factor(myDesign$patient))-1) %% 26,cex=0.7,
                           legend.txt=c(names(sampleTypeColor)[names(sampleTypeColor) %in% unique(as.character(myDesign$sampleType))]),
                           legend.txt2=unique(as.character(myDesign$patient))
                           )
    ## output clustering result in txt
    out.df <- myDesign
    out.df$clusterBySIMLR <- sprintf("cluster%s",res.simlr[[as.character(kk)]]$y$cluster)
    write.table(out.df,sprintf("%s/%s.k%s.cluster.txt",out.dir,sample.id,kk),sep = "\t",quote = F,row.names = F)
}
loginfo("... end ...")
save(Y,res.simlr,file=sprintf("%s/%s.all.RData",out.dir,sample.id))
