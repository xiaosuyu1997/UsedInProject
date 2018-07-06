#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-i", "--inFile", type="character", required=TRUE, help="infile")
parser$add_argument("-a", "--sampleDescFile", type="character", required=TRUE, help="sample.desc.file")
parser$add_argument("-b", "--geneFile", type="character", required=TRUE, help="gene list file")
parser$add_argument("-t", "--trainLabelFile", type="character", required=TRUE, help="train label file")
parser$add_argument("-p", "--predLabelFile", type="character", required=TRUE, help="prediction label file")
parser$add_argument("-o", "--outputPrefix", type="character", required=TRUE, help="output prefix")
parser$add_argument("-s", "--sample", type="character", default="SAMPLE", help="sample id  [default %(default)s]")
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
train.label <- args$trainLabelFile
pred.label <- args$predLabelFile
sample.id <- args$sample
args.log <- args$log
args.center <- args$center
mode.verbose <- args$verbose

#in.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/cross.patient/OUT.centered/liver.TC.level1.v2/liver.all.RData"
#sample.desc.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/cross.patient/new.clusteringResult/liver.leaf.addMajor.CD8.txt"
#gene.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/cross.patient/new.clusteringResult/g.list.forSVM"
#out.prefix <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/cross.patient/new.clusteringResult/liver.leaf.addMajor.CD8.svm"
#train.label <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/cross.patient/new.clusteringResult/train.label.list"
#pred.label <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/cross.patient/new.clusteringResult/pred.label.list"
#sample.id <- "liver"
#args.log <- F
#args.center <- F
#mode.verbose <- TRUE

source("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")

loginfo("begin ...")

#### sample data
sample.desc <- read.table(sample.desc.file,header = T,sep = "\t",check.names = F,stringsAsFactors = F)
rownames(sample.desc) <- sample.desc$sample
print(dim(sample.desc))
print(head(sample.desc))
print(str(sample.desc))

train.cluster <- read.table(train.label,header = F,sep = "\t",check.names = F,stringsAsFactors = F)$V1
pred.cluster <- read.table(pred.label,header = F,sep = "\t",check.names = F,stringsAsFactors = F)$V1

#### exp data
suppressPackageStartupMessages(library("R.utils"))
if(grepl("\\.scran\\.RData$",in.file,perl=T)){
    lenv <- loadToEnv(in.file)
    Y <- exprs(lenv[["sce.norm"]])
    args.log <- F
    args.center <- F
}else if(grepl("RData$",in.file,perl=T)){
    lenv <- loadToEnv(in.file)
    Y <- lenv[["Y"]]
}else{
    in.table <- read.table(in.file,header = T,sep = "\t",stringsAsFactors = F,check.names = F)
    rownames(in.table) <- in.table[,1]
    Y <- in.table[,c(-1,-2)]
}
sname <- intersect(rownames(sample.desc),colnames(Y))
sample.desc <- sample.desc[sname,,drop=F]
Y <- Y[,sname,drop=F]

f <- apply(Y,1,function(x){ nE <- sum(x>0); return( nE > 5 & nE/length(x) > 0.01 )  })
Y <- Y[f,]
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

g.GNAME <- entrezToXXX(rownames(Y))
names(g.GNAME) <- rownames(Y)

#### gene data
Y.inGeneList <- Y
if(!is.null(gene.file) && file.exists(gene.file)){
    gene.desc <- read.table(gene.file,header = T,sep = "\t",check.names = F,stringsAsFactors = F)
    g.f <- intersect(gene.desc$geneID,rownames(Y))
    Y.inGeneList <- Y[g.f,,drop=F]
}

train.dat <- Y.inGeneList[, subset(sample.desc, majorCluster %in% train.cluster, select="sample",drop=T) ]
pred.dat  <- Y.inGeneList[, subset(sample.desc, majorCluster %in% pred.cluster, select="sample",drop=T) ]

library("e1071")

model <- tryCatch(svm( t(train.dat), factor(sample.desc[colnames(train.dat),"majorCluster"]), kernel = "linear"), error = function(e) e)
pred <- predict(model, t(pred.dat))

out.df <- sample.desc
out.df$old.majorCluster <- sample.desc$majorCluster
out.df[names(pred),"majorCluster"] <- as.character(pred)

write.table(out.df,file = sprintf("%s.txt",out.prefix),sep = "\t",row.names = F,quote = F)

loginfo("end.")
