#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-i", "--infile", type="character", help="infile")
parser$add_argument("-c", "--countFile", type="character", help="countFile")
parser$add_argument("-t", "--tpmFile", type="character", help="tpmFile")
parser$add_argument("-o", "--outprefix", type="character", required=TRUE, help="outprefix")
parser$add_argument("-s", "--sample", type="character", default="SAMPLE", help="sample id  [default %(default)s]")
parser$add_argument("-d", "--designFile", required=TRUE, type="character", help="sample desc file")
parser$add_argument("-l", "--log", action="store_true", default=FALSE, help="whether do log transform [default %(default)s]")
parser$add_argument("-r", "--center", action="store_true", default=FALSE, help="whether center genes by individual [default %(default)s]")
args <- parser$parse_args()
print(args)
in.file <- args$infile
out.prefix <- args$outprefix
sample.id <- args$sample
sample.desc.file <- args$designFile
args.log <- args$log
args.center <- args$center
countFile <- args$countFile
tpmFile <- args$tpmFile
if(is.null(in.file) && is.null(countFile)){ 
    cat("one of --infile and --countFile must be specified\n")
    q()
}

#in.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/cross.patient/sf.byIndividual/sfDat.liver.txt.gz"
#out.prefix <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/gsnap_out/quantification/scater/liver.P5"
#sample.id <- "liver"
#sample.desc.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/cross.patient/cellCycle/liver.leaf.addClonotype.addCyclone.txt"
#args.log <- T
#args.center <- T
#countFile <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/cross.patient/rc/countDat.liver.txt.gz"
#tpmFile <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/cross.patient/tpm/tpmDat.liver.txt.gz"

source("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")
suppressPackageStartupMessages(library("scran"))
suppressPackageStartupMessages(library("scater"))

loginfo("make.SCESet.R begin ...")

sample.desc <- read.delim(sample.desc.file,header = T,sep = "\t",check.names = F,stringsAsFactors = F)
rownames(sample.desc) <- sample.desc$sample
print(dim(sample.desc))
print(head(sample.desc))
print(str(sample.desc))
sname <- rownames(sample.desc)
f <- NULL
Y <- NULL
Y.isExpressed <- NULL
countData <- NULL
tpmData <- NULL

if(!is.null(in.file) && file.exists(in.file)){
    if(grepl("RData$",in.file,perl=T)){
        suppressPackageStartupMessages(library("R.utils"))
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
    #f <- apply(Y,1,function(x){ nE <- sum(x>0); return( nE > 5 & nE/length(x) > 0.01 )  })
    #Y <- as.matrix(Y[f,])
    Y.isExpressed <- Y>0
    if(args.log) { Y <- log2(Y+1) }
    if(args.center){
        Y.new <- c()
        for(pp in unique(sample.desc$patient)){
            Y.block <- t(scale(t(Y[,subset(sample.desc,patient==pp,select="sample",drop=T)]),center = T,scale = F))
            Y.new <- cbind(Y.new,Y.block)
            print(apply(Y.block[1:4,],1,mean))
            print(apply(Y.block[1:4,],1,sd))
        }
        Y <- Y.new
        Y <- Y[,sname,drop=F]
    }
    print(dim(Y))
    print(Y[1:4,1:6])
}
if(!is.null(countFile) && file.exists(countFile)){
    countData <- read.table(countFile,header = T,sep = "\t",stringsAsFactors = F,check.names = F)
    rownames(countData) <- countData[,1]
    countData <- countData[,c(-1,-2)]
    countData <- as.matrix(countData[,sname])
    if(!is.null(f)){ countData <- countData[f,] }
}
if(!is.null(tpmFile) && file.exists(tpmFile)){
    tpmData <- read.table(tpmFile,header = T,sep = "\t",stringsAsFactors = F,check.names = F)
    rownames(tpmData) <- tpmData[,1]
    tpmData <- tpmData[,c(-1,-2)]
    tpmData <- as.matrix(tpmData[,sname])
    if(!is.null(f)){ tpmData <- tpmData[f,] }
}
### add ENSGXXX 
if(!is.null(Y)){
    fData <- data.frame(ensemblID=entrezToXXX(rownames(Y),type = "ENSG"),
                        geneID=rownames(Y),
                        geneSymbol=entrezToXXX(rownames(Y)),
                        stringsAsFactors = F)
}else{
    fData <- data.frame(ensemblID=entrezToXXX(rownames(countData),type = "ENSG"),
                        geneID=rownames(countData),
                        geneSymbol=entrezToXXX(rownames(countData)),
                        stringsAsFactors = F)
}
rownames(fData) <- fData$geneID
##f <- !is.na(fData$ensemblID) & !duplicated(fData$ensemblID)
##Y <- Y[f,]
##fData <- fData[f,]

if(!is.null(Y)){
    sce <- newSCESet(exprsData=Y, is_exprsData= Y.isExpressed, phenoData = new("AnnotatedDataFrame", data = sample.desc),
                     featureData = new("AnnotatedDataFrame",data=fData),logged = T)
}else{
    sce <- newSCESet(countData=countData,phenoData = new("AnnotatedDataFrame", data = sample.desc),
                     featureData = new("AnnotatedDataFrame",data=fData),logged = F)
}
if(!is.null(countData) && is.null(counts(sce))){
    counts(sce) <- countData
}
if(!is.null(tpmData) && is.null(tpm(sce))){
    tpm(sce) <- tpmData
}

save(sce,file = sprintf("%s.RData",out.prefix))

loginfo("make.SCESet.R run successfully")
