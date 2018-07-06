#!/usr/bin/env Rscript

################################################################
suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-s", "--sname", type="character", required=TRUE, help="serial name, such as GSE67801")
parser$add_argument("-p", "--pname", type="character", required=TRUE, help="platform name, such as GPL6480")

args <- parser$parse_args()
print(args)

sname <- args$sname
pname <- args$pname

#sname <- "GSE67801"
#pname <- "GPL6480"

library(Biobase)
library(GEOquery)

collapseGSet <- function(X){
    X <- X[!is.na(fData(X)$Gene.ID) & fData(X)$Gene.ID!="",]
    datET <- exprs(X)
    .fdat <- fData(X)
    datET[1:4,1:6]
    ### simple but faster version
    datET.rowmean <- sort(rowMeans(datET),decreasing = T)
    datET.geneID <- .fdat[names(datET.rowmean),"Gene.ID"]
    f.probe <- names(datET.rowmean)[which(!duplicated(datET.geneID))]
    datET.collapsed.a <- datET[f.probe,]
    rownames(datET.collapsed.a) <- .fdat[f.probe,"Gene.ID"]
    ### robust but slow version
    #datET.collapsed <- collapseRows(datET,.fdat$Gene.ID,rownames(datET),method = "MaxMean")
    ### construct ExpressionSet obj
    f.fDat <- .fdat[f.probe,]
    rownames(f.fDat) <- f.fDat[,"Gene.ID"]
    f.fDat[,"ID"] <- f.fDat[,"Gene.ID"]
    X.ret <- ExpressionSet(datET.collapsed.a,phenoData=phenoData(X),featureData=AnnotatedDataFrame(f.fDat))
    X.ret <- X.ret[order(rownames(X.ret)),]
}

# load series and platform data from GEO
gset <- getGEO(sname, GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep(pname, attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))
# log2 transform
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) || (qx[6]-qx[1] > 50 && qx[2] > 0) || (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { 
    ex[which(ex <= 0)] <- NaN
    exprs(gset) <- log2(ex) 
}

gset.g <- collapseGSet(X = gset)

save(gset,gset.g,file = sprintf("%s.%s.RData",sname,pname))

