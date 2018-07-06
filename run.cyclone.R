#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-i", "--infile", type="character", required=TRUE, help="infile")
parser$add_argument("-o", "--outprefix", type="character", required=TRUE, help="outprefix")
parser$add_argument("-s", "--sample", type="character", default="SAMPLE", help="sample id  [default %(default)s]")
parser$add_argument("-d", "--designFile", required=TRUE, type="character", help="sample desc file")
args <- parser$parse_args()
print(args)

in.file <- args$infile
out.prefix <- args$outprefix
sample.id <- args$sample
sample.desc.file <- args$designFile

#in.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/gsnap_out/quantification/P1116.count.tab.gz"
#out.prefix <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/cross.patient/cellCycle/P1116.cellCycle"
#sample.id <- "P1116"
#sample.desc.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/exclude_outlier/P1116.design.xOutlier.sf.txt"

source("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")
suppressPackageStartupMessages(library("scran"))

sample.desc <- read.table(sample.desc.file,header = T,sep = "\t",check.names = F,stringsAsFactors = F)
rownames(sample.desc) <- sample.desc$sample
print(dim(sample.desc))
print(head(sample.desc))
print(str(sample.desc))

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

f <- apply(Y,1,function(x){ nE <- sum(x>0); return( nE > 5 & nE/length(x) > 0.01 )  })
Y <- as.matrix(Y[f,])

print(dim(Y))
print(Y[1:4,1:6])

#### 
human.pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))
### need ENSGXXX 
fData <- data.frame(ensemblID=entrezToXXX(rownames(Y),type = "ENSG"),geneID=rownames(Y),geneSymbol=entrezToXXX(rownames(Y)),stringsAsFactors = F)
f <- !is.na(fData$ensemblID) & !duplicated(fData$ensemblID)
Y <- Y[f,]
fData <- fData[f,]
rownames(Y) <- fData$ensemblID
rownames(fData) <- fData$ensemblID

#sce <- newSCESet(tpmData=Y,phenoData = new("AnnotatedDataFrame", data = sample.desc),
#                 featureData = new("AnnotatedDataFrame",data=fData))
#assigned <- cyclone(sce, human.pairs,assay="tpm")
sce <- newSCESet(countData=Y,phenoData = new("AnnotatedDataFrame", data = sample.desc),
                 featureData = new("AnnotatedDataFrame",data=fData))
assigned <- cyclone(sce, human.pairs,assay="counts")

###out.df <- data.frame(sample=sample.desc$sample,check.names = F,stringsAsFactors = F)
out.df <- sample.desc
out.df <- cbind(out.df, assigned$score)
out.df$cyclePhase <- "S"
out.df$cyclePhase[assigned$scores$G1 > 0.5 ] <- "G1"
out.df$cyclePhase[assigned$scores$G2M > 0.5] <- "G2M"
out.df$cyclePhase[assigned$scores$G1 > 0.5 & assigned$scores$G2M > 0.5] <- "unknown"
table(out.df$phase)
write.table(out.df,sprintf("%s.txt",out.prefix),row.names = F,quote = F,sep = "\t")

pdf(sprintf("%s.pdf",out.prefix),width=8,height=8)
par(mar=c(6,6,4,2),cex.lab=2.0,cex.main=2.0,cex.axis=2.0)
plot(assigned$score$G1, assigned$score$G2M, col="lightblue", pch=16,xlab="G1",ylab="G2M",main=sprintf("CellCycle (%s)",sample.id))
abline(h=0.5,lty=2,lwd=1.2)
abline(v=0.5,lty=2,lwd=1.2)
dev.off()

