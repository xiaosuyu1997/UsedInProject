#!/usr/bin/env Rscript

args <- commandArgs(T)
if(length(args)<4)
{
    cat("filter.by.marker.R <sample id> <design file> <output.prefix> <in.tab> [exp threshold]\n")
    q()
}
sampleID <- args[1]
design.file <- args[2]
out.prefix <- args[3]
tab.file <- args[4]
if(length(args)>=5)
{
    exp.threshold <- args[5]
}else
{
    exp.threshold <- 1
}

source("/Share/BP/zhenglt/02.pipeline/cancer/lib/myFunc.R")
####source("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")

## threshold for unexpressed genes
#exp.threshold <- 10
#unexp.threshold <- 1
exp.threshold <- 30
unexp.threshold <- 3

#sampleID <- "P1116"
#design.file <- "../sample.design.P1116"
#out.prefix <-  "test.P1116.filter.by.marker"
#tab.file <- "RPKM.filter0/P1116.marker.TCellType/Marker.n7.P1116.txt"


loginfo(paste0("process sample ",sampleID))

in.table <- read.table(tab.file,sep="\t",stringsAsFactors = F,header = T,row.names = 1,check.names = F)
in.table <- as.matrix(in.table[,-1])

####myDesign<-read.table(design.file,header=T,check.names=F,colClasses=c("factor","character","factor","character"))
myDesign<-read.table(design.file,header=T,check.names=F)
rownames(myDesign) <- myDesign$sample
myDesign$sampleType <- factor(myDesign$sampleType,levels=unique(myDesign$sampleType))
in.table <- in.table[,rownames(myDesign)]

print(head(myDesign))
print.head(in.table)
#920     CD4
#925     CD8A
#926     CD8B
#915     CD3D
#916     CD3E
#917     CD3G
#3559    IL2RA

f.typeTC <- (myDesign$sampleType=="PTC" | myDesign$sampleType=="TTC" | myDesign$sampleType=="NTC" | myDesign$sampleType=="JTC")
f.typeTHR <- (!f.typeTC)
in.table.exp.binary <- (in.table >= exp.threshold)
in.table.unexp.binary <- (in.table < unexp.threshold)
print.head(in.table.exp.binary)
print.head(in.table.unexp.binary)
## CD3D
toRemove.T <- which(in.table.unexp.binary["915",])
toRemove.TC <-  which( (in.table.unexp.binary["925",] | in.table.exp.binary["920",])  & f.typeTC)
toRemove.THR <- which( (in.table.unexp.binary["920",] | in.table.exp.binary["925",]) & f.typeTHR)

toRemove <- unique(c(toRemove.T,toRemove.TC,toRemove.THR))
#loginfo("remove T cell")
#print(toRemove.T)
#loginfo("remove TC cell")
#print(toRemove.TC)
#loginfo("remove TH/TR cell")
#print(toRemove.THR)
if(length(toRemove)>0){
    myDesign.filtered <- myDesign[-toRemove,]
}else{
    myDesign.filtered <- myDesign
}
loginfo("before removal")
print(dim(myDesign))
loginfo("after removal")
print(dim(myDesign.filtered))
write.table(myDesign.filtered,paste0(out.prefix,".design.txt"),sep="\t",quote = F,row.names = F)
myDesign.Extend <- myDesign
myDesign.Extend$CD3.Failed <- rep(FALSE,nrow(myDesign))
myDesign.Extend$CD8.Failed <- rep(FALSE,nrow(myDesign))
myDesign.Extend$CD4.Failed <- rep(FALSE,nrow(myDesign))
myDesign.Extend$CD3.Failed[toRemove.T] <- TRUE
myDesign.Extend$CD8.Failed[toRemove.TC] <- TRUE
myDesign.Extend$CD4.Failed[toRemove.THR] <- TRUE
write.table(myDesign.Extend,paste0(out.prefix,".design.Extend.txt"),sep="\t",quote = F,row.names = F)


