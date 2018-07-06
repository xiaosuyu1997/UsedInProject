#!/usr/bin/env Rscript

args <- commandArgs(T)
if(length(args)<5)
{
    cat("class.differentiation.R <sample id> <design file> <output.prefix> <exp.tab> <CD45.tab> [exp threshold]\n")
    q()
}
sampleID <- args[1]
design.file <- args[2]
out.prefix <- args[3]
exp.tab.file <- args[4]
CD45.tab.file <- args[5]
if(length(args)>=6)
{
    exp.threshold <- args[6]
}else
{
    exp.threshold <- 1
}

source("/Share/BP/zhenglt/02.pipeline/cancer/lib/myFunc.R")
suppressPackageStartupMessages(library("RColorBrewer"))

#sampleID <- "P0205"
#design.file <- "P0205.filter.by.marker.design.txt"
#out.prefix <-  "classDifferentiation"
#exp.tab.file <- "RPKM.filter1.fltByMarker/P0205.marker.naive-mature/Marker.n9.P0205.txt"
#CD45.tab.file <- "PTPRC/P0205.CD45.type.slim.txt"
#exp.threshold <- 1

loginfo(paste0("process sample ",sampleID))
dir.create(out.prefix,recursive = T,showWarnings = F)

in.table <- read.table(exp.tab.file,sep="\t",stringsAsFactors = F,header = T,row.names = 1,check.names = F)
in.table <- as.matrix(in.table[,-1])

myDesign<-read.table(design.file,header=T,check.names=F,colClasses=c("factor","character","factor","character"))
rownames(myDesign) <- myDesign$sample
myDesign$sampleType <- factor(myDesign$sampleType,levels=unique(myDesign$sampleType))
in.table <- in.table[,rownames(myDesign)]

CD45.tab <- read.table(CD45.tab.file,row.names = 1,header = T,check.names = F,stringsAsFactors = F)
CD45.tab <- CD45.tab[rownames(myDesign),,drop=FALSE]
print(head(myDesign))
print(in.table[1:4,1:8])
print(head(CD45.tab))

out <- CD45.tab
#1236    CCR7
#6402    SELL    (CD62L)
out$CCR7Status <- "Negative"
out$CCR7Status[in.table["1236",] >= exp.threshold] <- "Positive"
out$CD62LStatus <- "Negative"
out$CD62LStatus[in.table["6402",] >= exp.threshold] <- "Positive"
out$DifferentiationStatus <- "Other"
out$DifferentiationStatus[out$CD45Status=="CD45RA" & out$CCR7Status=="Positive" & out$CD62LStatus=="Positive"] <- "naive"
out$DifferentiationStatus[out$CD45Status=="CD45RO" & out$CCR7Status=="Positive" & out$CD62LStatus=="Positive"] <- "central memory"
out$DifferentiationStatus[out$CD45Status=="CD45RO" & out$CCR7Status=="Negative" & out$CD62LStatus=="Negative"] <- "effector memory"
out$DifferentiationStatus[out$CD45Status=="CD45RA" & out$CCR7Status=="Negative" & out$CD62LStatus=="Negative"] <- "effector memory"
out <- cbind(myDesign,out)
write.table(out,paste0(out.prefix,"/",sampleID,".design.Extend.txt"),sep="\t",quote = F,row.names = F)
diff.stat.summary <- t(with(out,table(sampleType, DifferentiationStatus)))
out.diff.stat.summary <- data.frame(DiffStat=rownames(diff.stat.summary))
out.diff.stat.summary <- cbind(out.diff.stat.summary,diff.stat.summary[,])
write.table(out.diff.stat.summary,paste0(out.prefix,"/",sampleID,".diffStat.summary.txt"),sep="\t",quote = F,row.names = F)
print(diff.stat.summary)

stat.type <- unique(out$DifferentiationStatus)
colSet <- brewer.pal(length(stat.type),"Blues")
names(colSet) <- stat.type
dat.plot <- apply(diff.stat.summary,2,function(x){ x/sum(x)  })
pdf(paste0(out.prefix,"/",sampleID,".diffStat.pdf"),width=12,height=8)
par(mar=c(5,6,4,12),cex.lab=2)
barplot(dat.plot[stat.type,],beside = F,xlab="Cell Type",ylab="Differentiation Status",main=sampleID,col=colSet,cex.main=2)
legend("right",legend=stat.type,fill=colSet,inset = -0.25,xpd = T,horiz = F,cex=1.3)
dev.off()


