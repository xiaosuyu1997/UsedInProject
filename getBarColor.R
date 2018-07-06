#!/usr/bin/env Rscript

args <- commandArgs(T)
if(length(args)<1)
{
    cat("getBarColor.R <design file> <out.prefix>\n")
    q()
}
designFile <- args[1]
out.prefix <- args[2]

suppressPackageStartupMessages(library("RColorBrewer"))

designM<-read.table(designFile,header=T,row.names="sample",check.names=F,colClasses=c("factor","character","factor","factor"))

if(length(levels(designM$sampleType))<=9)
{
    colSet <- brewer.pal(9,"Set1")
    patientcolors <- colSet[as.numeric(designM$sampleType)]
}else
{
    patientcolors <- as.numeric(designM$sampleType)
}
designM$color <- patientcolors
#print(head(designM))
#print(str(designM))
legendDat <- unique(designM[,c("sampleType","color")])
pdf(paste0(out.prefix,".colorScheme.pdf"),width=8,height=8)
plot(x=c(1:10),y=c(1:10),type = "n",xlab="",ylab="")
legend("center",legend=legendDat$sampleType,fill=legendDat$color,border=legendDat$color,cex=2)
dev.off()
