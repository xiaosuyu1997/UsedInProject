#!/usr/bin/env Rscript

args <- commandArgs(T)
if(length(args)<2)
{
    cat(sprintf("ERCC.QC.multi.analyze.R <in.file> <out.prefix>\n"))
    q()
}
in.file <- args[1]
out.prefix <- args[2]
##in.file <- "ERCC.QC.multi.txt"
##out.prefix <- "ERCC.QC.multi.analyze"
in.table <- read.table(in.file,header = T,check.names = F)
print(head(in.table))
EPSILON <- 1e-6

library(beeswarm)
pdf(paste0(out.prefix,".pdf"),width = 8,height = 8)
par(mar=c(5,6,4,2),cex.lab=1.5)
boxplot(nERCCRatio ~ dilution, data = data.frame(dilution=as.factor(in.table$dilution),nERCCRatio=log10(in.table[,"nERCC:nSample"]+EPSILON)),
        outline = FALSE,main = "nERCC:nSample ~ Dilution",xlab="Dilution",ylab=expression(log[10](nERCC:nSample)))
beeswarm(nERCCRatio ~ dilution, data = data.frame(dilution=as.factor(in.table$dilution),nERCCRatio=log10(in.table[,"nERCC:nSample"]+EPSILON)),col = 4, pch = 16, add = TRUE)

tryCatch({
boxplot(pearson ~ dilution,data=data.frame(pearson=in.table$pearson,dilution=as.factor(in.table$dilution)),outline = FALSE,main="Pearson ~ Dilution",xlab="Dilution",ylab="Pearson")
beeswarm(pearson ~ dilution,data=data.frame(pearson=in.table$pearson,dilution=as.factor(in.table$dilution)),outline = FALSE,main="Pearson ~ Dilution",col=4,pch=16,add=T)
},error=function(e){e})

tryCatch({
plot(pearson ~ nERCCRatio,data=data.frame(pearson=in.table$pearson,nERCCRatio=log10(in.table[,"nERCC:nSample"]+EPSILON)),
     ylim=c(0,1), xlab=expression(log[10](nERCC:nSample)),ylab="Pearson",main="Pearson ~ nERCC:nSample")
abline(h=0.8,lty=2)
abline(h=0.9,lty=2)
},error=function(e){e})

plot(detected ~ nERCCRatio,data=data.frame(detected=in.table$detected,nERCCRatio=log10(in.table[,"nERCC:nSample"]+EPSILON)),
     ylim=c(0,100), xlab=expression(log[10](nERCC:nSample)),ylab="Detected",main="Detected ~ nERCC:nSample")
abline(h=30,lty=2)
abline(h=40,lty=2)
dev.off()
