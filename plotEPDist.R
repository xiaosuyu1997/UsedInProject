#!/usr/bin/env Rscript

### mclust 5.1 work ###

args<-commandArgs(T)
if(length(args)<2)
{
  cat("plotEPDist.R <ifile> <prefix> <keep.ERCC>\n")
  q()
}
ifile <- args[1]
prefix <- args[2]

print(ifile)
print(prefix)
print(getwd())


#ifile <- "231N.35e6.counts_gene.tab"
#prefix <- "231N.35e6.counts_gene"
#ifile <- "./NTC101-0617.counts_gene_exonic.tab"
#prefix <- "./NTC101-0617.counts_gene_exonic.tab.out"

library("data.table")

inTable<-fread(ifile)
#inTable<-read.table(ifile,header=T)


if(length(args)<3)
{
    f <- !grepl("^ERCC-",inTable$name,perl=T)
    inTable<-inTable[f,]
}
inTable$TPM<-90*inTable$count/inTable$width
t<-sum(inTable$TPM)
inTable$TPM<-inTable$TPM*1e6/t
write.table(inTable,paste0(prefix,".TPM.tab"),sep="\t",col.names=T,row.names=F,quote=F)
cat(paste0("Average RPKM of ",prefix,": ",mean(inTable$rpkm),"\n"))
cat(paste0("Average TPM of ",prefix,": ",mean(inTable$TPM),"\n"))

fitMixModel <- function(x,ofile,...)
{
  require(mclust)
  f<-is.finite(x)
  x<-x[f]
  x_mix<-densityMclust(x)
  x_mix_summary<-summary(x_mix)
  #print(x_mix)
  print(x_mix_summary)
  
  png(ofile,width=800,height=600)
  old_par<-par(no.readonly=T)
  #pdf(ofile,width=8,height=6)
  layout(matrix(c(1,1,1,1,1,1,2,3,4), 3, 3, byrow = TRUE))
  a_par<-par(cex.axis=2,cex.lab=2,cex.main=1.8,mar=c(5,6,4,2)+0.1)
  plot(x_mix,what="density",data=x,breaks=50,col="darkgreen",lwd=2,main="",...)
  abline(v=x_mix_summary$mean,lty=2)
  
  for(i in 1:x_mix_summary$G)
  {
    i_mean<-x_mix_summary$mean[i]
    i_sd<-sqrt(x_mix_summary$variance[i])
    i_pro<-x_mix_summary$pro[i]
    #i_sd<-RC_mix_summary$variance[i]
    d<-qnorm(c(0.0013,0.9987),i_mean,i_sd)
    e<-i_pro*dnorm(i_mean,i_mean,i_sd)
    lines(seq(d[1],d[2],by=0.01),i_pro*dnorm(seq(d[1],d[2],by=0.01),i_mean,i_sd),col="orange",lwd=2)
    #rect(d[1],0,d[2],e+0.02,col=rgb(1,0,0,0.2),border=NA)
  }
  plot(x_mix,data=x,breaks=20,col="darkgreen",lwd=2,what="BIC")
  densityMclust.diagnostic(x_mix,type = "cdf",cex.lab=1.5)
  densityMclust.diagnostic(x_mix,type = "qq")
  dev.off()
  
  #par(old_par)
    
  x_mix
}

x_mix<-fitMixModel(log2(inTable$rpkm), paste0(prefix,".RPKM.png"), xlab="log2(RPKM)")
x_mix<-fitMixModel(log2(inTable$TPM), paste0(prefix,".TPM.png"), xlab="log2(TPM)")




