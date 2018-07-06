#!/usr/bin/env Rscript

### mclust 5.1 work ###

suppressPackageStartupMessages(library("argparse"))
parser <- ArgumentParser()
parser$add_argument("-i", "--inFile", type="character", required=TRUE, help="infile")
parser$add_argument("-o", "--outprefix", type="character",required=TRUE, help="output prefix")
parser$add_argument("-G", "--GMax", type="integer",default=9, help="G.max [default %(default)s]")
parser$add_argument("-s", "--sample", type="character",default="SAMPLE", help="sample id [default %(default)s]")
parser$add_argument("-l", "--geneList", type="character",default="", help="what genes whant to analysis [default %(default)s]")
parser$add_argument("-p", "--pattern", type="character",help="only consider cells with names match this pattern [default %(default)s]")
parser$add_argument("-k", "--keepERCC", action="store_true", default=FALSE, help="wheter keep ERCC [default %(default)s]")
args <- parser$parse_args()
print(args)

source("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")

ifile <- args$inFile
prefix <- args$outprefix
glist.file <- args$geneList
keepERCC <- args$keepERCC
g.max <- args$GMax
sample.id <- args$sample
args.pattern <- args$pattern
####source("/WPS1/zhenglt/work/TCR <- chunhong/integrated.20150707/insilicoFACS/test.cofig.r")

inTable<-read.table(ifile,header=T,row.names = 1,check.names = F,stringsAsFactors = F)
gname <- inTable[,1]
names(gname) <- rownames(inTable)
inTable <- inTable[,-1]
### filter out ERCC ?
if(!keepERCC)
{
    f <- !grepl("^ERCC-",rownames(inTable),perl=T)
    inTable<-inTable[f,]
}
if(!is.null(args.pattern)){
    f <- grepl(pattern = args.pattern,x = colnames(inTable),perl = T)
    inTable <- inTable[,f]
}
### filter out genes have all 0s
f <- apply(inTable,1,function(x){ sum(x>0)>=3})
inTable <- inTable[f,]

inTable <- t(inTable)

fitMixModel <- function(x,ofile,G.max=9,uncertainty.threshold=0.1,my.main="",...)
{
  require(mclust)
  x.binary <- structure(rep(NA,length(x)),.Names=names(x))
  f<-is.finite(x)
  x.binary[!f] <- 0
  x<-x[f]
  x_mix<-densityMclust(x,G=1:G.max)
  x_mix_summary<-summary(x_mix)
  #print(x_mix)
  print(x_mix_summary)
  n.positive <- names(x)[x_mix$classification > 1 & x_mix$uncertainty < uncertainty.threshold]
  n.negative <- names(x)[x_mix$classification == 1 & x_mix$uncertainty < uncertainty.threshold]
  n.na <- names(x)[x_mix$uncertainty >= uncertainty.threshold]
  x.binary[n.positive] <- 1
  x.binary[n.negative] <- 0
  
  ### plot for diagnosis
  pdf(ofile,width=12,height=8)
  old_par<-par(no.readonly=T)
  layout(matrix(c(1,1,1,1,1,1,2,3,4), 3, 3, byrow = TRUE))
  a_par<-par(cex.axis=2,cex.lab=2,cex.main=1.8,mar=c(5,6,4,2)+0.1)
  plot(x_mix,what="density",data=x,breaks=50,col="darkgreen",lwd=2,main="",...)
  abline(v=x_mix_summary$mean,lty=2)
  
  for(i in 1:x_mix_summary$G)
  {
    i_mean<-x_mix_summary$mean[i]
    i_sd<-sqrt(x_mix_summary$variance[i])
    i_pro<-x_mix_summary$pro[i]
    ##d<-qnorm(c(0.0013,0.9987),i_mean,i_sd)
    d<-qnorm(c(0.005,0.995),i_mean,i_sd)
    e<-i_pro*dnorm(i_mean,i_mean,i_sd)
    lines(seq(d[1],d[2],by=0.01),i_pro*dnorm(seq(d[1],d[2],by=0.01),i_mean,i_sd),col="orange",lwd=2)
    rect(d[1],0,d[2],e+0.02,col=NA,border=rgb(0.9,0,0,0.9),lty=2,lwd=2.0)
    mtext(text = my.main,cex = 2,line = 0.5)
    ##rect(d[1],0,d[2],e+0.02,col=rgb(1,0,0,0.2),border=NA)
  }
  plot(x_mix,data=x,breaks=20,col="darkgreen",lwd=2,what="BIC")
  densityMclust.diagnostic(x_mix,type = "cdf",cex.lab=1.5)
  tryCatch( densityMclust.diagnostic(x_mix,type = "qq") , error = function(e) e)
  par(old_par)
  dev.off()
    
  ##x_mix
  return(x.binary)
}

if(file.exists(glist.file)) { 
    goi <- read.table(glist.file,header = F,check.names = F,stringsAsFactors = F,colClasses = c("character"))$V1
    goi <- intersect(goi,colnames(inTable))
}else{
    goi <- colnames(inTable)
}
print(head(goi))
out.exp <- t(sapply(goi,function(x){ 
                     cat(sprintf("%s ...\n",gname[x]))
                     fitMixModel(log2(inTable[,x]),sprintf("%s.%s.pdf",prefix,gname[x]), 
                                 xlab="log2(Expression)",my.main = sprintf("%s (%s)",gname[x],sample.id),G.max = g.max)
                    }))
out.exp.df <- data.frame(geneID=rownames(out.exp),geneName=entrezToXXX(rownames(out.exp)))
out.exp.df <- cbind(out.exp.df,out.exp)
write.table(out.exp.df,sprintf("%s.binary.txt",prefix),sep="\t",quote = F,row.names = F,col.names = T)

