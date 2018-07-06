#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-i", "--inFile", type="character", required=TRUE, help="infile")
parser$add_argument("-o", "--outPrefix", type="character", required=TRUE, help="outprefix")
parser$add_argument("-m", "--main", default="TRB", type="character", help="main title")
args <- parser$parse_args()
print(args)

out.prefix <- args$outPrefix
in.file <- args$inFile
g.main <- args$main

#out.prefix <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/tracer/repertoire_characterization/S3.clonotype.TB"
#in.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/tracer/repertoire_characterization/S3.clonotype.TB.list"
#main <- "TRB"

in.table <- read.table(in.file,header = F,check.names = F,stringsAsFactors = F)
colnames(in.table) <- c("V","J","ID","CDR3","N")
in.table <- in.table[in.table$V!="-" & in.table$J!="-",]

my.plotBar <- function(dat.plot,out.prefix,pdf.width=10,pdf.height=8,nTop=0,...){
    require("plotrix")
    pvalues <- sapply(dat.plot,function(x){ prop.test(x,sum(dat.plot),1/length(dat.plot),alternative="greater")$p.value })
    ###pdf(sprintf("%s",out.prefix),width = pdf.width,height = pdf.height)
    ##par(mar=c(5,6,4,2),cex.lab=1.5,cex.main=1.5)
    dat.plot <- dat.plot/sum(dat.plot)
    eFreq <- 1/length(dat.plot)
    if(nTop>0){
        dat.plot <- tail(dat.plot,nTop)
        pvalues <- tail(pvalues,nTop)
    }
    xx <- barplot(dat.plot,xaxt="s",border=NA,horiz=T,axisnames=F,...)
    abline(v=eFreq,lty=2,lwd=2)
    staxlab(2,at = xx,labels=names(dat.plot),srt=0, cex=0.8,adj=0.5,top.line=0)
    if(length(dat.plot[pvalues<0.05])>0){
        #text(dat.plot[pvalues<0.05],xx[pvalues<0.05]-0.7,"*",offset = 0.5, pos = 4,cex=1.8,xpd=T)
        text(dat.plot[pvalues<0.05],xx[pvalues<0.05]-0.5,"*",adj=c(-0.5,0.5),cex=1.8,xpd=T)
    }
    ##dev.off()
}
pdf(sprintf("%s.GeneUsage.pdf",out.prefix),width = 6,height = if(g.main=="TRB") 12 else 16)
layout(matrix(c(1,2), 2, 1, byrow = TRUE),heights=if(g.main=="TRB") c(2.5,1.5) else c(1,1))
dat.plot <- table(in.table$V)
par(mar=c(2,6,4,2),cex.lab=2.0,cex.axis=1.5,cex.main=2.0)
my.plotBar(sort(dat.plot),sprintf("%s.VGeneUsage.pdf",out.prefix),ylab="V Gene",xlab="",main=g.main,nTop=25)
dat.plot <- table(in.table$J)
par(mar=c(2,6,4,2),cex.lab=2.0,cex.axis=1.5,cex.main=2.0)
my.plotBar(sort(dat.plot),sprintf("%s.JGeneUsage.pdf",out.prefix),ylab="J Gene",xlab="",main="",nTop=25)
dev.off()

### CDR3 
reg <- regexpr("^C(.+?)FG.G$",in.table$CDR3,perl=T)
reg.start <- attr(reg,"capture.start")
reg.len <- attr(reg,"capture.length")
in.table$CDR3.inner <- sapply(seq_along(reg),function(i){
    substring(in.table$CDR3[i],reg.start[i],reg.start[i]+reg.len[i]-1)
})
write.table(in.table,sprintf("%s.out.txt",out.prefix),row.names = F,quote = F,sep = "\t")

pdf(sprintf("%s.CDR3Len.pdf",out.prefix),width = 10,height = 8)
par(mar=c(5,6,4,2),cex.lab=1.5,cex.main=1.5)
dat.plot <- table(nchar(in.table$CDR3.inner))
dat.plot <- dat.plot[names(dat.plot)!="0"]
barplot(dat.plot/sum(dat.plot),xaxt="s",main=sprintf("CDR3 length(%s)",g.main))
#staxlab(2,at = xx,labels=names(dat.plot),srt=0, cex=0.8,adj=0.5,top.line=0)
##hist(nchar(in.table$CDR3.inner),xlab="",main=sprintf("CDR3 lenght(%s)",g.main))
dev.off()



