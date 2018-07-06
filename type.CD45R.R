#!/usr/bin/env Rscript

args <- commandArgs(T)
if(length(args)<3)
{
    cat("type.CD45R.R <exp.file> <l.file> <out.prefix>\n")
    q()
}
exp.file <- args[1]
l.file <- args[2]
out.prefix <- args[3]

#exp.file <- "P0205.CD45R.tab"
#l.file <- "marker.naive-mature.exon.list"
#out.prefix <- "P0205.CD45R"

source("/Share/BP/zhenglt/02.pipeline/cancer/lib/myFunc.R")

suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("RColorBrewer"))

exp.table <- read.table(exp.file,header = T,row.names = 1,check.names = F,stringsAsFactors = F)
exp.table <- exp.table[,-1]
print(exp.table[1:5,1:8])
l.table <- read.table(l.file,header = F,row.names = 1,check.names = F,stringsAsFactors = F)
names(l.table) <- c("Ann")
print(l.table[1:5,])

f.table <- exp.table[rownames(l.table),]
print(f.table[1:5,1:8])

CD45R.status <- apply(f.table,2,function(v){
      region.E4 <- rownames(l.table)[which(l.table$Ann=="PTPRC::Exon4")]
      region.E5 <- rownames(l.table)[which(l.table$Ann=="PTPRC::Exon5")]
      region.E6 <- rownames(l.table)[which(l.table$Ann=="PTPRC::Exon6")]
      if(v[region.E4] > 0 && v[region.E5] > 0 && v[region.E6] > 0) { "CD45R_ABC" }
      else if(v[region.E4] > 0 && v[region.E5] > 0 && v[region.E6] == 0) { "CD45R_AB" }
      else if(v[region.E4] > 0 && v[region.E5] == 0 && v[region.E6] > 0) { "CD45R_AC" }
      else if(v[region.E4] == 0 && v[region.E5] > 0 && v[region.E6] > 0) { "CD45R_BC" }
      else if(v[region.E4] > 0 && v[region.E5] == 0 && v[region.E6] == 0) { "CD45R_A" }
      else if(v[region.E4] == 0 && v[region.E5] > 0 && v[region.E6] == 0) { "CD45R_B" }
      else if(v[region.E4] == 0 && v[region.E5] == 0 && v[region.E6] > 0) { "CD45R_C" }
      else if(v[region.E4] == 0 && v[region.E5] == 0 && v[region.E6] == 0) { "CD45R_O" }
      else { "Unknown" }
})
out <- data.frame(sample=names(CD45R.status),CD45R.status=CD45R.status)
write.table(out,paste0(out.prefix,".status.txt"),row.names = F,sep="\t",quote = F)
out.slim <- out
out.slim$CD45R.status <- -1
out.slim$CD45R.status[out$CD45R.status=="CD45R_A"] <- 0
out.slim$CD45R.status[out$CD45R.status=="CD45R_O"] <- 1
write.table(out.slim,paste0(out.prefix,".status.slim.txt"),row.names = F,sep="\t",quote = F)

cell.type <- as.factor(str_match(rownames(out),"^(...)")[,2])
names(cell.type) <- rownames(out)
table(out$CD45R.status)
table(cell.type)

status.summary <- sapply(levels(cell.type),function(x){
       table(out[cell.type==x,"CD45R.status"])
})

stat.type <- levels(out$CD45R.status)
colSet <- brewer.pal(length(stat.type),"Blues")
names(colSet) <- stat.type

dat.plot  <- apply(status.summary,2,function(x){ x/sum(x) } )
pdf(paste0(out.prefix,".status.summary.pdf"),width=12,height=8)
par(mar=c(5,6,4,8),cex.lab=2)
barplot(dat.plot,beside = F,xlab="Cell Type",ylab="CD45R Isoform",main="",col=colSet)
legend("right",legend=stat.type,fill=colSet,inset = -0.15,xpd = T,horiz = F)
dev.off()

