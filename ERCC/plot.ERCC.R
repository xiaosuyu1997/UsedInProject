#!/usr/bin/env Rscript

args <- commandArgs(T)
if(length(args)<4)
{
    cat(sprintf("plot.ERCC.R <exp.file> <conc.file> <dilution.file> <out.prefix>\n"))
    q()
}
exp.file <- args[1]
conc.file <- args[2]
d.file <- args[3]
out.prefix <- args[4]
#exp.file <- "./OUT/P0617/P0617.ERCC.TPM.tab.gz"
#conc.file <-  "/DBS/DB_temp/zhangLab/genome/ERCC/cms_095046.txt"
#d.file <- "../../ERCC.dilution.purification.txt"
#out.prefix <- "./OUT/P0617/P0617.ERCC.QC"
EPSILON <- 1e-6

exp.table <- read.table(exp.file,header = T,row.names = 1,check.names = F,stringsAsFactors = F)
exp.table <- exp.table[,-1,drop=F]
#print(head(exp.table[,1:8]))
conc.table <- read.table(conc.file,header = T,sep = "\t",row.names = "ERCC ID",check.names = F)
print(head(conc.table))
##d.table <- read.table(d.file,header = T,sep = "\t",row.names=1)
d.table <- read.table(d.file,header = T,sep = "\t",row.names="sample")
d.table$AfterDilution <- as.numeric(d.table$AfterDilution)
print(head(d.table))

pdf(paste0(out.prefix,".pdf"),width = 8,height = 8)
par(mar=c(5,6,4,2),cex.lab=1.5)
out.table <- sapply(seq_len(ncol(exp.table)), function(i,A){
      dilution <- d.table[colnames(A)[i],"AfterDilution"]
      x <<- log10(conc.table[rownames(A),3]*dilution)
      y <<- log10(A[,i]+EPSILON)
      dat.test <<- data.frame(x=x[y>-6],y=y[y>-6])
      if(nrow(dat.test)<3) {
        cat(sprintf("Warning: ERCC of sample %s expressed lowly!\n",colnames(A)[i]))
        c(colnames(A)[i],dilution,"NA","NA",sum(y>-6),length(y))
      }else {
        lm.out <<- lm(y~x,dat.test)
        cor.test.out <<- cor.test(dat.test$x,dat.test$y)
        #plot(x,y,xlab=expression(log[10](Attomoles + epsilon)),ylab=expression(log[10](TPM + epsilon)),main=colnames(A)[i],ylim=c(-2,5),xlim=c(-3,2))
        plot(x,y,xlab=expression(log[10](Attomoles)),ylab=expression(log[10](TPM + epsilon)),main=colnames(A)[i])
        abline(h=0,lty=2,lwd=1)
        abline(v=0,lty=2,lwd=1)
        abline(lm.out,lty=1,lwd=1)
        #text(x = -6,y = -2,labels = paste0("p value: ",round(cor.test.out$p.value,4)),pos = 4)
        #text(x = -6,y = -2.3,labels = paste0("pearson: ",round(cor.test.out$estimate,4)),pos = 4)
        #text(x = -6,y = -2.6,labels = paste0("detection: ",sum(y>-6),"/",length(y)),pos = 4)
        legend("topleft",legend = c(paste0("p value: ",round(cor.test.out$p.value,4)),paste0("pearson: ",round(cor.test.out$estimate,4)),paste0("detection: ",sum(y>-6),"/",length(y))))
        c(colnames(A)[i],dilution,cor.test.out$p.value,cor.test.out$estimate,sum(y>-6),length(y))
      }
},exp.table)
dev.off()
out.table <- t(out.table)
colnames(out.table) <- c("sample","dilution","p.value","pearson","detected","total")
write.table(out.table,file = paste0(out.prefix,".out.txt"),quote = F,row.names = F,col.names = T,sep = "\t")
