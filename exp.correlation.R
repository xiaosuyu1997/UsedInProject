#!/usr/bin/env Rscript


suppressPackageStartupMessages(library("gplots"))
source("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")

args <- commandArgs(T)
if(length(args)<5)
{
    cat("usage: exp.correlation.R <in.file> <sample.id> <out.dir> <design.file> <cellType.color.file>\n")
    q()
}

in.file <- args[1]
sample.id <- args[2]
out.dir <- args[3]
design.file <- args[4]
cellType.color.file <- args[5]

#in.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/gsnap_out/quantification/P1202.TPM.tab.gz"
#sample.id <- "P1202"
#out.dir <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/correlation/OUT/P1202"
#design.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/clustering/marker/P1202.filter.by.marker.design.txt"
#cellType.color.file  <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/clustering/CellType.color"

dir.create(out.dir,showWarnings = F,recursive = T)

my.design <- read.table(design.file,header = T,stringsAsFactors = F)
rownames(my.design) <- my.design$sample
in.table <- read.table(in.file,header=T,check.names = F,stringsAsFactors = F,row.names=1)

in.table <- in.table[,my.design$sample]
dim(in.table)
in.table[1:4,1:8]

my.func<-function (x, y, col ='blue', bg = NA, pch = par("pch"), cex = 0.03, col.line = "blue", lty = par("lty"),pty=par("pty"),...) 
{
	upar <- par(usr = c(0, max(log(x)+3), 0, max(log(y)+3) ))
	on.exit(par(upar))
	points(log(x), log(y), pch = ".", col = col, bg = bg, cex = cex)
	#smoothScatter(x, y, add = TRUE)
	#abline(lm(log(y)~log(x),data=data.frame(x=x,y=y)))
}


plotPairs2 <- function(dat,png.file,width=800,height=800,meth="pearson",xlim=c(0,15),ylim=c(0,15), TPM.CUTOFF=2)
{
	rowExpPercentage <- apply(dat,1,function(x){ sum(x > TPM.CUTOFF)/length(x) } )
	select <- rowExpPercentage > 0.5
	dat.plot <-dat[select,]
	head(dat.plot)
	dim(dat.plot)
    
    m=ncol(dat.plot)
	#png(png.file,width,height)
	pdf(png.file,width,height)
	pairs2(dat.plot,cex.labels=0.7,cex=1.1,gap=0.8,lower.panel=panel.cor,upper.panel=my.func,cex.axis=1.5,xlim=c(0,15),ylim=c(0,15),meth=meth)
	dev.off()
}

plotCorMatrix <- function(dat,out.prefix,designM,colSet=NULL,TPM.CUTOFF=2)
{
    if(is.null(colSet))
    {
        nLevel <- length(levels(designM$sampleType))
        if(nLevel<=9)
        {
            suppressPackageStartupMessages(require("RColorBrewer"))
            colSet <- brewer.pal(nLevel,"Set1")
            names(colSet) <- levels(designM$sampleType)
        }else
        {
            colSet <- rainbow(nLevel)
            names(colSet) <- levels(designM$sampleType)
        }
    }else
    {
      colSet <- colSet[names(colSet) %in% unique(designM$sampleType)]
    }

    selectGene <- function(X)
    {
        rowExpPercentage <- apply(X,1,function(v){ sum(v > TPM.CUTOFF)/length(v) } )
        select <- rowExpPercentage > 0.5
        select
    }
    designM <- designM[colnames(dat),]
	dat.plot <-dat[selectGene(dat),]
	#print(head(dat.plot))
	print(dim(dat.plot))
    m=ncol(dat.plot)
    ### correlation with all other samples
    dat.plot.cor.mat<<-cor(dat.plot,method = "pearson")
    out.df <- data.frame(sample=rownames(dat.plot.cor.mat),stringsAsFactors = F)
    out.df <- cbind(out.df,t(sapply(seq(nrow(dat.plot.cor.mat)),function(i,X){ c(mean(X[i,-i]),length(X[i,-i])) },dat.plot.cor.mat)))
    colnames(out.df) <- c("sample","mean.cor","N")
    rownames(out.df) <- out.df$sample
    ### correlation with all other samples with the same sampleType
    cor.list <<- list()
    cor.mean.stype <- c()
    for(stype in unique(designM$sampleType))
    {
        f <- designM[colnames(dat),"sampleType"] %in% stype
        dat.block <- dat[,f,drop=F]
        dat.block <- dat.block[selectGene(dat.block),,drop=F]
        print(stype)
	    print(dim(dat.block))
        cor.list[[stype]] <<- cor(dat.block,method = "pearson")
        o.df <- data.frame(sample=rownames(cor.list[[stype]]),stringsAsFactors = F)
        o.df <- cbind(o.df,t(sapply(seq(nrow(cor.list[[stype]])),function(i,X){ c(mean(X[i,-i]),length(X[i,-i])) },cor.list[[stype]])))
        colnames(o.df) <- c("sample","mean.cor.stype","N.stype")
        rownames(o.df) <- o.df$sample
        cor.mean.stype <- rbind(cor.mean.stype,o.df)
    }
    out.df <- cbind(out.df,cor.mean.stype[out.df$sample,c("mean.cor.stype","N.stype")])
    write.table(out.df,sprintf("%s.summary.txt",out.prefix),sep = "\t",row.names = F,quote = F)
    ### plot
    patientcolors <- colSet[as.character(designM$sampleType)]
    nn <- nrow(dat.plot.cor.mat)
    mm <- ncol(dat.plot.cor.mat)
    pdf(sprintf("%s.cor.matrix.pdf",out.prefix),width = 12,height = 10)
    heatmap.2(dat.plot.cor.mat,col=bluered(10), ColSideColors=patientcolors, RowSideColors=patientcolors, Rowv = F, Colv = F, scale="none", density.info="none", dendrogram="none", keysize=1.2, trace="none", margin=c(15, 20), main="",cexRow=min(1.8,8/nn,cexCol=min(1.8,8/mm)))
    legend("topright",legend=names(colSet),fill=colSet,border=colSet,cex=1.5)
    sapply(names(cor.list),function(stype,X){
                                            if(nrow(X[[stype]])<3 || ncol(X[[stype]])<3) { loginfo(sprintf("Too few samples/genes: %s",stype));return(NULL); } 
                                            heatmap.2(X[[stype]],col=bluered(10), 
                                                      ColSideColors=colSet[as.character(designM[colnames(X[[stype]]),"sampleType"])], 
                                                      RowSideColors=colSet[as.character(designM[rownames(X[[stype]]),"sampleType"])], 
                                                      Rowv = F, Colv = F, 
                                                      scale="none", density.info="none", 
                                                      dendrogram="none", keysize=1.2, trace="none", 
                                                      margin=c(15, 20), main="",
                                                      cexRow=min(1.8,8/nn,cexCol=min(1.8,8/mm)))
                                            legend("topright",legend=names(colSet),fill=colSet,border=colSet,cex=1.5)
                                        },cor.list)
    dev.off()
    out.df
}

sampleTypeColor <- read.SampleTypeColor(cellType.color.file)
out.df <- plotCorMatrix(in.table[,],sprintf("%s/exp.correlation.pearson",out.dir),my.design,colSet = sampleTypeColor)
out.df <- out.df[order(out.df$mean.cor,decreasing = F),]
plotPairs2(in.table[,unique(c(head(rownames(out.df),n=10),tail(rownames(out.df),n=10)))],sprintf("%s/exp.correlation.pearson.pdf",out.dir),width=15,height=15)
#plotPairs2(in.table[,unique(c(head(rownames(out.df),n=10),tail(rownames(out.df),n=10)))],sprintf("%s/exp.correlation.pearson.png",out.dir),width=1000,height=1000)
