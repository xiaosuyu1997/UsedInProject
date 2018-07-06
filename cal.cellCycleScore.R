#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-i", "--input", type="character", required=TRUE, help="input file, TPM in .tab or .tab.gz")
parser$add_argument("-d", "--designFile", type="character", required=TRUE, help="design matrix file")
parser$add_argument("-c", "--ccFile", type="character", 
                    default="/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/cellCycle/CellCycle.tableS2.Macosko_et_al_Cell_2015.gID.txt", help="cell cycle file [default %(default)s]")
parser$add_argument("-o", "--out", type="character", required=TRUE, help="output prefix")
parser$add_argument("-y", "--cellTypeColorFile", type="character", 
                    default="/WPS1/zhenglt/work/TCR_chunhong/data/CellType.color", help="cellTypeColorFile [default %(default)s]")
args <- parser$parse_args()
in.file <- args$input
d.file <- args$designFile
cellCycle.file <- args$ccFile
out.prefix <- args$out
cellTypeColorFile <- args$cellTypeColorFile

print(args)

#in.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/gsnap_out/quantification/P0205.TPM.tab.gz"
#d.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/OUT.scLVM/P0205/sfIgnoreERCC/P0205.designUsed.txt"
#cellCycle.file <- "CellCycle.tableS2.Macosko_et_al_Cell_2015.gID.txt"
#out.prefix <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/cellCycle/CellCycleScore.P0205"
#cellTypeColorFile <- "/WPS1/zhenglt/work/TCR_chunhong/data/CellType.color"

### another version of cell cycle regulated genes
##cc.file <- "CellCycle.tableS2.Macosko_et_al_Cell_2015.txt"
##cc.table <- read.table(cc.file,header = T,check.names = F,stringsAsFactors = F,sep = "\t")
##cc.table$geneID <- XXXToEntrez(cc.table$geneSymbol)
##print(head(cc.table))
##write.table(cc.table,sprintf("CellCycle.tableS2.Macosko_et_al_Cell_2015.gID.txt"),sep = "\t",row.names = F,quote = F)


myDesign <- read.table(d.file,header = T,check.names = F,stringsAsFactors = F)
rownames(myDesign) <- myDesign$sample

source("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")
in.table <- read.table(in.file,row.names = 1,header = T,sep = "\t",check.names = F,stringsAsFactors = F)
in.table <- in.table[,myDesign$sample]
in.table <- log2(in.table+1)
dim(in.table)
in.table[1:4,1:8]

all.geneID <- rownames(in.table)
cellCycle.table <- read.table(cellCycle.file,header = T,check.names = F,stringsAsFactors = F,sep = "\t")
cellCycle.table$geneID <- as.character(cellCycle.table$geneID)
#rownames(cellCycle.table) <- cellCycle.table$geneID
sampleTypeColor <- read.SampleTypeColor(cellTypeColorFile)

out.list <- list()
ccExp <- sapply(c("G1/S","S","G2/M","M","M/G1"),function(x){
    gID <- subset(cellCycle.table,phase==x)$geneID
    gID <- intersect(gID,all.geneID)
    phaseExp <- in.table[gID,]
    f <- apply(phaseExp,1,function(x){ sum(x>0)/length(x) >= 0.3 })
    phaseExp <- phaseExp[f,,drop=F]
    phaseMean <- apply(phaseExp,2,mean)
    phaseCor <- apply(phaseExp,1,function(x){ cor(x,phaseMean,method = "pearson",use="pairwise.complete.obs") })
    phaseCore <- phaseCor[phaseCor > 0.2 & !is.na(phaseCor)]
    gIDCoreSet <- names(phaseCore)
    phaseMeanCoreSet <- apply(phaseExp[gIDCoreSet,,drop=F],2,mean)
    out.list[[x]][['exp']] <<- phaseExp
    out.list[[x]][['mean']] <<- phaseMean
    out.list[[x]][['phaseCor']] <<- phaseCor
    out.list[[x]][['phaseCore']] <<- phaseCore
    out.list[[x]][['meanCore']] <<- phaseMeanCoreSet
    print(length(phaseCore))
    out.df <- data.frame(geneID=names(phaseCore))
    out.df$geneSymbol <- entrezToXXX(out.df$geneID)
    colnames(out.df) <- sprintf("%s.%s",colnames(out.df),x)
    write.table(out.df,sprintf("%s.%s.gene.list",out.prefix,gsub("[/ ]","_",x,perl = T)),sep = "\t",row.names = F,quote = F)
    phaseMeanCoreSet
})
out.list[[1]][['exp']][1:4,1:8]
out.list[[1]][['mean']][1:4]
out.list[[1]][['phaseCor']][1:4]
head(ccExp)
ccExp.norm <- scale(ccExp)
head(ccExp.norm)
ccExp.norm <- t(apply(ccExp.norm,1,scale))
colnames(ccExp.norm) <- colnames(ccExp)
head(ccExp.norm)

potential.pattern <- matrix(c(1,0,0,0,0,
                              1,1,0,0,0,
                              0,1,0,0,0,
                              0,0,1,0,0,
                              0,0,1,1,0,
                              0,0,0,1,0,
                              0,0,0,0,1,
                              1,0,0,0,1
                              ),ncol=5,byrow = T)
colnames(potential.pattern) <- c("G1/S","S","G2/M","M","M/G1")

which.pattern <- function(phaseScore,pp){
    pp.cor <- apply(pp,1,function(y){
          cor(phaseScore,y,method = "pearson")
    })
    i <- which.max(pp.cor)
    return(c(i,pp.cor[i],ifelse(i>1,pp.cor[i-1],1),ifelse(i<8,pp.cor[i+1],1)))
}
ccOrder <- as.data.frame(t(apply(ccExp.norm,1,function(x){
      which.pattern(x,potential.pattern)
})))
colnames(ccOrder) <- c("pattern","pattern.cor","preceding.cor","succeeding.cor")
ccOrder <- ccOrder[with(ccOrder, order(pattern,-preceding.cor,succeeding.cor) ),]
print(head(ccOrder))

plot.ccOrder <- function(dat,dat.order,out.prefix,colColor){
    suppressPackageStartupMessages(require("gplots"))
    suppressPackageStartupMessages(require("RColorBrewer"))
	colSet <- rev(brewer.pal(9,"RdBu"))
    pdf(sprintf("%s.ccOrder.pdf",out.prefix),width=10,height = 8)
    heatmap.2(t(dat[rownames(dat.order),]),Rowv=F,Colv=F,dendrogram="none",
              ColSideColors=colColor[rownames(dat.order)],
              cexRow=min(1.8,55/ncol(dat)),cexCol=min(1.8,55/nrow(dat)),
              key=T,keysize=1.2,trace="none",density.info="none",
              col=colorRampPalette(colSet)(100))
    dev.off()
}
plot.ccOrder(ccExp.norm,ccOrder,out.prefix,structure(sampleTypeColor[myDesign$sampleType],names=rownames(myDesign)))

##phaseName <- colnames(ccExp.norm)
samplePhase <- data.frame(sampleID=rownames(ccExp.norm))
samplePhase$phase <- (apply(ccExp.norm,1,function(x){ names(which.max(x)) } ))
samplePhase$sampleType <- myDesign[rownames(ccExp.norm),"sampleType"]
head(samplePhase)
write.table(samplePhase,sprintf("%s.ccPhase.txt",out.prefix),sep = "\t",row.names = F,quote = F)
table(samplePhase$phase)
table(samplePhase$phase,samplePhase$sampleType)

### output
out.df <- data.frame(sampleID=rownames(ccExp.norm))
out.df <- cbind(out.df,ccExp.norm)
write.table(out.df,sprintf("%s.ccScore.txt",out.prefix),sep = "\t",row.names = F,quote = F)




