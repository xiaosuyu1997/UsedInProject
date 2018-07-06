#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-i", "--inputFile", type="character", required=TRUE, help="input file")
parser$add_argument("-s", "--sample", type="character", default="SAMPLE", help="sample id  [default %(default)s]")
parser$add_argument("-o", "--outPrefix", type="character", required=TRUE, help="output prefix")
parser$add_argument("-d", "--designFile", type="character", required=TRUE, help="design file")
parser$add_argument("-c", "--cellTypeColorFile", type="character", default="/WPS1/zhenglt/work/TCR_chunhong/data/CellType.color", help="cell type color file")
args <- parser$parse_args()

print(args)
sample.id <- args$sample
out.prefix <- args$outPrefix
in.file <- args$inputFile
designFile <- args$designFile
cellType.color.file <- args$cellTypeColorFile

#sample.id <- "P0205"
#out.prefix <- "P0205.test.alpha1.beta1"
#in.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/tracer/OUT/P0205/filtered_TCR_summary/P0205.summary.cell.reassigneClonotype.methodChunhong.r.txt"
#designFile <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/OUT.scLVM/P0205/sfEndo/P0205.designUsed.txt"
#cellType.color.file <- "/WPS1/zhenglt/work/TCR_chunhong/data/CellType.color"

source("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))
library("plyr")

#myDesign<-read.table(designFile,header=T,row.names="sample",check.names=F,colClasses=c("factor","character","factor","factor"))
myDesign<-read.table(designFile,header=T,row.names="sample",check.names=F)
if(!is.null(cellType.color.file) && file.exists(cellType.color.file)){
    sampleTypeColor <- read.SampleTypeColor(cellType.color.file)
}else{
    nTypes <- length(unique(myDesign$sampleType))
    sampleTypeColor <- structure(colorRampPalette(brewer.pal(nTypes,"Paired"))(nTypes),
                              names=unique(myDesign$sampleType))
}

in.table <- read.delim(in.file,header = T,check.names = F,sep = "\t",stringsAsFactors = F)
in.table$sampleType <- myDesign[in.table$Cell_Name,"sampleType"]
##in.table$sampleType <- factor(myDesign[in.table$Cell_Name,"sampleType"],levels=c("PTC","TTC","PTH","TTH","PTR","TTR","NTC","NTH","NTR"),ordered = T)
in.table$sampleColor <- sampleTypeColor[in.table$sampleType]

in.table[,"ID_Alpha1"] <- in.table[,"Identifier(Alpha1)"]
in.table[,"ID_Beta1"] <- in.table[,"Identifier(Beta1)"]
in.table[,"ID_Alpha1Beta1"] <- sprintf("%s/\n%s",in.table[,"ID_Alpha1"],in.table[,"ID_Beta1"])

t.IDAlpha1 <- table(in.table[,"ID_Alpha1"])
t.IDBeta1 <- table(in.table[,"ID_Beta1"])
##in.table[,"Alpha1.occurence"] <- t.IDAlpha1[in.table[,"ID_Alpha1"]]
##in.table[,"Beta1.occurence"] <- t.IDBeta1[in.table[,"ID_Beta1"]]
in.table[,"Alpha1.occurence"] <- structure(as.integer(t.IDAlpha1[in.table[,"ID_Alpha1"]]), .Name=names(t.IDAlpha1[in.table[,"ID_Alpha1"]]))
in.table[,"Beta1.occurence"] <- structure(as.integer(t.IDBeta1[in.table[,"ID_Beta1"]]), .Name=names(t.IDBeta1[in.table[,"ID_Beta1"]]))
in.table[,"hold.var"] <- 1

id.color.mapping <- in.table[,c("ID_Alpha1Beta1","sampleType","sampleColor")]
id.color.mapping <- ddply(id.color.mapping, .(ID_Alpha1Beta1,sampleType,sampleColor), nrow)
id.color.mapping <- id.color.mapping[with(id.color.mapping,order(ID_Alpha1Beta1,-V1)),]
id.color.mapping <- with(id.color.mapping, id.color.mapping[!duplicated(ID_Alpha1Beta1),])
rownames(id.color.mapping) <- id.color.mapping$ID_Alpha1Beta1

dat.plot <- unique(in.table[,c("ID_Alpha1Beta1","Alpha1.occurence","Beta1.occurence")])
t.IDAlpha1Beta1 <- table(in.table[,"ID_Alpha1Beta1"])
dat.plot$combine.occurence <- structure(as.integer(t.IDAlpha1Beta1[dat.plot$ID_Alpha1Beta1]), .Name=names(t.IDAlpha1Beta1[dat.plot$ID_Alpha1Beta1]))
dat.plot$note <- with(dat.plot,sprintf("%s(%d)",ID_Alpha1Beta1,combine.occurence))
dat.plot.highlight <- with(dat.plot, dat.plot[combine.occurence>3,])
dat.plot.highlight$sampleType <- id.color.mapping[dat.plot.highlight$ID_Alpha1Beta1,"sampleType"]

pdf(sprintf("%s.dist.pdf",out.prefix),width = 8,height = 6)

gout <- ggplot(in.table, aes(Alpha1.occurence, Beta1.occurence)) + 
    geom_count(aes(size = ..prop.., group = hold.var)) + scale_size_area(max_size = 12) +
    labs(list(title = sprintf("%s",sample.id), x = "Alpha1 Occurence", y = "Beta1 Occurence")) +
    geom_label_repel(aes(label=note,x=Alpha1.occurence,y=Beta1.occurence,fill=sampleType),
               data=dat.plot.highlight,
               arrow = arrow(length = unit(0.01, 'npc')), force = 5, segment.color = '#555555',
               colour = "white",show.legend = NA, nudge_y = 0.0,size=1.5)
               ###colour = "white",show.legend = NA, nudge_y = 0.5,size=1.5) + scale_fill_manual(values = sampleTypeColor)
print(gout)
gout <- ggplot(in.table, aes(Alpha1.occurence, Beta1.occurence)) + 
    geom_count(aes(size = ..prop.., group = hold.var)) + scale_size_area(max_size = 12) +
    labs(list(title = sprintf("%s",sample.id), x = "Alpha1 Occurence", y = "Beta1 Occurence")) +
    geom_label_repel(aes(label=combine.occurence,x=Alpha1.occurence,y=Beta1.occurence,fill=sampleType),
               data=dat.plot.highlight,
               arrow = arrow(length = unit(0.01, 'npc')), force = 1, segment.color = '#555555',
               colour = "white",show.legend = NA, nudge_y = 0.0,size=1.5)
print(gout)

#in.table$Alpha1.occurence <- in.table$Alpha1.occurence + runif(nrow(in.table),-0.5,0.5)
#in.table$Beta1.occurence <- in.table$Beta1.occurence + runif(nrow(in.table),-0.5,0.5)
    ##geom_count(aes(size = ..prop.., group = hold.var)) + scale_size_area(max_size = 12) +
gout <- ggplot(in.table, aes(Alpha1.occurence, Beta1.occurence)) + 
    geom_point(size=0.2) + geom_jitter(size=0.2,width=0.5*max(c(in.table$Alpha1.occurence,in.table$Beta1.occurence))/16,height=0.5*max(c(in.table$Alpha1.occurence,in.table$Beta1.occurence))/16) +
    labs(list(title = sprintf("%s",sample.id), x = "Alpha1 Occurence", y = "Beta1 Occurence")) +
    geom_label_repel(aes(label=combine.occurence,x=Alpha1.occurence,y=Beta1.occurence,fill=sampleType),
               data=dat.plot.highlight,
               arrow = arrow(length = unit(0.01, 'npc')), force = 1, segment.color = '#555555',
               colour = "white",show.legend = NA, nudge_y = 0.0,size=1.5)
print(gout)


dev.off()

