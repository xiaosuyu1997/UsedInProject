#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-i", "--inputFile", type="character", required=TRUE, help="input file")
parser$add_argument("-c", "--cellTypeColorFile", type="character", required=TRUE, help="cellTypeColorFile")
parser$add_argument("-d", "--designFile", type="character", required=TRUE, help="design matrix file")
parser$add_argument("-s", "--sample", type="character", default="SAMPLE", help="sample id  [default %(default)s]")
parser$add_argument("-o", "--outDir", type="character", required=TRUE, help="output directory")
parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="whether verbose mode [default %(default)s]")
args <- parser$parse_args()

print(args)

suppressPackageStartupMessages(library(reshape2))

tracer.summary.file <- args$inputFile
cellType.color.file <- args$cellTypeColorFile
design.file <- args$designFile
sample.id <- args$sample
out.dir <- args$outDir

#tracer.summary.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/tracer/P1118.summary"
#cellType.color.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/clustering/CellType.color"
#design.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/P1118.P2/sample.design/P1118.design.ERCCFlt.txt"
#sample.id <- "P1118"
#out.dir <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/tracer/OUT.landscape"

dir.create(out.dir,showWarnings = F,recursive = T)

plot.tcr.landscape <- function(dat.plot,out.prefix,colSet,sampleType,sampleTypeOrder=c("TTC","NTC","PTC"),main="",byType=T)
{
    suppressPackageStartupMessages(require("ComplexHeatmap"))
    suppressPackageStartupMessages(require("circlize"))
    suppressPackageStartupMessages(require("gridBase"))
    suppressPackageStartupMessages(require("dendextend"))
    suppressPackageStartupMessages(require("gplots"))
    sampleType <- as.character(sampleType)
    #print(str(sampleType))
    #print(str(colSet))
    #print(str(dat.plot))
    names(sampleType) <- colnames(dat.plot)
    ### add code to check dat.plot
    n <- nrow(dat.plot)
    m <- ncol(dat.plot)
    if(n<3||m<3){ cat(sprintf("Too few genes/samples: #genes: %s, #samples: %s\n",n,m)); return(NULL); }

    ### reorder according sampleTypeOrder
    dat.plot.reorder <- c()
    for(iType in sampleTypeOrder) {
        dat.block <- dat.plot[,sampleType[colnames(dat.plot)]==iType]
        dat.plot.reorder <- cbind(dat.plot.reorder,dat.block)
    }
    dat.plot <- dat.plot.reorder
    sampleType <- sampleType[colnames(dat.plot)]
    n <- nrow(dat.plot)
    m <- ncol(dat.plot)
    print(dim(dat.plot))

    dat.plot[dat.plot==0] <- 0
    dat.plot[dat.plot>0] <- 1

    pdf(sprintf("%s.pdf",out.prefix),width=22,height=12)
    par(mar=c(5,4,4,20),cex.main=1.5)
    
    ### plot 0
    plot.new()
    legend("topright",legend=names(colSet)[names(colSet) %in% sampleType],
           fill=colSet[names(colSet) %in% sampleType],
           border=colSet[names(colSet) %in% sampleType],
           cex=1.5,inset=c(-0.08,-0.06),xpd=T)
    title(main)
    vps <- baseViewports()
    pushViewport(vps$inner, vps$figure, vps$plot)

    annDF <- data.frame(sampleType=sampleType)
    annColList <- list(sampleType=colSet)
    ha.col <- HeatmapAnnotation(df = annDF, col = annColList, show_legend = FALSE)
    top_annotation_height <- unit(1.5*ncol(annDF), "cm")

    if(byType)
    {
        oncoprint_row_order <- function(){
            order(rowSums(dat.plot), decreasing = T)
        }
        rrorder <- oncoprint_row_order()
        print(rrorder)

        oncoprint_column_order <- function(){
            offset <- 0
            corder <- c()
            for(iType in sampleTypeOrder) {
                dat.block <- dat.plot[,sampleType[colnames(dat.plot)]==iType]
                #print(dim(dat.block))
                scoreCol = function(x) { 
                    score = 0
                    for (i in 1:length(x)) { 
                        if (x[i]) { 
                            score = score + 2^(length(x) - i * 1/x[i])
                        }
                    }
                    return(score)
                }
                scores = apply(dat.block[rrorder, ], 2, scoreCol)
                corder <- append(corder,c(order(scores, decreasing = TRUE) + offset))
                print(corder)
                offset  <- offset + ncol(dat.block)
            }
            return(corder)
        }
        ccorder <- oncoprint_column_order()

        ht <- oncoPrint(list(TCR=dat.plot),
                        show_heatmap_legend = F,
                        show_row_names=T,show_column_names=T,
                        row_names_gp = gpar(fontsize = 3.2*55/max(nrow(dat.plot),25)),
                        column_names_gp = gpar(fontsize = 3.2*160/max(ncol(dat.plot),25)),
                        pct_gp = gpar(fontsize = 3.2*120/max(nrow(dat.plot),25)),
                        show_column_barplot=F,
                        show_row_barplot=F,
                        top_annotation = ha.col,top_annotation_height = top_annotation_height,
                        row_order=rrorder,
                        column_order = ccorder)
        draw(ht, newpage= FALSE)
    }else {
        ht <- oncoPrint(list(TCR=dat.plot),
                        show_heatmap_legend = F,
                        show_row_names=T,show_column_names=T,
                        row_names_gp = gpar(fontsize = 3.2*55/max(nrow(dat.plot),25)),
                        column_names_gp = gpar(fontsize = 3.2*160/max(ncol(dat.plot),25)),
                        pct_gp = gpar(fontsize = 3.2*120/max(nrow(dat.plot),25)),
                        show_column_barplot=F,
                        show_row_barplot=F,
                        top_annotation = ha.col,top_annotation_height = top_annotation_height)
        draw(ht, newpage= FALSE)

    }
    dev.off()
}

source("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")

### ------ tmp file ------

tmp.design.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/P1118.P2/sample.design/P0508.design.ERCCFlt.txt"
tmp.sample.id <- "P0508"
tmp.out.dir <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/tracer/OUT.landscape"
tmp.tracer.summary.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/tracer/OUT.landscape/P0508.recurrent.txt"
tmp.cellType.color.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/clustering/CellType.color"
tmp.sampleTypeColor <- read.SampleTypeColor(tmp.cellType.color.file)
tmp.myDesign<-read.table(tmp.design.file,header=T,row.names="sample",check.names=F,colClasses=c("factor","character","factor","factor"))
tmp.table <- read.table(tmp.tracer.summary.file,header = T,check.names = F,stringsAsFactors = F,sep = "\t")
tmp.table$TCR <- paste0(tmp.table$TRA,"/",tmp.table$TRB)
tmp.table$P <- 1
tmp.dat.plot <- dcast(data = tmp.table,TCR~sampleID,value.var="P")
rownames(tmp.dat.plot) <- tmp.dat.plot$TCR
tmp.dat.plot <- as.matrix(tmp.dat.plot[,-1])
tmp.dat.plot[is.na(tmp.dat.plot)] <- 0
tmp.dat.plot[1:4,1:8]
plot.tcr.landscape(tmp.dat.plot,sprintf("%s/%s.tmp.byTypeT",tmp.out.dir,tmp.sample.id),tmp.sampleTypeColor,tmp.myDesign[colnames(tmp.dat.plot),"sampleType"],main=tmp.sample.id)
plot.tcr.landscape(tmp.dat.plot,sprintf("%s/%s.tmp.byTypeF",tmp.out.dir,tmp.sample.id),tmp.sampleTypeColor,tmp.myDesign[colnames(tmp.dat.plot),"sampleType"],main=tmp.sample.id,byType=F)
### ------ tmp file ------

sampleTypeColor <- read.SampleTypeColor(cellType.color.file)
myDesign<-read.table(design.file,header=T,row.names="sample",check.names=F,colClasses=c("factor","character","factor","factor"))
in.table <- read.table(tracer.summary.file,header = T,check.names = F,stringsAsFactors = F,sep = "\t")
in.table.productive <- in.table[as.logical(in.table$Productive),]

dat.plot <- dcast(data = in.table.productive,ID~Sample,value.var="TPM")
###dat.plot <- in.table.productive[in.table.productive$Sample %in% rownames(myDesign),]
rownames(dat.plot) <- dat.plot$ID
dat.plot <- as.matrix(dat.plot[,-1])
dat.plot[is.na(dat.plot)] <- 0
dat.plot[dat.plot > 0] <- 1
dat.plot[1:4,1:8]
## use TC
dat.plot <- dat.plot[,grepl("^.TC",colnames(dat.plot),perl=T)]
dat.plot <- dat.plot[,colnames(dat.plot) %in% rownames(myDesign)]
myDesign <- myDesign[colnames(dat.plot),]

plot.tcr.landscape(dat.plot[rowSums(dat.plot)>0,],sprintf("%s/%s.freq1.byTypeT",out.dir,sample.id),sampleTypeColor,myDesign$sampleType,main=sample.id,byType=T)
plot.tcr.landscape(dat.plot[rowSums(dat.plot)>0,],sprintf("%s/%s.freq1.byTypeF",out.dir,sample.id),sampleTypeColor,myDesign$sampleType,main=sample.id,byType=F)
col.f <- colSums(dat.plot[rowSums(dat.plot)>0,])>0
plot.tcr.landscape(dat.plot[rowSums(dat.plot)>0,col.f],sprintf("%s/%s.freq1.ColSum1.byTypeT",out.dir,sample.id),sampleTypeColor,myDesign$sampleType[col.f],main=sample.id,byType=T)
plot.tcr.landscape(dat.plot[rowSums(dat.plot)>0,col.f],sprintf("%s/%s.freq1.ColSum1.byTypeF",out.dir,sample.id),sampleTypeColor,myDesign$sampleType[col.f],main=sample.id,byType=F)

plot.tcr.landscape(dat.plot[rowSums(dat.plot)>1,],sprintf("%s/%s.freq2.byTypeT",out.dir,sample.id),sampleTypeColor,myDesign$sampleType,main=sample.id,byType=T)
plot.tcr.landscape(dat.plot[rowSums(dat.plot)>1,],sprintf("%s/%s.freq2.byTypeF",out.dir,sample.id),sampleTypeColor,myDesign$sampleType,main=sample.id,byType=F)
col.f <- colSums(dat.plot[rowSums(dat.plot)>1,])>0
plot.tcr.landscape(dat.plot[rowSums(dat.plot)>1,col.f],sprintf("%s/%s.freq2.ColSum1.byTypeT",out.dir,sample.id),sampleTypeColor,myDesign$sampleType[col.f],main=sample.id,byType=T)
plot.tcr.landscape(dat.plot[rowSums(dat.plot)>1,col.f],sprintf("%s/%s.freq2.ColSum1.byTypeF",out.dir,sample.id),sampleTypeColor,myDesign$sampleType[col.f],main=sample.id,byType=F)
