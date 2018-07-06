#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-i", "--inputFile", type="character", required=TRUE, help="input file")
parser$add_argument("-s", "--sample", type="character", default="SAMPLE", help="sample id  [default %(default)s]")
parser$add_argument("-o", "--outPrefix", type="character", required=TRUE, help="output prefix")
args <- parser$parse_args()

print(args)
sample.id <- args$sample
out.prefix <- args$outPrefix
in.file <- args$inputFile

#sample.id <- "P0205"
#out.prefix <- "P0205.test.alpha1.beta1"
#in.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/tracer/OUT/P0205/filtered_TCR_summary/P0205.summary.cell.reassigneClonotype.methodChunhong.r.txt"

source("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")
suppressPackageStartupMessages(library(reshape2))

in.table <- read.delim(in.file,header = T,check.names = F,sep = "\t",stringsAsFactors = F)

plot.tcr.landscape <- function(dat.plot,out.prefix,main="")
{
    suppressPackageStartupMessages(require("ComplexHeatmap"))
    suppressPackageStartupMessages(require("circlize"))
    suppressPackageStartupMessages(require("gridBase"))
    suppressPackageStartupMessages(require("dendextend"))
    suppressPackageStartupMessages(require("gplots"))
    
    n <- nrow(dat.plot)
    m <- ncol(dat.plot)
    if(n<2||m<2){ cat(sprintf("Too few genes/samples: #genes: %s, #samples: %s\n",n,m)); return(NULL); }

    print(dim(dat.plot))

    pdf(sprintf("%s.pdf",out.prefix),width=22,height=12)
    par(mar=c(5,4,4,20),cex.main=1.5)
    
    ### plot 0
    plot.new()
    #legend("topright",legend=names(colSet)[names(colSet) %in% sampleType],
    #       fill=colSet[names(colSet) %in% sampleType],
    #       border=colSet[names(colSet) %in% sampleType],
    #       cex=1.5,inset=c(-0.08,-0.06),xpd=T)
    title(main)
    vps <- baseViewports()
    pushViewport(vps$inner, vps$figure, vps$plot)

    ht <- oncoPrint(list(TCR=dat.plot),
                    show_heatmap_legend = F,
                    show_row_names=T,show_column_names=T,
                    row_names_gp = gpar(fontsize = 3.2*55/max(nrow(dat.plot),25)),
                    column_names_gp = gpar(fontsize = 3.2*160/max(ncol(dat.plot),25)),
                    pct_gp = gpar(fontsize = 3.2*120/max(nrow(dat.plot),25)),
                    show_column_barplot=F,
                    show_row_barplot=F
                    )
    draw(ht, newpage= FALSE)
    
    dev.off()
}

plot.tcr.pair <- function(dat.plot,out.prefix,sample.id)
{
    require("ComplexHeatmap")
    require("circlize")
    require("gridBase")
    require("dendextend")
    pdf(sprintf("%s.heatmap.pdf",out.prefix),width=12,height=12)
    par(mar=c(5,4,4,4))
    plot.new()
    title(main=sample.id,cex.main=2)
    mtext("TCR beta Genes",side = 1,line = 1,cex=2)
    mtext("TCR alpha Genes",side = 4,line = 1,cex=2)
    vps <- baseViewports()
    pushViewport(vps$inner, vps$figure, vps$plot)

    oncoprint_row_order <- function(){
        order(rowSums(dat.plot), decreasing = T)
    }
    rrorder <- oncoprint_row_order()

    oncoprint_column_order <- function(){
        offset <- 0
        corder <- c()

        scoreCol = function(x) { 
            score = 0
            for (i in 1:length(x)) { 
                if (x[i]) { 
                    score = score + 2^(length(x) - i * 1/x[i])
                }
            }
            return(score)
        }
        scores = apply(dat.plot[rrorder, ,drop=F], 2, scoreCol)
        corder <- order(scores, decreasing = TRUE)
        return(corder)
    }
    ccorder <- oncoprint_column_order()

    ht <- Heatmap(dat.plot, col = colorRamp2(c(0, 1, 10), c("white", "#FFB4B4", "red")),
                    column_names_gp = gpar(fontsize = 12*45/ncol(dat.plot)),row_names_gp = gpar(fontsize = 10*45/max(nrow(dat.plot),25)),
                    row_order = rrorder,column_order = ccorder,
                    cluster_rows = FALSE, cluster_columns = FALSE,
                    show_heatmap_legend = FALSE,
                    heatmap_legend_param = list(grid_width = unit(0.8, "cm"), 
                                                    grid_height = unit(0.8, "cm"), 
                                                    title_gp = gpar(fontsize = 14, fontface = "bold"),
                                                    at = 0:10,
                                                    grid_border = "black",
                                                    title = "occurrence",
                                                    label_gp = gpar(fontsize = 12), color_bar = "discrete")
                  )
    draw(ht, newpage= FALSE)
    dev.off()
}


in.table[,"ID_Alpha1"] <- in.table[,"Identifier(Alpha1)"]
in.table[,"ID_Beta1"] <- in.table[,"Identifier(Beta1)"]
in.table[,"occurence"] <- 1

dat.plot <- dcast(data = in.table,ID_Alpha1~ID_Beta1,value.var="occurence",fill = 0)
rownames(dat.plot) <- dat.plot[,"ID_Alpha1"]
dat.plot <- as.matrix(dat.plot[,-1])
##plot.tcr.landscape(dat.plot,out.prefix,main="")
plot.tcr.pair(dat.plot,out.prefix,sample.id)


