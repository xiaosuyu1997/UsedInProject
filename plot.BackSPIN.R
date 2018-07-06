#!/usr/bin/env Rscript

args <- commandArgs(T)
if(length(args)<6)
{
    cat("usage: plot.BackSPIN.R <out.dir> <sample.id> <gene.grp.file> <sample.grp.file> <exp.mat.file> <cellType.color.file>\n")
    q()
}

out.dir <- args[1]
sample.id <- args[2]
gene.grp.file <- args[3]
sample.grp.file <- args[4]
exp.mat.file <- args[5]
cellType.color.file <- args[6]

#out.dir <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/BackSPIN/OUT/P1118.N50/"
#sample.id <- "P1118"
#gene.grp.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/BackSPIN/OUT/P1118.N50/P1118_clustered.gene.grp"
#sample.grp.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/BackSPIN/OUT/P1118.N50/P1118_clustered.sample.grp"
#exp.mat.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/BackSPIN/OUT/P1118.N50/P1118_clustered.exp.mat"
#cellType.color.file <- "/WPS1/zhenglt/work/TCR_chunhong/data/CellType.color"

source("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")

plotBackSPIN <- function(out.prefix,dat.plot,sampleTypeColor,sampleType,sample.grp,gene.grp)
{
    require("ComplexHeatmap")
    require("circlize")
    require("gridBase")
    ###require("dendextend")
    require("gplots")
    pdf(sprintf("%s.pdf",out.prefix),width=15,height=15)
    #png(sprintf("%s.png",out.prefix),width=1500,height=1500)
    plot.new()
    legend("topright",legend=names(sampleTypeColor),fill=sampleTypeColor,border=sampleTypeColor,cex=1.5)
    par(mar=c(5,4,4,6))
    ### Integrating Grid Graphics Output with Base Graphics Output
    vps <- baseViewports()
    pushViewport(vps$inner, vps$figure, vps$plot)

    annDF <- data.frame(sampleType=sampleType,SPINLevel=sample.grp$Level_4_group)
    annColList <- list(sampleType=sampleTypeColor)
    #if(!is.null(clonotype.col))
    #{
    #    annDF$clonotype=clonotype.col$ctype[colnames(dat.plot)]
    #    annColList$clonotype=clonotype.col$col
    #}
    #print(annColList$clonotype)
    ha.col <- HeatmapAnnotation(df = annDF, col = annColList, show_legend = FALSE)
    top_annotation_height <- unit(1.5*ncol(annDF), "cm")
    dat.plot <- as.matrix(dat.plot)
    dat.plot.unscale <- dat.plot
    #### scale by row
    rowM <- rowMeans(dat.plot, na.rm = T)
    rowSD <- apply(dat.plot, 1, sd, na.rm = T)
    dat.plot <- sweep(dat.plot, 1, rowM)
    dat.plot <- sweep(dat.plot, 1, rowSD, "/")
    m <- ncol(dat.plot)
    n <- nrow(dat.plot)
    dat.plot[dat.plot<-4] <- -4
    dat.plot[dat.plot>4] <- 4
    ####
              #column_dend_height = unit(6, "cm"), row_dend_width = unit(6, "cm"),
    bk.range <- quantile(abs(dat.plot),probs=c(0.01,1))
    ht <- Heatmap(dat.plot, name="ZScore",
              col = colorRamp2(seq(-bk.range[2],bk.range[2],length=5), bluered(5),space="LAB"),
              column_names_gp = gpar(fontsize = 12*55/m),row_names_gp = gpar(fontsize = 10*55/max(n,25)),
              show_heatmap_legend = T, row_names_max_width = unit(10,"cm"),
              top_annotation_height = top_annotation_height,
              cluster_columns = F,
              cluster_rows = F,
              split = gene.grp$Level_4_group,
              heatmap_legend_param = list(grid_width = unit(0.8, "cm"), grid_height = unit(0.8, "cm"), title_gp = gpar(fontsize = 14, fontface = "bold"),
                                          label_gp = gpar(fontsize = 12), color_bar = "continuous"),
              top_annotation = ha.col)
    draw(ht, newpage= FALSE)
    dev.off()
}

gene.grp <- read.table(gene.grp.file,header = T,stringsAsFactors = F)
sample.grp <- read.table(sample.grp.file,header = T,stringsAsFactors = F)
exp.mat <- read.table(exp.mat.file,header = F,stringsAsFactors = F)
sampleTypeColor <- read.SampleTypeColor(cellType.color.file)

rownames(exp.mat) <- gene.grp$geneID
colnames(exp.mat) <- sample.grp$sample
f <- !is.na(gene.grp$gname)
dat.plot <- as.matrix(exp.mat)[f,]
rownames(dat.plot) <- gene.grp$gname[f]
dir.create(out.dir,showWarnings = F,recursive = T)
plotBackSPIN(sprintf("%s/%s.BackSPIN",out.dir,sample.id),dat.plot,sampleTypeColor,sample.grp$cellType,sample.grp,gene.grp[f,])
