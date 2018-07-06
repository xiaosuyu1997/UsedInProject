#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-i", "--inFile", type="character", required=TRUE, help="infile")
parser$add_argument("-a", "--sampleDescFile", type="character", required=TRUE, help="sample.desc.file")
parser$add_argument("-b", "--geneDescFile", type="character", required=TRUE, help="gene.desc.file")
parser$add_argument("-o", "--outputPrefix", type="character", required=TRUE, help="output prefix")
parser$add_argument("-s", "--sample", type="character", default="SAMPLE", help="sample id  [default %(default)s]")
parser$add_argument("-y", "--cellTypeColorFile", type="character", help="cellTypeColorFile")
parser$add_argument("-l", "--clonotypeFile", type="character", help="clonotype file")
parser$add_argument("-u", "--u", action="store_true", default=FALSE, help="unique gene [default %(default)s]")
parser$add_argument("-e", "--scale", action="store_true", default=FALSE, help="do scale [default %(default)s]")
args <- parser$parse_args()
print(args)

in.file <- args$inFile
sample.desc.file <- args$sampleDescFile
gene.desc.file <- args$geneDescFile
out.prefix <- args$outputPrefix
sample.id <- args$sample
cellTypeColorFile <- args$cellTypeColorFile
clonotype.file <- args$clonotypeFile
opt.u <- args$u
do.scale <- args$scale

##opt.u <- FALSE
##do.scale <- F
##in.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/gsnap_out/quantification/sizeFactor/P0508/P0508.all.countData.sfNormalized.txt.gz"
##sample.desc.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/iterative/P0508.sample.txt"
##gene.desc.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/iterative/P0508.markerGene.slim"
##out.prefix <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/iterative/test.P0508.markerGene"
##sample.id <- "P0508"
##cellTypeColorFile <- "/WPS1/zhenglt/work/TCR_chunhong/data/CellType.color"
##clonotype.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/tracer/OUT/P0508/filtered_TCR_summary/P0508.summary.cell.reassigneClonotype.methodChunhong.r.txt"

source("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")

suppressPackageStartupMessages(require("ComplexHeatmap"))
suppressPackageStartupMessages(require("circlize"))
suppressPackageStartupMessages(require("gridBase"))
suppressPackageStartupMessages(require("dendextend"))
suppressPackageStartupMessages(require("RColorBrewer"))
suppressPackageStartupMessages(require("gplots"))

### clonotype data
clonotype.data <- read.clonotype(in.file = clonotype.file,ctype.col = "C_strict")
#### exp data
in.table <- read.table(in.file,header = T,sep = "\t",check.names = F,stringsAsFactors = F,row.names = 1)
cnames <- entrezToXXX(rownames(in.table))
cnames.na <- which(is.na(cnames))
cnames[cnames.na] <- rownames(in.table)[cnames.na]
rownames(in.table) <- cnames
in.table <- in.table[,-1]
in.table <- as.matrix(log2(in.table+1))

print(dim(in.table))
print(in.table[1:4,1:8])

#### sample data
sample.desc <- read.table(sample.desc.file,header = T,sep = "\t",check.names = F,stringsAsFactors = F)
rownames(sample.desc) <- sample.desc$sample
print(dim(sample.desc))
print(head(sample.desc))

#### gene data
gene.desc <- read.table(gene.desc.file,header = F,sep = "\t",check.names = F,stringsAsFactors = F)
colnames(gene.desc) <- c("patient","geneName","markerClass")
gene.desc <- subset(gene.desc,markerClass!="uncharacterized")
print(dim(gene.desc))
print(head(gene.desc))

#### cell type color
sampleTypeColor <- read.SampleTypeColor(cellTypeColorFile)
sampleTypeColor <- sampleTypeColor[intersect(names(sampleTypeColor),unique(sample.desc$sampleType))]
leafClusterColor <- structure(colorRampPalette(brewer.pal(12,"Paired"))(length(unique(sample.desc$leafCluster))),
                              names=unique(sample.desc$leafCluster))
leafClusterColor["uncharacterized"] <- "gray"


m.cls.plot <- unique(gene.desc$markerClass)

### make new data to be used for plot
grp.plot.one <- c()
dat.plot.one <- c()
if(opt.u){
    gene.desc <- data.frame(patient=gene.desc$patient[1],geneName=unique(gene.desc$geneName),markerClass="M",stringsAsFactors = F)
    dat.plot.one <- in.table[gene.desc$geneName,sample.desc$sample,drop=F]
    grp.plot.one <- gene.desc$markerClass
}else{
    for(m.cls in m.cls.plot)
    {
      m.gene.desc <- subset(gene.desc,markerClass==m.cls)
      dat.plot <- in.table[m.gene.desc$geneName,sample.desc$sample,drop=F]
      rownames(dat.plot) <- sprintf("%s::%s",m.cls,rownames(dat.plot))
      dat.plot.one <- rbind(dat.plot.one,dat.plot)
      grp.plot.one <- append(grp.plot.one,m.gene.desc$markerClass)
    }
}
grp.plot.one <- factor(c(grp.plot.one),levels=c(unique(grp.plot.one)))
#### scale by row
if(do.scale){
    Z.LIMIT.LO <- -3
    Z.LIMIT.HI <- 3
    rowM <- rowMeans(dat.plot.one, na.rm = T)
    rowSD <- apply(dat.plot.one, 1, sd, na.rm = T)
    dat.plot.one <- sweep(dat.plot.one, 1, rowM)
    dat.plot.one <- sweep(dat.plot.one, 1, rowSD, "/")
}else
{
    Z.LIMIT.LO <- 0
    Z.LIMIT.HI <- 15
}

dat.plot <- c()
for(s.cls in unique(sample.desc$leafCluster))
{
    s.desc <- subset(sample.desc,leafCluster==s.cls)
    dat.t <- dat.plot.one[,s.desc$sample,drop=F]
    hc.t <- hclust(dist(t(dat.t)), "complete")
    dat.t <- dat.t[,hc.t$order,drop=F]
    dat.plot <- cbind(dat.plot,dat.t)
}
dat.plot.one <- dat.plot

dat.plot.one[dat.plot.one < Z.LIMIT.LO] <- Z.LIMIT.LO
dat.plot.one[dat.plot.one > Z.LIMIT.HI] <- Z.LIMIT.HI


pdf(sprintf("%s.%s.%s.pdf",out.prefix,sample.id,"all"),width=30,height=20)

plot.new()
par(mar=c(5,4,4,4))
title(main=sprintf("%s",sample.id),cex.main=2)
#legend("top",horiz=T,legend=names(sampleTypeColor),fill=sampleTypeColor,
#       border=sampleTypeColor,cex=2.5,inset=c(0,-0.18),xpd=T)
### Integrating Grid Graphics Output with Base Graphics Output
vps <- baseViewports()
pushViewport(vps$inner, vps$figure, vps$plot)
ht.list = NULL

allSampleCls <- unique(sample.desc$leafCluster)
i <- 0
n <- nrow(dat.plot.one)
m <- ncol(dat.plot.one)
###bk.range <- quantile(dat.plot.one,probs=c(0.01,1),na.rm=T)
bk.range <- quantile(dat.plot.one,probs=c(0,1),na.rm=T)

### sample annotation
annDF <- data.frame(leafCluster=sample.desc$leafCluster,
                    sampleType=sample.desc$sampleType)
annColList <- list(leafCluster=leafClusterColor,sampleType=sampleTypeColor)
if(!is.null(clonotype.data))
{
    annDF$clonotype=clonotype.data$ctype[sample.desc$sample]
    annDF$clonotype[is.na(annDF$clonotype)] <- "NA"
    annColList$clonotype=clonotype.data$col
}
  
ha.col <- HeatmapAnnotation(df = annDF, col = annColList, show_legend = T)
top_annotation_height <- unit(1.5*ncol(annDF), "cm")
### construct heatmap
                             ###col = colorRamp2(seq(bk.range[1],bk.range[2],length=50), 
ht.list <- Heatmap(dat.plot.one,name=sprintf("ht%d",1),
                             na_col = "white",
                             col = if(!do.scale) colorRamp2( 
                                                 seq(bk.range[1],bk.range[2],length=100), 
                                                 colorRampPalette(rev(brewer.pal(n = 7, name =  "RdYlBu")))(100),
                                              space="LAB")
                                 else colorRamp2( 
                                                 c(seq(bk.range[1],-1.96,length=20),
                                                   seq(-1.96,1.96,length=60),
                                                   seq(1.96,bk.range[2],length=20)), 
                                                 colorRampPalette(rev(brewer.pal(n = 7, name =  "RdBu")))(100),
                                              space="LAB"),
                             column_names_gp = gpar(fontsize = 16*100/max(m,30)),
                             row_names_gp = gpar(fontsize = 14*30/max(n,25)),
                             cluster_columns = F,cluster_rows = T, clustering_distance_rows = "spearman",
                             show_column_dend = F,show_row_dend = F, 
                             column_dend_height = unit(20, "mm"),
                             split = if(opt.u) NULL else grp.plot.one,
                             row_title_gp = gpar(fontsize = 24), row_title_rot = 0,
                             gap = unit(1,"mm"),
                             top_annotation_height = top_annotation_height,
                             show_heatmap_legend = T,
                             show_row_names = T,
                             show_column_names = F,
                             heatmap_legend_param = list(title="",
                                                         at=seq(Z.LIMIT.LO,Z.LIMIT.HI,2),
                                                         grid_width = unit(2, "cm"), 
                                                         grid_height = unit(2, "cm"), 
                                                         title_gp = gpar(fontsize = 24, fontface = "bold"),
                                                         label_gp = gpar(fontsize = 24), color_bar = "continuous"),
                             top_annotation = ha.col
)

draw(ht.list, newpage= F)
####popViewport(n = 3)

dev.off()

