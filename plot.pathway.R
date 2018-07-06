#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-o", "--outputPrefix", type="character", required=TRUE, help="output prefix")
parser$add_argument("-i", "--plistFile", type="character", required=TRUE, help="plist.file")
parser$add_argument("-r", "--resListFile", type="character", required=TRUE, help="res.list.file")
parser$add_argument("-w", "--clusterColorFile", type="character",default="/WPS1/zhenglt/work/proj_xy/integrated/cross.patient.rename/curated.clusterColor.txt", 
                    help="clusterColorFile [default %(default)s]")
parser$add_argument("-a", "--filterDiag", action="store_true", default=FALSE, help="filter.diag [default %(default)s]")
parser$add_argument("-b", "--filterDiagN", type="integer", default=1, help="filter.diag.n [default %(default)s]")
parser$add_argument("-n", "--max", type="integer", default=100000, help="max number of gene sets [default %(default)s]")
parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="verbose mode [default %(default)s]")
##parser$add_argument("-s", "--sample", type="character", default="SAMPLE", help="sample id  [default %(default)s]")
parser$add_argument("-d", "--pdfWidth", type="double", default=1.0, help="pdf width adjust [default %(default)s]")
###parser$add_argument("-u", "--pdfHeight", type="integer", default=8, help="pdf heihgt [default %(default)s]")
args <- parser$parse_args()
print(args)

out.prefix <- args$outputPrefix
plist.file <- args$plistFile
res.list.file <- args$resListFile
clusterColorFile <- args$clusterColorFile
pdfWidth.adj <- args$pdfWidth
args.filterDiag <- args$filterDiag
args.filterDiagN <- args$filterDiagN
args.max <- args$max

source("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")
suppressPackageStartupMessages(library("ComplexHeatmap"))
suppressPackageStartupMessages(library("circlize"))
suppressPackageStartupMessages(library("gridBase"))
suppressPackageStartupMessages(library("dendextend"))
suppressPackageStartupMessages(library("RColorBrewer"))
suppressPackageStartupMessages(library("gplots"))

#out.prefix <- "/WPS1/zhenglt/work/proj_xy/integrated/cross.patient.rename/test/test"
#plist.file <- "/WPS1/zhenglt/work/proj_xy/integrated/cross.patient.rename/pathway_list/CD8.GO.BP.plist"
#res.list.file <- "/WPS1/zhenglt/work/proj_xy/integrated/cross.patient.rename/pathway_list/CD8.GO.BP.res.list"
#clusterColorFile <- "/WPS1/zhenglt/work/proj_xy/integrated/cross.patient.rename/curated.clusterColor.txt"

###out.prefix <- "/WPS1/zhenglt/work/proj_xy/integrated/cross.patient.rename/pathway_list/CD8.KEGG"
###plist.file <- "/WPS1/zhenglt/work/proj_xy/integrated/cross.patient.rename/pathway_list/CD8.KEGG.plist"
###res.list.file <- "/WPS1/zhenglt/work/proj_xy/integrated/cross.patient.rename/pathway_list/CD8.KEGG.list"
###clusterColorFile <- "/WPS1/zhenglt/work/proj_xy/integrated/cross.patient.rename/curated.clusterColor.txt"

plist <- read.table(plist.file,header = F,stringsAsFactors = F,check.names = F)$V1
res.list <- read.table(res.list.file,header = F,stringsAsFactors = F,check.names = F)
rownames(res.list) <- res.list[,1]
    
majorClusterColor <- read.SampleTypeColor(clusterColorFile)
majorClusterColor <- majorClusterColor[names(majorClusterColor) %in% res.list$V1]

dat.plot <- (sapply(seq_len(nrow(res.list)),function(i){
    .i.table <- read.table(res.list[i,2],header = T,stringsAsFactors = F,check.names = F,sep="\t")
    rownames(.i.table) <- .i.table[,1]
    .ret <- .i.table[plist,"p.value"]
    #.ret <- data.frame(V=.i.table[plist,"p.value"])
    #colnames(.ret) <- res.list[i,1]
    return(.ret)
}))
colnames(dat.plot) <- res.list[,1]
rownames(dat.plot) <- plist

dat.plot <- -log10(dat.plot)

### sorting v1
#dat.plot.bin <- dat.plot > -log10(0.01)
#dat.plot.bin <- dat.plot.bin[order(apply(dat.plot.bin,1,function(x){ sum(x*(0.5^seq_along(x))) }),decreasing = T),]
#dat.plot <- dat.plot[rownames(dat.plot.bin),]
### sorting v2
    
#pdf(sprintf("%s.gene.clusters.clusMax.pdf",out.prefix),width = ifelse(ncol(dat.plot)>4,6,4),height = 10)
    #heatmap.2(dat.plot.tmp, col=colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu"))), 
    #          symbreak=TRUE, trace='none', dendrogram="none", 
    #          RowSideColors=clust.col[gene.desc[rownames(dat.plot),"Group"]],Rowv=F, 
    #          ColSideColors=column.color, Colv = F, scale="none", 
    #          margin=c(14, 6),cexRow=min(1.8,55/nrow(dat.plot)),cexCol=min(1.2,40/ncol(dat.plot)),
    #          breaks=seq(max(-2,dat.plot.range[1]), min(2,dat.plot.range[2]), length.out=21))
    #dev.off()
do.plot <- function(dat.plot,do.scale=F,filter.diag=F,filter.diag.n=2)
{
    #dat.plot.bin <- dat.plot > -log10(0.01)
    #dat.plot <- dat.plot[order(apply(dat.plot.bin,1,function(x){ sum(x*(0.01^seq_along(x))) }),decreasing = T),]
    dat.plot <- dat.plot[order(apply(dat.plot,1,function(x){ sum(x*(x>2)*(0.01^seq_along(x))) }),decreasing = T),]
    gset.name.vec <- c()
    if(filter.diag){
        .f.diag <- apply(dat.plot,1,function(x){ sum(x>2)<=filter.diag.n })
        dat.plot <- dat.plot[.f.diag,]
        for(i in seq_along(colnames(dat.plot))){
            #.out.df <- data.frame(geneSet=rownames(dat.plot)[dat.plot[,i]>2])
            #write.table(.out.df,file = sprintf("%s.dat.plot.DiagN%d.%s.txt",out.prefix,filter.diag.n,colnames(dat.plot)[i]),sep = "\t",row.names = F,col.names = F,quote = F)
            .i.table <- read.table(res.list[colnames(dat.plot)[i],2],header = T,stringsAsFactors = F,check.names = F,sep="\t")
            rownames(.i.table) <- .i.table[,1]
            .gset <- rownames(dat.plot)[dat.plot[,i]>2]
            gset.name.vec <- union(gset.name.vec,head(.gset,n=args.max))
            write.table(.i.table[.gset,],
                        file = sprintf("%s.dat.plot.DiagN%d.%s.txt",out.prefix,filter.diag.n,colnames(dat.plot)[i]),
                        sep = "\t",row.names = F,col.names = T,quote = F)
        }
    }
    dat.plot <- dat.plot[gset.name.vec,]

    dat.plot.df <- data.frame(geneSet=rownames(dat.plot))
    dat.plot.df <- cbind(dat.plot.df,dat.plot)
    write.table(dat.plot.df,file = sprintf("%s.dat.plot.txt",out.prefix),sep = "\t",row.names = F,col.names = T,quote = F)
    ###dat.plot <- dat.plot[order(apply(dat.plot,1,function(x){ sum(x*(0.5^seq_along(x))) }),decreasing = T),]
    n <- nrow(dat.plot)
    m <- ncol(dat.plot)
    pdf(sprintf("%s.scale%s.pdf",out.prefix,ifelse(do.scale,"T","F")),width = ifelse(ncol(dat.plot)>4,12*pdfWidth.adj,9*pdfWidth.adj),height = 10)
    par(mar=c(4,12,4,4))
    plot.new()
    title(main = "",cex.main=2)
    #legend("topright",legend=names(colSet),fill=colSet,border=colSet,cex=1.5,inset=c(-0.03,0),xpd=T)
    ### Integrating Grid Graphics Output with Base Graphics Output
    vps <- baseViewports()
    pushViewport(vps$inner, vps$figure, vps$plot)
    annDF <- data.frame(majorCluster=colnames(dat.plot))
    annColList <- list(majorCluster=majorClusterColor)
    g.show.legend <- T
    top_annotation_height <- unit(0.5 * ncol(annDF), "cm")
    ha.col <- HeatmapAnnotation(df = annDF, col = annColList, show_legend = T, annotation_legend_param = list())
    #ha.row <- HeatmapAnnotation(df = data.frame(Group=rowCluster), col = list(Group=rowClusterColor), 
    #                            show_legend = g.show.legend, annotation_legend_param = list(),which = "row")
    if(do.scale)
    {
        z.lo <- -3
        z.hi <- 3
        z.step <- 0.5
        rowM <- rowMeans(dat.plot, na.rm = T)
        rowSD <- apply(dat.plot, 1, sd, na.rm = T)
        dat.plot <- sweep(dat.plot, 1, rowM)
        dat.plot <- sweep(dat.plot, 1, rowSD, "/")
        dat.plot[dat.plot < z.lo] <- z.lo
        dat.plot[dat.plot > z.hi] <- z.hi
        ###print(dat.plot[1:4,1:8])
    }else{
        z.lo <- 0
        z.hi <- 5
        z.step <- 1
        dat.plot[dat.plot>z.hi] <- z.hi
        dat.plot[dat.plot<2] <- 0
        dat.plot <- floor(dat.plot)
        #tmp.var <- pretty(dat.plot,n=8)
        #z.lo <- tmp.var[1]
        #z.hi <- tmp.var[length(tmp.var)]
        #z.step <- tmp.var[2]-tmp.var[1]
    }
    ht <- Heatmap(dat.plot,"-log(p)",
                col = colorRamp2(seq(z.lo,z.hi,length=6), 
                                 if(do.scale) colorRampPalette(rev(brewer.pal(n = 6, name = "RdBu")))(6) else colorRampPalette((brewer.pal(n = 6, name = "Purples")))(6), 
                                 space="LAB" ),
                                 ##colorRampPalette(rev(brewer.pal(n = 7, name = ifelse(do.scale,"RdBu","Blues"))))(10), space="LAB"),
                                 ###colorRampPalette(rev(brewer.pal(n = 7, name = ifelse(do.scale,"RdBu","RdYlBu"))))(10), space="LAB"),
                column_dend_height = unit(6, "cm"), row_dend_width = unit(6, "cm"),
                column_names_gp = gpar(fontsize = 12*28/max(m,32)),row_names_gp = gpar(fontsize = 10*28/max(n,32)),
                show_row_names=T,
                show_heatmap_legend = T, row_names_max_width = unit(10,"cm"),
                top_annotation_height = top_annotation_height,
                cluster_columns = F,
                cluster_rows = F,
                heatmap_legend_param = list(grid_width = unit(0.8, "cm"), 
                                            grid_height = unit(0.8, "cm"), 
                                            at = seq(z.lo,z.hi,z.step),
                                            title_gp = gpar(fontsize = 14, fontface = "bold"),
                                            label_gp = gpar(fontsize = 12), color_bar = "discrete"),
                top_annotation = ha.col)
                                            ###label_gp = gpar(fontsize = 12), color_bar = "continuous"),
                ####cell_fun = function(j, i, x, y, w, h, col) { grid.text(sprintf("%4.1f",dat.plot[i, j]), x, y) },
                                 ###if(do.scale) colorRampPalette(rev(brewer.pal(n = 6, name = "RdBu")))(6) else c("white",colorRampPalette((brewer.pal(n = 6, name = "Blues")))(5)), 
    ComplexHeatmap::draw(ht, newpage= FALSE)
    for(i in seq_along(names(ha.col@anno_list))){
        decorate_annotation(names(ha.col@anno_list)[i], 
                            {grid.text(names(ha.col@anno_list)[i], unit(-4, "mm"),gp=gpar(fontsize=14),just = "right")})
    }
    dev.off()
}

##do.plot(dat.plot,F)
do.plot(dat.plot,F,filter.diag = args.filterDiag,filter.diag.n = args.filterDiagN)
