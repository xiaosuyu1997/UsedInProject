#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("gplots"))
suppressPackageStartupMessages(library("ComplexHeatmap"))
suppressPackageStartupMessages(library("circlize"))
suppressPackageStartupMessages(library("gridBase"))
suppressPackageStartupMessages(library("dendextend"))
suppressPackageStartupMessages(library("RColorBrewer"))
    
source("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")

#i.file <- "/WPSnew/zhenglt/work/proj_zl/integrate0924/tcr/clonalExpansion/mergeMKI67/colonC.byLoc.CD8.mcls.clone.txt"
i.file <- "/WPSnew/zhenglt/work/proj_zl/integrate0924/tcr/clonalExpansion/mergeMKI67.rm00/colonC.byLoc.CD8.mcls.clone.loc.txt"
out.prefix <- "./GZMK.cluster.tcr.sharing"
pdf.width <- 10
pdf.height <- 12
majorClusterColor.file <- "/WPSnew/zhenglt/work/proj_zl/integrate0924/cross.patient.b20170523/majorClusterColor.mod.txt"

majorClusterColor <- read.SampleTypeColor(majorClusterColor.file)
majorClusterColor <- majorClusterColor[grepl("^CD8",names(majorClusterColor))]

i.table <- read.table(i.file,header = T,check.names = F,stringsAsFactors = F)
rownames(i.table) <- i.table[,1]
i.table <- as.matrix(i.table[,-1])
print(colnames(i.table))
#i.table <- i.table[,c("P.CD8_C02-GPR183","P.CD8_C03-CX3CR1",
#                      "N.CD8_C02-GPR183","N.CD8_C03-CX3CR1","N.CD8_C04-CD160","N.CD8_C05-CD6","N.CD8_C06-GZMK","N.CD8_C07-HAVCR2",
#                      "T.CD8_C02-GPR183","T.CD8_C03-CX3CR1","T.CD8_C04-CD160","T.CD8_C05-CD6","T.CD8_C06-GZMK","T.CD8_C07-HAVCR2")]
i.table <- i.table[,c("P.CD8_C02-GPR183","N.CD8_C02-GPR183","T.CD8_C02-GPR183",
                      "P.CD8_C03-CX3CR1","N.CD8_C03-CX3CR1","T.CD8_C03-CX3CR1",
                      "N.CD8_C04-CD160","T.CD8_C04-CD160",
                      "N.CD8_C05-CD6","T.CD8_C05-CD6",
                      "N.CD8_C06-GZMK","T.CD8_C06-GZMK",
                      "N.CD8_C07-HAVCR2","T.CD8_C07-HAVCR2")]
i.table.gzmk <- ((i.table[,"T.CD8_C06-GZMK"]>0) & apply(i.table,1,function(x){ sum(x>0)>1 })) %>% i.table[.,]


do.plot <- function(dat.plot)
{
    ###dat.plot <- i.table.gzmk

    dat.plot <- dat.plot[hclust(dist(dat.plot))$order,]

    for(pp in c("N.CD8_C06-GZMK","P.CD8_C02-GPR183","P.CD8_C03-CX3CR1","T.CD8_C07-HAVCR2","T.CD8_C06-GZMK"))
    {
        j.vec <- which(grepl(pp,colnames(dat.plot),perl = T))
        dat.plot <- dat.plot[names(sort(apply(dat.plot,1,
                                              function(x){ 
                                                  w <- 0.6^seq_along(j.vec)
                                                  sum( (x[j.vec]>0)*w) }),decreasing = T)),]
    }

    dat.plot[dat.plot>10] <- 10
    m <- ncol(dat.plot)
    n <- nrow(dat.plot)
    top_annotation_height <- 1.2
    branch.col <- F
    branch.row <- F
    z.lo <- 0
    z.hi <- 8 
    z.step <- 1
    pdf(sprintf("%s.pdf",out.prefix),width=pdf.width,height=pdf.height)
    #par(mar=c(5,14,4,4))
    par(mar=c(4,12,4,4))
    plot.new()
    title(main = "",cex.main=2)
    #legend("topright",legend=names(colSet),fill=colSet,border=colSet,cex=1.5,inset=c(-0.03,0),xpd=T)
    ### Integrating Grid Graphics Output with Base Graphics Output
    vps <- baseViewports()
    pushViewport(vps$inner, vps$figure, vps$plot)
    #ha.col <- NULL
    m.pattern <- regexpr("^(.)", colnames(dat.plot))
    annDF <- data.frame(majorCluster=gsub(pattern = "^..","",x = colnames(dat.plot)),
                        loc=regmatches(colnames(dat.plot), m.pattern))
    annColList <- list(majorCluster=majorClusterColor,
                       loc=c("P"="#E41A1C","N"="#A65628","T"="#377EB8"))
    print(head(annDF))

    annotation_legend_param <- list()
    ha.col <- HeatmapAnnotation(df = annDF, col = annColList, show_legend = T, annotation_legend_param = annotation_legend_param)

    print("all T.CD8_C06-GZMK >0 ")
    all(dat.plot[,"T.CD8_C06-GZMK"]>0) %>% print
    clone.grp <- rep("G0",nrow(dat.plot))
    clone.grp[ (dat.plot[,"P.CD8_C02-GPR183"]>0 | dat.plot[,"P.CD8_C03-CX3CR1"]>0) & dat.plot[,"T.CD8_C07-HAVCR2"]==0 ] <- "linkToTemTeff"
    clone.grp[ !(dat.plot[,"P.CD8_C02-GPR183"]>0 | dat.plot[,"P.CD8_C03-CX3CR1"]>0) & dat.plot[,"T.CD8_C07-HAVCR2"]>0 ] <- "linkToTex"
    row.annDF <- data.frame(clone.grp=clone.grp,stringsAsFactors = F)
    write.table(c(data.frame(clone.id=rownames(dat.plot),stringsAsFactors = F),
                  row.annDF),sprintf("%s.row.annDF.txt",out.prefix),row.names = F,quote = F,sep = "\t")
    clone.grop.colSet <- c("linkToTemTeff"="#A6CEE3","linkToTex"="#B15928","G0"="gray")
    ha.row <- HeatmapAnnotation(df = row.annDF, col = list(clone.grp=clone.grop.colSet), 
                                show_legend = T, annotation_legend_param = list(),which = "row")
    ht <- Heatmap(dat.plot, name="",
            ##col = colorRamp2(seq(-bk.range[2],bk.range[2],length=100), 
            col = colorRamp2(seq(z.lo,z.hi,length=9),
                             colorRampPalette((c("lightgrey",brewer.pal(n = 6, name = "YlOrRd"))))(9), space="LAB"),
                             ##colorRampPalette((c(brewer.pal(n = 6, name = "Blues"))))(9), space="LAB"),
                             #colorRampPalette((c(brewer.pal(n = 6, name = "Purples"))))(9), space="LAB"),
                             ###colorRampPalette((brewer.pal(n = 7, name = "Blues")))(11), space="LAB"),
                             ###colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(11), space="LAB"),
            column_dend_height = unit(6, "cm"), row_dend_width = unit(6, "cm"),
            column_names_gp = gpar(fontsize = 12*28/max(m,32)),row_names_gp = gpar(fontsize = 10*28/max(n,32)),
            show_heatmap_legend = T, row_names_max_width = unit(10,"cm"),
            #top_annotation_height = top_annotation_height,
            cluster_columns = branch.col,
            cluster_rows = branch.row,
            heatmap_legend_param = list(grid_width = unit(0.8, "cm"), 
                                        grid_height = unit(0.8, "cm"), 
                                        at = seq(z.lo,z.hi,z.step),
                                        title_gp = gpar(fontsize = 14, fontface = "bold"),
                                        label_gp = gpar(fontsize = 12), color_bar = "discrete"),
            top_annotation = ha.col)
    ComplexHeatmap::draw(ha.row+ht, newpage= FALSE)
    for(i in seq_along(names(ha.col@anno_list))){
        decorate_annotation(names(ha.col@anno_list)[i], 
                            {grid.text(names(ha.col@anno_list)[i], unit(-4, "mm"),gp=gpar(fontsize=14),just = "right")})
    }
    dev.off()
}

do.plot(i.table.gzmk)
