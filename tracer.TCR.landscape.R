#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-i", "--inputFile", type="character", required=TRUE, help="input file")
parser$add_argument("-c", "--cellTypeColorFile", type="character", help="cellTypeColorFile")
parser$add_argument("-y", "--sampleTypeColorFile", type="character", required=TRUE, help="sampleTypeColorFile")
parser$add_argument("-d", "--designFile", type="character", required=TRUE, help="design matrix file")
parser$add_argument("-s", "--sample", type="character", default="SAMPLE", help="sample id  [default %(default)s]")
###parser$add_argument("-p", "--pattern", type="character", default="^.TC", help="regular pattern used to select subset cells [default %(default)s]")
parser$add_argument("-r", "--cellTypeOrder", type="character", default="TTC,NTC,PTC,TTR,NTR,PTR,TTH,NTH,PTH", help="cell type order [default %(default)s]")
parser$add_argument("-o", "--outDir", type="character", required=TRUE, help="output directory")
parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="whether verbose mode [default %(default)s]")
args <- parser$parse_args()

print(args)

suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("RColorBrewer"))

tracer.summary.file <- args$inputFile
cellType.color.file <- args$cellTypeColorFile
sampleType.color.file <- args$sampleTypeColorFile
design.file <- args$designFile
sample.id <- args$sample
out.dir <- args$outDir
##grp.pattern <- args$pattern
cellTypeOrder <- args$cellTypeOrder

groupBy <- "majorCluster"
#groupBy <- "sampleType"
sampleType <- "sampleType"

#tracer.summary.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/tracer/OUT.landscape/P1118/^.TC/P1118.input.txt"
#cellType.color.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/clustering/CellType.color"
#design.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/OUT.scLVM/P1118/sfEndo/P1118.designUsed.txt"
#sample.id <- "P1118"
#out.dir <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/tracer/OUT.landscape/test"
#grp.pattern <- "^.TC"
#cellType.color.file <- "/WPS1/zhenglt/work/proj_xy/integrated/cross.patient/clusteringResult/lung.tSNE.colorByCluster.useAOV.n3000.lung.majorClusterColor.txt"

dir.create(out.dir,showWarnings = F,recursive = T)

plot.tcr.landscape <- function(dat.plot,out.prefix,
                               colSet,sampleType,
                               colSetII=NULL,sampleTypeII=NULL,
                               sampleTypeOrder=c("TTC","NTC","PTC","TTR","NTR","PTR","TTH","NTH","PTH","CD8.T","CD8.N","CD8.P","CD4.T","CD4.N","CD4.P"),
                               main="",byType=T,tfirst=FALSE)
{
    suppressPackageStartupMessages(library("ComplexHeatmap"))
    suppressPackageStartupMessages(library("circlize"))
    suppressPackageStartupMessages(library("gridBase"))
    suppressPackageStartupMessages(library("dendextend"))
    suppressPackageStartupMessages(library("gplots"))
    sampleType <- as.character(sampleType)
    names(sampleType) <- colnames(dat.plot)
    if(!is.null(sampleTypeII)){
        sampleTypeII <- as.character(sampleTypeII)
        names(sampleTypeII) <- colnames(dat.plot)
    }
    n <- nrow(dat.plot)
    m <- ncol(dat.plot)
    if(n<2||m<2){ cat(sprintf("Too few genes/samples: #genes: %s, #samples: %s\n",n,m)); return(NULL); }

    ### reorder according sampleTypeOrder
    dat.plot.reorder <- c()
    for(iType in sampleTypeOrder) {
        dat.block <- dat.plot[,sampleType[colnames(dat.plot)]==iType,drop=F]
        dat.plot.reorder <- cbind(dat.plot.reorder,dat.block)
    }
    dat.plot <- dat.plot.reorder
    sampleType <- sampleType[colnames(dat.plot)]
    if(!is.null(sampleTypeII)){
        sampleTypeII <- sampleTypeII[colnames(dat.plot)]
    }

    n <- nrow(dat.plot)
    m <- ncol(dat.plot)
    print(dim(dat.plot))

    #dat.plot[dat.plot==0] <- 0
    #dat.plot[dat.plot>0] <- 1

    #pdf(sprintf("%s.pdf",out.prefix),width=22,height=12)
    pdf(sprintf("%s.pdf",out.prefix),width=14,height=8)
    par(mar=c(5,4,4,3),cex.main=1.5)
    
    ### plot 0
    plot.new()
    #legend("topright",legend=names(colSet)[names(colSet) %in% sampleType],
    #       fill=colSet[names(colSet) %in% sampleType],
    #       border=colSet[names(colSet) %in% sampleType],
    #       cex=1.5,inset=c(-0.23,-0.06),xpd=T)
    #legend("right",legend=names(colSetII)[names(colSetII) %in% sampleTypeII],
    #       fill=colSetII[names(colSetII) %in% sampleTypeII],
    #       border=colSetII[names(colSetII) %in% sampleTypeII],
    #       cex=1.5,inset=c(-0.23,-0.06),xpd=T)
    title(main)
    vps <- baseViewports()
    pushViewport(vps$inner, vps$figure, vps$plot)

    if(!is.null(sampleTypeII)){
        annDF <- data.frame(sampleType=sampleTypeII,cellType=sampleType)
        annColList <- list(sampleType=colSetII,cellType=colSet)
    }else{
        annDF <- data.frame(cellType=sampleType)
        annColList <- list(cellType=colSet)
    }
    print("TTTT:")
    print(annColList)
    ha.col <- HeatmapAnnotation(df = annDF, col = annColList, show_legend = T)
    top_annotation_height <- unit(0.8*ncol(annDF), "cm")

    if(byType)
    {

        oncoprint_row_order <- function(){
            order(rowSums(dat.plot), decreasing = T)
        }
        rrorder <- oncoprint_row_order()

        oncoprint_column_order <- function(){
            offset <- 0
            corder <- c()
            for(iType in sampleTypeOrder) {
                dat.block <- dat.plot[,sampleType[colnames(dat.plot)]==iType,drop=F]
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
                if(ncol(dat.block)==0) { next }
                scores = apply(dat.block[rrorder, ,drop=F], 2, scoreCol)
                corder <- append(corder,c(order(scores, decreasing = TRUE) + offset))
                offset  <- offset + ncol(dat.block)
            }
            return(corder)
        }
        ccorder <- oncoprint_column_order()
        ## put tumor specific in first part
        if(tfirst) {
            t.chose <- rep(FALSE,length(rrorder))
            t.rrorder <- c()
            for(iType in sampleTypeOrder) {
                #tType <- sampleType[colnames(dat.plot)[1]]
                dat.block <- dat.plot[rrorder,sampleType[colnames(dat.plot)]==iType,drop=F]
                tType.f <- apply(dat.block,1,function(x){ sum(x)>0 })
                t.chose.block <- !t.chose & tType.f 
                t.rrorder <- append(t.rrorder,rrorder[t.chose.block])
                t.chose <- t.chose | t.chose.block
                #print(iType)
                #print(tType.f)
                #print(t.chose)
            }
            t.rrorder <- append(t.rrorder,rrorder[!t.chose])
            #cat("final t.chose:\n")
            #print(t.chose)
            rrorder <- t.rrorder
        }

        #ht <- oncoPrint(list(TCR=dat.plot),
        #save.image(sprintf("%s.tmp.RData",out.prefix))
        save(dat.plot,rrorder,ccorder,ha.col,top_annotation_height,file=sprintf("%s.tmp.RData",out.prefix))

        #.tdim <- dim(dat.plot)
        #.t.dat.plot <- matrix(as.character(dat.plot),nrow=.tdim[1],byrow = T)
        #rownames(.t.dat.plot) <- rownames(dat.plot)
        #colnames(.t.dat.plot) <- colnames(dat.plot)
        #ht <- oncoPrint((.t.dat.plot),
        .t.dat.plot <- dat.plot[rrorder,]
        .t.dat.plot <- .t.dat.plot[,ccorder]
        #ht <- Heatmap(dat.plot, name="",
        print(table(.t.dat.plot))
        ht <- Heatmap(.t.dat.plot, name="",
                    #col = colorRamp2(seq(-bk.range[2],bk.range[2],length=100), 
                    col = c("0"="gray","1"="black","2"="red"),
                    #col = colorRamp2(seq(0,2,length=100), colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100), space="LAB"),
                    row_names_gp = gpar(fontsize = 3.5*55/max(nrow(dat.plot),25)),
                    column_names_gp = gpar(fontsize = 2.5*160/max(ncol(dat.plot),25)),
                    show_heatmap_legend = F, row_names_max_width = unit(10,"cm"),
                    top_annotation_height = top_annotation_height,
                    cluster_columns = FALSE, cluster_rows = FALSE,
                    top_annotation = ha.col)
         ComplexHeatmap::draw(ht, newpage= FALSE)
       
#        ht <- oncoPrint(list(TCR=dat.plot),
#                        show_heatmap_legend = F,
#                        show_row_names=T,show_column_names=T,
#                        row_names_gp = gpar(fontsize = 3.5*55/max(nrow(dat.plot),25)),
#                        column_names_gp = gpar(fontsize = 2.5*160/max(ncol(dat.plot),25)),
#                        pct_gp = gpar(fontsize = 2.2*120/max(nrow(dat.plot),25)),
#                        show_row_barplot=F,
#                        top_annotation = ha.col,top_annotation_height = top_annotation_height,
#                        row_order=rrorder,
#                        column_order = ccorder,col = c("0"="gray","1"="black","2"="red"))

        draw(ht, newpage= FALSE)
    }else {
        ht <- oncoPrint(list(TCR=dat.plot),
                        show_heatmap_legend = F,
                        show_row_names=T,show_column_names=T,
                        row_names_gp = gpar(fontsize = 3.5*55/max(nrow(dat.plot),25)),
                        column_names_gp = gpar(fontsize = 2.5*160/max(ncol(dat.plot),25)),
                        pct_gp = gpar(fontsize = 2.2*120/max(nrow(dat.plot),25)),
                        show_column_barplot=F,
                        show_row_barplot=F,
                        top_annotation = ha.col,top_annotation_height = top_annotation_height)
        draw(ht, newpage= FALSE)

    }
    dev.off()
}

source("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")

myDesign<-read.delim(design.file,header=T,row.names="sample",check.names=F,stringsAsFactors=F)
###cellTypeColor <- read.SampleTypeColor(cellType.color.file)
###sampleTypeColor <- read.SampleTypeColor(cellType.color.file)
in.table <- read.table(tracer.summary.file,header = T,check.names = F,stringsAsFactors = F,sep = "\t")
in.table.productive <- in.table[as.logical(in.table$Productive),]

dat.plot <- dcast(data = in.table.productive,ID~Sample,value.var="TPM")
rownames(dat.plot) <- dat.plot$ID
dat.plot <- as.matrix(dat.plot[,-1])
dat.plot[is.na(dat.plot)] <- 0
dat.plot[dat.plot > 0] <- 1
#### clonal & non-clonal
f.clonal <- rowSums(dat.plot)>=2
tmp.dat.blcok <- dat.plot[f.clonal,]
tmp.dat.blcok[tmp.dat.blcok==1] <- 2
dat.plot[f.clonal,] <- tmp.dat.blcok
print(table(dat.plot))
print(dim(dat.plot))
#dat.plot[1:4,1:8]
s.f <- intersect(rownames(myDesign),colnames(dat.plot))
dat.plot <- dat.plot[,s.f,drop=F]
myDesign <- myDesign[s.f,,drop=F]

print(str(myDesign))
print(head(myDesign))
hadGroupBy <- groupBy %in% colnames(myDesign)
cat(sprintf("had %s column ? %s\n",groupBy,hadGroupBy))
if(hadGroupBy){
    print(head(myDesign[s.f,groupBy]))
    nMajor <- length(unique(myDesign[,groupBy]))
    if(!is.null(cellType.color.file) && file.exists(cellType.color.file)){
        cellTypeColor <- read.SampleTypeColor(cellType.color.file)
        cellTypeColor <- cellTypeColor[ names(cellTypeColor) %in% unique(myDesign[,groupBy]) ]
    }else{
        cellTypeColor <- structure(colorRampPalette(brewer.pal(nMajor,"Paired"))(nMajor),
                                  names=unique(myDesign[,groupBy]))
    }
    print(cellTypeColor)
}

mMajor <- length(unique(myDesign[,sampleType]))
if(!is.null(sampleType.color.file) && file.exists(sampleType.color.file)){
    sampleTypeColor <- read.SampleTypeColor(sampleType.color.file)
    sampleTypeColor <- sampleTypeColor[ names(sampleTypeColor) %in% unique(myDesign[,sampleType]) ]
}else{
    sampleTypeColor <- structure(colorRampPalette(brewer.pal(mMajor,"Paired"))(mMajor),
                              names=unique(myDesign[,sampleType]))
}
print(sampleTypeColor)

#plot.tcr.landscape(dat.plot[rowSums(dat.plot)>0,,drop=F],sprintf("%s/%s.freq1.byTypeT",out.dir,sample.id),sampleTypeColor,myDesign$sampleType,sampleTypeOrder=unlist(strsplit(cellTypeOrder,split=",")),main=sample.id,byType=T,tfirst=T)
###col.f <- colSums(dat.plot[rowSums(dat.plot)>0,,drop=F])>0
###plot.tcr.landscape(dat.plot[rowSums(dat.plot)>0,col.f],sprintf("%s/%s.freq1.ColSum1.byTypeT",out.dir,sample.id),sampleTypeColor,myDesign$sampleType[col.f],main=sample.id,byType=T,tfirst=T)

loginfo(sprintf("begin"))

if(hadGroupBy){
    plot.tcr.landscape(dat.plot[rowSums(dat.plot)>1,,drop=F],
                       out.prefix = sprintf("%s/%s.freq2.byType.%s",out.dir,sample.id,groupBy),
                       colSet = cellTypeColor,
                       sampleType = myDesign[,groupBy],
                       colSetII = sampleTypeColor,
                       sampleTypeII = myDesign[,sampleType],
                       sampleTypeOrder=unique(myDesign[,groupBy]),
                       main=sample.id,byType=T,tfirst=T)

    plot.tcr.landscape(dat.plot[rowSums(dat.plot)>1,,drop=F],
                       out.prefix = sprintf("%s/%s.freq2.byType.%s.NULL",out.dir,sample.id,groupBy),
                       colSet = cellTypeColor,
                       sampleType = myDesign[,groupBy],
                       colSetII = NULL,
                       sampleTypeII = NULL,
                       sampleTypeOrder=unique(myDesign[,groupBy]),
                       main=sample.id,byType=T,tfirst=T)
}

#cat("by sampleType (2)\n")
#plot.tcr.landscape(dat.plot[rowSums(dat.plot)>1,,drop=F],
#                   out.prefix = sprintf("%s/%s.freq2.byType.%s",out.dir,sample.id,sampleType),
#                   colSetII = if(hadGroupBy) cellTypeColor else NULL,
#                   sampleTypeII = if(hadGroupBy) myDesign[,groupBy] else NULL,
#                   colSet = sampleTypeColor,
#                   sampleType = myDesign[,sampleType],
#                   sampleTypeOrder=unique(myDesign[,sampleType]),
#                   main=sample.id,byType=T,tfirst=T)

cat("by sampleType (2)\n")
##plot.tcr.landscape(dat.plot[rowSums(dat.plot)>1,,drop=F],
plot.tcr.landscape(dat.plot[apply(dat.plot,1,function(x){ sum(x>0)>=2 }),,drop=F],
                   out.prefix = sprintf("%s/%s.freq2.byType.%s.NULL",out.dir,sample.id,sampleType),
                   colSetII = NULL,
                   sampleTypeII = NULL,
                   colSet = sampleTypeColor,
                   sampleType = myDesign[,sampleType],
                   #sampleTypeOrder=unique(myDesign[,sampleType]),
                   main=sample.id,byType=T,tfirst=T)

cat("by sampleType (1)\n")
##plot.tcr.landscape(dat.plot[rowSums(dat.plot)>0,,drop=F],
plot.tcr.landscape(dat.plot[apply(dat.plot,1,function(x){ sum(x>0)>=1 }),,drop=F],
                   out.prefix = sprintf("%s/%s.freq1.byType.%s.NULL",out.dir,sample.id,sampleType),
                   colSetII = NULL,
                   sampleTypeII = NULL,
                   colSet = sampleTypeColor,
                   sampleType = myDesign[,sampleType],
                   #sampleTypeOrder=unique(myDesign[,sampleType]),
                   main=sample.id,byType=T,tfirst=T)


loginfo(sprintf("end"))
