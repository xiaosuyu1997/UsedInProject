#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-i", "--infile", type="character", required=TRUE, help="infile list, comma seperated")
parser$add_argument("-o", "--outprefix", type="character", required=TRUE, help="outprefix")
parser$add_argument("-s", "--sampleID", type="character", default="SAMPLE", help="title [default %(default)s]")
parser$add_argument("-c", "--cmp", type="character", required=TRUE, help="comparison")
parser$add_argument("-n", "--ngenes", type="integer", default=20, help="top n genes [default %(default)s]")
parser$add_argument("-r", "--reverse", action="store_true", default=FALSE, help="reverse the logFC [default %(default)s]")
parser$add_argument("-l", "--geneFile", type="character", help="for genes to highlight")
##parser$add_argument("-a", "--col", type="integer", default=1, help="column of infile [default %(default)s]")
parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="whether verbose mode [default %(default)s]")
args <- parser$parse_args()
print(args)

in.file <- args$infile
out.prefix <- args$outprefix
cmp <- args$cmp
args.r <- args$reverse
gene.file <- args$geneFile
args.verbose <- args$verbose
sample.id <- args$sampleID
args.ngenes <- args$ngenes
### Text mode
###in.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/cross.patient/tmp.patch20170209.someFigure/data_foreign/TTR_diff.txt"
###out.prefix <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/cross.patient/tmp.patch20170209.someFigure/data_foreign/volcano"
###cmp <- "Treg"
###args.r <- F
###gene.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/cross.patient/tmp.patch20170209.someFigure/data_foreign/OUT/Treg.cmp.v5.shared.txt"
###sample.id <- "Treg"
###args.ngenes <- 23
### RData mode
#in.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/cross.patient/geneClustering/OUT.rmIsCycle.kmeans/Treg.byPatientF/Treg.gene.clusters.all.RData"
#out.prefix <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/cross.patient/tmp.patch20170209.someFigure/volcano/Treg"
#cmp <- "CD4_Cluster3-CTLA4-CD4_Cluster2-IL2RA"
#args.r <- F
#gene.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/cross.patient/tmp.patch20170209.someFigure/OUT/Treg.cmp.v5.last.specific.txt"
#sample.id <- "Treg_C3_.vs._C2"
#args.ngenes <- 24


suppressPackageStartupMessages(library("R.utils"))
suppressPackageStartupMessages(require("ComplexHeatmap"))
library("ggplot2")
library("ggrepel")
###source("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")
if(grepl(pattern = "\\.RData$",in.file,perl = T)){
    lenv <- loadToEnv(in.file)
    aov.res <- lenv[["aov.res"]]
    g.GNAME <- aov.res[["aov.out"]][,"geneSymbol"]
    names(g.GNAME) <- rownames(aov.res[["aov.out"]])
    g.GNAME.f <- which(is.na(g.GNAME))
    g.GNAME[g.GNAME.f] <- names(g.GNAME)[g.GNAME.f]
    dat.plot <- data.frame(x=aov.res[["aov.out"]][,sprintf("HSD.diff.%s",cmp)],y=(aov.res[["aov.out"]][,sprintf("HSD.padj.%s",cmp)]))
    rownames(dat.plot) <- rownames(aov.res[["aov.out"]])
}else if(grepl(pattern = "\\.txt$",in.file,perl = T)){
    in.table <- read.table(in.file,header = T,sep = "\t",stringsAsFactors = F,check.names = F,row.names = 1)
    g.GNAME <- in.table[,"geneSymbol"]
    names(g.GNAME) <- rownames(in.table)
    g.GNAME.f <- which(is.na(g.GNAME))
    g.GNAME[g.GNAME.f] <- names(g.GNAME)[g.GNAME.f]
    dat.plot <- data.frame(x=in.table[,"logFC"],y=in.table[,"adj.P.Val"])
    rownames(dat.plot) <- rownames(in.table)
}
if(args.r){ dat.plot$x <- -dat.plot$x }
dat.plot$Significant <- abs(dat.plot$x)>1 & dat.plot$y<0.01
colnames(dat.plot) <- c("log2FoldChange","pvalue","Significant")
dat.plot$pvalue[is.infinite(-log10(dat.plot$pvalue)) | dat.plot$pvalue < 1e-300] <- 1e-300
y.range <- range(-log10(dat.plot$pvalue))

if(!is.null(gene.file) && file.exists(gene.file)){
    .t.table <- read.table(gene.file,header = T,sep = "\t",check.names = F,stringsAsFactors = F)
    rownames(.t.table) <- .t.table[,1]
    #print(head(.t.table))
    g.f <- as.character(.t.table[,1])
    g.n <- as.character(.t.table[,2])
    if(args.ngenes>0) { 
        g.f <- g.f[!grepl(pattern = "-AS1$",x = g.n,perl=T)]
        g.f <- head(g.f,n=args.ngenes)
        #print(g.f)
    }
}else{
    g.f <- rownames(subset(dat.plot,Significant==T))
    if(args.ngenes>0) { 
        g.f <- g.f[!grepl(pattern = "-AS1$",x = g.f,perl=T)]
        g.f <- head(g.f,n=args.ngenes) 
    }
}
#print(head(dat.plot))

    ##geom_point(aes(color = Significant)) + scale_color_manual(values = c("grey", "red")) + 
### scale_color_discrete(name = 'NStudies')
#tmpDat <- cbind(dat.plot[g.f,],data.frame(NStudies=paste("N",.t.table[g.f,"NStudies"],sep="")))
#print(head(tmpDat))
my.plot <- ggplot(dat.plot, aes(x = log2FoldChange, y = -log10(pvalue))) + 
    geom_point(aes(color = Significant),size=0.5,data = subset(dat.plot,Significant==F)) + 
    geom_point(aes(color = Significant),size=0.5,data = subset(dat.plot,Significant==T)) + 
    scale_color_manual(values = c("grey", "red","#1B9E77", "#D95F02", "#7570B3")) +  
    labs(title = sample.id) +
    coord_cartesian(ylim = c(0, y.range[2]+5),xlim=c(-7,7)) +
    theme_bw(base_size = 12) + 
    theme(legend.position = "bottom", 
          title=element_text(size=18),
          axis.text=element_text(size=14), 
          axis.title=element_text(size=16)) + 
    if("NStudies" %in% colnames(.t.table)){
        geom_text_repel(data = cbind(dat.plot[g.f,],data.frame(NStudies=paste("X",.t.table[g.f,"NStudies"],sep=""))) , 
                        aes(label = g.GNAME[rownames(dat.plot[g.f,])],color=factor(NStudies)), 
                        size = 4, box.padding = unit(0.35, "lines"),fontface=3, 
                        point.padding = unit(0.3, "lines"))
    }else{
        geom_text_repel(data = dat.plot[g.f,], 
                        aes(label = g.GNAME[rownames(dat.plot[g.f,])]), 
                        size = 4, box.padding = unit(0.35, "lines"),fontface=3, 
                        point.padding = unit(0.3, "lines"))
    }
pdf(sprintf("%s.volcano.pdf",out.prefix),width=6,height=6)
par(mar=c(6,5,4,2),cex.lab=1.5)
print(my.plot)
dev.off()

