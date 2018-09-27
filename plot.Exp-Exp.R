#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-i", "--infile", type="character", required=TRUE, help="infile list, comma seperated")
parser$add_argument("-o", "--outprefix", type="character", required=TRUE, help="outprefix")
parser$add_argument("-c", "--cmpColumn", type="character", required=TRUE, help="comparison columns, seperated by ','")
parser$add_argument("-p", "--pvalueColumn", type="character", required=TRUE, help="pvalue column")
parser$add_argument("-f", "--fcColumn", type="character", required=TRUE, help="log2foldChange column")
parser$add_argument("-r", "--fdr", type="double", default=0.05, help="threshold for p value [default %(default)s]")
parser$add_argument("-l", "--geneFile", type="character", help="for genes to highlight")
parser$add_argument("-s", "--sampleID", type="character", default="SAMPLE", help="title [default %(default)s]")
parser$add_argument("-n", "--ngenes", type="integer", default=20, help="top n genes [default %(default)s]")
parser$add_argument("-x", "--xlim", type="integer", default=7, help="xlim [default %(default)s]")
parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="whether verbose mode [default %(default)s]")
args <- parser$parse_args()
print(args)

in.file <- args$infile
out.prefix <- args$outprefix
cmp.column.str <- args$cmpColumn
pvalue.column.str <- args$pvalueColumn
fc.column.str <- args$fcColumn
args.fdr <- args$fdr
gene.file <- args$geneFile
sample.id <- args$sampleID
args.ngenes <- args$ngenes
args.verbose <- args$verbose
args.xlim <- args$xlim

### Text mode
#in.file <- "/WPSnew/zhenglt/work/proj_zl/integrate0924/cross.patient.b20170523/functional/tcrBase/version_merge_GPR183_CX3CR1_rm00/OUT.DE/gzmk.yellow.orange.red/gzmk.yellow.orange.red.de.aov.geneTable.all.txt"
#out.prefix <- "/WPSnew/zhenglt/work/proj_zl/integrate0924/cross.patient.b20170523/functional/tcrBase/version_merge_GPR183_CX3CR1_rm00/OUT.DE/gzmk.yellow.orange.red/gzmk.yellow.orange.red.de.aov.geneTable.all"
#cmp.column.str <- "avg.linkToTemTeff,avg.linkToTex"
####pvalue.column.str <- "HSD.padj.linkToTex-linkToTemTeff"
#pvalue.column.str <- "F.adjp"
#fc.column.str <- "HSD.diff.linkToTex-linkToTemTeff"
#args.fdr <- 0.1
#gene.file <- NULL
#sample.id <- "merge_GPR183_CX3CR1_rm00"
#args.ngenes <- 0
####

cmp.column <- unlist(strsplit(cmp.column.str,split = ",",perl = T))

library("ggplot2")
library("ggrepel")
if(file.exists("/lustre1/zeminz_pkuhpc/02.pipeline/cancer/lib/scRNAToolKit.R")){
	source("/lustre1/zeminz_pkuhpc/02.pipeline/cancer/lib/scRNAToolKit.R")
}else{
	source("/WPSnew/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")
}

in.table <- read.table(in.file,header = T,sep = "\t",stringsAsFactors = F,check.names = F)
rownames(in.table) <- in.table[,1]
g.GNAME <- in.table[,"geneSymbol"]
names(g.GNAME) <- rownames(in.table)
g.GNAME.f <- which(is.na(g.GNAME))
g.GNAME[g.GNAME.f] <- names(g.GNAME)[g.GNAME.f]
    
dat.plot <- data.frame(x=in.table[,cmp.column[1]],y=in.table[,cmp.column[2]],
                       p.value=in.table[,pvalue.column.str],
                       fc=in.table[,fc.column.str],
                       geneID=in.table[,"geneID"],
                       geneSymbol=in.table[,"geneSymbol"],
                       stringsAsFactors = F)
rownames(dat.plot) <- rownames(in.table)
dat.plot$Significant <- abs(dat.plot$fc)>1 & dat.plot$p.value<args.fdr

#dat.plot$pvalue[is.infinite(-log10(dat.plot$pvalue)) | dat.plot$pvalue < 1e-300] <- 1e-300
#y.range <- range(-log10(dat.plot$pvalue))

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

ngene.up <- sum(dat.plot$fc>1 & dat.plot$p.value<args.fdr)
print(ngene.up)
ann.df <- data.frame(xpos=-Inf,ypos=Inf,hjustvar=-0.5,vjustvar=2.0,annotateText=sprintf("UP:%d",ngene.up))
print(ann.df)
my.plot.ExpExp <- ggplot(dat.plot, aes(x = x, y = y)) + 
    geom_point(aes(color = Significant),size=0.5,data = subset(dat.plot,Significant==F)) + 
    geom_point(aes(color = Significant),size=1.5,data = subset(dat.plot,Significant==T)) + 
    scale_color_manual(values = c("grey", "red","#1B9E77", "#D95F02", "#7570B3")) +  
    geom_abline(intercept = 0, slope = 1) +
    geom_abline(intercept = -1, slope = 1,linetype="dashed") +
    geom_abline(intercept = 1, slope = 1,linetype="dashed") +
    labs(title = sample.id) + xlab(cmp.column[1]) + ylab(cmp.column[2]) + 
    #coord_cartesian(ylim = c(0, y.range[2]+5),xlim=c(-7,7)) +
    theme_bw(base_size = 12) + 
    theme(legend.position = "bottom", 
          title=element_text(size=18),
          axis.text=element_text(size=14), 
          axis.title=element_text(size=16)) +
    geom_text_repel(data = dat.plot[g.f,], 
                    aes(label = g.GNAME[rownames(dat.plot[g.f,])]), 
                    size = 4, box.padding = unit(0.35, "lines"),fontface=3, 
                    point.padding = unit(0.3, "lines"))
ggsave(sprintf("%s.Exp-Exp.pdf",out.prefix),width=6,height=6)


my.plot.volcano(dat.plot=dat.plot,out.prefix=sprintf("%s",out.prefix),
                my.xlim=c(-args.xlim,args.xlim),
                g.f=g.f,
                col.x="fc",col.y="p.value", sample.id=sample.id,g.GNAME=g.GNAME)

#dat.plot$p.value[is.infinite(-log10(dat.plot$p.value)) | dat.plot$p.value < 1e-300] <- 1e-300
#y.pretty <- pretty(-log10(dat.plot$p.value))
#y.range <- y.pretty[c(1,length(y.pretty))]
#my.plot.volcano <- ggplot(dat.plot, aes(x = fc, y = -log10(p.value))) + 
#    geom_point(aes(color = Significant),size=0.5,data = subset(dat.plot,Significant==F)) + 
#    geom_point(aes(color = Significant),size=1.0,data = subset(dat.plot,Significant==T)) + 
#    scale_color_manual(values = c("grey", "red","#1B9E77", "#D95F02", "#7570B3")) +  
#    labs(title = sample.id) + xlab("log2FoldChange") + ylab("-log10(pvalue)") + 
#    geom_text(data=ann.df, aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText)) +
#    coord_cartesian(ylim = c(0, y.range[2]),xlim=c(-args.xlim,args.xlim)) +
#    theme_bw(base_size = 12) + 
#    theme(legend.position = "bottom", 
#          title=element_text(size=18),
#          axis.text=element_text(size=14), 
#          axis.title=element_text(size=16)) + 
#    geom_text_repel(data = dat.plot[g.f,], 
#                    aes(label = g.GNAME[rownames(dat.plot[g.f,])]), 
#                    size = 4, box.padding = unit(0.35, "lines"),fontface=3, 
#                    point.padding = unit(0.3, "lines"))
#ggsave(sprintf("%s.volcano.pdf",out.prefix),width=6,height=6)
#
