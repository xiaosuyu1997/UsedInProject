#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-i", "--infile", type="character", required=TRUE, help="infile list, comma seperated")
parser$add_argument("-s", "--sample", type="character", required=TRUE, help="id list, comma seperated")
parser$add_argument("-l", "--linkFile", type="character", help="for link annotation")
parser$add_argument("-a", "--col", type="integer", default=1, help="column of infile [default %(default)s]")
parser$add_argument("-b", "--background", type="character", help="file used for background, column 1 should be the identifier")
parser$add_argument("-o", "--outprefix", type="character", required=TRUE, help="outprefix")
parser$add_argument("-m", "--method", type="character", default="venn", help="plot type: venn or heatmap [default %(default)s]")
parser$add_argument("-t", "--test", action="store_true", default=FALSE, help="whether test mode [default %(default)s]")
parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="whether verbose mode [default %(default)s]")
args <- parser$parse_args()
print(args)

filename.list <- args$infile
samplename.list <- args$sample
args.col <- args$col
outprefix <- args$outprefix
background.file <- args$background
args.verbose <- args$verbose
args.method <- args$method
link.file <- args$linkFile
args.test <- args$test

suppressPackageStartupMessages(require("ComplexHeatmap"))

#### TEST venn5
#filename.list <- "/WPS1/zhenglt/work/proj_xy/integrated/database/Plitas2016.mmc2.UP.txt,/WPS1/zhenglt/work/proj_xy/integrated/database/Plitas2016.mmc3.UP.txt,/WPS1/zhenglt/work/proj_xy/integrated/database/DeSimone_2016.Treg.txt,/WPS1/zhenglt/work/proj_xy/integrated/database/Tirosh2016.Treg.txt,/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/cross.patient/new.clusteringResult/sigGene/liver.CD4.G6.CTLA4.txt"
#samplename.list <- "Plitas2016.mmc2.UP,Plitas2016.mmc3.UP,DeSimone_2016.Treg,Tirosh2016.Treg,CD4.CTLA4"
#args.col <- 2
#background.file <-  "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/cross.patient/tmp.patch20170209.someFigure/expressedGene.HCC.list"
#outprefix <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/cross.patient/tmp.patch20170209.someFigure/OUT/Treg.cmp.v5" 
#args.method <- "heatmap"
#link.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/cross.patient/tmp.patch20170209.someFigure/OUT/Treg.cmp.v5.linkAnn.txt"
####
#filename.list <- "/WPS1/zhenglt/work/proj_xy/integrated/database/Plitas2016.mmc2.txt,/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/cross.patient/new.clusteringResult/sigGene/liver.CD4.G6.CTLA4.txt"
#samplename.list <- "Plitas2016.mmc2,CD4.CTLA4"
#args.col <- 1
#outprefix <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/cross.patient/tmp.patch20170209.someFigure/test"
#background.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/cross.patient/tmp.patch20170209.someFigure/expressedGene.HCC.list"
#args.verbose <- T

dir.create(dirname(outprefix),showWarnings = F,recursive = T)

link.list <- NULL
if(!is.null(link.file) && file.exists(link.file)){
    link.list <- as.character(read.table(link.file,header = T,sep = "\t",check.names = F,stringsAsFactors = F)[,args.col])
    link.list <- link.list[!is.na(link.list)]
}
file.list <- unlist(strsplit(filename.list,split = ",",perl = T))
sample.list <- strsplit(samplename.list,split = ",",perl = T)
v.list <- lapply(seq_len(length(file.list)),function(i){
           .t.v <- as.character(read.table(file.list[i],header = T,sep = "\t",check.names = F,stringsAsFactors = F)[,args.col])
           .t.v[!is.na(.t.v)]
})
names(v.list) <- unlist(sample.list)
background.list <- NULL
if(!is.null(background.file) && file.exists(background.file)){
    background.list <- as.character(read.table(background.file,header = T,sep = "\t",check.names = F,stringsAsFactors = F)[,args.col])
    background.list <- background.list[!is.na(background.list)]
    for(i in seq_along(v.list)){
        v.list[[i]] <- intersect(v.list[[i]],background.list)
    }
}

###source("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")
n.sample <- NULL
if(length(v.list)==2 && !is.null(background.list)){
    list.and <- intersect(v.list[[1]],v.list[[2]])
    list.diff <- setdiff(v.list[[2]],v.list[[1]])
    x <- length(list.and)
    m <- length(v.list[[1]])
    n <- length(background.list)-m
    k <- length(v.list[[2]])
    n.sample <- m+n
    p.value <- phyper(x-1, m, n, k, lower.tail = F)
    ###
    if(args.verbose){
        out.df <- data.frame(V1=list.and)
        colnames(out.df) <- sprintf("%s-%s",names(v.list)[1],names(v.list)[2])
        write.table(out.df,sprintf("%s.and.txt",outprefix),row.names = F,quote = F,sep = "\t")
        out.df <- data.frame(V1=list.diff)
        colnames(out.df) <- sprintf("%s-%s",names(v.list)[1],names(v.list)[2])
        write.table(out.df,sprintf("%s.diff.txt",outprefix),row.names = F,quote = F,sep = "\t")
    }
}
library("VennDiagram")
library("gridBase")

if(args.method=="venn"){
    pdf(sprintf("%s.venn.pdf",outprefix),width=6,height=6)
    if(length(v.list)==2){
        opar <- par(mar=c(6,6,6,6),cex.lab=1.5,cex.main=1.5,xpd=T)
    }else{
        opar <- par(mar=c(2,2,2,2),cex.lab=1.5,cex.main=1.5,xpd=T)
    }
    plot.new()
    #title(main="venn",sub=sprintf("p=%s",p.value),cex=1.2)
    vps <- baseViewports()
    pushViewport(vps$inner, vps$figure, vps$plot)
    ##venn.plot <- venn.diagram(v.list,filename = NULL, cat.cex=1.5,cex=1.5,
    venn.plot <- venn.diagram(v.list,filename = NULL, cat.cex=1.5,cex=1.5,
                              hyper.test=T,total.population=n.sample,lower.tail=F,sub.pos=c(0.5,0.25),sub.cex=1.5,
                              margin=0.2,cat.dist=0.15,na = "remove")
    grid.draw(venn.plot)
    par(opar)
    dev.off()
}else if(args.method=="heatmap"){
    group.downsample <- 50
    g.list <- Reduce(union,v.list)
    nSet <- length(v.list)
    dat.to.plot <- data.frame(matrix(0, nrow = length(g.list), ncol = length(v.list)))
    colnames(dat.to.plot) <- names(v.list)
    rownames(dat.to.plot) <- g.list
    for(vname in names(v.list)){
        dat.to.plot[v.list[[vname]],vname] <- 1
    }

    f.last.specific <- apply(dat.to.plot,1,function(x){ sum(x==c(rep(0,nSet-1),1))==nSet  })
    out.df <- data.frame(geneID=rownames(dat.to.plot)[f.last.specific])
    write.table(out.df,file = sprintf("%s.last.specific.list",outprefix),quote = F,sep = "\t",row.names = F)
    f.shared <- apply(dat.to.plot,1,function(x){ sum(x==rep(1,nSet))==nSet  })
    out.df <- data.frame(geneID=rownames(dat.to.plot)[f.shared])
    write.table(out.df,file = sprintf("%s.shared.list",outprefix),quote = F,sep = "\t",row.names = F)
    f.shared.m <- apply(dat.to.plot,1,function(x){ sum(head(x,nSet-1)==rep(1,nSet-1))==(nSet-1)  })
    out.df <- data.frame(geneID=rownames(dat.to.plot)[f.shared.m])
    write.table(out.df,file = sprintf("%s.shared.m.list",outprefix),quote = F,sep = "\t",row.names = F)
    out.df <- data.frame(sampleID=rownames(dat.to.plot))
    out.df <- cbind(out.df,dat.to.plot)
    write.table(out.df,file = sprintf("%s.cmpByHeatmap.txt",outprefix),quote = F,sep = "\t",row.names = F)
    if(args.test){
        f.sharedButNotLast <- apply(dat.to.plot,1,function(x){ sum(x[1:(nSet-1)]==c(rep(1,nSet-1)))==(nSet-1)  })
        dat.to.plot <- dat.to.plot[f.sharedButNotLast,]
    }

    n <- nrow(dat.to.plot)
    m <- ncol(dat.to.plot)
    hc.t <- hclust(dist((dat.to.plot)),"complete")
    dat.to.plot <- dat.to.plot[hc.t$order,,drop=F]
    pdf(sprintf("%s.heatmap.pdf",outprefix),width=5,height=8)
    par(mar=c(2,4,4,2))
    vps <- baseViewports()
    pushViewport(vps$inner, vps$figure, vps$plot)
    ha.link <- NULL
    ###if(n>group.downsample){
    if(!is.null(link.list)){
        ann.link.at <- match(link.list,rownames(dat.to.plot))
        ann.link.label <- link.list
        ha.link <- rowAnnotation(link = row_anno_link(at = ann.link.at, labels = ann.link.label,labels_gp=gpar(fontsize=8)), 
                                 width = unit(1, "cm") + max_text_width(ann.link.label))
    }
    ht <- Heatmap(dat.to.plot,"Called",
                col=structure(c("gray","brown1"), names = c("0", "1")),  
                column_dend_height = unit(6, "cm"), row_dend_width = unit(6, "cm"),
                column_names_gp = gpar(fontsize = 12*28/max(m,32)),row_names_gp = gpar(fontsize = 10*28/max(n,32)),
                show_row_names=ifelse(args.test,T,F),
                show_heatmap_legend = T, row_names_max_width = unit(10,"cm"),
                cluster_columns = F,
                cluster_rows = F)
    if(args.test){
        ComplexHeatmap::draw(ht, newpage= FALSE)
    }else{
        if(!is.null(ha.link)){
            ComplexHeatmap::draw(ht+ha.link, newpage= FALSE)
        }else{
            ComplexHeatmap::draw(ht, newpage= FALSE)
        }
    }
    dev.off()
}

