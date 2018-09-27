#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-i", "--inFile", type="character", required=TRUE, help="infile")
parser$add_argument("-o", "--outputPrefix", type="character", required=TRUE, help="output prefix")
parser$add_argument("-s", "--sample", type="character", default="SAMPLE", help="sample [default %(default)s]")
parser$add_argument("-c", "--clone", type="character", default="C_strict", help="column name of clone data [default %(default)s]")
parser$add_argument("-a", "--clusterColorFile", type="character", 
                    default="/WPSnew/zhenglt/work/proj_zl/integrate0924/cross.patient.b20180409/majorClusterColor.mod.txt",
                    help="clusterColorFile [default %(default)s]")
parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="verbose mode [default %(default)s]")
args <- parser$parse_args()
print(args)

out.prefix <- args$outputPrefix
in.file <- args$inFile
clone.colname <- args$clone
clusterColorFile <- args$clusterColorFile

#out.prefix <- "OUT.index/colonC.test"
#in.file <- "all.summary.cellV2.reassigneClonotype.methodChunhongV2.r.flt-zl.addInfo.mergeMKI67.txt"
#clone.colname <- "C_strict"
####source("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")
source("/WPSnew/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")
#clusterColorFile <- "/WPSnew/zhenglt/work/proj_zl/integrate0924/cross.patient.b20170523/majorClusterColor.mod.txt"
majorClusterColor <- read.SampleTypeColor(in.file = clusterColorFile)

dir.create(dirname(out.prefix),showWarnings = F,recursive = T)

library("ggplot2")
library("ggpubr")
library("gplots")
library("RColorBrewer")
library("plotrix")
library("reshape2")
library("igraph")
library("entropy")
library("dplyr")

##c.data <- read.clonotype(in.file,"C_strict")
c.data <- read.clonotype(in.file,clone.colname)

#dat.plot <- c.data$itable[,c("patient","sampleType","majorCluster","C_strict")]
dat.plot <- c.data$itable[!is.na(c.data$itable$majorCluster) & c.data$itable$majorCluster!="uncharacterized",
                          c("patient","sampleType","stype","majorCluster",clone.colname)]
clone.freq <- sort(table(dat.plot[,clone.colname]),decreasing = T)
dat.plot$ctype <- c.data$ctypeII[rownames(dat.plot)]
m <- regexpr("^(.)", dat.plot$sampleType,perl = T)
dat.plot$loc <- regmatches(dat.plot$sampleType, m)
dat.plot$loc <- factor(dat.plot$loc,levels = c("P","N","T"))
dat.plot$majorCluster <- gsub(pattern = "Cluster",replacement = "C",x = dat.plot$majorCluster)
dat.plot$loc.majorCluster <- sprintf("%s:%s",dat.plot$loc,dat.plot$majorCluster)

calIndex <- function(dat.plot,out.prefix)
{
    ret <- list()
    clust.size <- unclass(table(dat.plot$majorCluster))
    patient.size <- unclass(table(dat.plot$patient))
    clone.size <- unclass(sort(table(dat.plot[,clone.colname]),decreasing = T))
    clone2patient <- unique(dat.plot[,c("patient",clone.colname)])
    clone2patient <- structure(clone2patient$patient,names=clone2patient[,clone.colname])

    clone.index.df <- data.frame(clone.id=names(clone.size),
                                 clone.size=clone.size,
                                 patient=clone2patient[names(clone.size)],
                                 stringsAsFactors = F)
    clone.index.df$patient.size <- patient.size[clone.index.df$patient]
    clone.index.df$clone.global.p <- clone.index.df$clone.size/clone.index.df$patient.size
    #clone.global.p.star <- aggregate(clone.global.p[names(clone.size)],by=list(clone2patient[names(clone.size)]),function(x){ log2(x/min(x)) })
    #clone.global.p.star.vec <- c()
    #for(i in seq_along(clone.global.p.star$x)){
    #    clone.global.p.star.vec <- append(clone.global.p.star.vec,clone.global.p.star$x[[i]])
    #}
    clone.index.df$clone.index.exp <- log2(clone.index.df$clone.size)
    clone.index.df$clone.index.exp.v2 <- log2(clone.index.df$patient.size)*log2(clone.index.df$clone.global.p)
    clone.index.df$clone.index.exp.v3 <- (clone.index.df$patient.size/(mean(patient.size)))*log2(clone.index.df$clone.global.p)
    .t.cor.V2 <- cor(clone.index.df$clone.index.exp,clone.index.df$clone.index.exp.v2,method = "pearson")
    .t.cor.V3 <- cor(clone.index.df$clone.index.exp,clone.index.df$clone.index.exp.v3,method = "pearson")
    #pdf(sprintf("%s.clone.index.exp.cor.pdf",out.prefix),width = 4,height = 4)
    #plot(clone.index.df$clone.index.exp,clone.index.df$clone.index.exp.v2,xlab="log2(n)",ylab="log2(m)*log2(n/m)")
    #mtext(text = sprintf("Cor.pearson=%4.4f",.t.cor.V2),line = -1,adj = 0.05)
    #plot(clone.index.df$clone.index.exp,clone.index.df$clone.index.exp.v3,xlab="log2(n)",ylab="(m/mean(m))*log2(n/m)")
    #mtext(text = sprintf("Cor.pearson=%4.4f",.t.cor.V3),line = -1,adj = 0.05)
    #dev.off()

    clone.dist.loc <- table(dat.plot[,c(clone.colname,"loc")])
    clone.index.df$clone.index.mig <- apply(clone.dist.loc[rownames(clone.index.df),,drop=F],1,entropy.empirical,unit="log2")
    clone.index.df$clone.index.migNorm <- apply(clone.dist.loc[rownames(clone.index.df),,drop=F],1,
                                                function(x){ x[x>1] <- 1; entropy.empirical(x,unit="log2") })
    clone.dist.clust <- table(dat.plot[,c(clone.colname,"majorCluster")])
    clone.index.df$clone.index.dev <- apply(clone.dist.clust[rownames(clone.index.df),,drop=F],1,entropy.empirical,unit="log2")
    #f.stype <- grepl(pattern = sprintf("CD4"),colnames(clone.dist.clust),perl = T)
    #clone.index.dev.CD4 <- apply(clone.dist.clust[,f.stype],1,entropy.empirical,unit="log2")
    #clone.index.dev.CD4[is.na(clone.index.dev.CD4)] <- 0
    #clone.index.dev.CD8 <- apply(clone.dist.clust[,!f.stype],1,entropy.empirical,unit="log2")
    #clone.index.dev.CD8[is.na(clone.index.dev.CD8)] <- 0
    #clone.index.dev <- clone.index.dev.CD4+clone.index.dev.CD8
    #clone.index.df <- data.frame(clone.id=names(clone.size),
    #                             clone.size=clone.size,
    #                             global.p=clone.global.p[names(clone.size)],
    #                             global.p.star=clone.global.p.star.vec[names(clone.size)],
    #                             index.exp=clone.index.exp[names(clone.size)],
    #                             index.mig=clone.index.mig[names(clone.size)],
    #                             index.dev=clone.index.dev[names(clone.size)],
    #                             stringsAsFactors = F)
    cluster.index.df <- data.frame(majorCluster=colnames(clone.dist.clust),stringsAsFactors = F)
    cluster.index.df$avg.global.p <- apply(clone.dist.clust,2,function(x){ sum(x*clone.index.df[names(x),"clone.global.p"])/sum(x) })
    #cluster.index.df$index.byP <- apply(clone.dist.clust,2,function(x){ sum(x*clone.index.df[names(x),"global.p.star"])/sum(x) })
    cluster.index.df$index.exp <- apply(clone.dist.clust,2,function(x){ sum(x*clone.index.df[names(x),"clone.index.exp"])/sum(x) })
    cluster.index.df$index.exp.V2 <- apply(clone.dist.clust,2,function(x){ sum(x*clone.index.df[names(x),"clone.index.exp.v2"])/sum(x) })
    cluster.index.df$index.exp.V3 <- apply(clone.dist.clust,2,function(x){ sum(x*clone.index.df[names(x),"clone.index.exp.v3"])/sum(x) })

    cluster.index.df$index.mig <- apply(clone.dist.clust,2,function(x){ sum(x*clone.index.df[names(x),"clone.index.mig"])/sum(x) })
    ###
    cluster.index.df$index.migNorm <- apply(clone.dist.clust,2,function(x){ sum(x*clone.index.df[names(x),"clone.index.migNorm"])/sum(x) })
    cluster.index.df$index.dev <- apply(clone.dist.clust,2,function(x){ sum(x*clone.index.df[names(x),"clone.index.dev"])/sum(x) })
    cluster.index.df$shannon.entropy <- apply(clone.dist.clust,2,function(x){ entropy.empirical(x,unit="log2") })
    cluster.index.df$shannon.entropy.max <- apply(clone.dist.clust,2,function(x){ log2(sum(x>0)) })
    cluster.index.df$shannon.entropy.scale <- cluster.index.df$shannon.entropy/cluster.index.df$shannon.entropy.max
    cluster.index.df$index.exp.V4 <- 1-cluster.index.df$shannon.entropy.scale
    cluster.index.df$n <- apply(clone.dist.clust,2,sum)
    ##cluster.index.df$n <- clust.size[cluster.index.df$majorCluster]
    cluster.index.df$index.exp.max <- log2(cluster.index.df$n)
    cluster.index.df$index.mig.max <- log2(3)
    cluster.index.df$index.dev.max <- 0
    f.stype <- grepl(pattern = sprintf("CD4"),colnames(clone.dist.clust),perl = T)
    cluster.index.df[grepl("CD4",cluster.index.df$majorCluster,perl = T),"index.dev.max"] <- log2(sum(f.stype))
    cluster.index.df[!grepl("CD4",cluster.index.df$majorCluster,perl = T),"index.dev.max"] <- log2(sum(!f.stype))
    cluster.index.df$index.exp.scale <- cluster.index.df$index.exp/cluster.index.df$index.exp.max
    cluster.index.df$index.mig.scale <- cluster.index.df$index.mig/cluster.index.df$index.mig.max
    cluster.index.df$index.dev.scale <- cluster.index.df$index.dev/cluster.index.df$index.dev.max
    write.table(clone.index.df,sprintf("%s.clone.index.txt",out.prefix),row.names = F,quote = F,sep = "\t")
    write.table(cluster.index.df,sprintf("%s.cluster.index.txt",out.prefix),row.names = F,quote = F,sep = "\t")
    ret[["clone.index"]] <- clone.index.df
    ret[["cluster.index"]] <- cluster.index.df

    ####### index given two cluster or loc

    ## mig 
    clone.dist.loc.majorCluster <- table(dat.plot[,c("majorCluster",clone.colname,"loc")])
    cls.mig.index.df <- data.frame()
    for(i in seq_len(dim(clone.dist.loc.majorCluster)[1])){
        dat.cls <- clone.dist.loc.majorCluster[i,,]
        i.name <- dimnames(clone.dist.loc.majorCluster)[["majorCluster"]][i]
        givenT.loc <- as.data.frame(t(combn(colnames(dat.cls),2)),stringsAsFactors=F)
        .index.mig <- t(apply(givenT.loc,1,function(x){
             dat.block <- dat.cls[,x]     
             dat.block.clone.index.mig <- apply(dat.block,1,entropy.empirical,unit="log2")
             dat.block.clone.index.mig[is.na(dat.block.clone.index.mig)] <- 0
             .index.a <- sum(dat.block.clone.index.mig*rowSums(dat.block)/sum(dat.block))
             ##dat.block[dat.block>1] <- 1
             dat.block.clone.index.migNorm <- apply(dat.block,1,function(x){ 
                                                        x[x>1] <- 1
                                                        entropy.empirical(x,unit="log2") })
             dat.block.clone.index.migNorm[is.na(dat.block.clone.index.migNorm)] <- 0
             .index.b <- sum(dat.block.clone.index.migNorm*rowSums(dat.block)/sum(dat.block))
             c(.index.a,.index.b)
                                     }))
        colnames(.index.mig) <- c("index.mig","index.migNorm")
        givenT.loc <- cbind(data.frame(majorCluster=rep(i.name,nrow(givenT.loc)),stringsAsFactors = F),
                            givenT.loc,.index.mig)
        cls.mig.index.df <- rbind(cls.mig.index.df,givenT.loc)
    }
    cls.mig.index.df$crossLoc <- sprintf("%s-%s",cls.mig.index.df$V1,cls.mig.index.df$V2)
    out.cls.mig.index.df <- dcast(cls.mig.index.df,majorCluster ~ crossLoc,value.var = "index.mig")
    write.table(out.cls.mig.index.df,sprintf("%s.cluster.2loc.index.txt",out.prefix),row.names = F,quote = F,sep = "\t")
    ret[["cls2.mig.index"]] <- out.cls.mig.index.df
    out.cls.migNorm.index.df <- dcast(cls.mig.index.df,majorCluster ~ crossLoc,value.var = "index.migNorm")
    write.table(out.cls.migNorm.index.df,sprintf("%s.cluster.2locNorm.index.txt",out.prefix),row.names = F,quote = F,sep = "\t")
    ret[["cls2.migNorm.index"]] <- out.cls.migNorm.index.df

    ## dev
    clone.dist.clust %>% head
    ##givenT.cls <- as.data.frame(t(combn(colnames(clone.dist.clust),2)),stringsAsFactors=F)
    givenT.cls <- expand.grid(colnames(clone.dist.clust),colnames(clone.dist.clust),stringsAsFactors = F)
    givenT.cls <- givenT.cls[givenT.cls[,1]!=givenT.cls[,2],]
    givenT.cls$index.dev <- apply(givenT.cls,1,function(x){
             dat.block <- clone.dist.clust[,x]     
             dat.block.clone.index <- apply(dat.block,1,entropy.empirical,unit="log2")
             dat.block.clone.index[is.na(dat.block.clone.index)] <- 0
             sum(dat.block.clone.index*rowSums(dat.block)/sum(dat.block))
                                     })
    out.cls.dev.index.df <- dcast(givenT.cls,Var2~Var1,value.var = "index.dev")
    write.table(out.cls.dev.index.df,sprintf("%s.cluster.2cls.index.txt",out.prefix),row.names = F,quote = F,sep = "\t")
    ret[["cls2.dev.index"]] <- out.cls.dev.index.df
    return(ret)
}

res.all <- calIndex(dat.plot,out.prefix)
patient.size <- dat.plot$patient %>% table
res.byPatient <- list()
for(pid in (dat.plot$patient %>% unique)){
    if(patient.size[pid]>3){
        res.byPatient[[pid]] <- calIndex(subset(dat.plot,patient==pid),sprintf("%s.%s",out.prefix,pid))
    }
}

### cluster.index
#for(.cname in c("avg.global.p","shannon.entropy",
#                "index.exp.scale","index.exp","index.exp.V2","index.exp.V3","index.exp.V4",
#                "index.mig.scale", "index.mig",
#                "index.dev.scale", "index.dev")){
for(.cname in c("index.exp.V4",
                "index.mig.scale", "index.mig", "index.migNorm",
                "index.dev.scale", "index.dev"))
{
    .dplot <- NULL
    for(i in seq_along(res.byPatient)){
        .pid <- names(res.byPatient)[i]
        .dblock <- res.byPatient[[.pid]][["cluster.index"]][,c("majorCluster",.cname)]
        colnames(.dblock) <- c("majorCluster",.pid)
        if(is.null(.dplot)){ 
            .dplot <- .dblock 
        }else{
            .dplot <- full_join(.dplot,.dblock)
        }
    }
    .dplot <- melt(.dplot)
    colnames(.dplot) <- c("majorCluster","patient",.cname)
    .dplot <- .dplot[!is.na(.dplot[,.cname]),]
    .dplot <- .dplot[order(.dplot$majorCluster),]
    print(.cname)
    print(summary(.dplot[,3]))
    print(which.max(.dplot[,3]) %>% .dplot[.,])
    my_comparisons <- list( c("CD8_C07-LAYN","CD8_C03-CX3CR1"), c("CD8_C07-LAYN","CD8_C04-GZMK"),
                           c("CD8_C03-CX3CR1", "CD8_C04-GZMK") )
    p <- ggboxplot(subset(.dplot,!grepl("^CD4",majorCluster,perl = T)),
                   x = "majorCluster", y = .cname, color = "steelblue", add = "point", outlier.colour=NULL) +
                   ##x = "majorCluster", y = .cname, color = "steelblue", add = "jitter", outlier.colour=NULL) +
                   ####x = "majorCluster", y = .cname, color = "steelblue", add = "dotplot", outlier.colour=NULL) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        theme(axis.title.y = element_text(size = rel(1.5)))
    if(.cname=="index.exp.V4"){
        p <- p + stat_compare_means(ref.group = "CD8_C01-LEF1", aes(label = ..p.signif..),method = "wilcox.test",hide.ns = T,label.y=0.61)
        p <- p + stat_compare_means(label.y = 0.5,label.x=1.5)
    }else{
        p <- p + stat_compare_means(comparisons = my_comparisons, method="wilcox.test")
        p <- p + stat_compare_means(label.y = NULL,label.x=1.5)
    }
    ggsave(sprintf("%s.cluster.%s.CD8.pdf",out.prefix,.cname),width = 8,height = 5)

    my_comparisons <- list( c("CD4_C08-IL23R", "CD4_C07-GZMK"),c("CD4_C09-CXCL13", "CD4_C07-GZMK"),
                           c("CD4_C12-CTLA4", "CD4_C07-GZMK"),c("CD4_C03-GNLY","CD4_C07-GZMK") )
    p <- ggboxplot(subset(.dplot,grepl("^CD4",majorCluster,perl = T)),
                   x = "majorCluster", y = .cname, color = "steelblue", add = "point", outlier.colour=NULL) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        theme(axis.title.y = element_text(size = rel(1.5)))
    if(.cname=="index.exp.V4"){
        p <- p + stat_compare_means(ref.group = "CD4_C01-CCR7", aes(label = ..p.signif..),method = "wilcox.test",hide.ns = T,label.y=0.10)
        p <- p + stat_compare_means(label.y = 0.15,label.x=1.5)
    }else{
        p <- p + stat_compare_means(comparisons = my_comparisons, method="wilcox.test")
        p <- p + stat_compare_means(label.y = NULL,label.x=1.5)
    }
    ggsave(sprintf("%s.cluster.%s.CD4.pdf",out.prefix,.cname),width = 8,height = 5)
}


### cls2.mig.index & cls2.dev.index
for(.cname in c("cls2.mig.index","cls2.dev.index"))
{

    .dplot <- NULL
    for(i in seq_along(res.byPatient)){
        .pid <- names(res.byPatient)[i]
        .dblock <- res.byPatient[[.pid]][[.cname]]
        .dblock <- melt(.dblock)
        colnames(.dblock) <- c("majorCluster","cls2",.cname)
        .dblock$patient <- .pid
        if(is.null(.dplot)){ 
            .dplot <- .dblock 
        }else{
            .dplot <- rbind(.dplot,.dblock)
        }
    }
    #.dplot[is.na(.dplot[,.cname]),.cname] <- 0
    .dplot <- .dplot[!is.na(.dplot[,.cname]),]
    .dplot <- .dplot[order(.dplot$majorCluster),]
    .dplot$cls2 <- factor(.dplot$cls2,levels=sort(unique(as.character(.dplot$cls2))))
    if(.cname=="cls2.mig.index"){
        .dplot$cls2 <- factor(.dplot$cls2,levels=c("P-N","P-T","N-T"))
    }

    ######################### CD8 ##########################
    #my_comparisons <- list( c("CD8_C04-CD160","CD8_C06-GZMK"), c("CD8_C07-HAVCR2","CD8_C06-GZMK"),
    #                       c("CD8_C03-CX3CR1", "CD8_C06-GZMK") )
    my_comparisons <- list( c("CD8_C03-CX3CR1", "CD8_C04-GZMK") )
    p <- ggboxplot(subset(.dplot,!grepl("^CD4",majorCluster,perl = T) & !grepl("^CD4",cls2,perl=T) ),
                   x = "majorCluster", facet.by ="cls2", color="cls2", y = .cname, add = "point",
                   outlier.colour=NULL,ylim=c(0,if(.cname=="cls2.dev.index" && max(.dplot[,.cname])<0.5) 0.5 else 1)) +
        theme(axis.text.x = element_text(angle = if(.cname=="cls2.dev.index") 90 else 45, hjust = 1)) +
        theme(axis.title.y = element_text(size = rel(1.5))) + 
        #stat_compare_means(comparisons = my_comparisons, method="wilcox.test")+ 
        stat_compare_means(label.y = if(.cname=="cls2.dev.index") 0.4 else 0.9,label.x=1.5)
    ggsave(sprintf("%s.cluster.%s.CD8.pdf",out.prefix,.cname),width = 10,height = if(.cname=="cls2.dev.index") 10 else 5)

    if(.cname=="cls2.dev.index"){
        for(.mcls in unique(.dplot$majorCluster)){
            if(!grepl("^CD4",.mcls,perl = T)){
                p <- ggboxplot(subset(.dplot,majorCluster==.mcls & !grepl("^CD4",cls2,perl=T) ),
                               x = "cls2", color=majorClusterColor[.mcls], y = .cname, add = "point",
                               outlier.colour=NULL,ylim=c(0,if(.cname=="cls2.dev.index" && max(.dplot[,.cname])<0.5) 0.5 else 1)) +
                    theme(axis.text.x = element_text(angle = if(.cname=="cls2.dev.index") 90 else 45, hjust = 1)) +
                    theme(axis.title.y = element_text(size = rel(1.5))) +
                    stat_compare_means(label.y = if(.cname=="cls2.dev.index") 0.4 else 0.9,label.x=1.5)
                ggsave(sprintf("%s.cluster.%s.CD8.byCluster.%s.pdf",out.prefix,.cname,.mcls),
                       width = 5,height = if(.cname=="cls2.dev.index") 5 else 3)
            }
        }
    }

    if(.cname=="cls2.mig.index"){
        p <- ggboxplot(subset(.dplot,!grepl("^CD4",majorCluster,perl = T) & !grepl("^CD4",cls2,perl=T) ),
                       x = "cls2", facet.by ="majorCluster", color="majorCluster", y = .cname, add = "point",
                       outlier.colour=NULL,ylim=c(0,if(.cname=="cls2.dev.index") 0.5 else 1)) +
            theme(axis.text.x = element_text(angle = if(.cname=="cls2.dev.index") 90 else 45, hjust = 1)) +
            theme(axis.title.y = element_text(size = rel(1.5))) + 
            #stat_compare_means(comparisons = my_comparisons, method="wilcox.test")+ 
            stat_compare_means(label.y = if(.cname=="cls2.dev.index") 0.4 else 0.9,label.x=1.5,aes(label=paste0("p = ", ..p.format..)))
        ggsave(sprintf("%s.cluster.%s.CD8.facetCluster.pdf",out.prefix,.cname),width = 6,height = if(.cname=="cls2.dev.index") 10 else 8)

        for(.mcls in unique(.dplot$majorCluster)){
            if(!grepl("^CD4",.mcls,perl = T)){
                p <- ggboxplot(subset(.dplot,majorCluster==.mcls), x = "cls2", color=majorClusterColor[.mcls], y = .cname, add = "point",
                               outlier.colour=NULL,ylim=c(0,if(.cname=="cls2.dev.index") 0.5 else 1)) +
                    theme(axis.text.x = element_text(angle = if(.cname=="cls2.dev.index") 90 else 45, hjust = 1)) +
                    theme(axis.title.y = element_text(size = rel(1.5))) + 
                    stat_compare_means(label.y = if(.cname=="cls2.dev.index") 0.4 else 0.9,label.x=1.5,aes(label=paste0("p = ", ..p.format..)))
                ggsave(sprintf("%s.cluster.%s.CD8.byCluster.%s.pdf",out.prefix,.cname,.mcls),
                       width = 3,height = if(.cname=="cls2.dev.index") 5 else 3)
            }
        }

    }

    ######################### CD4 ##########################
    
    #my_comparisons <- list( c("CD4_C08-IL23R", "CD4_C06-IL7R"),c("CD4_C07-CXCL13", "CD4_C06-IL7R"),
    #                       c("CD4_C11-CCR8", "CD4_C06-IL7R"),c("CD4_C03-CX3CR1","CD4_C06-IL7R") )
    my_comparisons <- NULL
    p <- ggboxplot(subset(.dplot,grepl("^CD4",majorCluster,perl = T) & !grepl("^(CD8|MAIT)",cls2,perl=T)),
                   x = "majorCluster", facet.by ="cls2", color="cls2", y = .cname, add = "point",
                   outlier.colour=NULL,ylim=c(0,if(.cname=="cls2.dev.index") 0.3 else 0.5)) +
        theme(axis.text.x = element_text(angle = if(.cname=="cls2.dev.index") 90 else 45, hjust = 1)) +
        theme(axis.title.y = element_text(size = rel(1.5))) + 
        stat_compare_means(label.y = if(.cname=="cls2.dev.index") 0.2 else 0.4,label.x=2.5) +
        stat_compare_means(ref.group = "CD4_C01-CCR7", aes(label = ..p.signif..),method = "wilcox.test",hide.ns = T)
    ggsave(sprintf("%s.cluster.%s.CD4.pdf",out.prefix,.cname),width = 10,height = if(.cname=="cls2.dev.index") 10 else 5)

    if(.cname=="cls2.dev.index"){
        for(.mcls in unique(.dplot$majorCluster)){
            #print(str(.dplot))
            #print(head(.dplot))
            #print(.dplot$cls2)
            if(grepl("^CD4",.mcls,perl = T)){
                p <- ggboxplot(subset(.dplot,majorCluster==.mcls & !grepl("^(CD8|MAIT)",cls2,perl=T)),
                               x = "cls2", color=majorClusterColor[.mcls], y = .cname, add = "point",
                               outlier.colour=NULL,ylim=c(0,if(.cname=="cls2.dev.index") 0.3 else 0.5)) +
                    theme(axis.text.x = element_text(angle = if(.cname=="cls2.dev.index") 90 else 45, hjust = 1)) +
                    theme(axis.title.y = element_text(size = rel(1.5))) +
                    stat_compare_means(label.y = if(.cname=="cls2.dev.index") 0.2 else 0.4,label.x=2.5) +
                    stat_compare_means(ref.group = "CD4_C01-CCR7", aes(label = ..p.signif..),method = "wilcox.test",hide.ns = T)
                ggsave(sprintf("%s.cluster.%s.CD4.byCluster.%s.pdf",out.prefix,.cname,.mcls),
                       width = 5,height = if(.cname=="cls2.dev.index") 5 else 3)
            }
        }
    }



    if(.cname=="cls2.mig.index"){

        p <- ggboxplot(subset(.dplot,grepl("^CD4",majorCluster,perl = T) & !grepl("^(CD8|MAIT)",cls2,perl=T)),
                       x = "cls2", facet.by ="majorCluster", color="majorCluster", y = .cname, add = "point",
                       outlier.colour=NULL,ylim=c(0,if(.cname=="cls2.dev.index") 0.3 else 0.5),ncol=3) +
            theme(axis.text.x = element_text(angle = if(.cname=="cls2.dev.index") 90 else 45, hjust = 1)) +
            theme(axis.title.y = element_text(size = rel(1.5))) + 
            stat_compare_means(label.y = if(.cname=="cls2.dev.index") 0.2 else 0.4,label.x=2.5,
                               aes(label=paste0("p = ", ..p.format..))) 
            #stat_compare_means(ref.group = "CD4_C01-CCR7", aes(label = ..p.signif..),method = "wilcox.test",hide.ns = T)
        ggsave(sprintf("%s.cluster.%s.CD4.facetCluster.pdf",out.prefix,.cname),
               width = 6,height = if(.cname=="cls2.dev.index") 10 else 10)

        for(.mcls in unique(.dplot$majorCluster)){
            if(grepl("^CD4",.mcls,perl = T)){
                p <- ggboxplot(subset(.dplot,majorCluster==.mcls & !grepl("^(CD8|MAIT)",cls2,perl=T)),
                       x = "cls2", color=majorClusterColor[.mcls], y = .cname, add = "point",
                       outlier.colour=NULL,ylim=c(0,if(.cname=="cls2.dev.index") 0.3 else 0.5),ncol=3) +
                    theme(axis.text.x = element_text(angle = if(.cname=="cls2.dev.index") 90 else 45, hjust = 1)) +
                    theme(axis.title.y = element_text(size = rel(1.5))) +
                    stat_compare_means(label.y = if(.cname=="cls2.dev.index") 0.2 else 0.4,label.x=2.5,
                               aes(label=paste0("p = ", ..p.format..))) 
                ggsave(sprintf("%s.cluster.%s.CD4.byCluster.%s.pdf",out.prefix,.cname,.mcls),
                       width = 3,height = if(.cname=="cls2.dev.index") 10 else 3)
            }
        }

    }
}

save.image(sprintf("%s.all.RData",out.prefix))

