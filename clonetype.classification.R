#!/usr/bin/env Rscript

#suppressPackageStartupMessages(library("argparse"))
#
#parser <- ArgumentParser()
#parser$add_argument("-i", "--inFile", type="character", required=TRUE, help="infile")
#parser$add_argument("-w", "--clusterColorFile", type="character", help="cluster color file [default %(default)s]")
#parser$add_argument("-a", "--cellTypeColorFile", default="/WPS1/zhenglt/work/TCR_chunhong/data/CellType.color", 
#                    type="character", help="cellTypeColorFile [default %(default)s]")
#parser$add_argument("-o", "--outputPrefix", type="character", required=TRUE, help="output prefix")
#parser$add_argument("-s", "--sample", type="character", default="SAMPLE", help="sample [default %(default)s]")
#parser$add_argument("-c", "--clone", type="character", default="C_strict", help="column name of clone data [default %(default)s]")
#parser$add_argument("-f", "--facet", type="character", default="loc", help="facet by [default %(default)s]")
#parser$add_argument("-d", "--pdfWidth", type="double", default=10, help="pdf width [default %(default)s]")
#parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="verbose mode [default %(default)s]")
#args <- parser$parse_args()
#print(args)
#
#
#in.file <- args$inFile
#out.prefix <- args$outputPrefix
#pdf.width <- args$pdfWidth
#sample.id <- args$sample
#clone.colname <- args$clone
#args.facet <- args$facet
#clusterColorFile <- args$clusterColorFile
#cellTypeColorFile <- args$cellTypeColorFile

in.file <- "/WPSnew/zhenglt/work/proj_zl/integrate0924/tcr/all.summary.cellV2.reassigneClonotype.methodChunhongV2.r.flt-zl.addInfo.mergeMKI67.txt"
out.prefix <- "/WPSnew/zhenglt/work/proj_zl/integrate0924/tcr/clonalExpansion/mig.dev.example.mergeMKI67/colonC.mig.dev.example.mergeMKI67"
pdf.width <- 10
sample.id <- "colonC"
clone.colname <- "C_strict"
args.facet <- "loc"
clusterColorFile <- "/WPSnew/zhenglt/work/proj_zl/integrate0924/cross.patient.b20170523/majorClusterColor.txt"
cellTypeColorFile <- "/WPS1/zhenglt/work/TCR_chunhong/data/CellType.color"

dir.create(dirname(out.prefix),showWarnings = F,recursive = T)

library("ggplot2")
library("gplots")
library("RColorBrewer")
library("plotrix")
library("reshape2")
library("igraph")

source("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")

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

#### cluster color
nMajor <- length(unique(c.data$itable$majorCluster))
myDesign <- c.data$itable
if(!is.null(clusterColorFile) && file.exists(clusterColorFile)){
    majorClusterColor <- read.SampleTypeColor(clusterColorFile)
}else{
    majorClusterColor <- structure(colorRampPalette(brewer.pal(nMajor,"Paired"))(nMajor),
                              names=unique(myDesign$majorCluster))
}
names(majorClusterColor) <- gsub(pattern = "Cluster",replacement = "C",x = names(majorClusterColor))

if(!is.null(cellTypeColorFile) && file.exists(cellTypeColorFile)){
    sampleTypeColor <- read.SampleTypeColor(cellTypeColorFile)
    #sampleTypeColor <- sampleTypeColor[names(sampleTypeColor) %in% unique(as.character(myDesign[,"sampleType"]))]
}else{
    sampleTypeColor <- auto.colSet(n=length(unique(myDesign[,"sampleType"])),name="Paired")
    names(sampleTypeColor) <- unique(myDesign[,"sampleType"])
}

locColor <- structure(c("#E41A1C","#A65628","#377EB8"),names=c("P","N","T"))


dat.clone.mcl.dist <- table(dat.plot[,c(clone.colname,"loc.majorCluster")])

dat.clone.mcl.dist.CD8.CX3CR1 <- dat.clone.mcl.dist[,grepl("CD8_C03-CX3CR1",colnames(dat.clone.mcl.dist),perl = T)]
f.clone <- apply(dat.clone.mcl.dist.CD8.CX3CR1,1,function(x){ sum(x>0)>=2 })
dat.clone.mcl.dist.CD8.CX3CR1 <- dat.clone.mcl.dist.CD8.CX3CR1[f.clone,]
dat.plot.CD8.CX3CR1 <- subset(dat.plot, C_strict %in% rownames(dat.clone.mcl.dist.CD8.CX3CR1) & majorCluster=="CD8_C03-CX3CR1")

dat.clone.mcl.dist.T <- dat.clone.mcl.dist[,grepl("^T:",colnames(dat.clone.mcl.dist),perl = T)]
dat.clone.mcl.dist.T.C4 <- dat.clone.mcl.dist[,grepl("(T:CD8_C03-CX3CR1|T:CD8_C06-GZMK|T:CD8_C07-HAVCR2|T:CD8_C05-CD6)",colnames(dat.clone.mcl.dist),perl = T)]
f.clone <- apply(dat.clone.mcl.dist.T.C4,1,function(x){ sum(x>0)>=2 })
dat.clone.mcl.dist.T.C4 <- dat.clone.mcl.dist.T.C4[f.clone,]
dat.plot.T.C4 <- subset(dat.plot, loc=="T" & 
                        C_strict %in% rownames(dat.clone.mcl.dist.T.C4) & 
                        grepl("CD8_C03-CX3CR1|CD8_C06-GZMK|CD8_C07-HAVCR2|CD8_C05-CD6",majorCluster,perl = T))

## for tcr sharing 
my.plot.sharing <- function(dat.plot,out.prefix,byLoc=T,TH.n=1,use.diag=F)
{
    clone.size <- sort(table(dat.plot[,clone.colname]),decreasing = T)
    ##f.clonal <- clone.size[dat.plot[,clone.colname]]>1
    if(byLoc){
        dat.plot.clonal.mcl <- table(dat.plot[,c(clone.colname,"loc.majorCluster")])
        vnames <- table(dat.plot[,"loc.majorCluster"])
    }else{
        dat.plot.clonal.mcl <- table(dat.plot[,c(clone.colname,"majorCluster")])
        vnames <- table(dat.plot[,"majorCluster"])
    }

    vnames <- data.frame(name=names(vnames),ncells=unclass(vnames),stringsAsFactors = F)
    vnames$size <- 2*(vnames$ncells+1)^0.5
    vnames$id <- sprintf("%02d",seq_len(nrow(vnames)))
    if(byLoc){
        vnames$majorCluster <- gsub("^..","",vnames$name)
        vnames$loc <- gsub(":.+","",vnames$name)
        ### if "majorCluster"==MAIT, CD8 MAIT
        m <- regexpr("^(...)", vnames$majorCluster,perl = T)
        .t <- regmatches(vnames$majorCluster, m)
        .t[.t=="MAI"] <- "CD8"
        vnames$locColor <- sampleTypeColor[sprintf("%s.%s",.t,vnames$loc)]
        #vnames$locColor <- locColor[vnames$loc]
    }else{
        vnames$majorCluster <- vnames$name
    }
    vnames$display.label <- sprintf("%s(%s,%s)",vnames$majorCluster,vnames$id,vnames$ncells)
    vnames$majorClusterColor <- majorClusterColor[vnames$majorCluster]
    write.table(vnames,sprintf("%s.n%d.useDiag%s.node.attribute.txt",out.prefix,TH.n,if(use.diag) "T" else "F"),
                row.names = F,quote = F,sep = "\t")
    ##### edges
    f.dat.plot.clonal.mcl.multiple <- apply(dat.plot.clonal.mcl,1,function(x){ sum(x>0)>1 })
    f.dat.plot.clonal.mcl.single <- apply(dat.plot.clonal.mcl,1,function(x){ sum(x>0)==1 })
    dat.plot.clonal.mcl.multiple <- dat.plot.clonal.mcl[f.dat.plot.clonal.mcl.multiple,]
    dat.plot.clonal.mcl.single <- dat.plot.clonal.mcl[f.dat.plot.clonal.mcl.single,]

    ###pnames <- expand.grid(colnames(dat.plot.clonal.mcl.multiple),colnames(dat.plot.clonal.mcl.multiple))
    pnames <- as.data.frame(t(combn(colnames(dat.plot.clonal.mcl.multiple),2)),stringsAsFactors=F)
    pnames <- rbind(pnames,
                    cbind(colnames(dat.plot.clonal.mcl.multiple),
                          colnames(dat.plot.clonal.mcl.multiple)))
    colnames(pnames) <- c("Var1","Var2")
    ##pnames <- pnames[apply(pnames,1,function(x){ x[1]!=x[2] }),]
    pair.info <- t(apply(pnames,1,function(x){
                             if(x[1]==x[2]){
                                 f <- dat.plot.clonal.mcl.single[,x[1]]>=TH.n
                                 n.clone <- sum(f)
                                 n.left <- sum(dat.plot.clonal.mcl.single[f,x[1]])
                                 #n.right <- sum(dat.plot.clonal.mcl.single[f,x[2]])
                                 n.right <- 0
                             }else{
                                 f <- (dat.plot.clonal.mcl.multiple[,x[1]]>=TH.n & dat.plot.clonal.mcl.multiple[,x[2]]>=TH.n)
                                 n.clone <- sum(f)
                                 n.left <- sum(dat.plot.clonal.mcl.multiple[f,x[1]])
                                 n.right <- sum(dat.plot.clonal.mcl.multiple[f,x[2]])
                             }
                             structure(c(n.clone,n.left,n.right),names=c("n.clone","n.left","n.right"))
                        }))
    pnames <- cbind(pnames,pair.info)
    pnames$label.cex <- 0.4
    pnames$width <- log2(pnames$n.clone+1)
    pnames$label <- sprintf("%s(%s:%s,%s:%s)",pnames$n.clone,
                            vnames[pnames$Var1,"id"],pnames$n.left,
                            vnames[pnames$Var2,"id"],pnames$n.right)
    pnames$label.slim <- pnames$label
    pnames$label.slim[pnames$n.clone<=1] <- ""
    pnames$ID <- sprintf("%s-%s",pnames$Var1,pnames$Var2)

    ### matrix of shared clone (numClone)
#    for(.stype in c("CD4","CD8|MAIT")){
#        f.stype <- grepl(pattern = sprintf("%s",.stype),pnames$Var1,perl = T) & grepl(pattern = sprintf("%s",.stype),pnames$Var2,perl = T)
#        pnames.stype <- pnames[f.stype,]
#        pnames.mtx <- dcast(pnames.stype[,c("Var1","Var2","n.clone")],Var2~Var1,value.var="n.clone")
#        rownames(pnames.mtx) <- pnames.mtx[,"Var2"]
#        pnames.mtx <- as.matrix(pnames.mtx[,-1])
#        diag.v <- apply(dat.plot.clonal.mcl.single[,grepl(.stype,colnames(dat.plot.clonal.mcl.single),perl = T)],2,function(x){ sum(x>0)  })
#        print("check colnames and rownames:")
#        print(all(colnames(pnames.mtx)==rownames(pnames.mtx)))
#        print(all(colnames(pnames.mtx)==names(diag.v)))
#        diag(pnames.mtx) <- diag.v
#        pnames.mtx.df <- data.frame(ID=rownames(pnames.mtx),stringsAsFactors = F)
#        pnames.mtx.df <- cbind(pnames.mtx.df,pnames.mtx)
#        write.table(pnames.mtx.df,sprintf("%s.n%d.useDiag%s.matrix.%s.txt",out.prefix,TH.n,if(use.diag) "T" else "F",gsub("\\|","-",.stype)),
#                    row.names = F,quote = F,sep = "\t")
#
#        ### matrix of shared clone (numCell)
#        pnames.numCell.mtx <- dcast(pnames.stype[,c("Var1","Var2","n.right")],Var2~Var1,value.var="n.right")
#        rownames(pnames.numCell.mtx) <- pnames.numCell.mtx[,"Var2"]
#        pnames.numCell.mtx <- as.matrix(pnames.numCell.mtx[,-1])
#        for(i in seq_len(nrow(pnames.stype))){
#            if(is.na(pnames.numCell.mtx[pnames.stype[i,"Var1"], pnames.stype[i,"Var2"]])){
#                pnames.numCell.mtx[pnames.stype[i,"Var1"], pnames.stype[i,"Var2"]] <- pnames.stype[i,"n.left"]
#            }
#        }
#        diag.v.numCell <- apply(dat.plot.clonal.mcl.single[,grepl(.stype,colnames(dat.plot.clonal.mcl.single),perl = T)],2,sum)
#        print("check colnames and rownames:")
#        print(all(colnames(pnames.numCell.mtx)==rownames(pnames.numCell.mtx)))
#        print(all(colnames(pnames.numCell.mtx)==names(diag.v.numCell)))
#        diag(pnames.numCell.mtx) <- diag.v.numCell
#        pnames.numCell.mtx.df <- data.frame(ID=rownames(pnames.numCell.mtx),stringsAsFactors = F)
#        pnames.numCell.mtx.df <- cbind(pnames.numCell.mtx.df,pnames.numCell.mtx)
#        write.table(pnames.numCell.mtx.df,sprintf("%s.n%d.useDiag%s.matrix.%s.cellNum.txt",out.prefix,TH.n,if(use.diag) "T" else "F",gsub("\\|","-",.stype)),
#                    row.names = F,quote = F,sep = "\t")
#    }

    if(use.diag){
        f.pnames <- pnames$n.clone>0
        v.in.pnames <- unique(c(pnames[f.pnames,"Var1"],pnames[f.pnames,"Var2"]))
        dat.plot.clonal.mcl.sharingGraph <- graph_from_data_frame(d = pnames[f.pnames,,drop=F],
                                                                  vertices = vnames[v.in.pnames,,drop=F], directed = F)
    }else{
        f.pnames <- pnames$n.clone>0 & pnames$Var1!=pnames$Var2
        v.in.pnames <- unique(c(pnames[f.pnames,"Var1"],pnames[f.pnames,"Var2"]))
        dat.plot.clonal.mcl.sharingGraph <- graph_from_data_frame(d=pnames[f.pnames,,drop=F],
                                                                  vertices = vnames[v.in.pnames,,drop=F], directed = F)
    }
    write.table(pnames,sprintf("%s.n%d.useDiag%s.edge.attribute.txt",out.prefix,TH.n,if(use.diag) "T" else "F"),
                row.names = F,quote = F,sep = "\t")

    pdf(file = sprintf("%s.n%d.useDiag%s.igraph.pdf",out.prefix,TH.n,if(use.diag) "T" else "F"),width = 10,height = 10)
    plot(dat.plot.clonal.mcl.sharingGraph,margin=c(0.1,0.1,0.1,0.2),layout=layout_with_kk)
    dev.off()
    write_graph(dat.plot.clonal.mcl.sharingGraph,
                file = sprintf("%s.n%d.useDiag%s.igraph.graphml",
                               out.prefix,TH.n,if(use.diag) "T" else "F"),
                format = "graphml")

}

my.plot.sharing(dat.plot.CD8.CX3CR1,out.prefix=sprintf("%s.mig.CD8.CX3CR1",out.prefix),byLoc=T,TH.n=1,use.diag=F)
my.plot.sharing(dat.plot.CD8.CX3CR1,out.prefix=sprintf("%s.mig.CD8.CX3CR1",out.prefix),byLoc=T,TH.n=2,use.diag=F)

my.plot.sharing(dat.plot.T.C4,out.prefix=sprintf("%s.mig.T.C4",out.prefix),byLoc=T,TH.n=1,use.diag=F)
my.plot.sharing(dat.plot.T.C4,out.prefix=sprintf("%s.mig.T.C4",out.prefix),byLoc=T,TH.n=2,use.diag=F)


make.clone.cluster.network <- function(dat.clone.cluster.dist,out.prefix,byLoc=T,Th.link=1)
{
    dat.clone.size <- apply(dat.clone.cluster.dist,1,sum)
    dat.cluster.size <- apply(dat.clone.cluster.dist,2,sum)
    node.clone.df <- data.frame(name=names(dat.clone.size),num.cell=dat.clone.size,ntype="clone",
                                majorCluster="NA",loc="NA",locColor="NA",majorClusterColor="NA",
                                display.label="")

    node.cluster.df <- data.frame(name=names(dat.cluster.size),num.cell=dat.cluster.size,ntype="cluster")
    if(byLoc){
        node.cluster.df$majorCluster <- gsub("^..","",node.cluster.df$name)
        node.cluster.df$loc <- gsub("\\..+","",node.cluster.df$name)
        ### if "majorCluster"==MAIT, CD8 MAIT
        m <- regexpr("^(...)", node.cluster.df$majorCluster,perl = T)
        .t <- regmatches(node.cluster.df$majorCluster, m)
        .t[.t=="MAI"] <- "CD8"
        node.cluster.df$locColor <- sampleTypeColor[sprintf("%s.%s",.t,node.cluster.df$loc)]
    }else{
        node.cluster.df$majorCluster <- node.cluster.df$name
    }
    node.cluster.df$majorClusterColor <- majorClusterColor[node.cluster.df$majorCluster]
    node.cluster.df$display.label=node.cluster.df$majorCluster

    tmp.dat.clone.cluster.dist <- dat.clone.cluster.dist
    tmp.dat.clone.cluster.dist[tmp.dat.clone.cluster.dist<Th.link] <- 0
    f.row <- apply(tmp.dat.clone.cluster.dist,1,function(x){ sum(x>0)>=2 })
    num.link.cluster <- apply(tmp.dat.clone.cluster.dist,2,function(x){ sum(x>0) })
    num.link.clone <- apply(tmp.dat.clone.cluster.dist,1,function(x){ sum(x>0) })
    node.cluster.df$num.link <- num.link.cluster[rownames(node.cluster.df)]
    node.clone.df$num.link <- num.link.clone[rownames(node.clone.df)]

    node.df <- rbind(node.clone.df,node.cluster.df)
    node.df$size <- 2*(node.df$num.cell+1)^0.5
    ##node.df$size <- log2(node.df$num.cell+1)
    write.table(node.df,sprintf("%s.clone.cluster.node.attribute.ThLink%d.txt",out.prefix,Th.link), row.names = F,quote = F,sep = "\t")

    edge.df <- melt(tmp.dat.clone.cluster.dist[f.row,,drop=F])
    colnames(edge.df) <- c("clone","cluster","num.cell")
    edge.df[,1] <- as.character(edge.df[,1])
    edge.df[,2] <- as.character(edge.df[,2])
    edge.df$ID <- sprintf("%s(-)%s",edge.df[,1],edge.df[,2])
    ##edge.df <- subset(edge.df,num.cell>=Th.link)
    edge.df$label.cex <- 0.4
    edge.df$width <- log2(edge.df$num.cell+1)
    u.names <- unique(c(edge.df[,1],edge.df[,2]))
    write.table(edge.df,sprintf("%s.clone.cluster.edge.attribute.ThLink%d.txt",out.prefix,Th.link), row.names = F,quote = F,sep = "\t")
    edge.df.iGraph <- graph_from_data_frame(d = edge.df,vertices = node.df[u.names,], directed = F)

    pdf(file = sprintf("%s.clone.cluster.igraph.ThLink%d.pdf",out.prefix,Th.link),width = 10,height = 10)
    plot(edge.df.iGraph,margin=c(0.1,0.1,0.1,0.2),layout=layout_with_kk)
    dev.off()
    write_graph(edge.df.iGraph, file = sprintf("%s.clone.cluster.igraph.ThLink%d.graphml", out.prefix,Th.link), format = "graphml")
}


make.clone.cluster.network(dat.clone.mcl.dist.CD8.CX3CR1,
                           out.prefix=sprintf("%s.clonal.clone-cluster.network.CD8.CX3CR1",out.prefix),byLoc=T,Th.link=1)
make.clone.cluster.network(dat.clone.mcl.dist.CD8.CX3CR1,
                           out.prefix=sprintf("%s.clonal.clone-cluster.network.CD8.CX3CR1",out.prefix),byLoc=T,Th.link=2)
make.clone.cluster.network(dat.clone.mcl.dist.T.C4,
                           out.prefix=sprintf("%s.clonal.clone-cluster.network.T.C4",out.prefix),byLoc=T,Th.link=1)
make.clone.cluster.network(dat.clone.mcl.dist.T.C4,
                           out.prefix=sprintf("%s.clonal.clone-cluster.network.T.C4",out.prefix),byLoc=T,Th.link=2)



############################# END ###############


