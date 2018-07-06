#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-i", "--inFile", type="character", required=TRUE, help="infile")
parser$add_argument("-w", "--clusterColorFile", type="character", help="cluster color file [default %(default)s]")
parser$add_argument("-a", "--cellTypeColorFile", default="/WPS1/zhenglt/work/TCR_chunhong/data/CellType.color", 
                    type="character", help="cellTypeColorFile [default %(default)s]")
parser$add_argument("-o", "--outputPrefix", type="character", required=TRUE, help="output prefix")
parser$add_argument("-s", "--sample", type="character", default="SAMPLE", help="sample [default %(default)s]")
parser$add_argument("-c", "--clone", type="character", default="C_strict", help="column name of clone data [default %(default)s]")
parser$add_argument("-f", "--facet", type="character", default="loc", help="facet by [default %(default)s]")
parser$add_argument("-d", "--pdfWidth", type="double", default=10, help="pdf width [default %(default)s]")
parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="verbose mode [default %(default)s]")
args <- parser$parse_args()
print(args)


in.file <- args$inFile
out.prefix <- args$outputPrefix
pdf.width <- args$pdfWidth
sample.id <- args$sample
clone.colname <- args$clone
args.facet <- args$facet
clusterColorFile <- args$clusterColorFile
cellTypeColorFile <- args$cellTypeColorFile

#in.file <- "/WPSnew/zhenglt/work/panC/data/merge.liverC.lungC.colonC/tracer/merge.v20180423.hclust.onlyPNT/all.summary.cellV2.reassigneClonotype.methodChunhongV2.r.flt.addMajorCluster.txt"
#out.prefix <- "/WPSnew/zhenglt/work/panC/ana/merge.liverC.lungC.colonC/tcr.sharing/clonalExpansion/all/panC.byLoc"
#pdf.width <- 10
#sample.id <- "panC"
#clone.colname <- "C_strict"
#args.facet <- "loc"
#clusterColorFile <- NULL
#cellTypeColorFile <- "/WPS1/zhenglt/work/TCR_chunhong/data/CellType.color"
#clusterColorFile <- "/WPSnew/zhenglt/work/proj_zl/integrate0924/cross.patient.b20170523/majorClusterColor.txt"

dir.create(dirname(out.prefix),showWarnings = F,recursive = T)

library("ggplot2")
library("gplots")
library("RColorBrewer")
library("plotrix")
library("reshape2")
library("igraph")

if(file.exists("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")){
    source("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")
}else{
    source("/WPSnew/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")
}

##c.data <- read.clonotype(in.file,"C_strict")
c.data <- read.clonotype(in.file,clone.colname)

#dat.plot <- c.data$itable[,c("patient","sampleType","majorCluster","C_strict")]
dat.plot <- c.data$itable[!is.na(c.data$itable$majorCluster) & c.data$itable$majorCluster!="uncharacterized",
                          c("patient","sampleType","stype","majorCluster",clone.colname)]
clone.freq <- sort(table(dat.plot[,clone.colname]),decreasing = T)
dat.plot$ctype <- c.data$ctypeII[rownames(dat.plot)]
m <- regexpr("^(.)", dat.plot$sampleType,perl = T)
dat.plot$loc <- regmatches(dat.plot$sampleType, m)
dat.plot$loc <- factor(dat.plot$loc,levels = c("P","N","T","L"))
dat.plot$majorCluster <- gsub(pattern = "Cluster",replacement = "C",x = dat.plot$majorCluster)
dat.plot$loc.majorCluster <- sprintf("%s:%s",dat.plot$loc,dat.plot$majorCluster)

#ftable(dat.plot[,c("patient","ctype","sampleType")])
#table(dat.plot[,c("ctype","sampleType","patient")])
#("count")
#("propotion")
sink(file = sprintf("%s.stat.txt",out.prefix))
cat("sampleType\n")
addmargins(table(dat.plot[,c("ctype","sampleType","patient")]))
cat("majorCluster\n")
addmargins(table(dat.plot[,c("ctype","majorCluster","patient")]))
#cat("freq\n")
#addmargins(prop.table(table(dat.plot[,c("ctype","sampleType","patient")])),margin = c(1))
sink()

#### TCR diversity
dat.plot.for.diversity <- table(dat.plot[,c(clone.colname,"majorCluster")])
### culumulative
###dat.plot.diversity.stat <- t(apply(dat.plot.for.diversity,2,function(x){ c(sum(x>0),sum(x>1),sum(x>2),sum(x>3),
###                                                                           sum(x),sum(x[x>1]),sum(x[x>2]),sum(x[x>3]) )  }))
### point density
dat.plot.diversity.stat <- t(apply(dat.plot.for.diversity,2,function(x){ c(sum(x==1),sum(x==2),sum(x==3),sum(x>=4),
                                                                           sum(x[x==1]),sum(x[x==2]),sum(x[x==3]),sum(x[x>=4]) )  }))
colnames(dat.plot.diversity.stat) <- c("nclones.k1","nclones.k2","nclones.k3","nclones.k4","ncells.k1","ncells.k2","ncells.k3","ncells.k4")
dat.plot.diversity.stat <- rbind(dat.plot.diversity.stat,apply(dat.plot.diversity.stat,2,sum))
rownames(dat.plot.diversity.stat) <- c(rownames(dat.plot.diversity.stat)[1:(nrow(dat.plot.diversity.stat)-1)],"Sum")
dat.plot.diversity.stat.df <- data.frame(majorCluster=rownames(dat.plot.diversity.stat),stringsAsFactors = F)
dat.plot.diversity.stat.df <- cbind(dat.plot.diversity.stat.df,dat.plot.diversity.stat)
write.table(dat.plot.diversity.stat.df,sprintf("%s.diversity.txt",out.prefix), row.names = F,quote = F,sep = "\t")

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

## for tcr sharing 
my.plot.sharing <- function(dat.plot,byLoc=T,TH.n=1,use.diag=F){
    clone.size <- sort(table(dat.plot[,clone.colname]),decreasing = T)
    f.clonal <- clone.size[dat.plot[,clone.colname]]>1
    if(byLoc){
        dat.plot.clonal.mcl <- table(dat.plot[f.clonal,c(clone.colname,"loc.majorCluster")])
        vnames <- table(dat.plot[,"loc.majorCluster"])
    }else{
        dat.plot.clonal.mcl <- table(dat.plot[f.clonal,c(clone.colname,"majorCluster")])
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
    for(.stype in c("CD4","CD8|MAIT")){
        f.stype <- grepl(pattern = sprintf("%s",.stype),pnames$Var1,perl = T) & grepl(pattern = sprintf("%s",.stype),pnames$Var2,perl = T)
        pnames.stype <- pnames[f.stype,]
        pnames.mtx <- dcast(pnames.stype[,c("Var1","Var2","n.clone")],Var2~Var1,value.var="n.clone")
        rownames(pnames.mtx) <- pnames.mtx[,"Var2"]
        pnames.mtx <- as.matrix(pnames.mtx[,-1])
        diag.v <- apply(dat.plot.clonal.mcl.single[,grepl(.stype,colnames(dat.plot.clonal.mcl.single),perl = T)],2,function(x){ sum(x>0)  })
        print("check colnames and rownames:")
        print(all(colnames(pnames.mtx)==rownames(pnames.mtx)))
        print(all(colnames(pnames.mtx)==names(diag.v)))
        diag(pnames.mtx) <- diag.v
        pnames.mtx.df <- data.frame(ID=rownames(pnames.mtx),stringsAsFactors = F)
        pnames.mtx.df <- cbind(pnames.mtx.df,pnames.mtx)
        write.table(pnames.mtx.df,sprintf("%s.n%d.useDiag%s.matrix.%s.txt",out.prefix,TH.n,if(use.diag) "T" else "F",gsub("\\|","-",.stype)),
                    row.names = F,quote = F,sep = "\t")

        ### matrix of shared clone (numCell)
        pnames.numCell.mtx <- dcast(pnames.stype[,c("Var1","Var2","n.right")],Var2~Var1,value.var="n.right")
        rownames(pnames.numCell.mtx) <- pnames.numCell.mtx[,"Var2"]
        pnames.numCell.mtx <- as.matrix(pnames.numCell.mtx[,-1])
        for(i in seq_len(nrow(pnames.stype))){
            if(is.na(pnames.numCell.mtx[pnames.stype[i,"Var1"], pnames.stype[i,"Var2"]])){
                pnames.numCell.mtx[pnames.stype[i,"Var1"], pnames.stype[i,"Var2"]] <- pnames.stype[i,"n.left"]
            }
        }
        diag.v.numCell <- apply(dat.plot.clonal.mcl.single[,grepl(.stype,colnames(dat.plot.clonal.mcl.single),perl = T)],2,sum)
        print("check colnames and rownames:")
        print(all(colnames(pnames.numCell.mtx)==rownames(pnames.numCell.mtx)))
        print(all(colnames(pnames.numCell.mtx)==names(diag.v.numCell)))
        diag(pnames.numCell.mtx) <- diag.v.numCell
        pnames.numCell.mtx.df <- data.frame(ID=rownames(pnames.numCell.mtx),stringsAsFactors = F)
        pnames.numCell.mtx.df <- cbind(pnames.numCell.mtx.df,pnames.numCell.mtx)
        write.table(pnames.numCell.mtx.df,sprintf("%s.n%d.useDiag%s.matrix.%s.cellNum.txt",out.prefix,TH.n,if(use.diag) "T" else "F",gsub("\\|","-",.stype)),
                    row.names = F,quote = F,sep = "\t")
    }

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
#my.plot.sharing(dat.plot,use.diag = F,TH.n = 2)
#my.plot.sharing(dat.plot,use.diag = T,TH.n = 2)
#my.plot.sharing(dat.plot,use.diag = F,TH.n = 1)
#my.plot.sharing(dat.plot,use.diag = T,TH.n = 1)
###my.plot.sharing(dat.plot,useClone=F)

#q()

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
        node.cluster.df$loc <- NA
        node.cluster.df$locColor <- NA
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

#####
## select CD8 data
f.CD8 <- grepl(pattern = "^(CD8|MAIT)",x = dat.plot$majorCluster,perl = T)
#f.CD8 <- grepl(pattern = "^(CD8)",x = dat.plot$majorCluster,perl = T)
dat.plot.CD8 <- dat.plot[f.CD8,]
f.CD4 <- grepl(pattern = "^(CD4)",x = dat.plot$majorCluster,perl = T)
dat.plot.CD4 <- dat.plot[f.CD4,]

#f.clonal <- !grepl(":1$",dat.plot[,clone.colname],perl=T)
f.clonal  <- dat.plot$ctype %in% c("N2","N3","N4")
dat.plot.clonal.CD8 <- dat.plot[f.clonal & f.CD8,]
if(args.facet=="loc"){
    dat.plot.clonal.CD8$loc.mcl <- paste(dat.plot.clonal.CD8$loc,dat.plot.clonal.CD8$majorCluster,sep = ".")
}else{
    dat.plot.clonal.CD8$loc.mcl <- paste(dat.plot.clonal.CD8$majorCluster,dat.plot.clonal.CD8$loc,sep = ".")
}
dat.plot.clonal.CD4 <- dat.plot[f.clonal & f.CD4,]
if(args.facet=="loc"){
    dat.plot.clonal.CD4$loc.mcl <- paste(dat.plot.clonal.CD4$loc,dat.plot.clonal.CD4$majorCluster,sep = ".")
}else{
    dat.plot.clonal.CD4$loc.mcl <- paste(dat.plot.clonal.CD4$majorCluster,dat.plot.clonal.CD4$loc,sep = ".")
}

##### clones distribute across cluster or loc.cluster
## CD8
dat.plot.clonal.CD8.cloneDist <- table(dat.plot.clonal.CD8[,c(clone.colname,"loc.mcl")])
if(args.facet=="loc"){
}else{
    print(colnames(dat.plot.clonal.CD8.cloneDist))
    dat.plot.clonal.CD8.cloneDist <- dat.plot.clonal.CD8.cloneDist[,sort(colnames(dat.plot.clonal.CD8.cloneDist))]
}
dat.plot.clonal.CD8.cloneDist.df <- data.frame(clone.id=rownames(dat.plot.clonal.CD8.cloneDist),stringsAsFactors = F)
dat.plot.clonal.CD8.cloneDist.df <- cbind(dat.plot.clonal.CD8.cloneDist.df,as.data.frame.matrix(dat.plot.clonal.CD8.cloneDist))
write.table(dat.plot.clonal.CD8.cloneDist.df,sprintf("%s.CD8.mcls.clone.loc.txt",out.prefix),row.names = F,sep = "\t",quote = F)

make.clone.cluster.network(dat.plot.clonal.CD8.cloneDist,
                           out.prefix=sprintf("%s.CD8.clonal.clone-cluster.network",out.prefix),byLoc=T,Th.link=1)
make.clone.cluster.network(dat.plot.clonal.CD8.cloneDist,
                           out.prefix=sprintf("%s.CD8.clonal.clone-cluster.network",out.prefix),byLoc=T,Th.link=2)

######### not by location
dat.plot.clonal.CD8.cloneDist.notByLoc <- table(dat.plot.clonal.CD8[,c(clone.colname,"majorCluster")])
make.clone.cluster.network(dat.plot.clonal.CD8.cloneDist.notByLoc,
                           out.prefix=sprintf("%s.CD8.clonal.clone-cluster.network.notByLoc",out.prefix),byLoc=F,Th.link=1)
make.clone.cluster.network(dat.plot.clonal.CD8.cloneDist.notByLoc,
                           out.prefix=sprintf("%s.CD8.clonal.clone-cluster.network.notByLoc",out.prefix),byLoc=F,Th.link=2)
#########


dat.plot.clonal.CD8.mcls.clone <- table(dat.plot.clonal.CD8[,c(clone.colname,"majorCluster")])
dat.plot.clonal.CD8.mcls.clone.df <- data.frame(clone.id=rownames(dat.plot.clonal.CD8.mcls.clone),stringsAsFactors = F)                           
dat.plot.clonal.CD8.mcls.clone.df <- cbind(dat.plot.clonal.CD8.mcls.clone.df,as.data.frame.matrix(dat.plot.clonal.CD8.mcls.clone))
dat.plot.clonal.CD8.mcls.clone.df$num.cell <- apply(dat.plot.clonal.CD8.mcls.clone,1,sum)
dat.plot.clonal.CD8.mcls.clone.df$num.cluster <- apply(dat.plot.clonal.CD8.mcls.clone,1,function(x){ sum(x>0) })
write.table(dat.plot.clonal.CD8.mcls.clone.df,sprintf("%s.CD8.mcls.clone.txt",out.prefix),row.names = F,sep = "\t",quote = F)

## CD4
dat.plot.clonal.CD4.cloneDist <- table(dat.plot.clonal.CD4[,c(clone.colname,"loc.mcl")])
if(args.facet=="loc"){
}else{
    print(colnames(dat.plot.clonal.CD4.cloneDist))
    dat.plot.clonal.CD4.cloneDist <- dat.plot.clonal.CD4.cloneDist[,sort(colnames(dat.plot.clonal.CD4.cloneDist))]
}
dat.plot.clonal.CD4.cloneDist.df <- data.frame(clone.id=rownames(dat.plot.clonal.CD4.cloneDist),stringsAsFactors = F)
dat.plot.clonal.CD4.cloneDist.df <- cbind(dat.plot.clonal.CD4.cloneDist.df,as.data.frame.matrix(dat.plot.clonal.CD4.cloneDist))
write.table(dat.plot.clonal.CD4.cloneDist.df,sprintf("%s.CD4.mcls.clone.loc.txt",out.prefix),row.names = F,sep = "\t",quote = F)
dat.plot.clonal.CD4.mcls.clone <- table(dat.plot.clonal.CD4[,c(clone.colname,"majorCluster")])
dat.plot.clonal.CD4.mcls.clone.df <- data.frame(clone.id=rownames(dat.plot.clonal.CD4.mcls.clone),stringsAsFactors = F)                           
dat.plot.clonal.CD4.mcls.clone.df <- cbind(dat.plot.clonal.CD4.mcls.clone.df,as.data.frame.matrix(dat.plot.clonal.CD4.mcls.clone))
dat.plot.clonal.CD4.mcls.clone.df$num.cell <- apply(dat.plot.clonal.CD4.mcls.clone,1,sum)
dat.plot.clonal.CD4.mcls.clone.df$num.cluster <- apply(dat.plot.clonal.CD4.mcls.clone,1,function(x){ sum(x>0) })
write.table(dat.plot.clonal.CD4.mcls.clone.df,sprintf("%s.CD4.mcls.clone.txt",out.prefix),row.names = F,sep = "\t",quote = F)

make.clone.cluster.network(dat.plot.clonal.CD4.cloneDist,
                           out.prefix=sprintf("%s.CD4.clonal.clone-cluster.network",out.prefix),byLoc=T,Th.link=1)
make.clone.cluster.network(dat.plot.clonal.CD4.cloneDist,
                           out.prefix=sprintf("%s.CD4.clonal.clone-cluster.network",out.prefix),byLoc=T,Th.link=2)

######### not by location
dat.plot.clonal.CD4.cloneDist.notByLoc <- table(dat.plot.clonal.CD4[,c(clone.colname,"majorCluster")])
make.clone.cluster.network(dat.plot.clonal.CD4.cloneDist.notByLoc,
                           out.prefix=sprintf("%s.CD4.clonal.clone-cluster.network.notByLoc",out.prefix),byLoc=F,Th.link=1)
make.clone.cluster.network(dat.plot.clonal.CD4.cloneDist.notByLoc,
                           out.prefix=sprintf("%s.CD4.clonal.clone-cluster.network.notByLoc",out.prefix),byLoc=F,Th.link=2)
#########




## cross multiple clusters or only one cluster
f.clonal.multiComp <- apply(dat.plot.clonal.CD8.cloneDist,1,function(x){ sum(x>0)>1 })
dat.plot.clonal.CD8.cloneDist.multi <- dat.plot.clonal.CD8.cloneDist[f.clonal.multiComp,]
f.clonal.singleComp <- apply(dat.plot.clonal.CD8.cloneDist,1,function(x){ sum(x>0)==1 })
dat.plot.clonal.CD8.cloneDist.single <- dat.plot.clonal.CD8.cloneDist[f.clonal.singleComp,]
### reorder
if(args.facet=="loc"){
    dat.plot.clonal.CD8.cloneDist.multi <- 
        dat.plot.clonal.CD8.cloneDist.multi[,c(which(grepl("^P",colnames(dat.plot.clonal.CD8.cloneDist.multi),perl = T)),
                                               which(grepl("^N",colnames(dat.plot.clonal.CD8.cloneDist.multi),perl = T)),
                                               which(grepl("^T",colnames(dat.plot.clonal.CD8.cloneDist.multi),perl = T)))]
}else{
    #dat.plot.clonal.CD8.cloneDist.multi <- 
    #    dat.plot.clonal.CD8.cloneDist.multi[,c(which(grepl(".P$",colnames(dat.plot.clonal.CD8.cloneDist.multi),perl = T)),
    #                                           which(grepl(".N$",colnames(dat.plot.clonal.CD8.cloneDist.multi),perl = T)),
    #                                           which(grepl(".T$",colnames(dat.plot.clonal.CD8.cloneDist.multi),perl = T)))]
}
### plot TCR sharing matrix
tmp.dat.plot <- dat.plot.clonal.CD8.cloneDist.multi
tmp.dat.plot[tmp.dat.plot>4] <- 4 
if(nrow(tmp.dat.plot)>2 && ncol(tmp.dat.plot)>2){
    tmp.dat.plot <- tmp.dat.plot[hclust(dist(tmp.dat.plot))$order,]
    if(sample.id=="lungC"){
        #-ZNF683
        j.ZNF683 <- which(grepl("-ZNF683",colnames(tmp.dat.plot),perl = T))
        tmp.dat.plot <- tmp.dat.plot[names(sort(apply(tmp.dat.plot,1,
                                                      function(x){ w <- rev(0.5^seq_along(j.ZNF683))
                                                                   sum( (x[j.ZNF683]>0)*w) }),decreasing = T)),]
        # P.CD8_C2-CCL5, T.CD8_C4-GZMK
        ##j.GZMK <- which(grepl("(P.CD8_C2-CCL5|T.CD8_C4-GZMK)",colnames(tmp.dat.plot),perl = T))
        ##j.GZMK <- which(grepl("^P.+CD8_C2-CCL5",colnames(tmp.dat.plot),perl = T))
        j.GZMK <- which(grepl("CD8_C2-CCL5",colnames(tmp.dat.plot),perl = T))
        tmp.dat.plot <- tmp.dat.plot[names(sort(apply(tmp.dat.plot,1,
                                                      function(x){ w <- 0.6^seq_along(j.GZMK)
                                                                   sum( (x[j.GZMK]>0)*w) }),decreasing = T)),]
        #-CXCL13
        j.CXCL13 <- which(grepl("-CXCL13",colnames(tmp.dat.plot),perl = T))
        tmp.dat.plot <- tmp.dat.plot[names(sort(apply(tmp.dat.plot,1,
                                                      function(x){ w <- rev(0.9^seq_along(j.CXCL13))
                                                                   sum( (x[j.CXCL13]>0)*w) }),decreasing = T)),]
    }else if(sample.id=="colonC"){
        j.vec <- which(grepl("^P.+-CCL5",colnames(tmp.dat.plot),perl = T))
        tmp.dat.plot <- tmp.dat.plot[names(sort(apply(tmp.dat.plot,1,
                                                      function(x){ w <- 0.6^seq_along(j.vec)
                                                                   sum( (x[j.vec]>0)*w) }),decreasing = T)),]
        #-CXCL13
        j.vec <- which(grepl("-HAVCR2",colnames(tmp.dat.plot),perl = T))
        tmp.dat.plot <- tmp.dat.plot[names(sort(apply(tmp.dat.plot,1,
                                                      function(x){ w <- rev(0.9^seq_along(j.vec))
                                                                   sum( (x[j.vec]>0)*w) }),decreasing = T)),]
    }else if(sample.id=="liverC"){
        j.vec <- which(grepl("SLC4A10",colnames(tmp.dat.plot),perl = T))
        tmp.dat.plot <- tmp.dat.plot[names(sort(apply(tmp.dat.plot,1,
                                                      function(x){ w <- 0.6^seq_along(j.vec)
                                                                   sum( (x[j.vec]>0)*w) }),decreasing = T)),]
        j.vec <- which(grepl("-LAYN",colnames(tmp.dat.plot),perl = T))
        tmp.dat.plot <- tmp.dat.plot[names(sort(apply(tmp.dat.plot,1,
                                                      function(x){ w <- rev(0.9^seq_along(j.vec))
                                                                   sum( (x[j.vec]>0)*w) }),decreasing = T)),]
    }
    if(T){
        ## CX3CR1
        j.CX3CR1 <- which(grepl("-CX3CR1",colnames(tmp.dat.plot),perl = T))
        tmp.dat.plot <- tmp.dat.plot[names(sort(apply(tmp.dat.plot,1,
                                                      function(x){ w <- 0.3^seq_along(j.CX3CR1)
                                                                   sum( (x[j.CX3CR1]>0)*w) }),decreasing = T)),]
    }
    pnt.color <- structure(c("#E41A1C","#A65628","#377EB8"),names=c("P","N","T"))

	cc <- colnames(tmp.dat.plot)
    if(args.facet=="loc"){
        cc <- pnt.color[regmatches(cc,regexpr(pattern = "^(.)",text = cc,perl = T))]
    }else{
        cc <- majorClusterColor[gsub("..$","",cc)]
    }
    print(cc)
    pdf(sprintf("%s.CD8.TCRSharing.pdf",out.prefix),width=7,height=10)
    par(fig = c(0.05, 0.95, 0, 1.0), mar = c(0, 0, 0, 0), xpd = NA) 
          ###col=colorRampPalette(c("#2f98ce", "#e0f3db", "#f47104")),
    heatmap.2(tmp.dat.plot, na.rm=F, na.col="grey50",
          col=colorRampPalette(c(brewer.pal(5,"Blues"))),
          cexRow=min(1.0,50/nrow(tmp.dat.plot)),
          Colv=F,Rowv=F,dendrogram ="none",ColSideColors=cc,breaks=6,
          density.info="none", trace="none",scale="none",margins=c(14,10),
          key.par=list(mar=c(12,3,4,3)))
    dev.off()
}

####### plot clonity distribution
### 
dd <- list("CD8"=dat.plot.CD8,"CD4"=dat.plot.CD4)
print(str(dd))
sink(file = sprintf("%s.clone.cluster.txt",out.prefix))
for(nn in names(dd)){
    this.dd <- subset(dd[[nn]],loc %in% c("P","N","T"))
	if(nrow(this.dd)<3){next}
	#pdf(sprintf("%s.%s.bar.count.pdf",out.prefix,nn),width=pdf.width,height=4)
	p <- ggplot(this.dd[,c("majorCluster","ctype","loc")]) +
		geom_bar(aes(majorCluster,fill=ctype)) + 
		theme_bw(base_size = 12) + 
		facet_wrap( ~ loc, ncol=3,scales = "fixed") + 
		scale_fill_manual(values = c.data$colII) +
		theme(axis.text.x = element_text(angle = 45, hjust = 1),strip.placement = "inside")
	ggsave(sprintf("%s.%s.bar.count.pdf",out.prefix,nn),width=pdf.width,height=4)
	#pdf(sprintf("%s.%s.bar.freq.pdf",out.prefix,nn),width=pdf.width,height=4)
	p <- ggplot(this.dd[,c("majorCluster","ctype","loc")]) +
		geom_bar(aes(majorCluster,fill=ctype), position = "fill") + 
		labs(list(y = "freq")) +
		theme_bw(base_size = 12) + 
		facet_wrap( ~ loc, ncol=3,scales = "fixed") + 
		scale_fill_manual(values = c.data$colII) +
		theme(axis.text.x = element_text(angle = 45, hjust = 1),strip.placement = "inside")
	ggsave(sprintf("%s.%s.bar.freq.pdf",out.prefix,nn),width=pdf.width,height=4)

	for(ll in unique(this.dd$loc)){
		cat(sprintf("%s\t%s\n",nn,ll))
		with(subset(this.dd,loc==ll), table(ctype,majorCluster)) %>% print
	}
}
sink()






