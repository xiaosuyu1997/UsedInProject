#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-i", "--inFile", type="character", required=TRUE, help="infile")
parser$add_argument("-a", "--sampleDescFile", type="character", required=TRUE, help="sample.desc.file (contain stype, sampleType etc)")
parser$add_argument("-b", "--geneFile", type="character", help="gene list file")
parser$add_argument("-o", "--outputPrefix", type="character", required=TRUE, help="output prefix")
parser$add_argument("-k", "--groupBy", type="character",default="majorCluster,clonotype,sampleType", help="groupBy [default %(default)s]")
parser$add_argument("-q", "--groupHave", type="character",default="majorCluster,clonotype,sampleType", help="groupHave [default %(default)s]")
parser$add_argument("-s", "--sample", type="character", default="SAMPLE", help="sample id  [default %(default)s]")
parser$add_argument("-c", "--cloneColName", type="character", default="C_strict", help="colname of clone data  [default %(default)s]")
parser$add_argument("-x", "--rename", action="store_true", default=FALSE, help="rename sampleType (need stype field in designFile) [default %(default)s]")
parser$add_argument("-y", "--cellTypeColorFile", type="character", default="/WPS1/zhenglt/work/TCR_chunhong/data/CellType.color", 
                    help="cellTypeColorFile [default %(default)s]")
parser$add_argument("-w", "--clusterColorFile", type="character", 
                    help="clusterColorFile [default %(default)s]")
parser$add_argument("-f", "--clonotypeFile", type="character", help="clonotype file")
parser$add_argument("-t", "--rowSort", action="store_true", default=FALSE, help="row sort [default %(default)s]")
parser$add_argument("-l", "--log", action="store_true", default=FALSE, help="whether do log transform [default %(default)s]")
parser$add_argument("-r", "--center", action="store_true", default=FALSE, 
                    help="whether center genes by individual [default %(default)s]")
parser$add_argument("-n", "--geneOnTopAnnotation", type="character", help="gene list file specify the genes showed on top annotation(contian geneID, geneSymbol)")
parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="verbose mode [default %(default)s]")
##parser$add_argument("-m", "--useTPM", action="store_true", default=FALSE, help="useTPM [default %(default)s]")
parser$add_argument("-z", "--measure", type="character",
                    help="measurement to be used (exprs,norm_exprs,tpm,counts) [default %(default)s]")
parser$add_argument("-e", "--disableGSymbol", action="store_true", default=FALSE, help="disable gene symbol transform [default %(default)s]")
parser$add_argument("-j", "--disableFilterGene", action="store_true", default=FALSE, help="disable filter gene [default %(default)s]")
parser$add_argument("-p", "--scale", action="store_true", default=FALSE, help="do scale [default %(default)s]")
parser$add_argument("-d", "--pdfWidth", type="integer", default=10, help="pdf width [default %(default)s]")
parser$add_argument("-u", "--pdfHeight", type="integer", default=8, help="pdf heihgt [default %(default)s]")
args <- parser$parse_args()
print(args)

in.file <- args$inFile
sample.desc.file <- args$sampleDescFile
gene.file <- args$geneFile
out.prefix <- args$outputPrefix
sample.id <- args$sample
cellTypeColorFile <- args$cellTypeColorFile
clusterColorFile <- args$clusterColorFile
clonotype.file <- args$clonotypeFile
args.rowSort <- args$rowSort
args.log <- args$log
args.center <- args$center
###scoreMatrix.file <- args$scoreMatrix
mode.verbose <- args$verbose
args.scale <- args$scale
args.clone.colname <- args$cloneColName
args.rename <- args$rename
args.measure <- args$measure

#args.rowSort <- FALSE
pdf.width <- args$pdfWidth
pdf.height <- args$pdfHeight
#bin.exp.file <- args$binFile
args.groupBy <- args$groupBy
args.groupHave <- args$groupHave
args.disableGSymbol <-  args$disableGSymbol
args.disableFilterGene <- args$disableFilterGene
geneOnTopAnnotation.file <- args$geneOnTopAnnotation
args.useTPM <- args$useTPM
geneOnTopAnnotation.discrete <- F

if(file.exists("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")){
    source("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")
}else{
    source("/WPSnew/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")
}

suppressPackageStartupMessages(library("ComplexHeatmap"))
suppressPackageStartupMessages(library("circlize"))
suppressPackageStartupMessages(library("gridBase"))
suppressPackageStartupMessages(library("dendextend"))
suppressPackageStartupMessages(library("RColorBrewer"))
suppressPackageStartupMessages(library("gplots"))
suppressPackageStartupMessages(library("vioplot"))
suppressPackageStartupMessages(library("plotrix"))
suppressPackageStartupMessages(library("scater"))

loginfo("begin ...")
dir.create(dirname(out.prefix),showWarnings = F,recursive = T)

#### sample data
sample.desc <- read.delim(sample.desc.file,header = T,sep = "\t",check.names = F,stringsAsFactors = F)
rownames(sample.desc) <- sample.desc$sample
print(dim(sample.desc))
print(head(sample.desc))
print(str(sample.desc))
##q()

### clonotype data
clonotype.data <- read.clonotype(in.file = clonotype.file,ctype.col = args.clone.colname)
#### exp data
g.inputD <- processInput(sample.desc.file,cellTypeColorFile,in.file,args.disableFilterGene,geneFile=NULL,
                         args.center=args.center,args.log=args.log,args.measure=args.measure)
myDesign <- g.inputD$myDesign
sample.desc <- myDesign
#sampleTypeColor <- g.inputD$sampleTypeColor
patient.col.list <- g.inputD$patient.col.list
Y <- g.inputD$Y
g.GNAME <- g.inputD$g.GNAME
##### 

#### rename rownames
Y.newRowname <- Y
newName <- g.GNAME[rownames(Y)]
names(newName) <- rownames(Y)
f.gname.na <- is.na(newName)
newName[f.gname.na] <- names(newName)[f.gname.na]
rownames(Y.newRowname) <- unname(newName)
#if(args.geneIDType=="ensemblID"){ Y <- Y.newRowname }

print(dim(Y))
print(Y[1:4,1:6])

#### gene data
if(!is.null(gene.file)){
    gene.desc <- read.delim(gene.file,header = T,sep = "\t",check.names = F,stringsAsFactors = F)
    gene.desc$geneID <- as.character(gene.desc$geneID)
    rownames(gene.desc) <- gene.desc$geneID
}else{
    gene.desc <- data.frame(geneID=rownames(Y),geneSymbol=g.GNAME[rownames(Y)],stringsAsFactors = F)
}
print(str(gene.desc))

geneOnTopAnnotation <- NULL
if(!is.null(geneOnTopAnnotation.file) && file.exists(geneOnTopAnnotation.file)){
    geneOnTopAnnotation <- read.delim(geneOnTopAnnotation.file,row.names="geneID",header = T,sep = "\t",check.names=F,stringsAsFactors=F)
    print(str(geneOnTopAnnotation))
}

g.f <- intersect(gene.desc$geneID,rownames(Y))
Y.inGeneList <- Y[g.f,,drop=F]
gene.desc <- gene.desc[g.f,,drop=F]
### cell type color
if(args.rename){
    sample.desc$sampleType <- paste(sample.desc$stype,sapply(strsplit(sample.desc$sampleType,""),function(x){x[1]}),sep=".")
}


if(!is.null(cellTypeColorFile) && file.exists(cellTypeColorFile)){
    sampleTypeColor <- read.SampleTypeColor(cellTypeColorFile)
    sampleTypeColor <- sampleTypeColor[ names(sampleTypeColor) %in% unique(sample.desc[colnames(Y),"sampleType"]) ]
}else{
    sampleTypeColor <- auto.colSet(n=length(unique(myDesign[,"sampleType"])),name="Paired")
    names(sampleTypeColor) <- unique(myDesign[,"sampleType"])
}
print("sampleTypeColor:")
print(sampleTypeColor)

### major cluster color
if("majorCluster" %in% colnames(sample.desc)){
    sample.desc <- sample.desc[order(sample.desc$majorCluster),]
    nMajor <- length(unique(sample.desc$majorCluster))
    if(!is.null(clusterColorFile) && file.exists(clusterColorFile)){
        majorClusterColor <- read.SampleTypeColor(clusterColorFile)
        print(majorClusterColor)
        majorClusterColor <- majorClusterColor[ names(majorClusterColor) %in% unique(sample.desc[colnames(Y),"majorCluster"]) ]
    }else{
        majorClusterColor <- structure(colorRampPalette(brewer.pal(12,"Paired"))(nMajor),
                                  names=unique(sample.desc$majorCluster))
    }
    majorClusterColor <- majorClusterColor[order(names(majorClusterColor))]
    #print("majorClusterColor")
    #print(nMajor)
    #print(names(majorClusterColor))
    #print(majorClusterColor)
    sample.desc$majorCluster <- factor(sample.desc$majorCluster,levels = unique(sample.desc$majorCluster))
}
print("ZZZZZ")
dat.plot <- data.frame()
### row order
if(args.rowSort){
    loginfo("sort row by hclust")
    dat.plot <- Y.inGeneList
}else{
    if("Group" %in% colnames(gene.desc)){
        loginfo("sort row by Group then hclust")
        for(g.cls in unique(gene.desc$Group))
        {
            g.desc <- subset(gene.desc,Group==g.cls)
            dat.t <- Y.inGeneList[g.desc$geneID,,drop=F]
            # tmp disable it
            if(nrow(dat.t)>2){
                tryCatch({
                    hc.t <- hclust(dist(dat.t),"complete")
                    dat.t <- dat.t[hc.t$order,,drop=F]
                },error=function(e) e)
            }
            if(nrow(dat.plot)==0){
                dat.plot <- dat.t
            }else{
                dat.plot <- rbind(dat.plot,dat.t)
            }
        }
    }else{
        dat.plot <- Y.inGeneList
    }
}
print("XXXXXXXXXXX")
### column order
if(ncol(dat.plot)>2){
        tryCatch({
            #hc.t <- hclust(dist(t(dat.plot)), "complete")
            hc.t <- hclust(as.dist(1-cor(dat.plot)),"complete")
            dat.plot <- dat.plot[,hc.t$order,drop=F]
        },error = function(e){
            hc.t <- hclust(dist(t(dat.plot)), "complete")
            dat.plot <- dat.plot[,hc.t$order,drop=F]
            e
        })
}

## add more properties to the sample.desc
if(!is.null(clonotype.data)){
    clonotype.d.t <- clonotype.data$ctypeII[rownames(sample.desc)]
    clonotype.d.t[is.na(clonotype.d.t)] <- "N0"
    sample.desc$clonotype <- clonotype.d.t
}
#sample.desc$majorCluster <- factor(sample.desc$majorCluster,levels = unique(sample.desc$majorCluster))
if(!is.null(geneOnTopAnnotation)){
    d.t <- t(Y[rownames(geneOnTopAnnotation),rownames(sample.desc),drop=F])
    colnames(d.t) <- geneOnTopAnnotation$geneSymbol
    sample.desc<-cbind(sample.desc[,!colnames(sample.desc) %in% colnames(d.t)],d.t)
    print("summary of d.t:")
    print(summary(d.t))
    print("sample.desc:")
    print(summary(sample.desc))
    rm(d.t)
}
### sort
cat(sprintf("with(sample.desc,order(%s))\n",args.groupBy))
print(head(sample.desc))
sample.desc <- sample.desc[eval(parse(text=sprintf("with(sample.desc,order(%s))",args.groupBy))),]
#print("after ordering...")
#print(str(sample.desc))
#print(dim(dat.plot))
#print(dat.plot[1:4,1:8])
if(!is.null(clonotype.data)){ sample.desc$clonotype[sample.desc$clonotype=="N0"] <- "NA" }
dat.plot <- dat.plot[,rownames(sample.desc),drop=F]
##
#print(dat.plot[1:4,1:8])

topAnnotation.legend.para.list <- NULL
grpl <- unlist(strsplit(args.groupHave,split=","))
grpl.sel <- sapply(grpl,function(x){ unlist(strsplit(x,":"))[1] })
grpl.isSimple <- sapply(grpl,function(x){ length(unlist(strsplit(x,":")))==1 })
annDF <- sample.desc[,grpl.sel[grpl.isSimple],drop=F]
print("annDF:")
print(summary(annDF))
#print(head(annDF))
annCol <- list()
annFunc <- list()
for(i in seq_along(grpl)){
    if(grpl[i]=="majorCluster")
    {
        annCol[[grpl[i]]]=majorClusterColor
    }else if(grpl[i]=="sampleType"){
        annCol[[grpl[i]]]=sampleTypeColor
    }else if(grpl[i]=="clonotype" && !is.null(clonotype.data)){
        annCol[[grpl[i]]]=clonotype.data$colII
    }else{
        grpl.i <- unlist(strsplit(grpl[i],":"))
        if(length(grpl.i)==1){
            annDF[,grpl.i[1]] <- sample.desc[,grpl.i[1]]
            ###annDF[is.na(annDF[,grpl.i[1]]),grpl.i[1]] <- "NA"
            if(class(sample.desc[,grpl.i[1]])=="numeric"){
                if(all(annDF[,grpl.i[1]]<=1) && all(annDF[,grpl.i[1]]>=0)){
                    Y.level <- c(0,1) 
                }else{
                    Y.level <- pretty(annDF[,grpl.i[1]],n=8)
                }
                print(sprintf("Y.level(%s):",grpl.i[1]))
                print(Y.level)
                #print(str(annDF[,grpl.i[1]]))
                #print(summary(annDF[,grpl.i[1]]))
                # discreate version
                if(geneOnTopAnnotation.discrete){
                    annDF[ annDF[,grpl.i[1]] > Y.level[length(Y.level)], grpl.i[1]] <- Y.level[length(Y.level)]
                    ann.t <- cut(annDF[,grpl.i[1]], breaks=Y.level,include.lowest=T,ordered_result = T)
                    annDF[,grpl.i[1]] <- as.character(ann.t)
                    ##nLevel <- length(unique(annDF[,grpl[i]]))
                    nLevel <- nlevels(ann.t)
                    ###annCol[[grpl.i[1]]] <- rev(structure(colorRampPalette(rev(brewer.pal(nLevel,ifelse(args.scale,"RdBu","RdYlBu"))))(nLevel),
                    annCol[[grpl.i[1]]] <- rev(structure(colorRampPalette(rev(brewer.pal(nLevel,"RdYlBu")))(nLevel),
                                                   names=levels(ann.t)))
                    if(is.null(topAnnotation.legend.para.list)){
                        topAnnotation.legend.para.list <- list()
                    }
                    topAnnotation.legend.para.list[[grpl.i[1]]] <- list(color_bar="discrete",legend_direction="horizontal")
                }else{
                    # continious version
                    annCol[[grpl.i[1]]]=colorRamp2(seq(Y.level[1],Y.level[length(Y.level)],length=7),
                                                 rev(brewer.pal(n = 7, name = "RdYlBu")), 
                                                 ###rev(brewer.pal(n = 7, name = ifelse(args.scale,"RdBu","RdYlBu"))), 
                                                 space="LAB")
                    if(is.null(topAnnotation.legend.para.list)){
                        topAnnotation.legend.para.list <- list()
                    }
                    topAnnotation.legend.para.list[[grpl.i[1]]] <- list(color_bar="continuous",legend_direction="horizontal",legend_width=unit(4, "cm"),legend_height=unit(2, "cm"))
                }
            }else{
                vGrpI <- unique(sample.desc[,grpl[i]])
                ###annCol[[grpl[i]]]=structure(colorRampPalette(brewer.pal(length(vGrpI),"Paired"))(length(vGrpI)), names=vGrpI)
                annCol[[grpl[i]]]=structure(colorRampPalette(brewer.pal(12,"Paired"))(length(vGrpI)), names=vGrpI)
            }
        }else if(grpl.i[2]=="b"){
            annFunc[[grpl.i[1]]] <- SingleAnnotation(name = grpl.i[1], 
                                                     fun = anno_barplot(sample.desc[,grpl.i[1]], axis = F,border = F,
                                                                        gp=gpar(fill="#CCCCCC",col="#CCCCCC",lwd=0)), 
                                                     which = "column")
        }
    }
}

ha.col <- HeatmapAnnotation(df = annDF, col = annCol, show_legend = T, annotation_legend_param = topAnnotation.legend.para.list)
#for(x in names(annFunc)){ ha.col@anno_list[[x]] <- annFunc[[x]] }
ha.col@anno_list <- c(rev(annFunc),ha.col@anno_list)
n_anno <- length(ha.col@anno_list)
anno_size <- c(rep(0.8,ncol(annDF)),rep(0.8,length(annFunc)))
ha.col@gap <- rep(unit(0, "mm"),n_anno)
if(ncol(annDF)>=2){
    anno_size <- anno_size/sum(anno_size) * (unit(1, "npc") - sum(ha.col@gap))
}else{
    #anno_size <- anno_size/sum(anno_size) * (unit(2.5, "npc") - sum(ha.col@gap))
    anno_size <- anno_size/sum(anno_size) * (unit(1.0, "npc") - sum(ha.col@gap))
}
ha.col@anno_size = anno_size[rev(seq_len(n_anno))]

#print("annCol")
#print(annCol)
#print(str(ha.col))
#print(head(dat.plot[,1:6]))
#print(str(sampleTypeColor))
runHierarchicalClusteringAnalysis(dat.plot,mytitle = sample.id, pdf.width=pdf.width,pdf.height=pdf.height, 
                                  sprintf("%s.%s.selectedExpressionHeatmap",out.prefix,sample.id), 
                                  do.clustering.col=F, 
                                  do.clustering.row=args.rowSort, 
                                  sampleType=NULL, 
                                  colSet=sampleTypeColor,
                                  ha.col=ha.col,
                                  ann.extra.df = annDF, ann.extra.df.col = annCol, ann.bar.height = min(1,6/length(ha.col@anno_list)), 
                                  k.row=1,gid.mapping=g.GNAME,
                                  clonotype.col=NULL,ntop=NULL, 
                                  row.names.original=args.disableGSymbol, 
                                  ##annotation_legend_param=topAnnotation.legend.para.list,z.lo=-2.5,z.hi=2.5,z.step=0.5,z.title="ZScore", 
                                  annotation_legend_param=topAnnotation.legend.para.list,z.lo=-2.5,z.hi=2.5,z.step=0.5,z.title="Exp", 
                                  complexHeatmap.use=TRUE,verbose=FALSE,do.scale=args.scale)

### output gene list to text file
if(mode.verbose){
    dat.plot.df <- data.frame(geneID=rownames(dat.plot),stringsAsFactors = F)
    dat.plot.df$geneSymbol=g.GNAME[dat.plot.df$geneID]
    dat.plot.df <- cbind(dat.plot.df,dat.plot)
    write.table(dat.plot.df,file = sprintf("%s.%s.selectedExpressionHeatmap.txt",out.prefix,sample.id), 
                quote = F,sep = "\t",row.names = F,col.names = T)
}

loginfo("end.")
