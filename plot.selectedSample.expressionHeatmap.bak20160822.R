#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-i", "--inFile", type="character", required=TRUE, help="infile")
parser$add_argument("-a", "--sampleDescFile", type="character", required=TRUE, help="sample.desc.file")
parser$add_argument("-b", "--geneFile", type="character", required=TRUE, help="gene list file")
parser$add_argument("-o", "--outputPrefix", type="character", required=TRUE, help="output prefix")
parser$add_argument("-z", "--binFile", type="character", help="binarized expression file")
parser$add_argument("-s", "--sample", type="character", default="SAMPLE", help="sample id  [default %(default)s]")
parser$add_argument("-y", "--cellTypeColorFile", type="character", default="/WPS1/zhenglt/work/TCR_chunhong/data/CellType.color", help="cellTypeColorFile [default %(default)s]")
parser$add_argument("-l", "--clonotypeFile", type="character", help="clonotype file")
parser$add_argument("-w", "--disableLog", action="store_true", default=FALSE, help="disable log transform [default %(default)s]")
parser$add_argument("-r", "--rowSort", action="store_true", default=FALSE, help="row sort [default %(default)s]")
parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="verbose mode [default %(default)s]")
#parser$add_argument("-b", "--disableGSymbol", action="store_true", default=FALSE, help="disable gene symbol transform [default %(default)s]")
#parser$add_argument("-e", "--scale", action="store_true", default=FALSE, help="do scale [default %(default)s]")
#parser$add_argument("-d", "--pdfWidth", type="integer", default=18, help="pdf width [default %(default)s]")
args <- parser$parse_args()
print(args)

in.file <- args$inFile
sample.desc.file <- args$sampleDescFile
gene.file <- args$geneFile
out.prefix <- args$outputPrefix
sample.id <- args$sample
###sample.type <- args$sampleType
cellTypeColorFile <- args$cellTypeColorFile
clonotype.file <- args$clonotypeFile
disable.log <- args$disableLog
mode.verbose <- args$verbose
#pdf.width <- args$pdfWidth
bin.exp.file <- args$binFile
args.rowSort <- args$rowSort

##in.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/gsnap_out/quantification/sizeFactor/P0322/P0322.all.countData.sfNormalized.txt.gz"
##sample.desc.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/byTCellType/selectedSamples/P0322/P0322.FOXP3.expression.TTC.txt"
##gene.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/byTCellType/geneList/cytotoxic.FOXP3"
##out.prefix <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/byTCellType/selectedSamples/P0322/P0322.TTC"
##sample.id <- "P0322"
##cellTypeColorFile <- "/WPS1/zhenglt/work/TCR_chunhong/data/CellType.color"
##clonotype.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/tracer/OUT/P0322/filtered_TCR_summary/P0322.summary.cell.reassigneClonotype.methodChunhong.r.txt"
##disable.log <- FALSE
##mode.verbose <- FALSE
##bin.exp.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/byTCellType/selectedSamples/P0322/P0322.goi.exp.mmodel.out.binary.txt"
#args.rowSort <- F

source("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")

suppressPackageStartupMessages(require("ComplexHeatmap"))
suppressPackageStartupMessages(require("circlize"))
suppressPackageStartupMessages(require("gridBase"))
suppressPackageStartupMessages(require("dendextend"))
suppressPackageStartupMessages(require("RColorBrewer"))
suppressPackageStartupMessages(require("gplots"))
suppressPackageStartupMessages(library("vioplot"))
suppressPackageStartupMessages(library("plotrix"))

loginfo("begin ...")

### clonotype data
clonotype.data <- read.clonotype(in.file = clonotype.file,ctype.col = "C_strict")
print(str(clonotype.data))
#### exp data
in.table <- read.table(in.file,header = T,sep = "\t",check.names = F,stringsAsFactors = F,row.names = 1)
g.GNAME <- in.table[,1]
names(g.GNAME) <- rownames(in.table)
in.table <- as.matrix(in.table[,-1])
if(!disable.log){ in.table <- as.matrix(log2(in.table+1)) }

print(dim(in.table))
print(in.table[1:4,1:8])
### binarized exp data
##if(!is.null(bin.exp.file) && file.exists(bin.exp.file)){
##    bin.exp.table <- read.table(bin.exp.file,header = T,sep = "\t",check.names = F,stringsAsFactors = F,row.names = 1)
##    bin.exp.table <- bin.exp.table[,-1]
##}

#### gene data
gene.desc <- read.table(gene.file,header = T,sep = "\t",check.names = F,stringsAsFactors = F,colClasses = c("character","character"))
#### sample data
sample.desc <- read.table(sample.desc.file,header = T,sep = "\t",check.names = F,stringsAsFactors = F)
rownames(sample.desc) <- sample.desc$sample
print(dim(sample.desc))
print(head(sample.desc))

s.f <- intersect(sample.desc$sample,colnames(in.table))
sample.desc <- sample.desc[s.f,,drop=F]
in.table <- in.table[,s.f,drop=F]

g.f <- intersect(gene.desc$geneID,rownames(in.table))
in.table <- in.table[g.f,,drop=F]
### cell type color
sampleTypeColor <- read.SampleTypeColor(cellTypeColorFile)
sampleTypeColor <- sampleTypeColor[ names(sampleTypeColor) %in% unique(sample.desc[colnames(in.table),"sampleType"]) ]

#dat.plot <- in.table
### reorder samples by leaf cluster than hclust
dat.plot <- c()
for(s.cls in unique(sample.desc$leafCluster))
{
    s.desc <- subset(sample.desc,leafCluster==s.cls)
    dat.t <- in.table[,rownames(s.desc),drop=F]
    if(ncol(dat.t)>2){
        tryCatch({
            hc.t <- hclust(dist(t(dat.t)), "complete")
            dat.t <- dat.t[,hc.t$order,drop=F]
        },error = function(e) e)
    }
    dat.plot <- cbind(dat.plot,dat.t)
}
### leaf cluster color
nLeaf <- length(unique(sample.desc$leafCluster))
leafClusterColor <- structure(colorRampPalette(brewer.pal(nLeaf,"Paired"))(nLeaf),
                              names=unique(sample.desc$leafCluster))
#leafClusterColor["uncharacterized"] <- "gray"

annDF <- data.frame(clonotype=clonotype.data$ctypeII[colnames(dat.plot)],stringsAsFactors = F,leafCluster=sample.desc[colnames(dat.plot),"leafCluster"])
annDF[is.na(annDF)] <- "NA"
annCol <- list(clonotype=clonotype.data$colII,leafCluster=leafClusterColor)
#if(!is.null(scoreMatrix.file)){
#    annDF <- cbind(annDF,t(scoreMatrix[c("TReg","exhaustion_G5","cytotoxic_G8","naive_G4"),colnames(dat.1)]))
#    annCol <- c(annCol,list(TReg=scoreM.col.fun,
#                            exhaustion_G5=scoreM.col.fun,
#                            cytotoxic_G8=scoreM.col.fun,
#                            naive_G4=scoreM.col.fun))
#}else{
#}
runHierarchicalClusteringAnalysis(dat.plot,mytitle = sample.id, pdf.width=18,pdf.height=8, 
                                  sprintf("%s.%s.selectedExpressionHeatmap",out.prefix,sample.id), 
                                  do.clustering.col=F, 
                                  do.clustering.row=args.rowSort, 
                                  sampleType=sample.desc[colnames(dat.plot),"sampleType"], 
                                  colSet=sampleTypeColor,
                                  ann.extra.df = annDF, ann.extra.df.col = annCol, ann.bar.height = 1.1, 
                                  k.row=2, 
                                  clonotype.col=NULL,ntop=NULL, 
                                  row.names.original=F, 
                                  annotation_legend_param=NULL, 
                                  complexHeatmap.use=TRUE,verbose=FALSE,do.scale=F)
### output gene list to text file
if(mode.verbose){
    dat.plot.df <- data.frame(geneID=rownames(dat.plot),stringsAsFactors = F)
    dat.plot.df$geneSymbol=g.GNAME[dat.plot.df$geneID]
    dat.plot.df <- cbind(dat.plot.df,dat.plot)
    write.table(dat.plot.df,file = sprintf("%s.%s.selectedExpressionHeatmap.txt",out.prefix,sample.id), 
                quote = F,sep = "\t",row.names = F,col.names = T)
}

loginfo("end.")
