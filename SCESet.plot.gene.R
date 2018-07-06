#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-i", "--inFile", type="character", required=TRUE, help="infile")
#parser$add_argument("-a", "--sampleDescFile", type="character", required=TRUE, help="sample.desc.file")
#parser$add_argument("-b", "--geneFile", type="character", required=TRUE, help="gene list file")
parser$add_argument("-o", "--outputPrefix", type="character", required=TRUE, help="output prefix")
parser$add_argument("-l", "--geneFile", type="character", required=TRUE, help="gene list file")
parser$add_argument("-a", "--groupPlot", type="character", default="majorCluster:all", help="group [default %(default)s]")
parser$add_argument("-e", "--exprs", type="character", default="tpm", help="which measurement to use [default %(default)s]")
#parser$add_argument("-s", "--sample", type="character", default="SAMPLE", help="sample id  [default %(default)s]")
#parser$add_argument("-y", "--cellTypeColorFile", type="character", default="/WPS1/zhenglt/work/TCR_chunhong/data/CellType.color", 
#                    help="cellTypeColorFile [default %(default)s]")
#parser$add_argument("-w", "--clusterColorFile", type="character", 
#                    help="clusterColorFile [default %(default)s]")
parser$add_argument("-f", "--figurePerPage", type="integer",default=10, help="figure per page [default %(default)s]")
parser$add_argument("-k", "--logV", action="store_true", default=FALSE, help="log value [default %(default)s]")
parser$add_argument("-n", "--noColorBy", action="store_true", default=FALSE, help="no color by [default %(default)s]")
parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="verbose mode [default %(default)s]")
#parser$add_argument("-j", "--disableFilterGene", action="store_true", default=FALSE, help="disable filter gene [default %(default)s]")
#parser$add_argument("-d", "--pdfWidth", type="integer", default=25, help="pdf width [default %(default)s]")
args <- parser$parse_args()
print(args)

in.file <- args$inFile
out.prefix <- args$outputPrefix
gene.list.file <- args$geneFile
strGroupPlot <- args$groupPlot
fpp <- args$figurePerPage
exprs_values <- args$exprs
logV <- args$logV
noColorBy <- args$noColorBy

## figures per page

#in.file <- "/WPS1/zhenglt/work/proj_xy/integrated/knownFactors/gsva.out.MSigDB.c2c5c7.crossPatient/lung.P6/lung.P6.CD4.MSigDB.c2c5c7-c-a.gsva.RData"
#out.prefix <- "/WPS1/zhenglt/work/proj_xy/integrated/knownFactors/gsva.out.MSigDB.c2c5c7.crossPatient/lung.P6/plot/CD4/lung.P6.CD4.MSigDB.c2c5c7-c-a"
#gene.list.file <- "/WPS1/zhenglt/work/proj_xy/integrated/knownFactors/gsva.out.MSigDB.c2c5c7.crossPatient/lung.P6/diffGene/CD4/lung.CD4.withMAIT.G7.gene.desc.all.txt"
#exprs_values <- "exprs"
#logV <- F
#noColorBy <- T

#in.file <- "/WPS1/zhenglt/work/proj_xy/integrated/cross.patient/exp/lung.P6.SCESet.afterClustering.RData"
#out.prefix <- "/WPS1/zhenglt/work/proj_xy/integrated/cross.patient/CXCL13.CD4/test"
#gene.list.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/cross.patient/checkExp/liver.CD4.CXCL13.gene.desc.txt"
#gene.list.file <- "/WPS1/zhenglt/work/proj_xy/integrated/cross.patient/clusteringResult/lung.CD4.withMAIT.G7.gene.desc.txt"
#strGroupPlot <- "majorCluster:CD4_P_helper,CD4_P_FOXP3,CD4_NKG7,CD4_RGS1_GZMA,CD4_RGS1_CXCL13,CD4_RGS1_FOXP3,CD4_RGS1_FOXP3_CCR8"
#fpp <- 10
source("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")
suppressPackageStartupMessages(library("scran"))
suppressPackageStartupMessages(library("scater"))

loginfo(" begin ...")

lname <- load(in.file)
#rm(is_exprs,envir=sce@assayData)
if(noColorBy){
    sce[["score"]] <- T
}

groupPlot <- list()
gg <- strsplit(strGroupPlot,":")
for(i in length(gg)){
    if(gg[[i]][2]=="all"){
        groupPlot[[ gg[[i]][1] ]] <- unique(sce[[ gg[[i]][1] ]])
    }else{
        groupPlot[[ gg[[i]][1] ]] <- unlist(strsplit(gg[[i]][2],","))
    }
}

gene.desc <- read.delim(gene.list.file,row.names="geneID",header = T,sep = "\t",stringsAsFactors=F)
g.list <- intersect(featureNames(sce),rownames(gene.desc))
##g.list <- g.list[1:25]
for(nn in names(groupPlot)){
    
    for(i in seq_len(ceiling(length(g.list)/fpp))){
        gidx <- intersect( ( ((i-1)*fpp+1):(i*fpp) ),seq_along(g.list))
        print(gidx)       
        sce.plot <- sce[g.list[gidx],sce[[nn]] %in% groupPlot[[nn]]]
        ##featureNames(sce.plot)<-featureData(sce.plot)$geneSymbol
        featureNames(sce.plot)<-gene.desc[featureNames(sce.plot),"geneSymbol"]
        pdf(sprintf("%s.sce.%s.Page%04d%s.pdf",out.prefix,nn,i,ifelse(fpp==1,sprintf(".%s",featureNames(sce.plot)[1]),"")),width=10,height=2*ceiling(dim(sce.plot)[1]/2)+2)
        p2 <- plotExpression(sce.plot, featureNames(sce.plot), x = nn, show_median=T, exprs_values = exprs_values,log2_values = logV,colour_by=ifelse(noColorBy,"score",NULL)) +
            theme(axis.text.x  = element_text(angle=45, vjust=0.5)) +
            labs(list(y = ifelse(exprs_values=="tpm","log2(TPM+1)","Score")))
        multiplot(p2,cols = 1)
        dev.off()
    }
}


