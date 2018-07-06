#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-d", "--dfile", type="character", required=TRUE, help="sample info file")
parser$add_argument("-o", "--outprefix", type="character", required=TRUE, help="outprefix")
parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="whether verbose mode [default %(default)s]")
args <- parser$parse_args()
print(args)

out.prefix <- args$outprefix
designFile <- args$dfile

#out.prefix <- "OUT.DE/Prol.hi.lo.CD8.Tex/Prol.hi.lo.CD8.Tex.de.limma"
#designFile <- "d/colonC.MKI67.addG.all.refine.CD8.Tex.txt"

cellTypeColorFile <- "/WPS1/zhenglt/work/TCR_chunhong/data/CellType.color"
inputFile <- "/WPSnew/zhenglt/work/proj_zl/integrate0924/cross.patient.b20180409/exp/colonC.P12.scran.RData"
args.notFilter <- F
geneFile <- NULL
args.center <- F
args.log <- F
args.norm.exprs <- F
args.verbose <- T
args.geneIDType <- ""
args.measure <- F

if(file.exists("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")){
    source("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")
}else{
    source("/WPSnew/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")
}

suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("RColorBrewer"))
suppressPackageStartupMessages(library("ROCR"))

dir.create(dirname(out.prefix),showWarnings = F,recursive = T)

g.inputD <- processInput(designFile,cellTypeColorFile,inputFile,args.notFilter,geneFile,
                         args.center,args.log,args.norm.exprs = args.norm.exprs,args.measure=args.measure)
myDesign <- g.inputD$myDesign
sampleTypeColor <- g.inputD$sampleTypeColor
patient.col.list <- g.inputD$patient.col.list
Y <- g.inputD$Y
if(!args.notFilter){
    f <- apply(Y,1,function(x){ nE <- sum(x>0); return( nE > 5 & nE/length(x) > 0.01 )  })
    Y <- Y[f,]
}
g.GNAME <- g.inputD$g.GNAME
##### 

#### rename rownames
Y.newRowname <- Y
newName <- g.GNAME[rownames(Y)]
names(newName) <- rownames(Y)
f.gname.na <- is.na(newName)
newName[f.gname.na] <- names(newName)[f.gname.na]
rownames(Y.newRowname) <- unname(newName)
if(args.geneIDType=="ensemblID"){ Y <- Y.newRowname }

#.clust <- as.numeric(myDesign[,sprintf("%s.bin",gene.beTested)])+1
.clust <- myDesign[,"majorCluster"]
names(.clust) <- rownames(myDesign)
clusterColor <- structure(colorRampPalette(brewer.pal(12,"Paired"))(length(unique(.clust))),
                          names=sort(as.character(unique(.clust))))


d.forLimma <- myDesign
d.forLimma$sampleType <- d.forLimma$majorCluster
res.limma <- run.limma.from.matrixFile(Y,d.forLimma,T.fdr=0.01,T.logFC=1,
                                       withPatient=F,verbose=args.verbose,
                                       gid.mapping=g.GNAME,do.voom=F)

res.limma$all %>% head
res.limma$sig %>% head


write.table(res.limma$all, file = sprintf("%s.%s.txt",out.prefix,args.measure), row.names = F,quote = F,sep = "\t")
write.table(res.limma$sig, file = sprintf("%s.%s.sig.txt",out.prefix,args.measure), row.names = F,quote = F,sep = "\t")


