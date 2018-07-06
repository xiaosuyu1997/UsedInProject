#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-d", "--designFile", type="character", required=TRUE, help="design matrix file")
parser$add_argument("-i", "--inputFile", type="character", required=TRUE, help="input file")
parser$add_argument("-c", "--cellTypeColorFile", default="/WPS1/zhenglt/work/TCR_chunhong/data/CellType.color", 
                    type="character", help="cellTypeColorFile")
parser$add_argument("-e", "--clonotypeFile", type="character", help="clonotype file")
parser$add_argument("-o", "--outDir", type="character", required=TRUE, help="output directory")
parser$add_argument("-s", "--sample", type="character", default="SAMPLE", help="sample id  [default %(default)s]")
###parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="whether verbose mode [default %(default)s]")
parser$add_argument("-l", "--log", action="store_true", default=FALSE, help="whether do log transform [default %(default)s]")
parser$add_argument("-n", "--nKeep", type="integer", default=1500, help="use top nKeep genes [default %(default)s]")
args <- parser$parse_args()

print(args)

designFile <- args$designFile
inputFile <- args$inputFile
cellTypeColorFile <- args$cellTypeColorFile
clonotype.file <- args$clonotypeFile
out.dir <- args$outDir
sample.id <- args$sample
args.log <- args$log
nKeep <- args$nKeep

##args <- commandArgs(T)
##if(length(args)<1){
##    cat("w.patch.P1116.R <nKeep>\n")
##    q()
##}
##nKeep=args[1]

#### TEST DATA 
##designFile <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/byTCellType/OUT.scLVM.byTCellType/P1116/sfIgnoreERCC/TC/P1116.designUsed.txt"
##inputFile <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/gsnap_out/quantification/sizeFactor/P1116/P1116.all.countData.sfNormalized.txt.gz"
##cellTypeColorFile <- "/WPS1/zhenglt/work/TCR_chunhong/data/CellType.color"
##clonotype.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/tracer/OUT/P1116/filtered_TCR_summary/P1116.summary.cell.reassigneClonotype.txt"
##out.dir <- sprintf("/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/byTCellType/patch.P1116/OUT.tryTopVar.%s",nKeep)
##sample.id <- "P1116"
##args.log <- TRUE
##nKeep <- 1500

dir.create(out.dir,recursive = T,showWarnings = F)
source("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")

myDesign<-read.table(designFile,header=T,row.names="sample",check.names=F,colClasses=c("factor","character","factor","factor"))
sampleTypeColor <- read.SampleTypeColor(cellTypeColorFile)

clonotype.strict.data <- read.clonotype(in.file = clonotype.file,ctype.col = "C_strict")
patient.col.list <- patientColorListFromMyDesign(myDesign)

if(grepl("RData$",inputFile,perl=T)){
    suppressPackageStartupMessages(library("R.utils"))
    lenv <- loadToEnv(inputFile)
    obj.scdn <- lenv[["obj.scdn"]]
    ##Y <- rbind(obj.scdn@normalized.endo,obj.scdn@normalized.ERCC)
    Y <- obj.scdn@normalized.endo
}else{
    in.table <- read.table(inputFile,header = T,sep = "\t",stringsAsFactors = F,check.names = F)
    rownames(in.table) <- in.table[,1]
    Y <- in.table[,c(-1,-2)]
}

sname <- intersect(rownames(myDesign),colnames(Y))
myDesign <- myDesign[sname,,drop=F]
Y <- Y[,sname,drop=F]

f <- apply(Y,1,function(x){ nE <- sum(x>0); return( nE > 5 & nE/length(x) > 0.01 )  })
Y <- Y[f,]
if(args.log) { Y <- log2(Y+1) }

loginfo("... all samples.")

doit <- function(dat.plot,extra=""){
    pca.res <- runPCAAnalysis(dat.plot,sprintf("%s/%s%s.het.PCA",out.dir,sample.id,extra),
                   myDesign[,"sampleType"],sampleTypeColor[names(sampleTypeColor) %in% unique(as.character(myDesign$sampleType))],
                   ntop=NULL,main=sample.id)

    tsne.res <- runTSNEAnalysis(dat.plot,sprintf("%s/%s%s.het.tSNE",out.dir,sample.id,extra),
                    col.points = sampleTypeColor[as.character(myDesign$sampleType)],
                    legend=c(names(sampleTypeColor)[names(sampleTypeColor) %in% unique(as.character(myDesign$sampleType))],levels(myDesign$libType)),
                    col.legend=c(sampleTypeColor[names(sampleTypeColor) %in% unique(as.character(myDesign$sampleType))],rep("black",length(levels(myDesign$libType)))),
                    pch=(as.numeric(myDesign$libType)-1) %% 26,cex=0.7,
                    pch.legend=c(rep(16,sum(names(sampleTypeColor) %in% unique(as.character(myDesign$sampleType)))),(seq_along(levels(myDesign$libType))-1) %% 26)
                    )

    for(nn in c(50,100,150,200,250,300,350,400,450,500,1000,2000,100000))
    {
    runHierarchicalClusteringAnalysis(dat.plot,
                    sprintf("%s/%s%s.het.hclustering.n%s",out.dir,sample.id,extra,ifelse(nn>99999,"All",nn)),
                    myDesign[sname,"sampleType"],sampleTypeColor[names(sampleTypeColor) %in% unique(as.character(myDesign$sampleType))],
                    clonotype.col=clonotype.strict.data,
                    patient.col.list = patient.col.list,
                    ntop=nn,
                    complexHeatmap.use=TRUE,do.cuttree=T,
                    verbose=TRUE,main="Variable Genes")
    }

    sc3.res <- runSC3Analysis(dat.plot,sprintf("%s/%s%s.SC3",out.dir,sample.id,extra),
                   myDesign[sname,"sampleType"],sampleTypeColor[names(sampleTypeColor) %in% unique(as.character(myDesign[sname,"sampleType"]))],
                   do.log.scale=FALSE,n.cores=10)
    ###save(sc3.res,file=sprintf("%s/%s.SC3.RData",out.dir,sample.id))
}

Y.sd <- apply(Y, 1, sd)
Y.sd.sort <- sort(Y.sd,decreasing=TRUE)
Y.sd.df <- data.frame(geneID=names(Y.sd.sort),geneSymbol=entrezToXXX(names(Y.sd.sort)))
Y.sd.df <- cbind(Y.sd.df,Y.sd.sort)
write.table(Y.sd.df,sprintf("%s/%s.Y.sd.sort.txt",out.dir,sample.id),row.names = F,sep = "\t",quote = F)
doit(Y[names(Y.sd.sort)[1:nKeep],],extra=".uncorrected")

save.image(sprintf("%s/%s.all.RData",out.dir,sample.id))

