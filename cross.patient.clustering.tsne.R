#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-d", "--designFile", type="character", required=TRUE, help="design matrix file")
parser$add_argument("-i", "--inputFile", type="character", required=TRUE, help="input file")
parser$add_argument("-f", "--geneFile", type="character", help="gene list file")
parser$add_argument("-c", "--cellTypeColorFile", default="/WPS1/zhenglt/work/TCR_chunhong/data/CellType.color", 
                    type="character", help="cellTypeColorFile [default %(default)s]")
##parser$add_argument("-e", "--clonotypeFile", type="character", help="clonotype file")
parser$add_argument("-o", "--outDir", type="character", required=TRUE, help="output directory")
parser$add_argument("-s", "--sample", type="character", default="SAMPLE", help="sample id  [default %(default)s]")
###parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="whether verbose mode [default %(default)s]")
parser$add_argument("-l", "--log", action="store_true", default=FALSE, help="whether do log transform [default %(default)s]")
parser$add_argument("-r", "--center", action="store_true", default=FALSE, help="whether center genes by individual [default %(default)s]")
parser$add_argument("-y", "--dbscan", action="store_true", default=FALSE, help="do dbscan clustering [default %(default)s]")
parser$add_argument("-k", "--kNN", type="integer", default=5, help="kNN [default %(default)s]")
parser$add_argument("-e", "--eps", type="character", default="0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.5,3,3.5,4,5", help="eps.clus [default %(default)s]")
#parser$add_argument("-k", "--notSC3", action="store_true", default=FALSE, help="don't run SC3 pipeline [default %(default)s]")
parser$add_argument("-m", "--notFilter", action="store_true", default=FALSE, help="don't filtering genes [default %(default)s]")
parser$add_argument("-n", "--nKeep", type="integer", default=1500, help="use top nKeep genes [default %(default)s]")
parser$add_argument("-t", "--myseed", type="character", help="seed [default %(default)s]")
parser$add_argument("-w", "--geneIDType", type="character",default="entrezID", help="geneID type (entrezID, ensemblID) [default %(default)s]")
args <- parser$parse_args()

print(args)

designFile <- args$designFile
inputFile <- args$inputFile
geneFile <- args$geneFile
cellTypeColorFile <- args$cellTypeColorFile
clonotype.file <- args$clonotypeFile
out.dir <- args$outDir
sample.id <- args$sample
args.log <- args$log
args.center <- args$center
nKeep <- args$nKeep
#args.notSC3 <- args$notSC3
args.notFilter <- args$notFilter

args.dbscan <- args$dbscan
args.k <- args$kNN
args.eps <- args$eps
if(!is.null(args.eps)) { args.eps <- as.numeric(unlist(strsplit(args.eps,","))) }
args.myseed <- args$myseed
args.geneIDType <- args$geneIDType

#clonotype.strict.data <- read.clonotype(in.file = clonotype.file,ctype.col = "C_strict")

##### frequently used code
dir.create(out.dir,recursive = T,showWarnings = F)
source("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")
g.inputD <- processInput(designFile,cellTypeColorFile,inputFile,args.notFilter,geneFile,args.center,args.log)
myDesign <- g.inputD$myDesign
sampleTypeColor <- g.inputD$sampleTypeColor
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
if(args.geneIDType=="ensemblID"){ Y <- Y.newRowname }


doit.tsne <- function(Y,g.f,extra="",myseed=NULL,data.forDE=NULL,do.dbscan=args.dbscan,do.refine=T)
{
    tsne.res <- runTSNEAnalysis(Y[g.f,],sprintf("%s/%s%s.het.tSNE.run00",out.dir,sample.id,extra),
                    col.points = sampleTypeColor[as.character(myDesign$sampleType)],
                    legend=c(names(sampleTypeColor)[names(sampleTypeColor) %in% unique(as.character(myDesign$sampleType))],
                             unique(myDesign$libType)),
                    col.legend=c(sampleTypeColor[names(sampleTypeColor) %in% unique(as.character(myDesign$sampleType))],
                                 rep("black",length(unique(myDesign$libType)))),
                    pch=(as.numeric(as.factor(myDesign$patient))+19) %% 26,
                    pch.legend=c(rep(16,sum(names(sampleTypeColor) %in% unique(as.character(myDesign$sampleType)))),
                                 (seq_along(unique(myDesign$patient))+19) %% 26),
                    myseed=myseed,do.dbscan = args.dbscan,k = args.k,eps.clus=args.eps,do.clustering = F,data.forDE=data.forDE,do.scale = F)
   
   if(args.dbscan){
       ### chose the best solution
       #val.df <- as.data.frame.matrix(tsne.res$`30`$dbscan.clustering.validation.df)
       #val.df$avg.silwidth.rank <- rank(val.df[,"avg.silwidth"],na.last = F)
       #val.df$dunn.rank <- rank(val.df[,"dunn"],na.last = F)
       #val.df$score <- apply(val.df[,c("avg.silwidth.rank","dunn.rank")],1,mean)
       #val.df.best <- subset(val.df[order(val.df$score,decreasing = T),],clustered.sample/total.sample>0.85)
       #cat("INFO best solution (first-run)\n")
       #print(val.df.best)
       #write.table(val.df.best,sprintf("%s/%s%s.het.tSNE.run00.val.best.txt",out.dir,sample.id,extra),row.names = T,sep = "\t",quote = F)

       ###
       val.df.best <- tsne.res$`30`$dbscan.clustering.validation.df.best
       print(tsne.res$`30` %>% names)
       #print("TEST: YYY")
       #print(str(val.df.best))
       if(nrow(val.df.best)>0 && do.refine){
           eps.best <- val.df.best[1,"eps"]
           iter.g.f <- XXXToEntrez(tsne.res$`30`$dbscan.diffGene.list[[as.character(eps.best)]]$geneSymbol)
           iter.g.f <- head(iter.g.f[!is.na(iter.g.f)],n=nKeep)

           iter.tsne.res <<- runTSNEAnalysis(Y[iter.g.f,],sprintf("%s/%s%s.het.tSNE.refine",out.dir,sample.id,extra),
                            col.points = sampleTypeColor[as.character(myDesign$sampleType)],
                            legend=c(names(sampleTypeColor)[names(sampleTypeColor) %in% unique(as.character(myDesign$sampleType))],
                                     unique(myDesign$patient)),
                            col.legend=c(sampleTypeColor[names(sampleTypeColor) %in% unique(as.character(myDesign$sampleType))],
                                         rep("black",length(unique(myDesign$patient)))),
                            pch=(as.numeric(as.factor(myDesign$patient))+19) %% 26,
                            pch.legend=c(rep(16,sum(names(sampleTypeColor) %in% unique(as.character(myDesign$sampleType)))),
                                         (seq_along(unique(myDesign$patient))+19) %% 26),
                            myseed=myseed,do.dbscan = T,k = args.k,eps.clus=args.eps,do.clustering = F,data.forDE=data.forDE,do.scale = F)
            
           #val.df <- as.data.frame.matrix(iter.tsne.res$`30`$dbscan.clustering.validation.df)
           #val.df$avg.silwidth.rank <- rank(val.df[,"avg.silwidth"],na.last = F)
           #val.df$dunn.rank <- rank(val.df[,"dunn"],na.last = F)
           #val.df$score <- apply(val.df[,c("avg.silwidth.rank","dunn.rank")],1,mean)
           #val.df.best <- subset(val.df[order(val.df$score,decreasing = T),],clustered.sample/total.sample>0.85)
           #cat("INFO best solution (refine-run)\n")
           #print(val.df.best)
           #write.table(val.df.best,sprintf("%s/%s%s.het.tSNE.refine.val.best.txt",out.dir,sample.id,extra),row.names = T,sep = "\t",quote = F)
        }
   }
}

#### select variable genes
Y.sd <- apply(Y, 1, sd)
Y.sd.sort <- sort(Y.sd,decreasing=TRUE)
Y.sd.df <- data.frame(geneID=names(Y.sd.sort),geneSymbol=g.GNAME[names(Y.sd.sort)])
Y.sd.df <- cbind(Y.sd.df,Y.sd.sort)
write.table(Y.sd.df,sprintf("%s/%s.Y.sd.sort.txt",out.dir,sample.id),row.names = F,sep = "\t",quote = F)

#### rename rownames
#Y.newRowname <- Y
#newName <- g.GNAME[rownames(Y)]
#names(newName) <- rownames(Y)
#f.gname.na <- is.na(newName)
#newName[f.gname.na] <- names(newName)[f.gname.na]
#rownames(Y.newRowname) <- unname(newName)

g.f <- head(names(Y.sd.sort),n=nKeep)
doit.tsne(Y,g.f,extra="",myseed=args.myseed,data.forDE=Y.newRowname)

