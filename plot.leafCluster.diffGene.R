#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-i", "--inFile", type="character", required=TRUE, help="infile")
parser$add_argument("-a", "--sampleDescFile", type="character", required=TRUE, help="sample.desc.file")
parser$add_argument("-c", "--contrastFile", type="character", required=TRUE, help="contrast file")
parser$add_argument("-o", "--outputPrefix", type="character", required=TRUE, help="output prefix")
parser$add_argument("-z", "--binFile", type="character", help="binarized expression file")
parser$add_argument("-m", "--scoreMatrix", type="character", help="score matrix used for annotating samples")
parser$add_argument("-s", "--sample", type="character", default="SAMPLE", help="sample id  [default %(default)s]")
###parser$add_argument("-t", "--subtype", type="character", default="all", help="consider subtype (TTC,PTC, et.al)  [default %(default)s]")
parser$add_argument("-y", "--cellTypeColorFile", type="character", default="/WPS1/zhenglt/work/TCR_chunhong/data/CellType.color", help="cellTypeColorFile [default %(default)s]")
parser$add_argument("-l", "--clonotypeFile", type="character", help="clonotype file")
parser$add_argument("-k", "--groupBy", type="character",default="leafCluster", help="group by")
parser$add_argument("-n", "--nCores", type="integer", help="number of cpu cores to use [default NULL (auto detect)]")
parser$add_argument("-q", "--FDR", type="double", default=0.05, help="FDR threshold [default %(default)s]")
parser$add_argument("-f", "--FC", type="double", default=1, help="FC threshold [default %(default)s]")
parser$add_argument("-u", "--u", action="store_true", default=FALSE, help="unique gene [default %(default)s]")
parser$add_argument("-w", "--disableLog", action="store_true", default=FALSE, help="disable log transform [default %(default)s]")
parser$add_argument("-x", "--filterGene", action="store_true", default=FALSE, help="filter lowly expressed genes [default %(default)s]")
parser$add_argument("-b", "--disableGSymbol", action="store_true", default=FALSE, help="disable gene symbol transform [default %(default)s]")
parser$add_argument("-e", "--scale", action="store_true", default=FALSE, help="do scale [default %(default)s]")
parser$add_argument("-p", "--ensembl", action="store_true", default=FALSE, help="ensembl id? [default %(default)s]")
parser$add_argument("-d", "--pdfWidth", type="integer", default=18, help="pdf width [default %(default)s]")
args <- parser$parse_args()
print(args)

in.file <- args$inFile
contrastFile <- args$contrastFile
sample.desc.file <- args$sampleDescFile
out.prefix <- args$outputPrefix
sample.id <- args$sample
###sample.type <- args$sampleType
cellTypeColorFile <- args$cellTypeColorFile
clonotype.file <- args$clonotypeFile
opt.u <- args$u
do.scale <- args$scale
n.cores <- args$nCores
args.FDR <- args$FDR
args.FC <- args$FC
scoreMatrix.file <- args$scoreMatrix
disable.log <- args$disableLog
disable.gsymbol <- args$disableGSymbol
pdf.width <- args$pdfWidth
bin.exp.file <- args$binFile
args.filterGene <- args$filterGene
args.groupBy <- args$groupBy
args.ensembl <- args$ensembl
            
###args.FDR <- 0.05
###args.FC <- 1
###disable.gsymbol <- F
###disable.log <- F
###pdf.width <- 18
###n.cores <- 8
###in.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/gsnap_out/quantification/sizeFactor/P0322/P0322.all.countData.sfNormalized.txt.gz"
###contrastFile <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/byTCellType/clusteringResult/P0322/contrast.P0322"
###sample.desc.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/byTCellType/clusteringResult/P0322/clustering.P0322.iterative.leaf.txt"
###out.prefix <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/byTCellType/OUT.contrast/test/contrast.P0322.geneLevel.out"
###sample.id <- "P0322"
###cellTypeColorFile <- "/WPS1/zhenglt/work/TCR_chunhong/data/CellType.color"
###clonotype.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/tracer/OUT/P0322/filtered_TCR_summary/P0322.summary.cell.reassigneClonotype.methodChunhong.r.txt"
###scoreMatrix.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/byTCellType/gsva.out/P0322/P0322.bindea.gsva.txt"
###bin.exp.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/gsnap_out/quantification/binarized/P0322/P0322.RPKM.tab.binarized.txt.gz"

source("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")

suppressPackageStartupMessages(library("ComplexHeatmap"))
suppressPackageStartupMessages(library("circlize"))
suppressPackageStartupMessages(library("gridBase"))
suppressPackageStartupMessages(library("dendextend"))
suppressPackageStartupMessages(library("RColorBrewer"))
suppressPackageStartupMessages(library("gplots"))
suppressPackageStartupMessages(library("vioplot"))
suppressPackageStartupMessages(library("plotrix"))
suppressPackageStartupMessages(library("R.utils"))

loginfo("begin ...")

### contrast
myContrast<-read.table(contrastFile,header=T) 
head(myContrast)
### clonotype data
clonotype.data <- read.clonotype(in.file = clonotype.file,ctype.col = "C_strict")
#print(str(clonotype.data))

g.inputD <- processInput(designFile=sample.desc.file,cellTypeColorFile,inputFile=in.file,args.notFilter=F,geneFile=NULL,
                         args.center=F,args.log=disable.log,args.norm.exprs = F)
myDesign <- g.inputD$myDesign
sampleTypeColor <- g.inputD$sampleTypeColor
patient.col.list <- g.inputD$patient.col.list
Y <- g.inputD$Y
g.GNAME <- g.inputD$g.GNAME
##### 
in.table <- Y

print(dim(in.table))
print(in.table[1:4,1:8])

if(!is.null(bin.exp.file) && file.exists(bin.exp.file)){
    if(grepl("RData$",bin.exp.file,perl = T)){
        lenv <- loadToEnv(bin.exp.file)
        in.bin.table <- lenv[["exp.bin"]]
    }else{
        in.bin.table <- read.table(bin.exp.file,header = T,sep = "\t",check.names = F,stringsAsFactors = F,row.names = 1)
        in.bin.table <- in.bin.table[,-1]
    }
}else{
    in.bin.table <- NULL
}

#### sample data
sample.desc <- myDesign
#sample.desc <- read.delim(sample.desc.file,header = T,sep = "\t",row.names = "sample",check.names = F,stringsAsFactors = F)
#rownames(sample.desc) <- sample.desc$sample

#### scoreMatrix data
if(!is.null(scoreMatrix.file))
{
    scoreMatrix <- read.table(scoreMatrix.file,header = T,row.names = 1,check.names = F,stringsAsFactors = F)
    scoreM.col.fun <- colorRamp2(seq(-1,1,length=100),colorRampPalette(rev(brewer.pal(n = 7, name =  "RdBu")))(100), space="LAB")
    scoreM.legend.para.list <- list(TReg=list(color_bar="continuous",legend_width=unit(2, "cm"),legend_height=unit(4, "cm")),
                                exhaustion_G5=list(color_bar="continuous",legend_width=unit(2, "cm"),legend_height=unit(4, "cm")),
                                cytotoxic_G8=list(color_bar="continuous",legend_width=unit(2, "cm"),legend_height=unit(4, "cm")),
                                naive_G4=list(color_bar="continuous",legend_width=unit(2, "cm"),legend_height=unit(4, "cm"))
                                )
}else
{
    scoreM.legend.para.list <- NULL
}

runOneContrast <- function(l.ref,l.alt,l.name,test.sample.type=c("TTC","TTR","TTH"),FDR.THRESHOLD=0.05,FC.THRESHOLD=1){
    out.prefix <- sprintf("%s.%s",out.prefix,l.name)
    for(sample.type in test.sample.type){
        loginfo(sprintf("begin test for %s .vs. %s (%s)",l.ref,l.alt,sample.type))
        if(l.alt=="-"){
            grps <- unlist(strsplit(x = l.ref,split = ";",perl = T))
            print(grps)
            print(str(grps))
            #if(sample.type=="all")
            #    dat.g <- in.table[,rownames(subset(sample.desc,leafCluster %in% grps)),drop=F]
            #else
            #    dat.g <<- in.table[,rownames(subset(sample.desc,leafCluster %in% grps & sampleType==sample.type)),drop=F]
            
            dat.g <- data.frame()
            cmp.grp <- c()
            grps.list <- c()
            if(sample.type=="all"){
                for(i in seq_along(grps)){
                    cls.ref <- unlist(strsplit(x = grps[i],split = ",",perl = T))
                    grps.list[[sprintf("Grp%02d",i)]] <- c()
                    for(j in seq_along(cls.ref)){ 
                        ###tmp.dat <- in.table[,rownames(subset(sample.desc,leafCluster==cls.ref[j])),drop=F]
                        tmp.dat <- in.table[,rownames(sample.desc[ sample.desc[,args.groupBy]==cls.ref[j], ,drop=F])]
                        if(ncol(dat.g)==0) { 
                            dat.g <- tmp.dat
                        }else{
                            dat.g <- cbind(dat.g, tmp.dat) 
                        }
                        cmp.grp <- append(cmp.grp,rep(sprintf("Grp%02d",i),ncol(tmp.dat)))
                        grps.list[[sprintf("Grp%02d",i)]] <- append(grps.list[[sprintf("Grp%02d",i)]],colnames(tmp.dat))
                    }
                }
            }else{
                for(i in seq_along(grps)){
                    cls.ref <- unlist(strsplit(x = grps[i],split = ",",perl = T))
                    for(j in seq_along(cls.ref)){ 
                        ###tmp.dat <- in.table[,rownames(subset(sample.desc,leafCluster==cls.ref[j] & sampleType==sample.type)),drop=F]
                        tmp.dat <- in.table[,rownames(sample.desc[ sample.desc[,args.groupBy]==cls.ref[j] & sample.desc[,"sampleType"]==sample.type , ,drop=F])]
                        if(ncol(dat.g)==0) { 
                            dat.g <- tmp.dat
                        }else{
                            dat.g <- cbind(dat.g, tmp.dat) 
                        }
                        cmp.grp <- append(cmp.grp,rep(sprintf("Grp%02d",i),ncol(tmp.dat)))
                        grps.list[[sprintf("Grp%02d",i)]] <- append(grps.list[[sprintf("Grp%02d",i)]],colnames(tmp.dat))
                    }
                }
            }
            if(ncol(dat.g)<3)
            {
                loginfo(sprintf("Two few samples: ncol(dat.g)=%d",ncol(dat.g)))
                next
            }
            ##mgeneTest.out <<- runMultiGroupSpecificGeneTest(dat.g,sample.desc[colnames(dat.g),"leafCluster"],
            mgeneTest.out <<- runMultiGroupSpecificGeneTest(dat.g,cmp.grp,
                                                           sprintf("%s.%s.FDR%g",out.prefix,
                                                                  sample.type,100*FDR.THRESHOLD),
                                                           FDR.THRESHOLD=FDR.THRESHOLD,FC.THRESHOLD=FC.THRESHOLD,
                                                           verbose=T,n.cores = n.cores,gid.mapping=g.GNAME)
            g.list <- as.character(rownames(mgeneTest.out$aov.out.sig))
            if(!is.null(in.bin.table)){
                for(j in seq_along(grps.list))
                {
                    mgeneTest.out$aov.out.sig[, sprintf("cluster.prevalence.%s",names(grps.list)[j])] <- apply(in.bin.table[ g.list, grps.list[[j]], drop=F ],1,function(x){ sum(x>0)/length(x) })
                }
                cls.pre.names <- names(mgeneTest.out$aov.out.sig)[grepl("^cluster.prevalence.",names(mgeneTest.out$aov.out.sig),perl=T)]
                mgeneTest.out$aov.out.sig$specificClusterPrevalence <- as.numeric(apply(mgeneTest.out$aov.out.sig,1,function(x){ if(x["cluster.lable"]=="NA") return(NA) 
                                                                                     else return(x[sprintf("cluster.prevalence.%s",x["cluster.lable"])]) 
                                                                                    }))
                mgeneTest.out$aov.out.sig$otherClusterPrevalence <- as.numeric(apply(mgeneTest.out$aov.out.sig,1,function(x){ 
                                                                                  if(x["cluster.direction"]=="UP"){
                                                                                    return( max(x[cls.pre.names[!cls.pre.names %in% sprintf("cluster.prevalence.%s",x["cluster.lable"])]]) ) 
                                                                                  }else if(x["cluster.direction"]=="DOWN"){
                                                                                    return( min(x[cls.pre.names[!cls.pre.names %in% sprintf("cluster.prevalence.%s",x["cluster.lable"])]]) ) 
                                                                                  }else return(NA)
                                                                                }))
                write.table(mgeneTest.out$aov.out.sig,file = sprintf("%s.%s.FDR%g.aov.sig.prevalence.txt",out.prefix,sample.type,100*FDR.THRESHOLD),
                            quote = F,row.names = F,col.names = T,sep = "\t")
            }
            ### per gene comparison
            geneDetailDir=sprintf("%s.%s.aov.sig.noClusteringCol.FDR%g.GeneDetail",out.prefix,sample.type,100*FDR.THRESHOLD)
            dir.create(geneDetailDir,showWarnings = F,recursive = T)
            for(j in seq_along(g.list)){
                loginfo(sprintf("plot gene: %s",g.GNAME[g.list[j]]))
                pdf(sprintf("%s/%s.pdf",geneDetailDir,g.GNAME[g.list[j]]),width = 8,height = 8)
                par(mar=c(5,5,4,2),cex.lab=1.8,cex.main=1.8)
                ylim <- pretty(range(dat.g[g.list[j],]),n=5)
                plot(0,0,type="n",xlim=c(0.5,length(grps.list)+0.5), ylim=ylim[c(1,length(ylim))],
                     xaxt = 'n', xlab ="",ylab="Expression",main=sprintf("%s",g.GNAME[g.list[j]]))
                print(str(dat.g))
                for(jj in seq_along(grps.list)){
                    #print(str(na.omit(dat.g[g.list[j],grps.list[[jj]],drop=T])))
                    #tryCatch(vioplot(na.omit(dat.g[g.list[j],grps.list[[jj]],drop=T]), 
                    tryCatch(vioplot(unlist(na.omit(dat.g[g.list[j],grps.list[[jj]],drop=T])), 
                                     at = jj, add = T, col = auto.colSet(length(grps.list),name = "Dark2")[jj]), 
                             error = function(e) { loginfo(sprintf("vioplot() catch exception: gene %s",g.GNAME[g.list[j]])); e })
                }
                staxlab(1,at = seq_along(grps.list),labels=names(grps.list),srt=NA, cex=1.3,adj=0.5,top.line=2)
                ##print("##### TEST #####")
                ##print(str(mgeneTest.out$aov.out.sig))
                mtext(sprintf("FDR: %4.4f",mgeneTest.out$aov.out.sig[g.list[j],"F.adjp"]),side = 1,line = -2.1,adj = 0.05,cex=1.2)
                dev.off()
                if(j>50) break
            }
            ### heatmap
            do.heatmap.plot <- function(g.list,extra="",kk=1,do.scale=T,verbose=T)
            {
                dat.g.sig <- dat.g[g.list,,drop=F]
                dat.tmp <- data.frame()
                for(i in seq_along(grps.list))
                {
                    dat.tmp.1 <- dat.g.sig[,grps.list[[i]],drop=F]
                    if(ncol(dat.tmp)==0){
                        dat.tmp <- dat.tmp.1[,hclust(dist(t(dat.tmp.1)))$order,drop=F]
                    }else{
                        dat.tmp <- cbind(dat.tmp,dat.tmp.1[,hclust(dist(t(dat.tmp.1)))$order,drop=F])
                    }
                }
                dat.g.sig <- dat.tmp
                annDF <- data.frame(cluster=structure(cmp.grp,names=colnames(dat.g.sig)),
                                    clonotype=clonotype.data$ctypeII[colnames(dat.g.sig)],stringsAsFactors = F)
                annDF[is.na(annDF)] <- "NA"
                annCol <- list(cluster=structure(auto.colSet(length(grps.list),name = "Dark2"),names=unique(cmp.grp)),
                               clonotype=clonotype.data$colII)
                print(head(annDF))
                if(!is.null(scoreMatrix.file)){
                    ##annDF <- cbind(annDF,t(scoreMatrix[c("TReg","exhaustion_G5","cytotoxic_G8","naive_G4","CXCL13_G1"),colnames(dat.g.sig)]))
                    annDF <- cbind(annDF,t(scoreMatrix[c("TReg","exhaustion_G5","cytotoxic_G8","naive_G4"),colnames(dat.g.sig)]))
                    annCol <- c(annCol,list(TReg=scoreM.col.fun,
                                   exhaustion_G5=scoreM.col.fun,
                                   cytotoxic_G8=scoreM.col.fun,
                                   naive_G4=scoreM.col.fun) )
                                   ####CXCL13_G1=scoreM.col.fun) )
                }else{
                    print("nothing to do")
                    #annDF <- NULL
                    #annCol <- NULL
                }
                runHierarchicalClusteringAnalysis(dat.g.sig,mytitle = sample.id,
                                              pdf.width=pdf.width,pdf.height=15,do.clustering.col=F,
                                              sprintf("%s.%s.aov.sig.noClusteringCol.FDR%g%s",
                                                      out.prefix,sample.type,100*FDR.THRESHOLD,extra), 
                                              sampleType=sample.desc[colnames(dat.g.sig),"sampleType"], 
                                              colSet=sampleTypeColor,
                                              ann.extra.df = annDF,
                                              ann.extra.df.col = annCol,
                                              ann.bar.height = 0.6, annotation_legend_param=scoreM.legend.para.list,
                                              k.row=kk,clonotype.col=NULL,ntop=NULL, 
                                              row.names.original=disable.gsymbol,
                                              complexHeatmap.use=TRUE,verbose=FALSE,do.scale=do.scale,save.obj=F,gid.mapping=g.GNAME)
                                              ###k.row=kk,clonotype.col=clonotype.data,ntop=NULL, 
                if(verbose){
                    dat.g.sig.df <- data.frame(geneID=rownames(dat.g.sig))
                    ##dat.g.sig.df$geneSymbol=entrezToXXX(dat.g.sig.df$geneID)
                    dat.g.sig.df$geneSymbol=g.GNAME[dat.g.sig.df$geneID]
                    dat.g.sig.df <- cbind(dat.g.sig.df, dat.g.sig)
                    write.table(dat.g.sig.df,file = sprintf("%s.%s.aov.sig.noClusteringCol.FDR%g%s.txt",
                                                          out.prefix,sample.type,100*FDR.THRESHOLD,extra),
                                    quote = F,sep = "\t",row.names = F,col.names = T)
                }

            }
            ## aov genes
            g.list <- as.character(rownames(mgeneTest.out$aov.out.sig))
            if(length(g.list)>=3){ 
                do.heatmap.plot(g.list,extra="") 
                if(length(g.list)>30){ do.heatmap.plot(head(g.list,n=30),extra=".top30") }
            }
            ## cluster specific genes
            g.list <- as.character(rownames(subset(mgeneTest.out$aov.out.sig,is.clusterSpecific==TRUE & cluster.direction!="INCONSISTANT" )))
            if(length(g.list)>=3){ do.heatmap.plot(g.list,extra=".cluster.specific", kk=1) }
            g.list <- c()
            for(t.grp in unique(cmp.grp)){
                g.list <- append(g.list, head(as.character(rownames(subset(mgeneTest.out$aov.out.sig,is.clusterSpecific==TRUE & cluster.direction!="INCONSISTANT" & cluster.lable==t.grp ))),n=10))
            }
            if(length(g.list)>=3){ do.heatmap.plot(g.list,extra=".cluster.specific.top10", kk=1) }
            ## cluster specific genes & up-regulated
            g.list <- as.character(rownames(subset(mgeneTest.out$aov.out.sig,is.clusterSpecific==TRUE & cluster.direction=="UP")))
            if(length(g.list)>=3){ do.heatmap.plot(g.list,extra=".cluster.specific.UP", kk=length(unique(subset(mgeneTest.out$aov.out.sig[g.list,],cluster.lable!="NA")$cluster.lable))) }
            g.list <- c()
            for(t.grp in unique(cmp.grp)){
                g.list <- append(g.list, head(as.character(rownames(subset(mgeneTest.out$aov.out.sig,is.clusterSpecific==TRUE & cluster.direction=="UP" & cluster.lable==t.grp))),n=10))
            }
            if(length(g.list)>=3){ do.heatmap.plot(g.list,extra=".cluster.specific.UP.top10", kk=length(unique(subset(mgeneTest.out$aov.out.sig[g.list,],cluster.lable!="NA")$cluster.lable))) }
            ## cluster specific genes & down-regulated
            g.list <- as.character(rownames(subset(mgeneTest.out$aov.out.sig,is.clusterSpecific==TRUE & cluster.direction=="DOWN")))
            if(length(g.list)>=3){ do.heatmap.plot(g.list,extra=".cluster.specific.DOWN", kk=length(unique(subset(mgeneTest.out$aov.out.sig[g.list,],cluster.lable!="NA")$cluster.lable))) }
            g.list <- c()
            for(t.grp in unique(cmp.grp)){
                g.list <- append(g.list, head(as.character(rownames(subset(mgeneTest.out$aov.out.sig,is.clusterSpecific==TRUE & cluster.direction=="DOWN" & cluster.lable==t.grp))),n=10))
            }
            if(length(g.list)>=3){ do.heatmap.plot(g.list,extra=".cluster.specific.DOWN.top10", kk=length(unique(subset(mgeneTest.out$aov.out.sig[g.list,],cluster.lable!="NA")$cluster.lable))) }

            if(!is.null(in.bin.table)){
                ## cluster specific genes & up-regulated, also with prevalence filter
                g.list <- as.character(rownames(subset(mgeneTest.out$aov.out.sig,
                                                       is.clusterSpecific==TRUE & cluster.direction=="UP" 
                                                       & specificClusterPrevalence > 0.4 & otherClusterPrevalence < 0.4 )))
                par.kk <- length(unique(subset(mgeneTest.out$aov.out.sig[g.list,],cluster.lable!="NA")$cluster.lable))
                if(length(g.list)>=3){ 
                    do.heatmap.plot(g.list,extra=".cluster.specific.UP.prevalence", kk=par.kk) 
                    do.heatmap.plot(g.list,extra=".cluster.specific.UP.prevalence.scaleF", kk=par.kk,do.scale = F,verbose = F)
                    if(length(g.list)>30){
                        do.heatmap.plot(head(g.list,n=30),extra=".cluster.specific.UP.prevalence.top30", kk=par.kk) 
                        do.heatmap.plot(head(g.list,n=30),extra=".cluster.specific.UP.prevalence.scaleF.top30", kk=par.kk,do.scale = F,verbose = F)
                    }
                }
                g.list <- c()
                for(t.grp in unique(cmp.grp)){
                    g.list <- append(g.list, head(as.character(rownames(subset(mgeneTest.out$aov.out.sig, 
                                                                               is.clusterSpecific==TRUE 
                                                                               & cluster.direction=="UP" 
                                                                               & specificClusterPrevalence > 0.4 
                                                                               & otherClusterPrevalence < 0.4 
                                                                               & cluster.lable==t.grp ))),n=10))
                }
                par.kk <- length(unique(subset(mgeneTest.out$aov.out.sig[g.list,],cluster.lable!="NA")$cluster.lable))
                if(length(g.list)>=3){ 
                    do.heatmap.plot(g.list,extra=".cluster.specific.UP.prevalence.eachClusterTop10", kk=par.kk) 
                    do.heatmap.plot(g.list,extra=".cluster.specific.UP.prevalence.scaleF.eachClusterTop10", 
                                    kk=par.kk,do.scale = F,verbose = F)
                }
                ## cluster specific genes & down-regulated, also with prevalence filter
                g.list <- as.character(rownames(subset(mgeneTest.out$aov.out.sig,
                                                       is.clusterSpecific==TRUE & cluster.direction=="DOWN" 
                                                       & specificClusterPrevalence < 0.15 & otherClusterPrevalence > 0.15 )))
                par.kk <- length(unique(subset(mgeneTest.out$aov.out.sig[g.list,],cluster.lable!="NA")$cluster.lable))
                if(length(g.list)>=3){ 
                       do.heatmap.plot(g.list,extra=".cluster.specific.DOWN.prevalence", kk=par.kk)
                       do.heatmap.plot(g.list,extra=".cluster.specific.DOWN.prevalence.scaleF", kk=par.kk,do.scale = F,verbose = F)
                   if(length(g.list)>30){
                       do.heatmap.plot(head(g.list,n=30),extra=".cluster.specific.DOWN.prevalence.top30", kk=par.kk)
                       do.heatmap.plot(head(g.list,n=30),extra=".cluster.specific.DOWN.prevalence.scaleF.top30", kk=par.kk,do.scale = F,verbose = F)
                   } 
                }
                g.list <- c()
                for(t.grp in unique(cmp.grp)){
                    g.list <- append(g.list, head(as.character(rownames(subset(mgeneTest.out$aov.out.sig, 
                                                                               is.clusterSpecific==TRUE 
                                                                               & cluster.direction=="DOWN" 
                                                                               & specificClusterPrevalence < 0.15 
                                                                               & otherClusterPrevalence > 0.15 
                                                                               & cluster.lable==t.grp ))),n=10))
                }
                par.kk <- length(unique(subset(mgeneTest.out$aov.out.sig[g.list,],cluster.lable!="NA")$cluster.lable))
                if(length(g.list)>=3){ 
                    do.heatmap.plot(g.list,extra=".cluster.specific.DOWN.prevalence.eachClusterTop10", kk=par.kk) 
                    do.heatmap.plot(g.list,extra=".cluster.specific.DOWN.prevalence.scaleF.eachClusterTop10", 
                                    kk=par.kk,do.scale = F,verbose = F)
                }
            }

        }else{
            cls.ref <- unlist(strsplit(x = l.ref,split = ";",perl = T))
            cls.alt <- unlist(strsplit(x = l.alt,split = ";",perl = T))
            if(sample.type=="all")
            {
                dat.g1 <- c()
                dat.g2 <- c()
                ##for(i in seq_along(cls.ref)) { dat.g1 <- cbind(dat.g1,in.table[,rownames(subset(sample.desc,leafCluster==cls.ref[i])),drop=F]) }
                ##for(i in seq_along(cls.alt)) { dat.g2 <- cbind(dat.g2,in.table[,rownames(subset(sample.desc,leafCluster==cls.alt[i])),drop=F]) }
                for(i in seq_along(cls.ref)) { dat.g1 <- cbind(dat.g1,in.table[,rownames(sample.desc[ sample.desc[,args.groupBy]==cls.ref[i], ,drop=F])]) }
                for(i in seq_along(cls.alt)) { dat.g2 <- cbind(dat.g2,in.table[,rownames(sample.desc[ sample.desc[,args.groupBy]==cls.alt[i], ,drop=F])]) }
            }else{
                dat.g1 <- c()
                dat.g2 <- c()
                ##for(i in seq_along(cls.ref)) { dat.g1 <- cbind(dat.g1,in.table[,rownames(subset(sample.desc,leafCluster==cls.ref[i] & sampleType==sample.type)),drop=F]) }
                ##for(i in seq_along(cls.alt)) { dat.g2 <- cbind(dat.g2,in.table[,rownames(subset(sample.desc,leafCluster==cls.alt[i] & sampleType==sample.type)),drop=F]) }
                for(i in seq_along(cls.ref)) { dat.g1 <- cbind(dat.g1,in.table[,rownames(sample.desc[ sample.desc[,args.groupBy]==cls.ref[i] & sample.desc[,"sampleType"]==sample.type, ,drop=F])]) }
                for(i in seq_along(cls.alt)) { dat.g2 <- cbind(dat.g2,in.table[,rownames(sample.desc[ sample.desc[,args.groupBy]==cls.alt[i] & sample.desc[,"sampleType"]==sample.type, ,drop=F])]) }
            }
            print(dim(dat.g1))
            print(dim(dat.g2))
            if(ncol(dat.g1)<3 || ncol(dat.g2)<3)
            {
                loginfo(sprintf("Two few samples: ncol(dat.g1)=%d, ncol(dat.g2)=%d",ncol(dat.g1),ncol(dat.g2)))
                next
            }
            ttest.out <<- runTTest(dat.g1,dat.g2,sprintf("%s.%s.%s.%s.FDR%g",out.prefix,"Group1","Group2",sample.type,100*FDR.THRESHOLD),FDR.THRESHOLD,FC.THRESHOLD,verbose=T,n.cores = n.cores,gid.mapping=g.GNAME)
            ## for genes diff significantly
            g.list.1 <- as.character(rownames(ttest.out$ttest.out.sig))
            ## add prevalence
            if(!is.null(in.bin.table)){
                ttest.out$ttest.out.sig$x.prevalence <- apply(in.bin.table[g.list.1,colnames(dat.g1),drop=F],1,function(x){ sum(x>0)/length(x) })
                ttest.out$ttest.out.sig$y.prevalence <- apply(in.bin.table[g.list.1,colnames(dat.g2),drop=F],1,function(x){ sum(x>0)/length(x) })
                write.table(ttest.out$ttest.out.sig,file = sprintf("%s.ttest.sig.prevalence.txt",out.prefix),quote = F,row.names = F,col.names = T,sep = "\t")
            }
            ### per gene comparison
            geneDetailDir=sprintf("%s.%s.%s.%s.FDR%g.GeneDetail",out.prefix,"Group1","Group2",sample.type,100*FDR.THRESHOLD)
            dir.create(geneDetailDir,showWarnings = F,recursive = T)
            for(j in seq_along(g.list.1)){
                loginfo(sprintf("plot gene: %s",g.GNAME[g.list.1[j]]))
                pdf(sprintf("%s/%s.pdf",geneDetailDir,g.GNAME[g.list.1[j]]),width = 8,height = 8)
                par(mar=c(5,5,4,2),cex.lab=1.8,cex.main=1.8)
                ylim <- pretty(range(c(dat.g1[g.list.1[j],],dat.g2[g.list.1[j],])),n=5)
                plot(0,0,type="n",xlim=c(0.5,2.5), ylim=ylim[c(1,length(ylim))],
                     xaxt = 'n', xlab ="",ylab="Expression",main=sprintf("%s",g.GNAME[g.list.1[j]]))
                tryCatch(vioplot(na.omit(dat.g1[g.list.1[j],]), at = 1, add = T, col = auto.colSet(2,name = "Dark2")[1]), 
                         error = function(e) { loginfo(sprintf("vioplot() catch exception: gene %s",g.GNAME[g.list.1[j]])); e })
                tryCatch(vioplot(na.omit(dat.g2[g.list.1[j],]), at = 2, add = T, col = auto.colSet(2,name = "Dark2")[2]), 
                         error = function(e) { loginfo(sprintf("vioplot() catch exception: gene %s",g.GNAME[g.list.1[j]])); e })
                staxlab(1,at = 1:2,labels=c("Group1","Group2"),srt=NA, cex=1.3,adj=0.5,top.line=2)
                mtext(sprintf("FDR: %4.4f",ttest.out$ttest.out.sig[g.list.1[j],"p.adj"]),side = 1,line = -2.1,adj = 0.05,cex=1.2)
                mtext(sprintf("FC: %4.4f",ttest.out$ttest.out.sig[g.list.1[j],"fc"]),side = 1,line = -1.1,adj = 0.05,cex=1.2)
                dev.off()
                if(j>50) break
            }
            ### heatmap
            do.ttest.heatmap.plot <- function(g.list.1,extra="",do.clustering.col=F,kk=1,do.scale=T,verbose=T){
                if(length(g.list.1)>=3){
                    dat.g1.1 <- dat.g1[g.list.1,,drop=F]
                    dat.g2.1 <- dat.g2[g.list.1,,drop=F]
                    dat.1 <- cbind(dat.g1.1[,hclust(dist(t(dat.g1.1)))$order,drop=F], dat.g2.1[,hclust(dist(t(dat.g2.1)))$order,drop=F])

                    annDF <- data.frame(cluster=structure( c(rep("Group1",ncol(dat.g1)),rep("Group2",ncol(dat.g2))),names=colnames(dat.1) ),
                                        clonotype=clonotype.data$ctypeII[colnames(dat.1)],stringsAsFactors = F)
                    annDF[is.na(annDF)] <- "NA"
                    annCol <- list(cluster=structure(auto.colSet(2,name = "Dark2"),names=c("Group1","Group2")),
                                   clonotype=clonotype.data$colII)
                    if(!is.null(scoreMatrix.file)){
                        ###annDF <- cbind(annDF,t(scoreMatrix[c("TReg","exhaustion_G5","cytotoxic_G8","naive_G4","CXCL13_G1"),c(colnames(dat.g1),colnames(dat.g2))]))
                        ###annDF <- cbind(annDF,t(scoreMatrix[c("TReg","exhaustion_G5","cytotoxic_G8","naive_G4","CXCL13_G1"),colnames(dat.1)]))
                        annDF <- cbind(annDF,t(scoreMatrix[c("TReg","exhaustion_G5","cytotoxic_G8","naive_G4"),colnames(dat.1)]))
                        annCol <- c(annCol,list(TReg=scoreM.col.fun,
                                       exhaustion_G5=scoreM.col.fun,
                                       cytotoxic_G8=scoreM.col.fun,
                                       naive_G4=scoreM.col.fun))
                                       ###CXCL13_G1=scoreM.col.fun) )
                    }else{
                        annDF <- NULL
                        annCol <- NULL
                    }
                    runHierarchicalClusteringAnalysis(dat.1,mytitle = sample.id,
                                                      pdf.width=pdf.width,pdf.height=15,
                                                      sprintf("%s.%s.%s.%s.sig.FDR%g%s%s",out.prefix,"Group1","Group2",sample.type,100*FDR.THRESHOLD,
                                                              extra,ifelse(do.clustering.col,".ColT",".ColF")), 
                                                      do.clustering.col=do.clustering.col,
                                                      sampleType=sample.desc[colnames(dat.1),"sampleType"], 
                                                      colSet=sampleTypeColor,
                                                      ann.extra.df = annDF,
                                                      ann.extra.df.col = annCol,
                                                      ann.bar.height = 0.6,
                                                      k.row=kk,
                                                      clonotype.col=NULL,ntop=NULL,
                                                      row.names.original=disable.gsymbol, annotation_legend_param=scoreM.legend.para.list,
                                                      complexHeatmap.use=TRUE,verbose=FALSE,do.scale=do.scale,gid.mapping=g.GNAME)
                    ### output gene list to text file
                    if(verbose){
                        dat.1.df <- data.frame(geneID=rownames(dat.1))
                        ##dat.1.df$geneSymbol=entrezToXXX(dat.1.df$geneID)
                        dat.1.df$geneSymbol=g.GNAME[dat.1.df$geneID]
                        dat.1.df <- cbind(dat.1.df,dat.1)
                        write.table(dat.1.df,file = sprintf("%s.%s.%s.%s.sig.FDR%g%s%s.txt",out.prefix,"Group1","Group2",sample.type,100*FDR.THRESHOLD,
                                                            extra,ifelse(do.clustering.col,".ColT",".ColF")),
                                    quote = F,sep = "\t",row.names = F,col.names = T)
                    }
                }
            }
            do.ttest.heatmap.plot(g.list.1,extra=".scaleT",do.clustering.col=F,kk=length(unique(ttest.out$ttest.out.sig$fc>0)),do.scale = T)
            do.ttest.heatmap.plot(g.list.1,extra=".scaleF",do.clustering.col=F,kk=length(unique(ttest.out$ttest.out.sig$fc>0)),do.scale = F,verbose=F)
            if(length(g.list.1)>30){
                do.ttest.heatmap.plot(head(g.list.1,n=30),extra=".scaleT.top30",do.clustering.col=F,kk=length(unique(ttest.out$ttest.out.sig$fc>0)),do.scale = T)
                do.ttest.heatmap.plot(head(g.list.1,n=30),extra=".scaleF.top30",do.clustering.col=F,kk=length(unique(ttest.out$ttest.out.sig$fc>0)),do.scale = F,verbose=F)
            }
            if(!is.null(in.bin.table)){
                tmp.var <- subset(ttest.out$ttest.out.sig, (x.prevalence > 0.4 & y.prevalence < 0.15) | (y.prevalence > 0.4 & x.prevalence < 0.15) )
                do.ttest.heatmap.plot(rownames(tmp.var),extra=".prevalence.scaleT",do.clustering.col=T,kk=1,do.scale = T)
                do.ttest.heatmap.plot(rownames(tmp.var),extra=".prevalence.scaleF",do.clustering.col=T,kk=1,do.scale = F,verbose=F)
                do.ttest.heatmap.plot(rownames(tmp.var),extra=".prevalence.scaleT",do.clustering.col=F,kk=1,do.scale = T,verbose=F)
                do.ttest.heatmap.plot(rownames(tmp.var),extra=".prevalence.scaleF",do.clustering.col=F,kk=1,do.scale = F,verbose=F)
                ###loginfo("#### TEST ####")
                ###print(head(tmp.var))
                ###print(dim(tmp.var))
                ###print(str(tmp.var))
                ###print(class(tmp.var))
                if(nrow(tmp.var)>30){
                    do.ttest.heatmap.plot(head(rownames(tmp.var),n=30),extra=".prevalence.scaleT.top30",do.clustering.col=F,kk=1,do.scale = T,verbose=F)
                }
                g.list.1 <- c()
                g.list.1 <- append(g.list.1,head((rownames(subset(ttest.out$ttest.out.sig, 
                                                                  (x.prevalence > 0.4 & y.prevalence < 0.15) & x.mean > y.mean ))),n=10))
                g.list.1 <- append(g.list.1,head((rownames(subset(ttest.out$ttest.out.sig, 
                                                                  (y.prevalence > 0.4 & x.prevalence < 0.15) & y.mean > x.mean ))),n=10))
                if(length(g.list.1>3)){
                    do.ttest.heatmap.plot(g.list.1,extra=".prevalence.scaleT.eachClusterTop10",do.clustering.col=F,kk=1,do.scale = T,verbose=F)
                }
            }
            
        }
        ### run functional enrichment
        loginfo(sprintf("run successful for %s .vs. %s (%s)",l.ref,l.alt,sample.type))
    }
}
###invisible(apply(myContrast,1,function(x){ runOneContrast(x[1],x[2],FDR.THRESHOLD = args.FDR,FC.THRESHOLD = args.FC) } ))
##invisible(apply(myContrast,1,function(x){ runOneContrast(x["ref"],x["alt"],x["comp_name"],test.sample.type = c("all"),FDR.THRESHOLD = args.FDR,FC.THRESHOLD = args.FC) } ))
invisible(apply(myContrast,1,function(x){ 
                    t.subtype <- unlist(strsplit(x = x["subtype"],split = ",",perl = T))
                    runOneContrast(x["ref"],x["alt"],x["comp_name"],test.sample.type = t.subtype,FDR.THRESHOLD = args.FDR,FC.THRESHOLD = args.FC) 
                } ))
loginfo("end.")
