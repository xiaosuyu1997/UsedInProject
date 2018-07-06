#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("plotrix"))
suppressPackageStartupMessages(library("survival"))
suppressPackageStartupMessages(library("beanplot"))
suppressPackageStartupMessages(require("circlize"))
suppressPackageStartupMessages(require("RColorBrewer"))
suppressPackageStartupMessages(library("vioplot"))

parser <- ArgumentParser()
parser$add_argument("-i", "--inFile", type="character", required=TRUE, help="infile")
parser$add_argument("-o", "--outPrefix", type="character", required=TRUE, help="outprefix")
parser$add_argument("-d", "--dFile", type="character", required=TRUE, help="design file")
parser$add_argument("-a", "--clinicalFile", type="character", help="clinical file")
parser$add_argument("-l", "--geneSetFile", type="character", required=TRUE, help="geneSetFile")
parser$add_argument("-s", "--sample", type="character",default="SAMPLE", required=TRUE, help="sample [default %(default)s]")
parser$add_argument("-u", "--cluster", action="store_true", default=FALSE, help="do clustering [default %(default)s]")
parser$add_argument("-b", "--scale", action="store_true", default=FALSE, help="do scale in heatmap plot [default %(default)s]")
parser$add_argument("-n", "--norm", action="store_true", default=FALSE, help="norma ssgsea as Barbie et al. (2009)  [default %(default)s]")
parser$add_argument("-c", "--colClustering", action="store_true", default=FALSE, help="do column clustering in heatmap plot [default %(default)s]")
parser$add_argument("-y", "--cellTypeColorFile", default="/WPS1/zhenglt/work/TCR_chunhong/data/CellType.color", type="character", help="cellTypeColorFile")
args <- parser$parse_args()
print(args)

in.file <- args$inFile
design.file <- args$dFile
geneSet.file <- args$geneSetFile
out.prefix <- args$outPrefix
sample.id <- args$sample
cellTypeColorFile <- args$cellTypeColorFile
do.clustering <- args$cluster
do.scale <- args$scale
do.colClustering <- args$colClustering
clinical.file <- args$clinicalFile
args.norm <- args$norm

source("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")

#args.norm <- T
#do.scale <- T
#do.colClustering <- T
#do.clustering <- T
#in.file <- "/WPSnew/zhenglt/work/TCR_chunhong/TCGA/quantification/mdata/TPM.KIRC.txt.gz"
#design.file <- "/WPSnew/zhenglt/work/TCR_chunhong/TCGA/quantification/mdata/TPM.KIRC.design.excludeNormal.txt"
#geneSet.file <- "/WPS1/zhenglt/work/TCR_chunhong/dataset/Bindea.markerGene.unique.immune.target.Angiogenesis.txt"
##geneSet.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/geneSignature/msigdb.h_cp_bp.v5.0.entrez.gmt"
#out.prefix <- "/WPSnew/zhenglt/work/TCR_chunhong/TCGA/gsva/test"
#sample.id <- "KIRC"
#cellTypeColorFile <- "/WPSnew/zhenglt/work/TCR_chunhong/TCGA/gsva/CellType.color"
#clinical.file <- "/WPSnew/zhenglt/work/TCR_chunhong/TCGA/clinical/TCGA.clinical.patient.slim.txt"

## clinical file
clinical.table <- read.table(clinical.file,header = T,check.names = F,stringsAsFactors = F,sep = "\t",row.names = "patient_barcode")
## cell type color 
sampleTypeColor <- read.SampleTypeColor(cellTypeColorFile)
## read in gene set data
if(grepl(pattern = "\\.gmt$",geneSet.file,perl = T)){
    print("gmt format")
    geneSet.list <- readGMT(file = geneSet.file)$gSet
}else{
    geneSet.table <- read.delim(geneSet.file,header = T,check.names = F,stringsAsFactors = F,sep = "\t")
    #geneSet.table <- geneSet.table[grepl("_G|Tcm|Tem|TFH|Th17_cells|Th1_cells|Th2_cells|TReg",geneSet.table$cellType,perl=T),]
    geneSet.list <- lapply(unique(geneSet.table$cellType),function(x){
                           as.character(subset(geneSet.table,cellType==x)$geneID)
    })
    names(geneSet.list) <- unique(geneSet.table$cellType)
}

## read in design file
myDesign <- read.table(design.file,header=T,row.names="sample",check.names=F,stringsAsFactors = F)

## read in exp data
in.table <- read.table(in.file,header = T,sep = "\t",check.names = F,stringsAsFactors = F,row.names = 1)
in.table <- in.table[,-1]

f <- intersect(rownames(myDesign),colnames(in.table))
myDesign <- myDesign[f,,drop=F]
in.table <- in.table[,f,drop=F]
in.table <- log2(in.table+1)
clinical.table.slim <- clinical.table[myDesign$patient,c("project","grade","ajcc_stage_pathologic_stage","ajcc_stage_tnm_categories","OS","vital_status")]


check.expressionProfile.geneSet <- function(dat,g.list,T.THRESHOLD=0.02)
{
    dat.check <- dat[g.list,,drop=F]
    apply(dat.check,1,function(v){ sum(v > 0)/length(v) >= T.THRESHOLD })
}
## f.bg: gene list used for GSVA
f.bg <- check.expressionProfile.geneSet(in.table,rownames(in.table))
## geneSet.validate.list: gene set used for GSVA
geneSet.validate.list <- lapply(names(geneSet.list),function(x){ 
                                f.gset <- check.expressionProfile.geneSet(in.table,geneSet.list[[x]])
                                geneSet.list[[x]][f.gset]
})
names(geneSet.validate.list) <- names(geneSet.list)

gsva.out <- run.GSVA(as.matrix(in.table[f.bg,,drop=F]),geneSet.validate.list,gsva.method = "ssgsea",bool.rnaseq=FALSE,ncores=4,ssgsea.norm=args.norm)
gsva.table <- gsva.out
#gsva.table <- gsva.out$es.obs
z.col <- colorRamp2(seq(-3,3,length=100),colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100), space="LAB")
inf.score.out <- calInfiltrationScore(gsva.table,use.scale=F)
###
clinical.table.toAnalyse <- cbind(myDesign,clinical.table.slim,inf.score.out$score)
clinical.table.toAnalyse$CD8_Treg_Ratio_Group <- "LOW"
clinical.table.toAnalyse$CD8_Treg_Ratio_Group[ clinical.table.toAnalyse$CD8_Treg_Ratio > median(clinical.table.toAnalyse$CD8_Treg_Ratio) ] <- "HIGHT"
clinical.table.toAnalyse$Th17_Th2_Ratio_Group <- "LOW"
clinical.table.toAnalyse$Th17_Th2_Ratio_Group[ clinical.table.toAnalyse$Th17_Th2_Ratio > median(clinical.table.toAnalyse$Th17_Th2_Ratio) ] <- "HIGHT"

cat(sprintf("median of CD8/Treg: %4.4f\n",median(clinical.table.toAnalyse$CD8_Treg_Ratio)))
cat(sprintf("median of Th17/Th2: %4.4f\n",median(clinical.table.toAnalyse$Th17_Th2_Ratio)))
out.df <- data.frame(sampleID=rownames(clinical.table.toAnalyse))
out.df <- cbind(out.df,clinical.table.toAnalyse)
write.table(out.df,file = sprintf("%s.gsva.clinical.txt",out.prefix),quote = F,row.names = F,col.names = T,sep = "\t")

vv <- data.frame(vname=c("CD8_Treg_Ratio","Th17_Th2_Ratio"),ann.name=c("CD8/Treg","Th17/Th2"),g.name=c("CD8_Treg_Ratio_Group","Th17_Th2_Ratio_Group"),stringsAsFactors = F)

### OS analysis
clinical.table.toAnalyse.noOSNA <- subset(clinical.table.toAnalyse, !is.na(OS))
attach(clinical.table.toAnalyse.noOSNA)
pdf(sprintf("%s.survival.pdf",out.prefix),width=8,height=8)
par(mar=c(5,5,4,2),cex.lab=1.5,cex.main=1.5)
for(jj in 1:nrow(vv))
{
    surv.obj <- Surv(OS, vital_status=="Dead")
    survdiff.res <- survdiff(as.formula(sprintf("surv.obj ~ %s",vv$g.name[jj])))
    survdiff.res.p <- pchisq(survdiff.res$chisq,1,lower.tail = F)
    survfit.res <- survfit(as.formula(sprintf("surv.obj ~ %s",vv$g.name[jj])))
    plot(survfit.res,col=c("red","blue"),
         xlab="Overall Survival (days)",ylab="Survival Function",main=vv$ann.name[jj])
    legend("topright",legend=gsub(sprintf("%s=",vv$g.name[jj]),"",names(survfit.res$strata)),col=c("red","blue"),lty=1)
    mtext(sprintf("n:%d (%s)\nn:%d (%s)\np value: %4.6f", survdiff.res$n[1], gsub(sprintf("^%s=",vv$g.name[jj]),"",names(survdiff.res$n[1]),perl=T),
                                     survdiff.res$n[2], gsub(sprintf("^%s=",vv$g.name[jj]),"",names(survdiff.res$n[2]),perl=T),
                                     survdiff.res.p),
          side = 1,line = -3,adj = 0.05)

}
dev.off()
detach(clinical.table.toAnalyse.noOSNA)

###  ajcc_pathologic_stage
if(sum(!is.na(clinical.table.toAnalyse$ajcc_stage_pathologic_stage)) > 3 )
{
    clinical.table.toAnalyse.noStageNA <- subset(clinical.table.toAnalyse,!is.na(ajcc_stage_pathologic_stage))
    pdf(sprintf("%s.stage.pdf",out.prefix),width=8,height=8)
    par(mar=c(5,5,4,2),cex.lab=1.5,cex.main=1.5)
    for(jj in 1:nrow(vv))
    {
        ### anova
        aov.out <- aov(as.formula(sprintf("%s ~ %s",vv$vname[jj],"ajcc_stage_pathologic_stage")),data=clinical.table.toAnalyse.noStageNA)
        aov.out.s <- summary(aov.out)
        aov.out.F <- unlist(aov.out.s[[1]]["ajcc_stage_pathologic_stage",c("F value","Pr(>F)")])
        aov.out.hsd <- TukeyHSD(aov.out)
        ### plot
        dat.plot <- split(clinical.table.toAnalyse.noStageNA[,vv$vname[jj]], clinical.table.toAnalyse.noStageNA$ajcc_stage_pathologic_stage)
        ylim <- pretty(range(clinical.table.toAnalyse.noStageNA[,vv$vname[jj]]),n=5)
        plot(0,0,type="n",xlim=c(0.5,length(dat.plot)+0.5), ylim=ylim[c(1,length(ylim))],  xaxt = 'n', xlab ="",ylab=vv$ann.name[jj],main=sprintf("%s",sample.id))
        for (i in 1:length(dat.plot)) { vioplot(na.omit(dat.plot[[i]]), at = i, add = T, col = "lightgray") }
        staxlab(1,at = 1:length(dat.plot),labels=names(dat.plot),srt=if(length(dat.plot) > 5) 45 else 0, cex=1.3,adj=0.5,top.line=2)
        mtext(sprintf("%s: %4.6f",names(aov.out.F)[1],aov.out.F[1]),side = 1,line = -2.1,adj = 0.05)
        mtext(sprintf("%s: %4.6f",names(aov.out.F)[2],aov.out.F[2]),side = 1,line = -1.1,adj = 0.05)
        
        print(aov.out.s)
        print(aov.out.hsd)
        psig=as.numeric(apply(aov.out.hsd$ajcc_stage_pathologic_stage[,2:3,drop=F],1,prod)>=0)+1
        op=par(mar=c(5,8,4,2))
        plot(aov.out.hsd,col=psig,yaxt="n")
        for(j in 1:length(psig)){
            axis(2,at=j,labels=rownames(aov.out.hsd$ajcc_stage_pathologic_stage)[length(psig)-j+1], las=1,cex.axis=.8,col.axis=psig[length(psig)-j+1])
        }
        par(op)
    }
    dev.off()
}
###
pdf(sprintf("%s.IS.pdf",out.prefix),width=10,height=8)
par(mar=c(5,5,4,2),cex.lab=1.5,cex.main=1.5)
plot(density(inf.score.out$score[,"TIS.raw"]),col="darkblue",lwd=2,xlab="TIS",main="")
plot(density(inf.score.out$score[,"OIS.raw"]),col="green",lwd=2,xlab="OIS",main="")
invisible(sapply(rownames(gsva.table),function(x){
                                plot(density(gsva.table[x,]),col="black",lwd=2,xlab=x,main="")
                            }))
dev.off()

out.df <- data.frame(gsetName=rownames(gsva.table),note=rownames(gsva.table))
out.df <- cbind(out.df,gsva.table)
write.table(out.df,file = sprintf("%s.gsva.txt",out.prefix),quote = F,row.names = F,col.names = T,sep = "\t")

if(do.clustering)
{
    grps <- unique(myDesign$leafCluster)
    #dat.tmp <- c()
    #for(g in grps) 
    #{
    #    dat.tmp.1 <- gsva.table[,rownames(subset(myDesign,leafCluster==g)),drop=F]
    #    if(ncol(dat.tmp.1)>2){ 
    #        dat.tmp <- cbind(dat.tmp,dat.tmp.1[,hclust(dist(t(dat.tmp.1)))$order,drop=F]) 
    #    }else{
    #        dat.tmp <- cbind(dat.tmp,dat.tmp.1)
    #    }
    #}
    #dat.plot <- dat.tmp
    #dat.plot <- gsva.table[!rownames(gsva.table) %in% c("Angiogenesis","PD1","CTLA4"),]
    dat.plot <- gsva.table
    annDF <- data.frame(cluster=paste0("Group",as.numeric(factor(myDesign[colnames(dat.plot),"leafCluster"]))))
    annDF <- cbind(annDF,inf.score.out$score[colnames(dat.plot),c("TIS","OIS")])
    annDF.col <- list(cluster=structure(auto.colSet(length(grps),name = "Dark2"),names=paste0("Group",seq_along(grps))))
    annDF.col <- c(annDF.col,inf.score.out$col)
    runHierarchicalClusteringAnalysis(dat.plot,mytitle = sample.id, 
                                      do.scale = do.scale,z.lo = -3,z.hi = 3,z.step = 0.5,
                                      pdf.width=26,pdf.height=12,do.clustering.col=do.colClustering,do.clustering.row=F,
                                      sprintf("%s.gsva",out.prefix),
                                      sampleType=myDesign[colnames(dat.plot),"sampleType"], 
                                      colSet=sampleTypeColor,ann.bar.height=0.8,
                                      ann.extra.df = annDF,
                                      ann.extra.df.col = annDF.col,
                                      clustering.distance = "euclidean", clustering.method = "ward.D2",
                                      k.row=1,ntop=NULL,complexHeatmap.use=TRUE,verbose=FALSE)

}
