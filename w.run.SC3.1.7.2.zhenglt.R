#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-i", "--inputFile", type="character", required=TRUE, help="input file")
parser$add_argument("-d", "--designFile", type="character", required=TRUE, help="design matrix file")
parser$add_argument("-o", "--outprefix", type="character", required=TRUE, help="output prefix")
parser$add_argument("-c", "--cellTypeColorFile", default="/WPS1/zhenglt/work/TCR_chunhong/data/CellType.color", 
                    type="character", help="cellTypeColorFile [default %(default)s]")
parser$add_argument("-k", "--kmax", type="integer", default=10, help="kmax [default %(default)s]")
parser$add_argument("-w", "--kOptimal", type="character", help="specified optimal k [default %(default)s]")
parser$add_argument("-m", "--measure", type="character", default="exprs", help="measurement to use [default %(default)s]")
parser$add_argument("-a", "--vgene", type="character", default="sd", help="vgene method [default %(default)s]")
#parser$add_argument("-b", "--geneIDMapping", type="character", help="gene ID mapping file [default %(default)s]")
parser$add_argument("-n", "--nKeep", type="integer", default=1500, help="use top nKeep genes [default %(default)s]")
parser$add_argument("-y", "--ncores", type="integer", default=8, help="num of cors to use [default %(default)s]")
parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="whether verbose mode [default %(default)s]")
parser$add_argument("-t", "--myseed", type="character", help="seed [default %(default)s]")
args <- parser$parse_args()

print(args)

library("sscClust")
library("data.table")
library("tictoc")
library("SC3")
library("scater")
library("scran")
library("ggplot2")
library("doParallel")
library("foreach")
library("doRNG")
source("/WPSnew/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")

exp.file <- args$inputFile
d.file <- args$designFile
out.prefix <- args$outprefix
sampleType.color.file <- args$cellTypeColorFile
k.max <- args$kmax
k.best  <- as.integer(args$kOptimal)
args.ncore <- args$ncores
nKeep <- args$nKeep
args.measure <- args$measure
args.vgene <- args$vgene


if(file.exists(sprintf("%s.sce.RData",out.prefix))){
    lname <- load(sprintf("%s.sce.RData",out.prefix))
    write.table(colData(sce),sprintf("%s.colData.txt",out.prefix),row.names = F,sep = "\t",quote = F)
    q()
}

sampleTypeColor <- read.SampleTypeColor(sampleType.color.file)

out.dir <- dirname(out.prefix)
dir.create(out.dir,showWarnings = F,recursive = T)

d.table <- fread(d.file,key="sample")

lname <- load(exp.file)
sce <- sce.norm
sce.norm <- NULL

### add feature_symbol or geneID will be used automatically
rowData(sce)[,"feature_symbol"] <- rowData(sce)[,"display.name"]

gid.mapping <- rowData(sce)[,"display.name"]

if("centered_norm_exprs" %in% assayNames(sce)){
    assay(sce,"exprs") <- assay(sce,"centered_norm_exprs")
}else if("rmB_norm_exprs" %in% assayNames(sce)){
    assay(sce,"exprs") <- assay(sce,"rmB_norm_exprs")
}
print(head(d.table))
print(assay(sce,"exprs")[1:4,1:8])
sce <- sce[,d.table[,sample]]
sce
gc(T)

sampleTypeColor <- sampleTypeColor[unique(sce$sampleType)]
print(sampleTypeColor)

k.max <- min(k.max,dim(sce)[2]-1)
k.batch <- 2:k.max
################## SC3 pipeline ############

### SC3 use logcounts ##
assay(sce,"logcounts") <- NULL
sce <- ssc.build(sce)

tic("clustering ")
sce <- ssc.run(sce,assay.name=args.measure,method.vgene = args.vgene,ncore = args.ncore,
                            method.reduction = "pca",pca.npc = 20,
                            method.clust = "SC3",SC3.biology=T,
                            k.batch = k.batch,out.prefix = out.prefix)
toc()

write.table(colData(sce),sprintf("%s.colData.txt",out.prefix),row.names = F,sep = "\t",quote = F)
save(sce,file=sprintf("%s.sce.RData",out.prefix))

##### plot using pca info
pdf(sprintf("%s.pca.scree.pdf",out.prefix),width = 8,height = 6)
ssc.plot.pca(sce)
dev.off()

##save.image(sprintf("%s.test.RData",out.prefix))

verbose <- T
n.cores <- 9
n.de <- 1000
for(k in k.batch){
    ssc.plot.tsne(sce, columns = c(sprintf("sc3_%d_clusters",k)),
                  out.prefix = sprintf("%s.sc3_%d_clusters.pca.tsne",out.prefix,k),
                  reduced.name = "pca.tsne",p.ncol = 2,base_aspect_ratio = 1.4)
    ######
    gene.table <- my.clusterMarkerGene(assay(sce,args.measure), clust=as.numeric(colData(sce)[,sprintf("sc3_%d_clusters",k)]),
                                            ann.col=sampleTypeColor[sce$sampleType],
                                            clust.col=NULL,
                                            out.prefix=sprintf("%s.k%d.clust", out.prefix,k),
                                            n.cores=n.cores,
                                            sampleType = sce$sampleType,sampleTypeColSet = sampleTypeColor,gid.mapping = gid.mapping)

    ###### tsne using the de genes
    de.out <- gene.table$aov.res
    #de.out <- sscClust:::findDEGenesByAOV(xdata = assay(sce,"exprs"),
    #                                      xlabel = colData(sce)[,sprintf("sc3_%d_clusters",k)], 
    #                                      gid.mapping = structure(rowData(sce)[,"display.name"],
    #                                                              names=rownames(sce)))
    metadata(sce)$ssc[["de.res"]][[sprintf("L1C1.k%d",k)]] <- de.out
    metadata(sce)$ssc[["variable.gene"]][[sprintf("de.L1C1.k%d",k)]] <- de.out$aov.out.sig$geneID
    metadata(sce)$ssc[["variable.gene"]][[sprintf("de.L1C1.k%d.n%d",k,n.de)]] <- head(de.out$aov.out.sig$geneID,n=n.de)
}

save(sce,file=sprintf("%s.sce.RData",out.prefix))

######### m1 (large memory consumption, no sure why)
####tic("tsne")
####cl <- makeCluster(n_cores, outfile = "")
####registerDoParallel(cl, cores = n_cores)
####de.pca.tsne.res <- foreach::foreach(k = k.batch,.packages = c("sscClust")) %dopar% {}
####stopCluster(cl)
####toc()

######### m2
tic("tsne")
registerDoParallel(cores = n.cores)
do.pca.tsne.res <- plyr::llply(k.batch,function(k){
    ret <- list()
    de.out <- metadata(sce)$ssc[["de.res"]][[sprintf("L1C1.k%d",k)]]
    if(!is.null(de.out) && nrow(de.out$aov.out.sig)>30){
        rdimname.k <- sprintf("pca.de.L1C1.k%d.n%d",k,n.de)
        rdimname.tsne.k <- sprintf("pca.de.L1C1.k%d.n%d.tsne",k,n.de)
        sce <- ssc.reduceDim(sce,assay.name = args.measure,method = "pca",
                             method.vgene = sprintf("de.L1C1.k%d.n%d",k,n.de),pca.npc = 20,
                      tSNE.usePCA = F,autoTSNE = T,dim.name = rdimname.k ,seed = NULL)
        ret[[rdimname.k]] <- reducedDim(sce,rdimname.k)
        ret[[rdimname.tsne.k]] <- reducedDim(sce,rdimname.tsne.k)
        pdf(sprintf("%s.pca.de.L1C1.k%d.n%d.scree.pdf",out.prefix,k,n.de),width = 8,height = 6)
        ssc.plot.pca(sce)
        dev.off()
        ssc.plot.tsne(sce, assay.name=args.measure, columns = c(sprintf("sc3_%d_clusters",k)),
                      out.prefix = sprintf("%s.pca.de.L1C1.k%d.n%d.tsne.sc3_%d_clusters",out.prefix,k,k,n.de),
                      reduced.name = rdimname.tsne.k,p.ncol = 2,base_aspect_ratio = 1.4)
        ssc.plot.tsne(sce, assay.name=args.measure, columns = NULL,plotDensity = T,
                      out.prefix = sprintf("%s.pca.de.L1C1.k%d.n%d.tsne.density",out.prefix,k,n.de),
                      reduced.name = rdimname.tsne.k,p.ncol = 2,base_aspect_ratio = 1.4)
        #ssc.plot.tsne(sce, assay.name=args.measure,gene = c("CCR7","GZMK","CCL5","HAVCR2","CXCL13","CX3CR1","FOXP3","S1PR1","RGS1"),
        #              reduced.name = rdimname.tsne.k,
        #              out.prefix = sprintf("%s.pca.de.L1C1.k%d.n%d.tsne.G00",out.prefix,k,n.de),width = 9,height = 8)
        #ssc.plot.tsne(sce, assay.name=args.measure,gene = c("CXCR6","ACTN1","MAL","SELL","ADAM19","CD27","EOMES","TCF7","CAPG"),
        #              reduced.name = rdimname.tsne.k,
        #              out.prefix = sprintf("%s.pca.de.L1C1.k%d.n%d.tsne.G01",out.prefix,k,n.de),width = 9,height = 8)
    }else{
        cat(sprintf("less than 30 de genes for k=%d\n",k))
    }
    return(ret)
},.progress = "none",.parallel=T)
toc()

for(i in seq_along(do.pca.tsne.res)){
    for(.dname in names(do.pca.tsne.res[[i]])){
        reducedDim(sce,.dname) <- do.pca.tsne.res[[i]][[.dname]]
    }
}

##########################################

save(sce,file=sprintf("%s.sce.RData",out.prefix))

q()

####################
###
####### DPClust refineGene=T,nIter=2
###tic("clustering ")
####sce.dp.iter2.refineG <- ssc.run(sce.dp.iter2.refineG,method.reduction = "pca",pca.npc = 15,method.clust = "dpclust",
###sce.dp.iter2.refineG <- ssc.run(sce,method.reduction = "pca",pca.npc = 15,method.clust = "dpclust",
###                               refineGene = T,nIter = 2,do.DE = F,reuse = F,seed = 9998,
###                               parfile = "/WPSnew/zhenglt/work/panC/data/merge.phase01/pca.dpclust.par.CD4.r",
###                               out.prefix = sprintf("%s.pca.tsne.dpclust.iter2.refineG",out.prefix))
###toc()
###
###pdf(sprintf("%s.pca.tsne.dpclust.iter1.PCA.pdf",out.prefix),width = 8,height = 6)
###ssc.plot.pca(sce.dp.iter2.refineG)
###dev.off()
###
###ssc.plot.tsne(sce.dp.iter2.refineG,gene = c("CCR7","GZMK","CCL5","HAVCR2","CXCL13","CX3CR1","FOXP3","S1PR1","RGS1"),
###              reduced.name = "pca.tsne",out.prefix = sprintf("%s.pca.tsne.dpclust.iter1.G00",out.prefix),width = 9,height = 8)
###ssc.plot.tsne(sce.dp.iter2.refineG,gene = c("IL2RA","IL17A","RORC","PDCD1","IKZF2","CXCR5","LAYN","GPR183","CXCR6"),
###              reduced.name = "pca.tsne",out.prefix = sprintf("%s.pca.tsne.dpclust.iter1.G01",out.prefix),width = 9,height = 8)
###ssc.plot.tsne(sce.dp.iter2.refineG,gene = c("BCL6","ACTN1","MAL","SELL","ADAM19","CD27","SELL","TCF7","PLAC8"),
###              reduced.name = "pca.tsne",out.prefix = sprintf("%s.pca.tsne.dpclust.iter1.G02",out.prefix),width = 9,height = 8)
###
###ssc.plot.tsne(sce.dp.iter2.refineG,gene = c("CCR7","GZMK","CCL5","HAVCR2","CXCL13","CX3CR1","FOXP3","S1PR1","RGS1"),
###              reduced.name = "de.pca.tsne",out.prefix = sprintf("%s.pca.tsne.dpclust.iter1.de.G00",out.prefix),width = 9,height = 8)
###ssc.plot.tsne(sce.dp.iter2.refineG,gene = c("IL2RA","IL17A","RORC","PDCD1","IKZF2","CXCR5","LAYN","GPR183","CXCR6"),
###              reduced.name = "de.pca.tsne",out.prefix = sprintf("%s.pca.tsne.dpclust.iter1.de.G01",out.prefix),width = 9,height = 8)
###ssc.plot.tsne(sce.dp.iter2.refineG,gene = c("BCL6","ACTN1","MAL","SELL","ADAM19","CD27","SELL","TCF7","PLAC8"),
###              reduced.name = "de.pca.tsne",out.prefix = sprintf("%s.pca.tsne.dpclust.iter1.de.G02",out.prefix),width = 9,height = 8)
###
###ssc.plot.tsne(sce.dp.iter2.refineG,columns = c("pca.tsne.dpclust.kauto"),
###              reduced.name = "de.pca.tsne",out.prefix = sprintf("%s.pca.tsne.dpclust.iter1.de.kauto",out.prefix),base_aspect_ratio = 1.4)
###
###
###
###### rerun using the de genes
###de.out <- sscClust:::findDEGenesByAOV(xdata = assay(sce.dp.iter2.refineG,"exprs"),
###                                      xlabel = colData(sce.dp.iter2.refineG)[,"pca.tsne.dpclust.kauto"], 
###                                      gid.mapping = structure(rowData(sce.dp.iter2.refineG)[,"display.name"],
###                                                              names=rownames(sce.dp.iter2.refineG)))
###metadata(sce.dp.iter2.refineG)$ssc[["de.res"]][["L1C1.final"]] <- de.out
###metadata(sce.dp.iter2.refineG)$ssc[["variable.gene"]][["de.L1C1.final.n1000"]] <- head(de.out$aov.out.sig$geneID,n=1500)
###
#####cols <- c("patient","sample","sampleType","stype","invariantTCR","total_counts","total_features","size_factor","cyclePhase","pca.tsne.dpclust.kauto","majorCluster")
#####write.table(colData(sce.dp.iter2.refineG)[,cols],sprintf("%s.cluster.txt",out.prefix),row.names = F,sep = "\t",quote = F)
###
######## second run
###tic("clustering ")
###sce.dp.run2.refineG <- ssc.run(sce.dp.iter2.refineG,method.vgene = "de.L1C1.final.n1000",
###                               method.reduction = "pca",pca.npc = 15,method.clust = "dpclust",
###                               refineGene = F,nIter = 1,do.DE = F,reuse = F,seed = 9998,
###                               parfile="/WPSnew/zhenglt/work/panC/data/merge.phase01/pca.dpclust.par.CD4.run2.r",
###                               out.prefix = sprintf("%s.pca.tsne.dpclust.run2.refineG",out.prefix))
###toc()
###
###
###ssc.plot.tsne(sce.dp.run2.refineG,gene = c("CCR7","GZMK","CCL5","HAVCR2","CXCL13","CX3CR1","FOXP3","S1PR1","RGS1"),
###              reduced.name = "pca.tsne",out.prefix = sprintf("%s.pca.tsne.pca.dpclust.run2.G00",out.prefix),width = 9,height = 8)
###ssc.plot.tsne(sce.dp.run2.refineG,gene = c("IL2RA","IL17A","RORC","PDCD1","IKZF2","CXCR5","LAYN","GPR183","CXCR6"),
###              reduced.name = "pca.tsne",out.prefix = sprintf("%s.pca.tsne.pca.dpclust.run2.G01",out.prefix),width = 9,height = 8)
###ssc.plot.tsne(sce.dp.run2.refineG,gene = c("BCL6","ACTN1","MAL","SELL","ADAM19","CD27","EOMES","TCF7","PLAC8"),
###              reduced.name = "pca.tsne",out.prefix = sprintf("%s.pca.tsne.pca.dpclust.run2.G02",out.prefix),width = 9,height = 8)
###
###
###save(sce.dp.iter2.refineG,sce.dp.run2.refineG,file=sprintf("%s.sce.RData",out.prefix))
###
####### Th17
###f <- sce.dp.run2.refineG$pca.tsne.dpclust.kauto %in% c("C6","C7","C11")
###tic("clustering ")
###sce.dp.G3 <- ssc.run(sce.dp.run2.refineG[,f],method.vgene = "sd",
###                               method.reduction = "pca",pca.npc = 15,method.clust = "dpclust",
###                               refineGene = T,nIter = 1,do.DE = F,reuse = F,seed = 9998,
###                               parfile="/WPSnew/zhenglt/work/panC/data/merge.phase01/pca.dpclust.par.CD4.run2.G3.r",
###                               out.prefix = sprintf("%s.pca.tsne.dpclust.G3",out.prefix))
###toc()
###ssc.plot.tsne(sce.dp.G3,gene = c("CCR7","GZMK","CCL5","HAVCR2","CXCL13","CX3CR1","FOXP3","S1PR1","RGS1"),
###              reduced.name = "de.pca.tsne",out.prefix = sprintf("%s.pca.tsne.pca.dpclust.G3.G00",out.prefix),width = 9,height = 8)
###ssc.plot.tsne(sce.dp.G3,gene = c("IL23R","IL17A","RORC","PDCD1","IKZF2","CXCR5","LAYN","GPR183","CXCR6"),
###              reduced.name = "de.pca.tsne",out.prefix = sprintf("%s.pca.tsne.pca.dpclust.G3.G01",out.prefix),width = 9,height = 8)
###ssc.plot.tsne(sce.dp.G3,gene = c("BCL6","ACTN1","MAL","SELL","ADAM19","CD27","EOMES","TCF7","PLAC8"),
###              reduced.name = "de.pca.tsne",out.prefix = sprintf("%s.pca.tsne.pca.dpclust.G3.G02",out.prefix),width = 9,height = 8)
###
###
###save(sce.dp.G3,file=sprintf("%s.sce.G3.RData",out.prefix))
###
###
####ssc.plot.tsne(sce.dp.run2.refineG,columns = c("majorCluster"),
####              reduced.name = "pca.tsne",out.prefix = sprintf("%s.pca.tsne.pca.dpclust.run2.majorCluster.old",out.prefix),base_aspect_ratio = 1.4)
###
###sce.dp.run2.refineG$pca.tsne.dpclust.kfinal  <- sce.dp.run2.refineG$pca.tsne.dpclust.kauto
###sce.dp.run2.refineG$pca.tsne.dpclust.kfinal[sce.dp.run2.refineG$pca.tsne.dpclust.kfinal=="C4"]  <- "CD8_C1-CCR7"
###sce.dp.run2.refineG$pca.tsne.dpclust.kfinal[sce.dp.run2.refineG$pca.tsne.dpclust.kfinal=="C13"]  <- "CD8_C2-GPR183"
###sce.dp.run2.refineG$pca.tsne.dpclust.kfinal[sce.dp.run2.refineG$pca.tsne.dpclust.kfinal=="C1"]  <- "CD8_C3-CX3CR1"
###sce.dp.run2.refineG$pca.tsne.dpclust.kfinal[sce.dp.run2.refineG$pca.tsne.dpclust.kfinal=="C10"]  <- "CD8_C3-CX3CR1"
###sce.dp.run2.refineG$pca.tsne.dpclust.kfinal[sce.dp.run2.refineG$pca.tsne.dpclust.kfinal=="C12"]  <- "CD8_C3-CX3CR1"
###sce.dp.run2.refineG$pca.tsne.dpclust.kfinal[sce.dp.run2.refineG$pca.tsne.dpclust.kfinal=="C3"]  <- "CD8_C4-GZMK"
###sce.dp.run2.refineG$pca.tsne.dpclust.kfinal[sce.dp.run2.refineG$pca.tsne.dpclust.kfinal=="C7"]  <- "CD8_C4-GZMK"
###sce.dp.run2.refineG$pca.tsne.dpclust.kfinal[sce.dp.run2.refineG$pca.tsne.dpclust.kfinal=="C6"]  <- "CD8_C5-ZNF683"
###sce.dp.run2.refineG$pca.tsne.dpclust.kfinal[sce.dp.run2.refineG$pca.tsne.dpclust.kfinal=="C11"]  <- "CD8_C5-ZNF683"
###sce.dp.run2.refineG$pca.tsne.dpclust.kfinal[sce.dp.run2.refineG$pca.tsne.dpclust.kfinal=="C5"]  <- "CD8_C6-HAVCR2"
###sce.dp.run2.refineG$pca.tsne.dpclust.kfinal[sce.dp.run2.refineG$pca.tsne.dpclust.kfinal=="C9"]  <- "CD8_C6-HAVCR2"
###sce.dp.run2.refineG$pca.tsne.dpclust.kfinal[sce.dp.run2.refineG$pca.tsne.dpclust.kfinal=="C2"]  <- "CD8_C7-CD6"
###sce.dp.run2.refineG$pca.tsne.dpclust.kfinal[sce.dp.run2.refineG$pca.tsne.dpclust.kfinal=="C8"]  <- "CD8_C8-CD160"
###
###sce.dp.run2.refineG$majorCluster  <- sce.dp.run2.refineG$pca.tsne.dpclust.kfinal
######sce.dp.run2.refineG$majorCluster[sce.dp.run2.refineG$majorCluster=="CD8_C4-GZMK.P"] <- "CD8_C4-GZMK"
######sce.dp.run2.refineG$majorCluster[sce.dp.run2.refineG$majorCluster=="CD8_C4-GZMK.NT"] <- "CD8_C4-GZMK"
###
###cols <- c("patient","sample","sampleType","stype","invariantTCR","total_counts","total_features","size_factor","cyclePhase","pca.tsne.dpclust.kfinal","majorCluster")
###write.table(colData(sce.dp.run2.refineG)[,cols],sprintf("%s.cluster.final.txt",out.prefix),row.names = F,sep = "\t",quote = F)
###ssc.plot.tsne(sce.dp.run2.refineG,columns = c("pca.tsne.dpclust.kfinal"),
###              reduced.name = "pca.tsne",out.prefix = sprintf("%s.pca.tsne.dpclust.run2.kfinal",out.prefix),base_aspect_ratio = 1.4)
###ssc.plot.tsne(sce.dp.run2.refineG,columns = c("majorCluster"),
###              reduced.name = "pca.tsne",out.prefix = sprintf("%s.pca.tsne.dpclust.run2.majorCluster",out.prefix),base_aspect_ratio = 1.4)
###
###print(table(sce.dp.iter2.refineG$sampleType,sce.dp.iter2.refineG$pca.tsne.dpclust.kfinal))
###print(table(sce.dp.run2.refineG$sampleType,sce.dp.run2.refineG$pca.tsne.dpclust.kfinal))
###
###### final DE genes
###de.out.run2 <- sscClust:::findDEGenesByAOV(xdata = assay(sce.dp.run2.refineG,"exprs"),
###                                      xlabel = colData(sce.dp.run2.refineG)[,"majorCluster"], 
###                                      gid.mapping = structure(rowData(sce.dp.run2.refineG)[,"display.name"],
###                                                              names=rownames(sce.dp.run2.refineG)))
###
###metadata(sce.dp.run2.refineG)$ssc[["de.res"]][["final.run2"]] <- de.out.run2
####metadata(sce.dp.run2.refineG)$ssc[["variable.gene"]][["de.final.n1000"]] <- head(de.out.run2$aov.out.sig$geneID,n=1000)
###
###
###
###save(sce.dp.iter2.refineG,sce.dp.run2.refineG,file=sprintf("%s.sce.RData",out.prefix))
###
###
###

#ssc.plot.tsne(sce, columns = c("sc3_6_clusters"),out.prefix = sprintf("%s.sc3_6_clusters.spearman_laplacian.D1D2",out.prefix),
#              reduced.name = "sc3.spearman_laplacian",p.ncol = 2,base_aspect_ratio = 1.4,reduced.dim = c(1,2))
#ssc.plot.tsne(sce, columns = c("sc3_6_clusters"),out.prefix = sprintf("%s.sc3_6_clusters.spearman_laplacian.D3D4",out.prefix),
#              reduced.name = "sc3.spearman_laplacian",p.ncol = 2,base_aspect_ratio = 1.4,reduced.dim = c(3,4))
#ssc.plot.tsne(sce, columns = c("sc3_6_clusters"),out.prefix = sprintf("%s.sc3_6_clusters.spearman_laplacian.D5D6",out.prefix),
#              reduced.name = "sc3.spearman_laplacian",p.ncol = 2,base_aspect_ratio = 1.4,reduced.dim = c(5,6))
#ssc.plot.tsne(sce, columns = c("sc3_6_clusters"),out.prefix = sprintf("%s.sc3_6_clusters.spearman_laplacian.tsne",out.prefix),
#              reduced.name = "sc3.spearman_laplacian.tsne",p.ncol = 2,base_aspect_ratio = 1.4)


