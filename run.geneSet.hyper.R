#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-i", "--geneAOVFile", type="character", required=TRUE, help="")
parser$add_argument("-a", "--bgFile", type="character", required=TRUE, help="")
parser$add_argument("-o", "--outputPrefix", type="character", required=TRUE, help="output prefix")
parser$add_argument("-c", "--cluster", type="character", required=TRUE, help="which cluster to test")
parser$add_argument("-r", "--reference", type="character", required=TRUE, help="reference cluster")
parser$add_argument("-t", "--invert", type="integer", default=0, help="invert the diff [default %(default)s]")
parser$add_argument("-n", "--ncores", type="integer", default=4, help="cpu cores [default %(default)s]")
parser$add_argument("-s", "--species", type="character", default="human", help="species (human, mouse) [default %(default)s]")
parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="verbose mode [default %(default)s]")

args <- parser$parse_args()
print(args)

if(file.exists("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")){
	source("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")
}else{
	source("/lustre1/zeminz_pkuhpc/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")
}
library("magrittr")
library("ggplot2")
library("RColorBrewer")

gene.file <- args$geneAOVFile
universe.gene.file <- args$bgFile
out.prefix <- args$outputPrefix
cluster.str <- args$cluster
reference.str <- args$reference
args.verbose <- args$verbose
args.invert <- args$invert
args.ncores <- args$ncores
args.species <- args$species

#gene.file <- "/WPSnew/zhenglt/work/proj_zl/integrate0924/cross.patient.b20170523/randomForest/OUT.AUC65.CD8NoMAIT/CD8NoMAIT.pred.predMCls.rm00.markerGene.q01.txt"
#universe.gene.file <- "/WPSnew/zhenglt/work/proj_zl/integrate0924/cross.patient.b20170523/exp/colonC.P8.expressedGene.list"
#out.prefix <- "/WPSnew/zhenglt/work/proj_zl/integrate0924/cross.patient.b20170523/functional/GO.CD8_C08-MKI67"
##cluster.str <- "CD8_C08-MKI67"
#cluster.str <- "CD8_C07-HAVCR2"
#reference.str <- "CD8_C06-GZMK"
#args.invert <- 1
#args.ncores <- 4

#cmp.group.str <- "HSD.diff.CD8_C08-MKI67-CD8_C06-GZMK"

### gene set database
if(args.species=="human"){
    if(file.exists("/DBS/DB_temp/zhangLab/MSigDB")){
        annfile.list <- c("GO.BP"="/DBS/DB_temp/zhangLab/MSigDB/msigdb_v5.2_files_to_download_locally/msigdb_v5.2_GMTs/c5.bp.v5.2.entrez.gmt",
                          "GO.MF"="/DBS/DB_temp/zhangLab/MSigDB/msigdb_v5.2_files_to_download_locally/msigdb_v5.2_GMTs/c5.mf.v5.2.entrez.gmt",
                          "immunoSig"="/DBS/DB_temp/zhangLab/MSigDB/msigdb_v5.2_files_to_download_locally/msigdb_v5.2_GMTs/c7.all.v5.2.entrez.gmt",
                          "canonicalPathways"="/DBS/DB_temp/zhangLab/MSigDB/msigdb_v5.2_files_to_download_locally/msigdb_v5.2_GMTs/c2.cp.v5.2.entrez.gmt",
                          "positional"="/DBS/DB_temp/zhangLab/MSigDB/msigdb_v5.2_files_to_download_locally/msigdb_v5.2_GMTs/c1.all.v5.2.entrez.gmt",
                          "motif"="/DBS/DB_temp/zhangLab/MSigDB/msigdb_v5.2_files_to_download_locally/msigdb_v5.2_GMTs/c3.all.v5.2.entrez.gmt",
                          "TRRUST"="/WPS1/zhenglt/work/proj_xy/integrated/database/TRRUST.entrez.gmt",
                          "ZhangLabCurated"="/WPS1/zhenglt/work/proj_xy/integrated/database/ZhangLabCurated.entrez.gmt")
    }else{
        annfile.list <- c("GO.BP"="/lustre1/zeminz_pkuhpc/00.database/MSigDB/msigdb_v6.0_files_to_download_locally/msigdb_v6.0_GMTs/c5.bp.v6.0.entrez.gmt",
                          "GO.MF"="/lustre1/zeminz_pkuhpc/00.database/MSigDB/msigdb_v6.0_files_to_download_locally/msigdb_v6.0_GMTs/c5.mf.v6.0.entrez.gmt",
                          "immunoSig"="/lustre1/zeminz_pkuhpc/00.database/MSigDB/msigdb_v6.0_files_to_download_locally/msigdb_v6.0_GMTs/c7.all.v6.0.entrez.gmt",
                          "canonicalPathways"="/lustre1/zeminz_pkuhpc/00.database/MSigDB/msigdb_v6.0_files_to_download_locally/msigdb_v6.0_GMTs/c2.cp.v6.0.entrez.gmt",
                          "motif"="/lustre1/zeminz_pkuhpc/00.database/MSigDB/msigdb_v6.0_files_to_download_locally/msigdb_v6.0_GMTs/c3.all.v6.0.entrez.gmt")
    }
}else{
    annfile.list <- c("GO.BP"="/DBS/DB_temp/zhangLab/MSigDB/gskb/MousePath_GO_BP_format.ensembl.gmt",
                      "GO.MF"="/DBS/DB_temp/zhangLab/MSigDB/gskb/MousePath_GO_MF_format.ensembl.gmt",
                      "TF.target"="/DBS/DB_temp/zhangLab/MSigDB/gskb/MousePath_TF_format.ensembl.gmt",
                      "Metabolic"="/DBS/DB_temp/zhangLab/MSigDB/gskb/MousePath_Metabolic_format.ensembl.gmt",
                      "Pathway"="/DBS/DB_temp/zhangLab/MSigDB/gskb/MousePath_Pathway_format.ensembl.gmt",
                      "miRNA.target"="/DBS/DB_temp/zhangLab/MSigDB/gskb/MousePath_miRNA_format.ensembl.gmt")
}

annDB <- list()
for(v in names(annfile.list))
{
    annDB[[v]] <- readGMT(annfile.list[[v]])
}

gene.table <- read.table(gene.file,header = T,check.names = F,stringsAsFactors = F,sep = "\t")
universe.gene.table <- read.table(universe.gene.file,header = T,check.names = F,stringsAsFactors = F,sep = "\t")
gene.table$geneID <- as.character(gene.table$geneID)
universe.gene.table$geneID <- as.character(universe.gene.table$geneID)
g.GNAME <- universe.gene.table$geneSymbol
names(g.GNAME) <- universe.gene.table$geneID
f.na <- which(is.na(g.GNAME))
g.GNAME[f.na] <- names(g.GNAME)[f.na]

if(cluster.str==reference.str){
    gene.table.cluster <- subset(gene.table,cluster==cluster.str)
    geneList.toTest <- structure(gene.table.cluster[,"F"],names=gene.table.cluster$geneID)
    loginfo("geneList.toTest got!")
    print(str(geneList.toTest))
    print(head(geneList.toTest))
}else{
    cmp.group.str <- (grepl(cluster.str,colnames(gene.table),perl=T) &
                      grepl(reference.str,colnames(gene.table),perl=T)) %>%
                        colnames(gene.table)[.]
    cmp.group.str.diff <- grepl("diff",cmp.group.str,perl=T) %>% cmp.group.str[.]
    cmp.group.str.padj <- grepl("padj",cmp.group.str,perl=T) %>% cmp.group.str[.]
    print(cmp.group.str.diff)
    print(cmp.group.str.padj)
    cmp.group.str.diff.idx <- which(grepl(cmp.group.str.diff,colnames(gene.table),perl=T))
    cmp.group.str.padj.idx <- which(grepl(cmp.group.str.padj,colnames(gene.table),perl=T))

    #gene.table.cluster <- subset(gene.table, cluster==cluster.str)
    gene.table.cluster <- gene.table
    if(args.invert==0)
    {
        g.f <- abs(gene.table.cluster[,cmp.group.str.diff]) >= 1 &
                    gene.table.cluster[,cmp.group.str.padj]<0.01
    }else if(args.invert==1){
        g.f <- gene.table.cluster[,cmp.group.str.diff] > 1 &
                    gene.table.cluster[,cmp.group.str.padj]<0.01
    }else if(args.invert==-1){
        g.f <- gene.table.cluster[,cmp.group.str.diff] < 1 &
                    gene.table.cluster[,cmp.group.str.padj]<0.01
    }
    gene.table.cluster <- gene.table.cluster[g.f,]
    rownames(gene.table.cluster) <- gene.table.cluster$geneID
    gene.table.cluster <- gene.table.cluster[order(gene.table.cluster[,cmp.group.str.padj],
                                               -gene.table.cluster[,cmp.group.str.diff]),]
    head(gene.table.cluster[,c(1:8,cmp.group.str.diff.idx,cmp.group.str.padj.idx)],n=4) %>% print
    geneList.toTest <- structure(gene.table.cluster[,cmp.group.str.diff],names=gene.table.cluster[,"geneID"])
}
           

#library(org.Hs.eg.db)
#library(clusterProfiler)

doit <- function(atype="GO",out.prefix,sigScore.max=10,padj.threshold=0.05)
{
    gset <- annDB[[atype]][["gSet"]]
    #### Hypergeometric test
    loginfo(sprintf("annotate using gset %s",atype))
    #print(str(gset))
    #print(str(names(geneList.toTest)))
    #print(str(universe.gene.table$geneID))
    #print(head(names(geneList.toTest)))
    #print(head(universe.gene.table$geneID))
    enrich.gset <- get.geneSet.hyper(gset,
                                     names(geneList.toTest),
                                     universe.gene.table$geneID,
                                     min.size=5,n.cores=args.ncores,verbose = T,IDMapping = g.GNAME)
    enrich.gset.sig <- subset(enrich.gset,p.adj<padj.threshold)
    write.table(enrich.gset,sprintf("%s.%s.hyper.txt", out.prefix,atype),
                quote = F,sep = "\t",row.names = F)
    write.table(enrich.gset.sig,sprintf("%s.%s.hyper.sig.txt", out.prefix,atype),
                quote = F,sep = "\t",row.names = F)

    #### plot
    dat.plot <- head(enrich.gset.sig[,1:9,drop=F],n=20)
    f.longname <- nchar(dat.plot$geneSet)>32
    dat.plot$geneSet <- strtrim(dat.plot$geneSet,32)
    dat.plot$geneSet[f.longname] <- sprintf("%s...",dat.plot$geneSet[f.longname])
    dat.plot$geneSet <- sprintf("%s(%d)",dat.plot$geneSet,seq_len(nrow(dat.plot)))
    idx <- order(dat.plot$ratio, decreasing = F)
    dat.plot$geneSet <- factor(dat.plot$geneSet, levels=dat.plot$geneSet[idx])
    dat.plot$sigScore <- -log10(dat.plot$p.adj)
    dat.plot$sigScore[dat.plot$sigScore>sigScore.max] <- sigScore.max
    p <- ggplot(dat.plot,aes(x=ratio, y=geneSet)) +
        geom_point(aes(size=observed,colour=sigScore)) +
        ylab("") + xlab("EnrichmentScore") +
        scale_color_gradientn(colours=brewer.pal(9,"YlOrRd"),limits=c(-log10(padj.threshold),sigScore.max)) +
        theme_bw() +
        theme(axis.text.x=element_text(size=12),
              axis.title.x=element_text(size=14,face="bold"))
        ##theme(axis.text.y = element_text(angle = 60, hjust = 1))
    ggsave(filename = sprintf("%s.%s.hyper.dotplot.pdf",out.prefix,atype),width = 7,height = 6)
    p <- ggbarplot(dat.plot, x="geneSet", y="sigScore",rotate=T,color="steelblue",fill="steelblue") +
            theme(axis.text.x=element_text(size=12), axis.title.x=element_text(size=14,face="bold")) +
            xlab("") + ylab("-log10(FDR)") + 
    ggsave(filename = sprintf("%s.%s.hyper.barplot.pdf",out.prefix,atype),width = 7,height = 6)

    #### genes not well annotated
    gene.list.annotated <- sapply(enrich.gset.sig$geneID,
                                  function(x){ unlist(strsplit(x,",")) }) %>%
                            unname %>% unlist %>% unique
    gene.list.notAnnotated <- setdiff(names(geneList.toTest),gene.list.annotated)
    gene.table.cluster.notAnnotated <- gene.table.cluster[gene.list.notAnnotated,]
    write.table(gene.table.cluster.notAnnotated,sprintf("%s.%s.hyper.notSig.list.txt", out.prefix,atype),
                quote = F,sep = "\t",row.names = F)
    #pdf(sprintf("%s.%s.hyper.cnet.pdf",out.prefix,atype),width = 8,height = 8)
    #g <- cnetplot(e.result, categorySize="pvalue",showCategory = 5,vertex.label.cex=0.3,
    #         foldChange = geneList.toTest,layout=igraph::layout.fruchterman.reingold)
    #dev.off()
    #igraph::write_graph(g,file = sprintf("%s.cnet.%s.graphml",out.prefix,atype),format = "graphml")
    
    return(list("enrich"=enrich.gset.sig,"ann"=gene.list.annotated,"notAnn"=gene.table.cluster.notAnnotated))
}
genrich.out <- list()
all.gene.annotated <- c()
for(nn in names(annDB)){
    genrich.out[[nn]] <- doit(atype=nn,out.prefix)
    if(args.species=="human"){
        if(nn %in% c("GO.BP","GO.MF","canonicalPathways","motif")){
            all.gene.annotated <- append(all.gene.annotated,genrich.out[[nn]][["ann"]])
        }
    }else{
        if(nn %in% c("GO.BP","GO.MF","TF.target","Metabolic","Pathway","miRNA.target")){
            all.gene.annotated <- append(all.gene.annotated,genrich.out[[nn]][["ann"]])
        }
    }
}
all.gene.annotated <- unique(all.gene.annotated)

write.table(gene.table.cluster[setdiff(rownames(gene.table.cluster),all.gene.annotated),],sprintf("%s.any.hyper.notSig.list.txt", out.prefix),
            quote = F,sep = "\t",row.names = F)

save(genrich.out,file=sprintf("%s.RData",out.prefix))

