#!/usr/bin/env Rscript

gene.file <- "/WPSnew/zhenglt/work/proj_zl/integrate0924/cross.patient.b20170523/randomForest/OUT.AUC65.CD8NoMAIT/CD8NoMAIT.pred.predMCls.rm00.markerGene.q01.txt"
universe.gene.file <- "/WPSnew/zhenglt/work/proj_zl/integrate0924/cross.patient.b20170523/exp/colonC.P8.expressedGene.list"
out.prefix <- "/WPSnew/zhenglt/work/proj_zl/integrate0924/cross.patient.b20170523/functional/GO.CD8_C08-MKI67"
cluster.str <- "CD8_C08-MKI67"
cmp.group.str <- "HSD.diff.CD8_C08-MKI67-CD8_C06-GZMK"

### gene set database
annfile.list <- c("GO.BP"="/DBS/DB_temp/zhangLab/MSigDB/msigdb_v5.2_files_to_download_locally/msigdb_v5.2_GMTs/c5.bp.v5.2.entrez.gmt",
                  "GO.MF"="/DBS/DB_temp/zhangLab/MSigDB/msigdb_v5.2_files_to_download_locally/msigdb_v5.2_GMTs/c5.mf.v5.2.entrez.gmt",
                  "immunoSig"="/DBS/DB_temp/zhangLab/MSigDB/msigdb_v5.2_files_to_download_locally/msigdb_v5.2_GMTs/c7.all.v5.2.entrez.gmt",
                  "canonicalPathways"="/DBS/DB_temp/zhangLab/MSigDB/msigdb_v5.2_files_to_download_locally/msigdb_v5.2_GMTs/c2.cp.v5.2.entrez.gmt",
                  "positional"="/DBS/DB_temp/zhangLab/MSigDB/msigdb_v5.2_files_to_download_locally/msigdb_v5.2_GMTs/c1.all.v5.2.entrez.gmt",
                  "motif"="/DBS/DB_temp/zhangLab/MSigDB/msigdb_v5.2_files_to_download_locally/msigdb_v5.2_GMTs/c3.all.v5.2.entrez.gmt",
                  "ZhangLabCurated"="/WPS1/zhenglt/work/proj_xy/integrated/database/ZhangLabCurated.entrez.gmt")
annDB <- list()
for(v in names(annfile.list))
{
    annDB[[v]] <- read.gmt(annfile.list[[v]])
}

gene.table <- read.table(gene.file,header = T,check.names = F,stringsAsFactors = F,sep = "\t")
universe.gene.table <- read.table(universe.gene.file,header = T,check.names = F,stringsAsFactors = F,sep = "\t")
gene.table$geneID <- as.character(gene.table$geneID)
universe.gene.table$geneID <- as.character(universe.gene.table$geneID)
g.GNAME <- universe.gene.table$geneSymbol
names(g.GNAME) <- universe.gene.table$geneID
f.na <- which(is.na(g.GNAME))
g.GNAME[f.na] <- names(g.GNAME)[f.na]

gene.table.cluster <- subset(gene.table, cluster==cluster.str)
gene.table.cluster[1:4,1:8]

library(org.Hs.eg.db)
library(clusterProfiler)

g.f <- gene.table.cluster[,cmp.group.str] >= 1
geneList.toTest <- structure(gene.table.cluster[g.f,cmp.group.str],names=gene.table.cluster[g.f,"geneID"])

doit <- function(atype="GO",out.prefix)
{
    if(atype=="GO"){
        e.result <- enrichGO(gene = names(geneList.toTest),
                        universe = universe.gene.table$geneID,
                        OrgDb = org.Hs.eg.db, ont  = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.01,
                        qvalueCutoff  = 0.05,
                        readable = TRUE)
        e.result.slim <- simplify(e.result, cutoff=0.7, by="p.adjust", select_fun=min)
    }else if(atype %in% names(annDB)){
        e.result <- enricher(names(geneList.toTest), TERM2GENE=annDB[[atype]])
        #head(e.result[,1:7])
    }
    
    pdf(sprintf("%s.map.%s.pdf",out.prefix,atype),width = 8,height = 8)
    enrichMap(e.result,vertex.label.font = 0.3, vertex.label.cex=0.3,n = 10)
    dev.off()
    
    pdf(sprintf("%s.dotplot.%s.pdf",out.prefix,atype),width = 8,height = 8)
    dotplot(e.result)
    dev.off()

    pdf(sprintf("%s.cnet.%s.pdf",out.prefix,atype),width = 8,height = 8)
    g <- cnetplot(e.result, categorySize="pvalue",showCategory = 5,vertex.label.cex=0.3,
             foldChange = geneList.toTest,layout=igraph::layout.fruchterman.reingold)
    dev.off()
    igraph::write_graph(g,file = sprintf("%s.cnet.%s.graphml",out.prefix,atype),format = "graphml")
    
    write.table(e.result[,1:7],file = sprintf("%s.%s.txt",out.prefix,atype),row.names = F,quote = F,sep = "\t")

    return(e.result)
}
genrich.out <- list()
genrich.out[["canonicalPathways"]] <- doit(atype="canonicalPathways",out.prefix)


