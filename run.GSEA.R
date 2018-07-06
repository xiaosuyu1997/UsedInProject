#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-i", "--inFile", type="character", required=TRUE, help="infile")
parser$add_argument("-o", "--outputPrefix", type="character", required=TRUE, help="output prefix")
parser$add_argument("-a", "--cmp", type="character", help="comparison, comma seperated")
parser$add_argument("-s", "--scoring_scheme", type="character", default="classic",
                    help="scoring_scheme of GSEA: classic,weighted,weighted_p2,weighted_p1.5 [default %(default)s]")
parser$add_argument("-x", "--max", type="integer",default=500, help="maximum gene set size [default %(default)s]")
parser$add_argument("-d", "--db", type="character",default="/DBS/DB_temp/zhangLab/MSigDB/msigdb_v6.0_files_to_download_locally/msigdb_v6.0_GMTs/c2.cp.v6.0.symbols.gmt", 
                    help="db file [default %(default)s]")
parser$add_argument("-r", "--reverse", action="store_true", default=FALSE, help="reverse the sign [default %(default)s]")
parser$add_argument("-b", "--useSD", action="store_true", default=FALSE, help="dived foldChange by sd [default %(default)s]")
parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="verbose mode [default %(default)s]")

args <- parser$parse_args()
print(args)

in.file <- args$inFile
out.prefix <- args$outputPrefix
args.verbose <- args$verbose
args.cmp <- args$cmp
args.r <- args$reverse
db.file <- args$db
args.max <- args$max
args.score <- args$scoring_scheme
args.useSD <- args$useSD

##in.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/cross.patient/geneClustering/OUT.kmeans/CD8.RGS1.C2.inTumor.byPatientF/CD8.RGS1.C2.inTumor.gene.clusters.all.RData"
##out.prefix <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/cross.patient/geneClustering/OUT.kmeans/CD8.RGS1.C2.inTumor.byPatientF/CD8.RGS1.C2.inTumor.Cluster5-GZMK-CD8_Cluster4-LAYN.FC"
#in.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/cross.patient/geneClustering/OUT.rmIsCycle.kmeans/CD8.RGS1.C2.inTumor.byPatientF/CD8.RGS1.C2.inTumor.gene.clusters.all.RData"
#out.prefix <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/cross.patient/geneClustering/OUT.rmIsCycle.kmeans/CD8.RGS1.C2.inTumor.byPatientF/CD8.RGS1.C2.inTumor.Cluster5-GZMK-CD8_Cluster4-LAYN.FC"
#in.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/cross.patient/geneClustering/OUT.rmIsCycle.kmeans/CD8.byPatientF/CD8.gene.clusters.all.RData"
#out.prefix <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/hetero/cross.patient/geneClustering/OUT.rmIsCycle.kmeans/CD8.byPatientF/CD8.CX3CR1.LEF1"
#args.cmp <- "CD8_Cluster1-LEF1,CD8_Cluster2-CX3CR1"
#args.r <- F
#db.file <- "/WPS1/zhenglt/work/proj_xy/integrated/database/ZhangLabCurated.symbol.h8.gmt"
#args.max <- 500

##source("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")

suppressPackageStartupMessages(library("R.utils"))
if(grepl("\\.RData$",in.file,perl=T)){
    lenv <- loadToEnv(in.file)
    print(names(lenv))
    if("aov.res" %in% names(lenv)){
        aov.res <- lenv[["aov.res"]]
    }else{
        if("gene.table" %in% names(lenv[["clusterMarker.res"]])){
            aov.res <- list("aov.out"=lenv[["clusterMarker.res"]][["gene.table"]])
        }else{
            aov.res <- lenv[["clusterMarker.res"]][["aov.res"]]
        }
    }
}else{
    .i.table <- read.delim(in.file,header = T,check.names = F,stringsAsFactors = F,sep = "\t")
    rownames(.i.table) <- .i.table[,1]
    aov.res <- list("aov.out"=.i.table)
}
sd.sum <- NULL
if(is.na(args.cmp)){
    idxFC <- which(grepl("^HSD.diff",names(aov.res[["aov.out"]]),perl=T))[1]
    cmp.id <- "CMP"
}else{
    cmp.v <- unlist(strsplit(x = args.cmp,split = ",",perl=T))
    cmp.id <- paste(cmp.v,collapse = "-")
    idxFC <- which(grepl("^HSD.diff",names(aov.res[["aov.out"]]),perl=T) & 
                   grepl(cmp.v[1],names(aov.res[["aov.out"]]),perl=T) & 
                   grepl(cmp.v[2],names(aov.res[["aov.out"]]),perl=T))
    if(args.useSD){
        print(names(aov.res[["aov.out"]]))
        mu.a <- aov.res[["aov.out"]][,sprintf("avg.%s",cmp.v[1])]
        mu.b <- aov.res[["aov.out"]][,sprintf("avg.%s",cmp.v[2])]
        sd.a <- aov.res[["aov.out"]][,sprintf("sd.%s",cmp.v[1])]
        sd.b <- aov.res[["aov.out"]][,sprintf("sd.%s",cmp.v[2])]
        #mu.a[mu.a==0] <- 1
        #mu.b[mu.b==0] <- 1
        sd.min.a <- 0.2*abs(mu.a)
        sd.min.b <- 0.2*abs(mu.b)
        f.tooSmall.a <- sd.a<sd.min.a
        f.tooSmall.b <- sd.b<sd.min.b
        print("tooSmall sd:")
        print(sum(f.tooSmall.a))
        print(sum(f.tooSmall.b))
        if(sum(f.tooSmall.a)>0){ print(str(sd.a[f.tooSmall.a])) }
        if(sum(f.tooSmall.b)>0){ print(str(sd.b[f.tooSmall.b])) }
        sd.a[f.tooSmall.a] <- sd.min.a[f.tooSmall.a]
        sd.b[f.tooSmall.b] <- sd.min.b[f.tooSmall.b]
        sd.a[sd.a==0] <- 0.2
        sd.b[sd.b==0] <- 0.2
        sd.sum <- sd.a+sd.b
    }
}
print("TEST")
print(names(aov.res[["aov.out"]])[idxFC])
print(idxFC)
out.df <- data.frame(geneID=aov.res[["aov.out"]][,"geneID"],geneSymbol=aov.res[["aov.out"]][,"geneSymbol"],stringsAsFactors = F)
out.df$log2FC <- aov.res[["aov.out"]][,idxFC]
if(args.r){ out.df$log2FC <- -out.df$log2FC }
if(!is.null(sd.sum)){
    out.df$log2FC <- out.df$log2FC/sd.sum
}
write.table(out.df,file = sprintf("%s.txt",out.prefix),quote = F,sep = "\t",row.names = F)
write.table(out.df[,c("geneID","log2FC")],file = sprintf("%s.geneID.rnk",out.prefix),quote = F,sep = "\t",row.names = F,col.names = F)
write.table(out.df[!is.na(out.df$geneSymbol),c("geneSymbol","log2FC")],file = sprintf("%s.geneSymbol.rnk",out.prefix),quote = F,sep = "\t",row.names = F,col.names = F)
dir.create(sprintf("%s.%s.GSEA",out.prefix,args.score),showWarnings = F,recursive = T)
cmd.str <- sprintf("java -cp /Share/BP/zhenglt/01.bin/GSEA/gsea2-2.2.4.jar -Xmx2048m xtools.gsea.GseaPreranked -gmx %s -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk %s -scoring_scheme %s -rpt_label %s -include_only_symbols true -make_sets true -plot_top_x 200 -rnd_seed timestamp -set_max %s -set_min 5 -zip_report false -out %s -gui false",db.file,sprintf("%s.geneSymbol.rnk",out.prefix),args.score,cmp.id,args.max,sprintf("%s.%s.GSEA",out.prefix,args.score))
###cmd.str <- sprintf("java -cp /lustre1/zeminz_pkuhpc/01.bin/GSEA/gsea2-2.2.4.jar -Xmx2048m xtools.gsea.GseaPreranked -gmx %s -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk %s -scoring_scheme %s -rpt_label %s -include_only_symbols true -make_sets true -plot_top_x 100 -rnd_seed timestamp -set_max %s -set_min 5 -zip_report false -out %s -gui false",db.file,sprintf("%s.geneSymbol.rnk",out.prefix),args.score,cmp.id,args.max,sprintf("%s.%s.GSEA",out.prefix,args.score))
print(cmd.str)
ret <- system(cmd.str)

