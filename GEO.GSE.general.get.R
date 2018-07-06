#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-i", "--inFile", type="character", required=TRUE, help="infile")
parser$add_argument("-p", "--proj", type="character", default="PROJ", help="GEO accession number  [default %(default)s]")
parser$add_argument("-q", "--platform", type="character", default="PLATFORM", help="platform (GPLXX)  [default %(default)s]")
parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="verbose mode [default %(default)s]")
args <- parser$parse_args()
print(args)

in.file <- args$inFile
proj.name <- args$proj
platform.name <- args$platform

library(Biobase)
library(GEOquery)
library(limma)
library(magrittr)
library(WGCNA)
library("survival")
library("vioplot")
library("plotrix")

# load series and platform data from GEO
# RMA measure in log base 2 scale
in.file <- "GSE29623-GPL570_series_matrix.txt.gz"
proj.name <- "GSE29623"
platform.name <- "GPL570"

if(!is.null(in.file) && file.exists(in.file)){
    gset <- getGEO(filename = in.file, GSEMatrix =TRUE, AnnotGPL=TRUE, destdir="/DBS/DB_temp/zhangLab/GEO/GPL")
}else{
    gset <- getGEO(proj.name, GSEMatrix =TRUE, AnnotGPL=TRUE)
    if (length(gset) > 1) idx <- grep(platform.name, attr(gset, "names")) else idx <- 1
    gset <- gset[[idx]]
}
### clean the sample info
# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))
sample.desc <- data.frame(os_event=pData(gset) %>% .[,"characteristics_ch1.10"] %>% 
                            match(table=c("os event: alive",
                                          "os event: dead"),nomatch=NA) %>%
                            subtract(1),
                          dfs_event=pData(gset) %>% .[,"characteristics_ch1.12"] %>% 
                            match(table = c("dfs event: no recurrence",
                                            "dfs event: recurrence"),nomatch=NA) %>%
                            subtract(1),
                          os_time=pData(gset) %>% .[,"characteristics_ch1.9"] %>% 
                            gsub(pattern="overall survival \\(os\\): ",replacement="",x=.,perl=T) %>% as.numeric,
                          dfs_time=pData(gset) %>% .[,"characteristics_ch1.11"] %>%
                            gsub(pattern="disease-free survival time \\(dfs\\): ",replacement="",x=.,perl=T) %>% as.numeric,
                          title=pData(gset) %>% .[,"title"] %>% as.character,
                          gender=pData(gset) %>% .[,"characteristics_ch1.1"] %>% 
                            gsub(pattern="gender: ",replacement="",x=.,perl=T) %>% factor,
                          stage=pData(gset) %>% .[,"characteristics_ch1.6"] %>%
                            gsub(pattern="ajcc staging: ",replacement="",x=.,perl=T) %>% factor,
                          grade=pData(gset) %>% .[,"characteristics_ch1.5"] %>%
                            match(table=c("histologic grade: WELL DIFF",
                                          "histologic grade: MOD DIFF",
                                          "histologic grade: POORLY DIFF")) %>% factor,
                          metastasis=pData(gset) %>% .[,"characteristics_ch1.4"] %>%
                            match(table=c("m stage (0: no; 1: metastasis): 0","m stage (0: no; 1: metastasis): 1")) %>%
                            subtract(1),
                        stringsAsFactors = F)
rownames(sample.desc) <- pData(gset) %>% rownames
### collapse expression dataset from proble level to gene level
f.geneID <- !(fData(gset) %>% .[,"Gene.ID"] %>% equals(""))
gset.withGeneID <- gset[f.geneID,]
exp.collapsed <- collapseRows(exprs(gset.withGeneID),
                              fData(gset.withGeneID) %>% .[,"Gene.ID"],
                              rownames(gset.withGeneID),  method="MaxMean")
new.fData <- fData(gset.withGeneID) %>% .[,c("Gene.ID","Gene.symbol","Gene.title")] %>% unique
rownames(new.fData) <- new.fData$Gene.ID
new.fData <- new.fData[rownames(exp.collapsed$datETcollapsed),]
new.fData$selectedRowID <- exp.collapsed$group2row[rownames(new.fData),"selectedRowID"]
head(new.fData)
f.gene.uncert <- !(rownames(new.fData) %>% grepl(pattern="\\/",x=.,perl=T))
new.fData[f.gene.uncert,] %>% head(n=4)
exp.collapsed$datETcollapsed[f.gene.uncert,] %>% .[1:4,1:8]

gset.new <- ExpressionSet(assayData = exp.collapsed$datETcollapsed[f.gene.uncert,],
                          phenoData = AnnotatedDataFrame(sample.desc),
                          featureData = AnnotatedDataFrame(new.fData[f.gene.uncert,])
                          )
### check new gset object
exprs(gset.new) %>% .[1:4,1:8]
pData(gset.new) %>% head
fData(gset.new) %>% head

save(gset.new,file = sprintf("%s.%s.gset.new.RData",proj.name,platform.name))
