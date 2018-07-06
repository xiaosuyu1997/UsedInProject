#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-d", "--designFile", type="character", required=TRUE, help="design matrix file")
parser$add_argument("-i", "--inputFile", type="character", required=TRUE, help="input file")
parser$add_argument("-f", "--geneFile", type="character", help="gene list file")
parser$add_argument("-c", "--cellTypeColorFile", default="/WPS1/zhenglt/work/TCR_chunhong/data/CellType.color", 
                    type="character", help="cellTypeColorFile [default %(default)s]")
parser$add_argument("-o", "--outDir", type="character", required=TRUE, help="output directory")
parser$add_argument("-s", "--sample", type="character", default="SAMPLE", help="sample id  [default %(default)s]")
parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="whether verbose mode [default %(default)s]")
parser$add_argument("-l", "--log", action="store_true", default=FALSE, help="whether do log transform [default %(default)s]")
parser$add_argument("-r", "--center", action="store_true", default=FALSE,
                    help="whether center genes by individual [default %(default)s]")
parser$add_argument("-m", "--notFilter", action="store_true", default=FALSE, help="don't filtering genes [default %(default)s]")
parser$add_argument("-n", "--nKeep", type="integer", default=1500, help="use top nKeep genes [default %(default)s]")
parser$add_argument("-b", "--binFile", type="character", help="binarized exp file [default %(default)s]")
parser$add_argument("-u", "--iterMax", type="integer", default=4, help="number of iter [default %(default)s]")
parser$add_argument("-z", "--kOptimal", type="character", help="specified optimal k [default %(default)s]")
parser$add_argument("-x", "--densityPar", type="character", help="specified densityPar file [default %(default)s]")
parser$add_argument("-e", "--geneIDType", type="character",default="entrezID",
                    help="geneID type (entrezID, ensemblID) [default %(default)s]")
parser$add_argument("-j", "--measure", type="character", help="measure (tpm, exprs, norm_exprs)[default %(default)s]")
#### method specific
parser$add_argument("-y", "--transform", type="character", default="none",
                    help="do transformation(none, pca, eigenmap, pca_tsne, iCor_tsne) [default %(default)s]")
parser$add_argument("-k", "--kbatch", type="character", default="2,3,4,5,6,7,8,9,10", help="kbatch [default %(default)s]")
parser$add_argument("-a", "--valSpace", type="character", default="tSNE",
                    help="calculate distance in valSpace (original, tSNE, pca) [default %(default)s]")
parser$add_argument("-w", "--method", type="character", default="kmeans", 
                    help="clustering method (kmeans, hclust, density) [default %(default)s]")
parser$add_argument("-q", "--updateG", action="store_true", default=FALSE,
                    help="update g.f then cluster once [default %(default)s]")
parser$add_argument("-p", "--selectNCP", action="store_true", default=FALSE,
                    help="select number of componets [default %(default)s]")
parser$add_argument("-t", "--myseed", type="character", help="seed [default %(default)s]")
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
args.measure <- args$measure

args.transform <- args$transform
args.kbatch <- args$kbatch
if(!is.null(args.kbatch)) { args.kbatch <- as.numeric(unlist(strsplit(args.kbatch,","))) }
args.method <- args$method
args.myseed <- args$myseed
args.geneIDType <- args$geneIDType
args.valSpace <- args$valSpace
binFile <- args$binFile
args.iterMax <- args$iterMax
args.kOptimal <- args$kOptimal
if(!is.null(args.kOptimal) && file.exists(args.kOptimal)){
    source(args.kOptimal)
    print(args.kOptimal.list)
}else{
    args.kOptimal.list <- NULL
}
args.densityPar.file <- args$densityPar
if(!is.null(args.densityPar.file) && file.exists(args.densityPar.file)){
    source(args.densityPar.file)
    print(args.densityPar)
}else{
    args.densityPar <- NULL
}
args.select.ncp <- args$selectNCP
args.updateG <- args$updateG


#clonotype.strict.data <- read.clonotype(in.file = clonotype.file,ctype.col = "C_strict")

##### frequently used code
dir.create(out.dir,recursive = T,showWarnings = F)
if(file.exists("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")){
	source("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")
}else{
	source("/lustre1/zeminz_pkuhpc/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")
}
suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("densityClust"))

g.inputD <- processInput(designFile,cellTypeColorFile,inputFile,args.notFilter,geneFile,args.center,args.log,
                         args.measure=args.measure)
myDesign <- g.inputD$myDesign
sampleTypeColor <- g.inputD$sampleTypeColor
patient.col.list <- g.inputD$patient.col.list
Y <- g.inputD$Y
g.GNAME <- g.inputD$g.GNAME
#### rename rownames
Y.newRowname <- Y
newName <- g.GNAME[rownames(Y)]
names(newName) <- rownames(Y)
f.gname.na <- is.na(newName)
newName[f.gname.na] <- names(newName)[f.gname.na]
rownames(Y.newRowname) <- unname(newName)
#if(args.geneIDType=="ensemblID"){ Y <- Y.newRowname }

##### 

#### select variable genes
Y.sd <- apply(Y, 1, sd)
Y.sd.sort <- sort(Y.sd,decreasing=TRUE)
Y.sd.df <- data.frame(geneID=names(Y.sd.sort),geneSymbol=g.GNAME[names(Y.sd.sort)])
Y.sd.df <- cbind(Y.sd.df,Y.sd.sort)
write.table(Y.sd.df,sprintf("%s/%s.Y.sd.sort.txt",out.dir,sample.id),row.names = F,sep = "\t",quote = F)

if(args.transform=="iCor_tsne"){
    g.f <- names(Y.sd.sort)
}else{
    g.f <- head(names(Y.sd.sort),n=nKeep)
#doit.kmeans(Y,g.f,extra="",myseed=NULL,data.forDE=Y.newRowname)
}

### pca or other transformation
#if(args.transform=="pca"){
#    g.out.prefix <- sprintf("%s/%s.%s.pca",out.dir,sample.id,args.method)
#    #g.do.pca <- T
#}else if(args.transform=="pca_tsne"){
#    g.out.prefix <- sprintf("%s/%s.%s.pca_tsne",out.dir,sample.id,args.method)
#    #g.do.pca <- F
#}else if(args.transform=="none"){
#    g.out.prefix <- sprintf("%s/%s.%s.none",out.dir,sample.id,args.method)
#    #g.do.pca <- F
#}
g.out.prefix <- sprintf("%s/%s.%s.%s",out.dir,sample.id,args.method,args.transform)

do.it <- function(g.f,out.prefix,sampleType,sampleTypeColor,iter.samples,level=1,clusterName="L1C1",ngenes=c(1000,1500,2000,2500))
{
    ##clustering.out <- runSimpleClusteringAnalysis(Y[g.f,],g.out.prefix,myDesign$sampleType,colSet = sampleTypeColor,B=100,
    ##                  legend=c(names(sampleTypeColor)[names(sampleTypeColor) %in% unique(as.character(myDesign$sampleType))],
    dir.create(dirname(out.prefix),showWarnings = F,recursive = T)
    clustering.out <- runSimpleClusteringAnalysis(Y[g.f,iter.samples],out.prefix,sampleType,colSet = sampleTypeColor,B=100,
                      legend=c(names(sampleTypeColor)),
                      col.legend=c(sampleTypeColor),
                      col.points = sampleTypeColor[sampleType],
                      k.batch = args.kbatch,do.scale=T,data.forDE = Y,method=args.method,
                      myseed = args.myseed,val.space = args.valSpace,gid.mapping = g.GNAME,
                      select.ncp = args.select.ncp,transform = args.transform,
                      opt.rho = args.densityPar[[sprintf("%s.rho",clusterName)]],
                      opt.delta = args.densityPar[[sprintf("%s.delta",clusterName)]])
                      ###k.batch = args.kbatch,do.scale=T,data.forDE = Y.newRowname,method=args.method,
    .odir <- dirname(out.prefix)
    if(args.method=="density" && args.updateG && 
            length(clustering.out$tsne.res.de.list)>0 &&
            !file.exists(sprintf("%s.old",.odir)) &&
            !is.null(clustering.out) && 
            !is.null(clustering.out$tsne.res.de.list[[1]])){
        g.f <- clustering.out$diffGene.list[[1]]$geneID
        system(sprintf("mv %s %s.old",.odir,.odir))
        dir.create(.odir,showWarnings = F,recursive = T)
        clustering.out <- runSimpleClusteringAnalysis(Y[g.f,iter.samples],out.prefix,sampleType,colSet = sampleTypeColor,B=100,
                      legend=c(names(sampleTypeColor)),
                      col.legend=c(sampleTypeColor),
                      col.points = sampleTypeColor[sampleType],
                      k.batch = args.kbatch,do.scale=T,data.forDE = Y,method=args.method,
                      myseed = args.myseed,val.space = args.valSpace,gid.mapping = g.GNAME,
                      select.ncp = args.select.ncp,transform = args.transform,
                      ##select.ncp = args.select.ncp,transform = "pca_tsne",
                      opt.rho = args.densityPar[[sprintf("%s.rho.update",clusterName)]],
                      opt.delta = args.densityPar[[sprintf("%s.delta.update",clusterName)]])

    }
    if(is.null(clustering.out)){return(NULL)}
    g.clust.res[[sprintf("%s",clusterName)]] <<- clustering.out
    ##### recursively 
    best.k <- NULL
    if((level<args.iterMax)){
        if(!is.null(args.kOptimal.list)){
            best.k <- args.kOptimal.list[[clusterName]]
        }
        if(is.null(best.k) && nrow(clustering.out$clustering.validation.df.best)>0){
            best.k <- clustering.out$clustering.validation.df.best[1,"k"]
        }
    }
    if(!is.null(best.k) && best.k>1){
        this.cluster <- clustering.out$res.list[[as.character(best.k)]]$cluster
        this.cluster <- this.cluster[this.cluster>0]
        print(table(this.cluster))
        for(cls in sort(unique(this.cluster))){
            cls.sample <- names(this.cluster)[this.cluster==cls]
            print(cls.sample)
            if(length(cls.sample)>10){
                {
                    cls.Y.sd <- apply(Y[,cls.sample,drop=F], 1, sd)
                    cls.Y.sd.sort <- sort(cls.Y.sd,decreasing=TRUE)
                    if(args.transform=="iCor_tsne"){
                        cls.g.f <- names(cls.Y.sd.sort)
                    }else{
                        cls.g.f <- head(names(cls.Y.sd.sort),n=nKeep)
                    }
                }
                cls.aid <- sprintf("%s.L%sC%s",clusterName,level+1,cls)
                #cls.out.prefix <- sprintf("%s/%s/%s.simple",out.dir,cls.aid,sample.id)
                cls.out.prefix <- sprintf("%s/%s/%s.simple",dirname(dirname(out.prefix)),cls.aid,sample.id)
                ####dir.create(dirname(cls.out.prefix))
                print(cls.out.prefix)
                if(length(cls.g.f)>10){
                    do.it(cls.g.f,
                          out.prefix=cls.out.prefix,
                          sampleType=myDesign[cls.sample,"sampleType"],
                          sampleTypeColor=sampleTypeColor,
                          iter.samples=cls.sample,
                          level=(level+1),
                          clusterName=cls.aid,
                          ngenes=c(1000,1500,2000,2500))
                }
            }
        }
    }
}

g.clust.res <- list()

### recurvely clustering
do.it(g.f,
      out.prefix=sprintf("%s/L1C1/%s.simple",out.dir,sample.id),
      sampleType=myDesign$sampleType,
      sampleTypeColor=sampleTypeColor[names(sampleTypeColor) %in% unique(as.character(myDesign[,"sampleType"]))],
      iter.samples=colnames(Y),
      level=1,
      clusterName="L1C1",
      ngenes=c(1000,1500,2000,2500))


### collect result
get.iter.cluster <- function(g.clust.res)
{
    iter.cluster.df <- NULL
    for(aid in names(g.clust.res)){
        print(aid)
        aid.clustering.validation.df.best <- g.clust.res[[aid]][["clustering.validation.df.best"]]

        aid.best.k <- NULL
        if(!is.null(args.kOptimal.list)){
            aid.best.k <- args.kOptimal.list[[aid]]
        }
        if(is.null(aid.best.k) && nrow(aid.clustering.validation.df.best)>0){
            aid.best.k <- aid.clustering.validation.df.best[1,"k"]
        }
        cat(sprintf("best: %d\n",aid.best.k))
        if(!is.null(aid.best.k) && aid.best.k>1){
            aid.cluster <- g.clust.res[[aid]][["res.list"]][[as.character(aid.best.k)]][["cluster"]]
            aid.cluster.df <- data.frame(sample=names(aid.cluster),stringsAsFactors = F)
            aid.cluster.df[,aid] <- aid.cluster
            if(is.null(iter.cluster.df)){
                iter.cluster.df <- aid.cluster.df
            }else{
                iter.cluster.df <- left_join(iter.cluster.df,aid.cluster.df)
            }
        }
    }
    if(is.null(iter.cluster.df)) { return(iter.cluster.df) }
    cluster.label <- colnames(iter.cluster.df)[2:(ncol(iter.cluster.df))]
    cluster.label.L <- str_match_all(cluster.label,"L(\\d+)")%>%
        sapply(.,function(x){x[nrow(x),2]}) %>% as.integer
    m.label <- apply(iter.cluster.df[,cluster.label,drop=F],1,function(x){
              .i <- max(which(!is.null(x) & x>0))
              if(is.infinite(.i)){
                  return("L1C1.L2C0")
              }else{
                  return(sprintf("%s.L%sC%s",cluster.label[.i],cluster.label.L[.i]+1,x[.i]))
              }
          })
    ### if one m.label is maximum-prefix of other two m.labels, mark it as 'uncharacterized'
    m.label.set <- unique(m.label)
    m.label.set.maxPrefix <- sapply(str_match_all(m.label.set,"(.+\\d+)\\."),function(x){ x[1,2] })
    m.label.set.to.outlier  <- m.label.set %>% sapply(.,function(x){ sum(x==m.label.set.maxPrefix)>0 })
    m.label[m.label.set.to.outlier[m.label]] <- "L1C1.L2C0"
    ### cluster size  smaller than 5 also marked as 'uncharacterized'
    m.label.set.to.size <- table(m.label)
    m.label[m.label.set.to.size[m.label]<5] <- "L1C1.L2C0"
    print(table(m.label))
    iter.cluster.df[,"m.label"] <- m.label
    write.table(right_join(myDesign,iter.cluster.df),sprintf("%s/%s.iter.cluster.txt",out.dir,sample.id),row.names = F,quote = F,sep="\t")
    return(iter.cluster.df)
}

iter.cluster.df <- get.iter.cluster(g.clust.res)

if(!is.null(iter.cluster.df))
{
    rownames(iter.cluster.df) <- iter.cluster.df$sample

    ### re-validation

    clusterColor <- structure(colorRampPalette(brewer.pal(12,"Paired"))(length(unique(iter.cluster.df$m.label))),
                              names=sort(as.character(unique(iter.cluster.df$m.label))))
    clusterColor["L1C1.L2C0"] <- "gray"

    iter.sample.noUnc <- iter.cluster.df$sample[iter.cluster.df$m.label!="L1C1.L2C0"]

    ### diffGene
    iter.gene.table <- my.clusterMarkerGene(Y[,iter.sample.noUnc],
                                            clust=iter.cluster.df[iter.sample.noUnc,"m.label"],
                                            ann.col=sampleTypeColor[myDesign[iter.sample.noUnc,"sampleType"]],
                                            clust.col=clusterColor,
                                            out.prefix=sprintf("%s/%s.iter.clust", out.dir,sample.id),
                                            n.cores=args.ncores,
                                            sampleType = myDesign[iter.sample.noUnc,"sampleType"],
                                            sampleTypeColSet = sampleTypeColor,original.labels=(args.geneIDType=="ensemblID"))

    if(!is.null(iter.gene.table) && nrow(iter.gene.table$gene.table)>10){
        iter.tsne.res.update <- runTSNEAnalysis(Y[rownames(iter.gene.table$gene.table),iter.cluster.df$sample,drop=F],
                                                sprintf("%s/%s.iter.tsne.update.00",out.dir,sample.id),
                                 col.points = sampleTypeColor[myDesign[iter.cluster.df$sample,"sampleType"]],
                                 legend=names(sampleTypeColor),
                                 col.legend=sampleTypeColor,
                                 myseed=args.myseed,do.dbscan = F,do.scale = T,distance.metric = NULL)
        iter.tsne.Y <- iter.tsne.res.update$`30`$Rtsne.res$Y
        plot.tsne.points(iter.tsne.Y,
                    sprintf("%s/%s.iter.tsne.update.pdf",out.dir,sample.id),
                    tsne.col.points=clusterColor[iter.cluster.df$m.label],
                    col.tsne.legend=clusterColor,
                    tsne.legend=names(clusterColor),
                    pch=20,nclusters=length(clusterColor))
    
        ### re-evaluation
        val.res <- my.clusteringValidation(t(iter.tsne.Y)[,iter.sample.noUnc,drop=F],
                                       clust=iter.cluster.df[iter.sample.noUnc,"m.label"],
                                       out.prefix=sprintf("%s/%s.iter.clust.val.00",out.dir,sample.id),
                                       n.cores=args.ncores,dist.obj=NULL)
        
        val.res.sil.df <- data.frame(sample=iter.sample.noUnc,stringsAsFactors = F,
                                cluster.old=val.res$sil[,1],
                                neighbor=val.res$sil[,2],
                                sil_width=val.res$sil[,3],
                                cluster=val.res$sil[,1],
                                m.label.old=iter.cluster.df[iter.sample.noUnc,"m.label"],
                                m.label=iter.cluster.df[iter.sample.noUnc,"m.label"])
        val.res.sil.df$cluster[val.res.sil.df$sil_width<0] <- 0
        val.res.sil.df$m.label[val.res.sil.df$sil_width<0] <- "L1C1.L2C0"
        rownames(val.res.sil.df) <- val.res.sil.df$sample
        val.res.sil.df.flt <- subset(val.res.sil.df,sil_width>0)
        write.table(right_join(myDesign,val.res.sil.df),sprintf("%s/%s.iter.sil.txt",out.dir,sample.id),row.names = F,quote = F,sep="\t")
        rm.merged.outlier <- T
        if(rm.merged.outlier){
            val.res <- my.clusteringValidation(t(iter.tsne.Y)[,rownames(val.res.sil.df.flt),drop=F],
                                           clust=iter.cluster.df[rownames(val.res.sil.df.flt),"m.label"],
                                           out.prefix=sprintf("%s/%s.iter.clust.val",out.dir,sample.id),
                                           n.cores=args.ncores,dist.obj=NULL)
            plot.tsne.points(iter.tsne.Y[rownames(val.res.sil.df.flt),],
                     sprintf("%s/%s.iter.tsne.update.noUnc.pdf",out.dir,sample.id),
                     tsne.col.points=clusterColor[val.res.sil.df.flt[,"m.label"]],
                     col.tsne.legend=clusterColor,
                     tsne.legend=names(clusterColor),
                     pch=20,nclusters=length(clusterColor))
            plot.tsne.points(iter.tsne.Y[rownames(val.res.sil.df),],
                     sprintf("%s/%s.iter.tsne.update.noUnc.00.pdf",out.dir,sample.id),
                     tsne.col.points=clusterColor[val.res.sil.df[,"m.label"]],
                     col.tsne.legend=clusterColor,
                     tsne.legend=names(clusterColor),
                     pch=20,nclusters=length(clusterColor))

            iter.gene.table.rmLowSil <- my.clusterMarkerGene(Y[,rownames(val.res.sil.df.flt)],
                                            clust=iter.cluster.df[rownames(val.res.sil.df.flt),"m.label"],
                                            ann.col=sampleTypeColor[myDesign[rownames(val.res.sil.df.flt),"sampleType"]],
                                            clust.col=clusterColor,
                                            out.prefix=sprintf("%s/%s.iter.clust.rmLowSil", out.dir,sample.id),
                                            n.cores=args.ncores,
                                            sampleType = myDesign[rownames(val.res.sil.df.flt),"sampleType"],
                                            sampleTypeColSet = sampleTypeColor,original.labels=(args.geneIDType=="ensemblID"))

        }

        clustering.validation.df <- c()
        clustering.validation.df <- rbind(clustering.validation.df, structure(c(NA,NA,
                                            ifelse(!is.null(val.res),val.res$c.stats$avg.silwidth,NA),
                                            ifelse(!is.null(val.res),val.res$c.stats$dunn,NA),
                                            if(!rm.merged.outlier) length(iter.sample.noUnc) else sum(val.res.sil.df.flt$cluster>0),
                                            ncol(Y),
                                            length(unique(iter.cluster.df[iter.sample.noUnc,"m.label"]))),
                                                                              names=c("k","note","avg.silwidth","dunn",
                                                                                      "clustered.sample","total.sample",
                                                                                      "num.clusters")))
        clustering.validation.df <- as.data.frame.matrix(clustering.validation.df)
        print(clustering.validation.df)
        write.table(clustering.validation.df, file=sprintf("%s/%s.iter.clust.val.txt",out.dir,sample.id),
                    row.names = F,sep = "\t",quote = F)

        ##### buid random forest model
    }
}



