#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-p", "--proj", type="character", default="PROJ", help="proj name [default %(default)s]")
parser$add_argument("-i", "--inFile", type="character", required=TRUE, help="infile")
parser$add_argument("-s", "--survivalType", type="character", required=TRUE, help="output prefix")
parser$add_argument("-o", "--outputPrefix", type="character", required=TRUE, help="output prefix")
parser$add_argument("-l", "--geneList", type="character", required=TRUE, help="gene name list file, head with \"geneSymbol\"")
parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="verbose mode [default %(default)s]")
args <- parser$parse_args()
print(args)

suppressPackageStartupMessages(library("Biobase"))
##suppressPackageStartupMessages(library(GEOquery))
##suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library("magrittr"))
##suppressPackageStartupMessages(library(WGCNA))
suppressPackageStartupMessages(library("survival"))
suppressPackageStartupMessages(library("vioplot"))
suppressPackageStartupMessages(library("plotrix"))

proj.name <- args$proj
rdata.file <- args$inFile
survival.type <- args$survivalType
out.prefix <- args$outputPrefix
gene.file <- args$geneList
geneNameList <- read.table(gene.file,header = T,check.names = F,stringsAsFactors = F,sep = "\t")$geneSymbol

print(str(geneNameList))

### formatted RData
#proj.name <- "GSE39582"
#rdata.file <- "OUT.RObject/GSE39582.GPL570.gset.new.RData"
#survival.type <- "os"
#out.prefix <- "./test"
#geneNameList <- c("CCR8","MAGEH1","LAYN")

dir.create(dirname(out.prefix),showWarnings = F,recursive = T)

### genes to analyse
control.gene.list <- c("CD3G","CD3D","CD3E")

# measure in log base 2 scale
lname <- load(rdata.file)
g.GNAME <- fData(gset.new) %>% .[,"Gene.symbol"] %>% as.character
names(g.GNAME) <- fData(gset.new) %>% rownames

### check new gset object
exprs(gset.new) %>% .[1:4,1:6]
pData(gset.new) %>% head
fData(gset.new) %>% head

### make corrected gene expression table
gene.to.fetch <- unique(c(geneNameList,control.gene.list))
f.toAna <- fData(gset.new) %>% .[,"Gene.symbol"] %>% match(x=gene.to.fetch,table=.)
print("genes not in this dataset:")
print(gene.to.fetch[is.na(f.toAna)])
f.toAna <- f.toAna[!is.na(f.toAna)]
edat.block <- as.data.frame(t(exprs(gset.new)[f.toAna,]))
colnames(edat.block) <- fData(gset.new)[f.toAna,"Gene.symbol"]
geneNameList <- intersect(geneNameList,colnames(edat.block))

edat.block[,"GSet"] <- apply(edat.block[,geneNameList,drop=F],1,mean)
edat.block[,"CD3"] <-  apply(edat.block[,control.gene.list,drop=F],1,mean)
for(g in geneNameList)
{
  ### ratio to CD3
  edat.block[,sprintf("%s_%s",g,"CD3")] <- edat.block[,g]-edat.block[,"CD3"]
}
edat.block[,"GSet_CD3"] <- edat.block[,"GSet"]-edat.block[,"CD3"]
colnames(edat.block) <- make.names(colnames(edat.block))
glist.toTest <- colnames(edat.block)
dat.toAna <- cbind(pData(gset.new),edat.block)
#head(edat.block)
#head(dat.toAna)

### helper function, used to group samples according high or low expression of specified genes
groupByMedian <- function(x,colname)
{
  t.median <- median(x[,colname],na.rm = T)
  y.grp <- rep(NA,nrow(x))
  y.grp[ x[,colname] > t.median ] <- "hi"
  y.grp[ x[,colname] < t.median ] <- "lo"
  return(y.grp)
}

groupByMedianRange <- function(x,colname)
{
  t.median <- median(x[,colname],na.rm = T)
  t.mad <- mad(x[,colname],na.rm = T)
  y.grp <- rep(NA,nrow(x))
  y.grp[ x[,colname] > t.median + t.mad/10 ] <- "hi"
  y.grp[ x[,colname] < t.median - t.mad/10 ] <- "lo"
  return(y.grp)
}

run.core <- function(dat.toAna,vname,survival.type="dss",xlab.txt="(months)")
{
  f.clin <- !is.na(dat.toAna[,sprintf("%s_time",survival.type)]) &
      dat.toAna[,sprintf("%s_time",survival.type)]!="" &
      dat.toAna[,sprintf("%s_time",survival.type)]>0 &
      !is.na(dat.toAna[,vname]) &
      !is.infinite(dat.toAna[,vname])
  dat.toAna <- dat.toAna[f.clin,]
  
  dat.toAna[,sprintf("%s.Grp",vname)] <- groupByMedianRange(dat.toAna,vname)
  #print(table(dat.toAna[,sprintf("%s.Grp",vname)],useNA = "always"))
  dat.toAna <- dat.toAna[!is.na(dat.toAna[,sprintf("%s.Grp",vname)]),]
  grp <- table(dat.toAna[,sprintf("%s.Grp",vname)])
  if(length(grp)<2){ return(0) }
  event <- dat.toAna[,sprintf("%s_event",survival.type)]
  surv.obj <- Surv(dat.toAna[,sprintf("%s_time",survival.type)],event)
  t.formula <- as.formula(sprintf("surv.obj ~ %s",sprintf("%s.Grp",vname)))
  survdiff.res <- survdiff(t.formula,data = dat.toAna)
  survfit.res <- survfit(t.formula,data = dat.toAna)
  survdiff.res.p <- pchisq(survdiff.res$chisq,1,lower.tail = F)
  ucox.res <- coxph(t.formula, data=dat.toAna)
  ucox.res.summary <- summary(ucox.res)
  HR.Grp <- rownames(ucox.res.summary$conf.int) %>% gsub(".Grp",":",.) %>% gsub("[_.]","/",.)
  HR <- ucox.res.summary$conf.int[1]
  HR.lower.95 <- ucox.res.summary$conf.int[3]
  HR.upper.95 <- ucox.res.summary$conf.int[4]
  HR.coef  <- ucox.res.summary$coefficients[1]
  HR.coef.se <- ucox.res.summary$coefficients[3]
  HR.pvalue <- ucox.res.summary$coefficients[5]
  print(ucox.res.summary)

  plot(survfit.res,col=c("red","blue"),xlim=c(0,60),xlab=sprintf("%s %s",survival.type,xlab.txt),
       ylab="Percent survival",main=gsub("[_.]","/",vname),mark.time=T)
  legend("topright",legend=names(survfit.res$strata) %>% gsub(".Grp","",.,perl = T) %>% gsub("[_.]","/",.,perl = T),
         col=c("red","blue"),lty=1,lwd=2,cex=1.2)
  mtext(sprintf("n:%d (%s)\nn:%d (%s)\np value: %4.6f\nHR: %4.4f%s (%4.4f, %4.4f) (%s)\n", 
                survdiff.res$n[1], gsub("[_.]","/",gsub(".Grp","",names(survdiff.res$n[1]))),
                survdiff.res$n[2], gsub("[_.]","/",gsub(".Grp","",names(survdiff.res$n[2]))),
                survdiff.res.p,
                HR,
                if(HR.pvalue<0.001) "***" else if(HR.pvalue<0.01) "**" else if(HR.pvalue<0.05) "*" else "",
                HR.lower.95,HR.upper.95,
                HR.Grp), 
        side = 1,line = -1,adj = 0.05,cex=1.2)
  return(data.frame(proj.name=proj.name,
                    survival.type=survival.type,
                    vname=vname,n1=survdiff.res$n[1],n2=survdiff.res$n[2],
                    n1.name=gsub("[_.]","/",gsub(".Grp","",names(survdiff.res$n[1]))),
                    n2.name=gsub(".Grp","",names(survdiff.res$n[2])),
                    logRank.p=survdiff.res.p,
                    HR.Grp=HR.Grp,
                    HR=HR,
                    HR.lower.95=HR.lower.95,
                    HR.upper.95=HR.upper.95,
                    HR.pvalue=HR.pvalue,
                    HR.coef=HR.coef,
                    HR.coef.se=HR.coef.se,
                    stringsAsFactors = F))
}

run.Cox <- function(dat.ana,vname,survival.type,v.adj=NULL)
{
  #vname <- "SLC4A10"
  #v.adj <- "CD3"
  #survival.type <- "OS"
  EPSILON <- 0.1
  f.clin <- !is.na(dat.ana[,sprintf("%s",survival.type)]) & dat.ana[,sprintf("%s",survival.type)]>0 &
    !is.na(dat.ana[,vname]) & !is.infinite(dat.ana[,vname]) & dat.ana[,"Stage"]!=""
  dat.ana <- dat.ana[f.clin,]
  t.dat <- dat.ana
  if(!is.null(v.adj)){   t.dat[,vname] <- log2(t.dat[,vname]+EPSILON)-log2(t.dat[,v.adj]+EPSILON) }
  # summary(t.dat[,vname])
  # hist(t.dat[,vname],breaks = 20)
  # ###summary(dat.ana[,c("SLC4A10","LAYN","LAYN","CTLA4")])
  # res.cut <- surv_cutpoint(t.dat, time = survival.type, event = sprintf("%s_STATUS",survival.type), variables = vname)
  # summary(res.cut)
  # print(plot(res.cut, vname, palette = "npg"))
  # res.cat <- surv_categorize(res.cut)
  # print(table(res.cat[,c(sprintf("%s_STATUS",survival.type),vname)]))
  # 
  event <- t.dat[,sprintf("%s_STATUS",survival.type)]
  surv.obj <- Surv(t.dat[,survival.type],event)
  #t.formula <- as.formula(sprintf("surv.obj ~ %s",vname))
  #fit <- survfit(t.formula, data = res.cat)
  # pp <- ggsurvplot(fit, risk.table = "abs_pct", risk.table.y.text=F,pval=T)
  # print(pp)
  # survdiff.res <- survdiff(t.formula,data = res.cat)
  # print(survdiff.res)
  # res.cat <- cbind(res.cat,t.dat[,c("Stage","Risk")])
  
  #Univariate Cox regression
  #ucox.res <- coxph(as.formula(sprintf("surv.obj ~ %s", vname)), data=res.cat)
  t.dat[,vname] <- groupByMedianRange(t.dat,vname)
  ucox.res <- coxph(as.formula(sprintf("surv.obj ~ %s", vname)), data=t.dat)
  ucox.res.summary <- summary(ucox.res)
  print(ucox.res.summary)
  
  #multivariate Cox regression
  mcox.res <- coxph(as.formula(sprintf("surv.obj ~ factor(Stage) + %s", vname)), data=t.dat)
  mcox.res.summary <- summary(mcox.res)
  print(mcox.res.summary)
  # ## Testing proportional Hazards assumption
  # test.ph <- cox.zph(mcox.res)
  # print(test.ph)
  # pp <- ggcoxzph(test.ph)
  # print(pp)
  # ## Testing influential observations
  # pp <- ggcoxdiagnostics(mcox.res, type = "dfbeta",linear.predictions = FALSE, ggtheme = theme_bw())
  # print(pp)
  # pp <- ggcoxdiagnostics(mcox.res, type = "deviance",linear.predictions = FALSE, ggtheme = theme_bw())
  # print(pp)
  
  ## Testing non linearity (mainly for continious covariate)
  #####ggcoxfunctional(Surv(time, status) ~ age + log(age) + sqrt(age), data = res.cat)
  
}

###run.core(dat.toAna,"LAYN_CD3",survival.type="dss",xlab.txt="(months)")
#run.core(dat.toAna,"ENTPD1-AS1",survival.type="dss",xlab.txt="(months)")
v.ann <- data.frame(ann.name=glist.toTest, 
                    grp.name=glist.toTest,
                    stringsAsFactors = F)

### survival plot
pdf(sprintf("%s.survival.%s.%s.pdf",out.prefix,proj.name,survival.type),width=6,height=6)
par(mar=c(5,5,4,2),cex.lab=1.5,cex.main=1.5)
out.df <- c()
for(i in seq_len(nrow(v.ann)))
{
    #print(v.ann$ann.name[i])
    .i.df <- run.core(dat.toAna,v.ann$ann.name[i],survival.type,xlab.txt="(months)")
    out.df <- rbind(out.df,.i.df)
    #print("finished")
}
dev.off()
write.table(out.df,sprintf("%s.survival.%s.%s.tab",out.prefix,proj.name,survival.type),quote = F,row.names = F,sep = "\t")
