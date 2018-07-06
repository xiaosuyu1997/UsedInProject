#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-d", "--designFile", type="character", required=TRUE, help="design matrix file")
parser$add_argument("-i", "--inputFile", type="character", required=TRUE, help="input file")
parser$add_argument("-o", "--outDir", type="character", required=TRUE, help="output directory")
#parser$add_argument("-c", "--cellTypeColorFile", default="/WPS1/zhenglt/work/TCR_chunhong/data/CellType.color", 
#                    type="character", help="cellTypeColorFile [default %(default)s]")
parser$add_argument("-a", "--geneVec", type="character",default="IL2RA,FOXP3", help="gene list  [default %(default)s]")
parser$add_argument("-b", "--expT", type="character",default="10,65", help="gene exp threshold  [default %(default)s]")
parser$add_argument("-n", "--withGName", action="store_true", default=FALSE, help="whether treat 2nd column as gname [default %(default)s]")
parser$add_argument("-s", "--sample", type="character", default="SAMPLE", help="sample id  [default %(default)s]")
parser$add_argument("-u", "--measure", type="character",default="tpm", help="measure  [default %(default)s]")
parser$add_argument("-l", "--log", action="store_true", default=FALSE, help="whether do log transform [default %(default)s]")
parser$add_argument("-r", "--center", action="store_true", default=FALSE, help="whether center genes by individual [default %(default)s]")
parser$add_argument("-m", "--notFilter", action="store_true", default=FALSE, help="don't filtering genes [default %(default)s]")
args <- parser$parse_args()

print(args)

designFile <- args$designFile
inputFile <- args$inputFile
out.dir <- args$outDir
#cellTypeColorFile <- args$cellTypeColorFile
geneVec <- args$geneVec
expT <- args$expT
sample.id <- args$sample
args.log <- args$log
args.center <- args$center
args.notFilter <- args$notFilter
args.withGName <- args$withGName
args.measure <- args$measure

#### TEST DATA 
#designFile <- "/WPSnew/zhenglt/work/proj_zl/integrate0924/cross.patient/OUT.Tr/colonC.Tr.divCD4.txt"
#inputFile <- "/WPSnew/zhenglt/work/proj_zl/integrate0924/cross.patient/exp/colonC.P5.scran.RData"
#out.dir <- "/WPSnew/zhenglt/work/proj_zl/integrate0924/cross.patient/test"
#cellTypeColorFile <- "/WPS1/zhenglt/work/TCR_chunhong/data/CellType.color"
#geneVec <- "IL2RA,FOXP3"
#expT <- "30,30"
#sample.id <- "colonC.gateTR"
#args.log <- F
#args.center <- F
#args.notFilter <- F
#args.withGName <- T

out.prefix <- sprintf("%s/%s.%s",out.dir,sample.id,gsub(",",".",geneVec))
print(out.prefix)

geneVec <- unlist(strsplit(x = geneVec,split = ","))
expT <- as.numeric(unlist(strsplit(expT,split=",")))

dir.create(out.dir,recursive = T,showWarnings = F)

if(file.exists("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")){
    source("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")
}else{
    source("/WPSnew/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")
}

suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("R.utils"))
suppressPackageStartupMessages(library("gridBase"))
suppressPackageStartupMessages(library("grid"))
suppressPackageStartupMessages(library("ks"))
suppressPackageStartupMessages(library("mclust"))
suppressPackageStartupMessages(library("scater"))
suppressPackageStartupMessages(require("fields"))

out.prefix <- sprintf("%s/%s",out.dir,sample.id)
g.inputD <- processInput(designFile,cellTypeColorFile = NULL,inputFile,args.notFilter,geneFile = NULL,
                         args.center,args.log,args.norm.exprs = F,args.measure = args.measure)
myDesign <- g.inputD$myDesign
sampleTypeColor <- g.inputD$sampleTypeColor
patient.col.list <- g.inputD$patient.col.list
Y <- g.inputD$Y
g.GNAME <- g.inputD$g.GNAME
##### 

loginfo("... all samples.")
###
geneVecID <- names(g.GNAME)[match(geneVec,g.GNAME)]
print(geneVecID)

###
print(range(Y))
print(Y[1:4,1:8])


dat.plot <- data.frame(x=Y[geneVecID[1],],y=Y[geneVecID[2],])
print(head(dat.plot))
rownames(dat.plot) <- colnames(Y)
dat.plot$x.bin <- dat.plot$x >= log2(expT[1]+1)
dat.plot$y.bin <- dat.plot$y >= log2(expT[2]+1)

out.df <- dat.plot[rownames(myDesign),]
colnames(out.df) <- c(geneVec,paste(geneVec,"bin",sep = "."))
out.df <- cbind(myDesign,out.df)
write.table(out.df,sprintf("%s.gated.txt",out.prefix),sep="\t",col.names=T,row.names=F,quote=F)
#pdf(file=sprintf("%s.scatter.pdf",out.prefix),width=6,height=6)
#p <- ggplot(dat.plot, aes(x = x, y = y)) +
#    geom_density_2d() +
#    theme_bw(base_size = 12)
#print(p)
#dev.off()



d_mix <- densityMclust(dat.plot[,c("x","y")])
d_mix_summary <- summary(d_mix)
#print(x_mix)
print(d_mix_summary)



##################### use basic plot
plot.2DDensity <- function(dat.plot,type="point",d_mix_summary=NULL){

    dat.plot.range <- round(range(dat.plot[,c("x","y")])+0.5)
    .density <- kde(dat.plot[,c("x","y")])

    cor.pcc.out <- cor.test(dat.plot$x,dat.plot$y,method="pearson")
    cor.spe.out <- cor.test(dat.plot$x,dat.plot$y,method="spearman")

    pdf(file=sprintf("%s.scatter.%s.%s.pdf",out.prefix,type,ifelse(is.null(d_mix_summary),"MixNULL","MixM")),width=6,height=6)
   
    zz <- c(5,10,20,30,40,50,60,70,80,90,95) 
    par(mar=c(6,6,6,6),cex.lab=1.5,cex.main=1.5)
    plot(.density,display="filled.contour2", cont=zz,xlab=geneVec[1], ylab=geneVec[2])
    # hijack
    ###plot(.density,display="filled.contour2", cont=zz,xlab=geneVec[1], ylab=geneVec[2],xlim=c(-5,15),ylim=c(-5,15),main="")
    image.plot(zlim=c(0,zz[length(zz)]),legend.only=TRUE, 
               col = c("transparent", rev(heat.colors(length(zz)))), 
               axis.args=list( at=zz, labels=sprintf("%s%%",100-zz)), legend.width=2.5,legend.mar=5.0)
    .addAnn <- function()
    {
        abline(v=log2(1+expT[1]),col="red4",lty=2,lwd=1.5)
        abline(h=log2(1+expT[2]),col="red4",lty=2,lwd=1.5)
        nn <- sum(dat.plot$x >= log2(1+expT[1]) & dat.plot$y >= log2(1+expT[2]))
        ##mtext(text = sprintf("%4.2f %%\n(%d/%d)",nn*100/nrow(dat.plot),nn,nrow(dat.plot)),side = 3,line=-2.0,adj = 0.95,xpd=T)
        text(par('usr')[2],par('usr')[4]-1,labels = sprintf("%4.2f %%\n(%d/%d)",nn*100/nrow(dat.plot),nn,nrow(dat.plot)),adj = 1.1,xpd=T)
        nn <- sum(dat.plot$x < log2(1+expT[1]) & dat.plot$y >= log2(1+expT[2]))
        ##mtext(text = sprintf("%4.2f %%\n(%d/%d)",nn*100/nrow(dat.plot),nn,nrow(dat.plot)),side = 3,line=-2.0,adj = 0.05,xpd=T)
        text(log2(1+expT[1]),par('usr')[4]-1,labels = sprintf("%4.2f %%\n(%d/%d)",nn*100/nrow(dat.plot),nn,nrow(dat.plot)),adj = 1.1,xpd=T)
        nn <- sum(dat.plot$x < log2(1+expT[1]) & dat.plot$y < log2(1+expT[2]))
        text(log2(1+expT[1]),log2(1+expT[2])-1,labels = sprintf("%4.2f %%\n(%d/%d)",nn*100/nrow(dat.plot),nn,nrow(dat.plot)),adj=1.1,xpd=T)
        nn <- sum(dat.plot$x >= log2(1+expT[1]) & dat.plot$y < log2(1+expT[2]))
        text(par('usr')[2],log2(1+expT[2])-1,labels = sprintf("%4.2f %%\n(%d/%d)",nn*100/nrow(dat.plot),nn,nrow(dat.plot)),adj=1.1,xpd=T)
        mtext(text=sprintf("Pearson Cor: %4.2f (p value: %4.2e)",cor.pcc.out$estimate,cor.pcc.out$p.value),side=3,line=2.5,adj=0.5)
        mtext(text=sprintf("Spearman Cor: %4.2f (p value: %4.2e)",cor.spe.out$estimate,cor.spe.out$p.value),side=3,line=1.5,adj=0.5)
        mtext(text=sprintf("%s",sample.id),side=3,line=0.5,adj=0.5)
    }
    .addAnn()

    ## just points
    plot(dat.plot$x,dat.plot$y,main="",pch=16, cex=min(1.0,0.8*1000/nrow(dat.plot)),
         xlab=geneVec[1], ylab=geneVec[2],xlim=dat.plot.range,ylim=dat.plot.range)
    .addAnn()

    #dat.plot.range <- NULL
    par(mar=c(5,6,0,0),cex.lab=1.5,cex.main=1.5,fig=c(0,0.8,0,0.8))
    if(type=="point"){
        plot(dat.plot$x,dat.plot$y,main="",pch=16, cex=min(1.0,0.8*1000/nrow(dat.plot)),
             xlab=geneVec[1], ylab=geneVec[2],xlim=dat.plot.range,ylim=dat.plot.range)
    #}else if(type=="density"){
        ###tsne.density <- kde2d(tsne.out[["30"]]$Rtsne.res$Y[f.points,1],tsne.out[["30"]]$Rtsne.res$Y[f.points,2],n=50)
        ###contour(tsne.density,xlab="",ylab="",add=TRUE,nlevels = 3,drawlabels=F)
        plot(.density, cont=c(50,60,70,80,90,95),lwd=1,col="blue",add=T,drawlabels=F)

    }else if(type=="dens"){
        #.density <- kde(dat.plot[,c("x","y")])
        plot(.density,display="filled.contour2", cont=c(5,10,20,30,40,50,60,70,80,90,95),xlab=geneVec[1], ylab=geneVec[2],xlim=dat.plot.range,ylim=dat.plot.range)
    }
    if(!is.null(d_mix_summary)){
        points(d_mix_summary$mean["x",],d_mix_summary$mean["y",],pch=3,cex=1.2,col="black")
    }
    #abline(v=quantile(e.x,probs=0.05),col="orange",lty=2,lwd=1.5) 
    #abline(h=quantile(e.y,probs=0.05),col="orange",lty=2,lwd=1.5)
    abline(v=log2(1+3),col="orange",lty=2,lwd=1.5) 
    abline(h=log2(1+3),col="orange",lty=2,lwd=1.5)
    abline(v=log2(1+5),col="red",lty=2,lwd=1.5)
    abline(h=log2(1+5),col="red",lty=2,lwd=1.5)
    abline(v=log2(1+10),col="red2",lty=2,lwd=1.5)
    abline(h=log2(1+10),col="red2",lty=2,lwd=1.5)
    abline(v=log2(1+expT[1]),col="red4",lty=2,lwd=1.5)
    abline(h=log2(1+expT[2]),col="red4",lty=2,lwd=1.5)
    #par(fig=c(0,1,0,1),new=F,mar=c(5,6,4,2))
    #title(main=sprintf("%s (%d/%d)",ttype,sum(f),nrow(in.table)),xpd=T)
            
    par(fig=c(0,0.8,0.8,1),new=T,mar=c(0.1,6,4,0))
    h.x <- hist(dat.plot$x,breaks=40,plot=F)
    plot(NULL, xlim=dat.plot.range, ylim=c(0,max(h.x$density)),xlab="",ylab="",xaxt="n",yaxt="n",main="",bty="n")
    ww <- h.x$mids[2]-h.x$mids[1]
    h.x$mids <- h.x$mids-ww/2
    rect(h.x$mids-ww/2,0,h.x$mids+ww/2,h.x$density,col="gray")
    ##barplot(h.x$density,xlab="",ylab="",xaxt="n",yaxt="n",main="",space=0,xlim=c(0,12))

    par(fig=c(0.8,1,0,0.8),new=T,mar=c(5,0.1,0,2))
    h.y <- hist(dat.plot$y,breaks=40,plot=F)
    plot(NULL, xlim=c(0,max(h.y$density)), ylim=dat.plot.range,xlab="",ylab="",xaxt="n",yaxt="n",main="",bty="n")
    ww <- h.y$mids[2]-h.y$mids[1]
    h.y$mids <- h.y$mids-ww/2
    rect(0,h.y$mids-ww/2,h.y$density,h.y$mids+ww/2,col="gray")
    #barplot(h.y$density,horiz = T,xlab="",ylab="",xaxt="n",yaxt="n",main="",space=0,ylim=c(0,12))
    print(h.y$mids)
    dev.off()
}

#plot.2DDensity(dat.plot,type="point")
plot.2DDensity(dat.plot,type="dens")
#plot.2DDensity(dat.plot,type="dens",d_mix_summary=d_mix_summary)



