do.plot.integratedFigure <- function(newick.file=NULL,perplexity=30,CD8CD4=F){
    ###if(is.null(newick.file)) {
    ###    pdf(file=sprintf("%s.integrated.tsne.pdf",out.prefix),width=8,height=8)
    ###}else{
    ###    pdf(file=sprintf("%s.integrated.tsne.newick.pdf",out.prefix),width=8,height=8)
    ###}
    pdf(file=sprintf("%s.integrated.perplexity%d.tsne.contour.pdf",out.prefix,perplexity),width=8,height=8)
    ###layout(matrix(c(1,2,3,4,5,5,5,5),ncol=4,byrow = T),widths = c(0.55,0.15,0.30,0.25),heights = c(0.85,0.15))
    par(mar=c(5,5,4,2),cex.lab=1.5,cex.main=1.5)
    ###plot(tsne.out[["30"]]$Rtsne.res$Y[,1],tsne.out[[as.character(perplexity)]]$Rtsne.res$Y[,2], 
    plot(tsne.out[[as.character(perplexity)]]$Rtsne.res$Y, 
         t='p',pch=16,col="lightgray", main="BarnesHutSNE",xlab="Dim1",ylab="Dim2",cex.lab=2.0,cex.axis=1.5,cex.main=2.0)
    for(m.cls in unique(sample.desc$majorCluster)){
        f.points <- sample.desc$majorCluster==m.cls
        ###tsne.density <- kde2d(tsne.out[["30"]]$Rtsne.res$Y[f.points,1],tsne.out[["30"]]$Rtsne.res$Y[f.points,2],n=50)
        ###contour(tsne.density,xlab="",ylab="",add=TRUE,nlevels = 3,drawlabels=F)
        tsne.density <- kde(tsne.out[[as.character(perplexity)]]$Rtsne.res$Y[f.points,])
        plot(tsne.density, cont=c(80),lwd=4,col=majorClusterColor[m.cls],add=T,drawlabels=F)

    }
    dev.off()

    pdf(file=sprintf("%s.integrated.perplexity%d.tsne.colPoints.pdf",out.prefix,perplexity),width=8,height=8)
    par(mar=c(5,5,4,2),cex.lab=1.5,cex.main=1.5)
    ###plot(tsne.out[["30"]]$Rtsne.res$Y[,1],tsne.out[[as.character(perplexity)]]$Rtsne.res$Y[,2], 
    plot(tsne.out[[as.character(perplexity)]]$Rtsne.res$Y, 
         t='p',pch=16,col="lightgray", main="BarnesHutSNE",xlab="Dim1",ylab="Dim2",cex.lab=2.0,cex.axis=1.5,cex.main=2.0)
    f.points <- !grepl("uncharacterized",sample.desc$majorCluster,perl=T)
    points(tsne.out[[as.character(perplexity)]]$Rtsne.res$Y[f.points,],
           col=majorClusterColor[sample.desc$majorCluster[f.points]],pch=16)
    dev.off()
    ####legend("bottom",horiz = T,legend=names(leafClusterColor),fill = NULL,inset = -0.01,xpd = NA,cex=1.5,pch=16,border =NA,col = leafClusterColor)
    #opar <- par(mar=c(5,11,4,8),cex.lab=1.5,cex.main=1.5)
    if(!is.null(newick.file) && file.exists(newick.file)){
        hc.merged <- ReadDendrogram(newick.file,convertBlanks = F)
    }else{
        ### hc on leafClusters
        leaf.cluster <- unique(sample.desc$leafCluster)
        leaf.cent <- NULL
        for(i in seq_along(leaf.cluster)){
            leaf.cent <- cbind(leaf.cent,rowMeans(Y[,sample.desc[colnames(Y),"leafCluster"]==leaf.cluster[i],drop=F]))
        }
        colnames(leaf.cent) <- leaf.cluster
        if(CD8CD4){
            f <- grepl("^CD8",leaf.cluster,perl = T)
            if(sum(f)>1){
                hc.CD8.leafCluster <- as.dendrogram(hclust(dist(t(leaf.cent[,f,drop=F]))^2, method = "complete", 
                                     members = table(sample.desc$leafCluster)[leaf.cluster[f]]))
            }else{
                hc.CD8.leafCluster <- NULL
            }
            f <- grepl("^CD4",leaf.cluster,perl = T)
            if(sum(f)>1){
                hc.CD4.leafCluster <- as.dendrogram(hclust(dist(t(leaf.cent[,f,drop=F]))^2, method = "complete", 
                                     members = table(sample.desc$leafCluster)[leaf.cluster[f]]))
            }else{
                hc.CD4.leafCluster <- NULL
            }
            if(!is.null(hc.CD8.leafCluster) && is.null(hc.CD4.leafCluster)){
                hc.merged <- hc.CD8.leafCluster
            }else if(is.null(hc.CD8.leafCluster) && !is.null(hc.CD4.leafCluster)){
                hc.merged <- hc.CD4.leafCluster
            }else if(!is.null(hc.CD8.leafCluster) && !is.null(hc.CD4.leafCluster)){
                hc.merged <- merge(hc.CD8.leafCluster,hc.CD4.leafCluster)
            }else{
                loginfo("both hc.CD8.leafCluster and hc.CD4.leafCluster are NULL! error will occur!")
            }
        }else{
            hc.merged <- as.dendrogram(hclust(dist(t(leaf.cent[,,drop=F]))^2, method = "complete", 
                                                       members = table(sample.desc$leafCluster)[leaf.cluster]))
        }
    }
    ###
    mapping.leafToMajor <- unique(sample.desc[,c("leafCluster","majorCluster")])
    mapping.leafToMajor.v <- mapping.leafToMajor$majorCluster
    names(mapping.leafToMajor.v) <- mapping.leafToMajor$leafCluster
    print(str(mapping.leafToMajor.v))

    #hc.merged <- set(hc.merged,"branches_k_col",leafClusterColor[leaf.cluster[!grepl("uncharacterized",leaf.cluster)][order.dendrogram(hc.merged)]])
    hc.merged <- set(hc.merged,"branches_k_col",majorClusterColor[mapping.leafToMajor.v[labels(hc.merged)]])
    hc.merged <- set(hc.merged,"branches_lwd",3)
   
    pdf(file=sprintf("%s.integrated.perplexity%d.tree.pdf",out.prefix,perplexity),width=8,height=9)
    layout(matrix(c(1,2),ncol=2,byrow = T),widths = c(0.10,0.30),heights = c(1))
    par(mar=c(6,1,5,0),cex.lab=1.5,cex.main=1.5)
    plot(hc.merged,horiz=T,yaxt="n",leaflab="none")
    ### barplot describing cell origin
    dat.plot <- table(sample.desc[,c("sampleType","leafCluster")])
    dat.plot <- dat.plot[,labels(hc.merged)]
    dat.plot <- apply(dat.plot,2,function(x){x/sum(x)})
    dat.plot.tree <- dat.plot
    ###par(mar=c(5,11,4,8),cex.lab=1.5,cex.main=1.5)
    par(mar=c(6,1,5,14),cex.lab=1.5,cex.main=1.5)
    xx <- barplot(dat.plot,horiz = T,beside=F,col=sampleTypeColor[rownames(dat.plot)],yaxt="n",
                  sub = "Cell Origin",cex.axis = 1.5,cex.sub = 2.0)
    ####staxlab(2,at = xx,labels=colnames(dat.plot),srt=0, cex=1.0,adj=1,top.line=0.5)
    text(rep(1.02,length(xx)),xx,labels=colnames(dat.plot), cex=1.05,adj=0,xpd=NA)
    points(x=rep(-0.05,length(xx)),y=xx,cex=3,col=majorClusterColor[mapping.leafToMajor.v[colnames(dat.plot)]],pch=16,xpd=NA)
    ##legend("bottomright",legend=rownames(dat.plot),fill=sampleTypeColor[rownames(dat.plot)],
    ##       border=sampleTypeColor[rownames(dat.plot)],cex=1.5,inset=c(-0.08,0.05),xpd=NA,horiz = T)
    legend("top",legend=rownames(dat.plot),fill=sampleTypeColor[rownames(dat.plot)],ncol=3,
           border=sampleTypeColor[rownames(dat.plot)],cex=1.2,inset=c(0.0,-0.12),xpd=NA,horiz = F)
    ###legend("topright",legend=rownames(dat.plot),fill=sampleTypeColor[rownames(dat.plot)],
    ###       border=sampleTypeColor[rownames(dat.plot)],cex=1.5,inset=c(-0.30,0),xpd=NA)
    dev.off()

    ### barplot describing patient
    pdf(file=sprintf("%s.integrated.perplexity%d.tree2.pdf",out.prefix,perplexity),width=8,height=8)
    dat.plot <- table(sample.desc[,c("patient","leafCluster")])
    dat.plot <- dat.plot[,labels(hc.merged)]
    par(mar=c(6,12,4,6),cex.lab=1.5,cex.main=1.5)
    patientColor <- auto.colSet(n=nrow(dat.plot),name="Accent")
    xx <- barplot(dat.plot,horiz = T,beside=F,col=patientColor,yaxt="s",las=2,
                  sub = "Patient",cex.axis = 1.5,cex.sub = 2.0)
    legend("bottomright",legend=rownames(dat.plot),fill=patientColor,
           border=patientColor,cex=1.5,inset=c(-0.13,0.05),xpd=NA)
    dev.off()
  
    ### scater: leafCluster ~ sampleType, by patient
    for(pp in unique(sample.desc$patient)){
        f.CD8 <- grepl("^CD8",sample.desc$leafCluster,perl = T)
        f.CD4 <- grepl("^CD4",sample.desc$leafCluster,perl = T)
        pdf(file=sprintf("%s.integrated.perplexity%d.%s.scatter.pdf",out.prefix,perplexity,pp),
            width=8,height=5*(as.numeric(sum(f.CD8)>0) + as.numeric(sum(f.CD4)>0)) )
        par(mar=c(5,5,4,2),cex.lab=1.5,cex.main=1.5)
        pushViewport(viewport(layout=grid.layout( as.numeric(sum(f.CD8)>0) + as.numeric(sum(f.CD4)>0),1)))
        if(sum(f.CD8)>0){
            dat.plot <- table(subset(sample.desc,patient==pp & f.CD8, select=c("sampleType","leafCluster")))
            dat.plot <- melt(dat.plot)
            dat.plot <- subset(dat.plot,value>0)
            print(ggplot(dat.plot, aes(sampleType, leafCluster)) + geom_point(aes(size = value)) + scale_size_area(max_size = 12) +
                  theme_minimal() + ggtitle(sprintf("%s (CD8)",pp)) + xlab("") + ylab(""), 
                  vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
        }
        if(sum(f.CD4)>0){
            dat.plot <- table(subset(sample.desc,patient==pp & f.CD4,select=c("sampleType","leafCluster")))
            dat.plot <- melt(dat.plot)
            dat.plot <- subset(dat.plot,value>0)
            print(ggplot(dat.plot, aes(sampleType, leafCluster)) + geom_point(aes(size = value)) + scale_size_area(max_size = 12) +
                  theme_minimal() + ggtitle(sprintf("%s (CD4)",pp)) + xlab("") + ylab(""), 
                  vp=viewport(layout.pos.row = 1+as.numeric(sum(f.CD8)>0), layout.pos.col = 1))
        }
        dev.off()
    }
    ### bar plot descrie subtype in TTC and TTH(R)
    s.table <- aggregate(sample~patient+sampleType+majorCluster,data=sample.desc,FUN=length)
    s.table$patient <- factor(s.table$patient,
                                    levels=c("P0407", "P0205", "P0508", "P1116", "P0322"))
    s.table$majorCluster <- factor(s.table$majorCluster,
                                    levels=c("CD8_RGS1.exhausted","CD8_RGS1.GZMK","CD8_SLC4A10","CD8_CX3CR1","CD8_naive",
                                             "CD4_GZMA.cytotoxic","CD4_GZMA.main","CD4_GZMA.exhausted","CD4_P.other","CD4_P.PTR","CD4_FOXP3-CTLA4"))
    dat.plot <- s.table[grepl("^TTC",s.table$sampleType,perl=T) & s.table$majorCluster!="CD8_naive",,drop=F]
    if(nrow(dat.plot)>0){
        dat.plot <- acast(dat.plot,majorCluster~patient,value.var="sample")
        dat.plot[is.na(dat.plot)] <- 0
        ##apply(dat.plot,2,sum)
        dat.plot <- apply(dat.plot,2,function(x){ x/sum(x) } )
        print(dat.plot)
        dat.plot.df <- data.frame(subtype=rownames(dat.plot))
        dat.plot.df <- cbind(dat.plot.df,dat.plot)
        write.table(dat.plot.df,sprintf("%s.integrated.perplexity%d.CD8.subtypeDist.txt",out.prefix,perplexity),quote = F,sep = "\t",row.names = F)

        pdf(sprintf("%s.integrated.perplexity%d.CD8.subtypeDist.pdf",out.prefix,perplexity),width=8,height=8)
        par(mar=c(5,6,4,11),cex.lab=2,cex.main=2,cex.axis=1.5)
        xx <- barplot(dat.plot,beside = F,xlab="",ylab="Cell Subtype",main="CD8+ T cell",col=majorClusterColor[rownames(dat.plot)],xaxt="n")
        #staxlab(1,at = xx,labels=names(dat.plot),srt=if(length(dat.plot) >= 5) 45 else 0, cex=1.3,adj=0.5,top.line=2)
        staxlab(1,at = xx,labels=colnames(dat.plot),srt=if(length(dat.plot) >= 5) 45 else 0, cex=1.3)
        legend("topright",legend=rownames(dat.plot),fill=majorClusterColor[rownames(dat.plot)],inset = c(-0.45,0),xpd = T,horiz = F)
        dev.off()
    }
    dat.plot <- s.table[grepl("^(TTH|TTR)",s.table$sampleType,perl=T) & !grepl("^CD4_P",s.table$majorCluster,perl=T),,drop=F]
    if(nrow(dat.plot)>0){
        dat.plot <- acast(dat.plot,majorCluster~patient,value.var="sample",fun.aggregate=sum)
        dat.plot[is.na(dat.plot)] <- 0
        ##apply(dat.plot,2,sum)
        dat.plot <- apply(dat.plot,2,function(x){ x/sum(x) } )
        print(dat.plot)
        dat.plot.df <- data.frame(subtype=rownames(dat.plot))
        dat.plot.df <- cbind(dat.plot.df,dat.plot)
        write.table(dat.plot.df,sprintf("%s.integrated.perplexity%d.CD8.subtypeDist.txt",out.prefix,perplexity),quote = F,sep = "\t",row.names = F)

        pdf(sprintf("%s.integrated.perplexity%d.CD4.subtypeDist.pdf",out.prefix,perplexity),width=8,height=8)
        par(mar=c(5,6,4,11),cex.lab=2,cex.main=2,cex.axis=1.5)
        xx <- barplot(dat.plot,beside = F,xlab="",ylab="Cell Subtype",main="CD4+ T cell",col=majorClusterColor[rownames(dat.plot)],xaxt="n")
        #staxlab(1,at = xx,labels=names(dat.plot),srt=if(length(dat.plot) >= 5) 45 else 0, cex=1.3,adj=0.5,top.line=2)
        staxlab(1,at = xx,labels=colnames(dat.plot),srt=if(length(dat.plot) >= 5) 45 else 0, cex=1.3)
        legend("topright",legend=rownames(dat.plot),fill=majorClusterColor[rownames(dat.plot)],inset = c(-0.45,0),xpd = T,horiz = F)
        dev.off()
    }


    ### legend
    pdf(file=sprintf("%s.integrated.perplexity%d.legend.pdf",out.prefix,perplexity),width=16,height=2)
    par(mar=c(1,1,2,1),cex.lab=1.5,cex.main=1.5)
    plot.new()
    ##idx.leaf <- ceiling(length(colnames(dat.plot))/2)
    #leafCluster.name <- colnames(dat.plot)
    #leafCluster.name <- paste0(leafCluster.name,strrep("x",max(nchar(leafCluster.name))-nchar(leafCluster.name)))
    majorCluster.name <- unique(mapping.leafToMajor.v[colnames(dat.plot.tree)])
    print(majorCluster.name)
    ii <- 0
    mm <- 3
    nn <- 4
    for(yy in head(1-seq(0,1,1/mm),n=mm)){ 
        for(xx in head(seq(0,1,1/nn),n=nn)) {
            ii <- ii+1
            text(x = xx+0.030,y = yy,labels = majorCluster.name[ii],cex=1.6,adj = c(0,1),offset = 4)
            points(xx+0.012,yy-0.050,cex=4,pch=16,xpd=NA,col=majorClusterColor[majorCluster.name[ii]])
        } 
    }
    dev.off()
}

