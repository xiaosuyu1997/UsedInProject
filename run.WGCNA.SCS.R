#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-i", "--inputFile", type="character", required=TRUE, help="input file in these two format: RData contain vstMat.blind obj; TPM in txt format")
parser$add_argument("-d", "--designFile", type="character", help="design matrix file;required if --inputFile is in txt format")
#parser$add_argument("-c", "--clinicalFile", type="character", required=TRUE, help="clinical trait file")
parser$add_argument("-s", "--excludeSample", type="character", help="samples to be excluded, \",\" seperated")
parser$add_argument("-p", "--excludePatient", type="character",help="patients to be excluded, \",\" seperated")
parser$add_argument("-t", "--idType", default="geneID", type="character",help="ID type: geneID, SYMBOL or ENSG [default: %(default)s]")
parser$add_argument("-o", "--outDir", type="character", required=TRUE, help="output directory")
parser$add_argument("-b", "--beta", type="integer", help="soft threshold: if NULL, program only run to get it than exit; if not NULL, use it to construct network",metavar="INT")
parser$add_argument("-n", "--nvar", type="integer", default=5000, help="top [nvar] most variable genes to be used for network construct",metavar="INT")
parser$add_argument("-l", "--log", type="character", default="TRUE", help="whether do log transform for TPM data [default %(default)s]")

sample.exclude <- c()
patient.exclude <- c()
args <- parser$parse_args()
in.file <- args$inputFile
design.file <- args$designFile
arg.log <- as.logical(args$log)
#in.clinical.txt <- args$clinicalFile
out.dir <- args$outDir
NVAR <- args$nvar
id.type <- args$idType
if(!is.null(args$excludeSample)) sample.exclude <- unlist(strsplit(args$excludeSample,split=",",perl=T))
if(!is.null(args$excludePatient)) patient.exclude <- unlist(strsplit(args$excludePatient,split=",",perl=T))

useDEGOnly <- FALSE
loadIndividualTOM <- TRUE

softPower <- args$beta
print(args)
#q()
#in.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/gsnap_out/quantification/P0205.TPM.tab.gz"
#design.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/clustering/marker/P0205.filter.by.marker.design.txt"
#in.file <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/clustering/DESeq2.OUT.phase01-29.multiCore.postTCellMarkerFilter/P0205/DESeq2.RData"
#out.dir <- "/WPS1/zhenglt/work/TCR_chunhong/integrated.20150707/WGCNA/test"
#NVAR <- 1000
#softPower <- NULL
#softPower <- 2

suppressPackageStartupMessages(library("WGCNA"))

options(stringsAsFactors = FALSE)
allowWGCNAThreads()

source("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")
####source("/lustre1/zeminz_pkuhpc/zhenglt/02.pipeline/cancer/lib/myFunc.R")
dir.create(out.dir,recursive = T,showWarnings = F)


## clinical #
##in.clinical.dat <- read.table(in.clinical.txt,header = T,sep = "\t",stringsAsFactors = T,colClasses=NA)
#class(in.clinical.dat$Sample.ID) <- "character"
#in.clinical.dat$T.Stage <- as.factor(in.clinical.dat$T.Stage)
#in.clinical.dat$N.Stage <- as.factor(in.clinical.dat$N.Stage)
#in.clinical.dat <- in.clinical.dat[order(in.clinical.dat$Sample.ID),]
#rownames(in.clinical.dat) <- in.clinical.dat$Sample.ID
#in.clinical.dat <- in.clinical.dat[,c(-1,-2,-4,-8,-14,-16)]
#sel <- !(rownames(in.clinical.dat) %in% patient.exclude)
#in.clinical.dat <- in.clinical.dat[sel,]
#in.clinical.dat

#nSets <- 2
#setLabels = c("Normal Tissue", "Cancer Tissue")
#shortLabels = c("Normal", "Cancer")
nSets <- 1
setLabels = c("SingleCell")
shortLabels = c("SingleCell")

multiExpr <- vector(mode = "list", length = nSets)

## read in data ###
LOW.THRESHOLD <- 0
datFormat <- ""
if(grepl(".RData$",in.file,perl=T))
{
    ### VST 
    suppressPackageStartupMessages(library("R.utils"))
    lenv <- loadToEnv(in.file)
    vstMat.blind <- lenv[["vstMat.blind"]]
    multiExpr[[1]]=list(data=(as.data.frame(t(vstMat.blind))))
    LOW.THRESHOLD <- range(multiExpr[[set]]$data)[1]
    ###multiExpr[[2]]=list(data=(as.data.frame(t(vstMat[,grepl("T$",colnames(vstMat),perl=T)]))))
    datFormat <- "RData"
}else if(grepl(".tab.gz$",in.file,perl=T))
{
    ### TPM
    in.table <- read.table(in.file,row.names="geneID",header = T,sep = "\t",check.names = F)
    in.table <- in.table[,-1]
    myDesign <- read.table(design.file,header=T,row.names="sample",check.names=F,colClasses=c("factor","character","factor","factor"))
    in.table <- in.table[,rownames(myDesign)]
    multiExpr[[1]]=list(data=t(in.table))
    LOW.THRESHOLD <- 1
    datFormat <- "tab"
}

for(set in 1:nSets)
{
    f.geneExpressed <- apply(multiExpr[[set]]$data,2,function(x){ sum(x>LOW.THRESHOLD)/length(x) > 0.1 })
	multiExpr[[set]]$data = multiExpr[[set]]$data[, f.geneExpressed]
    if(datFormat=="tab" && arg.log) {
        multiExpr[[set]]$data <- log2(1+multiExpr[[set]]$data)
    }
}
###if(useDEGOnly)
###{
###	for(set in 1:nSets)
###		multiExpr[[set]]$data = multiExpr[[set]]$data[, rownames(resSigStrict)]
###}

###for(set in 1:nSets)
###{
###	sel <- !(rownames(multiExpr[[set]]$data) %in% sample.exclude)
###	multiExpr[[set]]$data <- multiExpr[[set]]$data[sel,]
###}

if(!is.na(NVAR) && ncol(multiExpr[[1]]$data) > NVAR)
{
	col.var <- apply(multiExpr[[1]]$data,2,var)
	multiExpr[[1]]$data <- multiExpr[[1]]$data[,order(col.var,decreasing = T)]
	multiExpr[[1]]$data <- multiExpr[[1]]$data[,1:NVAR]
}

#for(set in 1:nSets)
#{
	###multiExpr[[set]]$data = scale(multiExpr[[set]]$data)
	######rownames(multiExpr[[set]]$data) <- gsub("[NT]$","",rownames(multiExpr[[set]]$data),perl = T)
#}

exprSize = checkSets(multiExpr)
exprSize
gsg = goodSamplesGenesMS(multiExpr, verbose = 3)
gsg$allOK

if (!gsg$allOK)
{
	# Print information about the removed genes:
	if (sum(!gsg$goodGenes) > 0)
		printFlush(paste("Removing genes:", paste(names(multiExpr[[1]]$data)[!gsg$goodGenes],collapse = ", ")))
	for (set in 1:exprSize$nSets)
	{
		if (sum(!gsg$goodSamples[[set]]))
		printFlush(paste("In set", setLabels[set], "removing samples", 
				 paste(rownames(multiExpr[[set]]$data)[!gsg$goodSamples[[set]]], collapse = ", ")))
		# Remove the offending genes and samples
		multiExpr[[set]]$data = multiExpr[[set]]$data[gsg$goodSamples[[set]], gsg$goodGenes];
	}
	# Update exprSize
	exprSize = checkSets(multiExpr)
}
exprSize
print(head(multiExpr[[1]]$data[,1:10]))

sampleTrees = list()
for (set in 1:nSets)
{
	sampleTrees[[set]] = hclust(dist(multiExpr[[set]]$data), method = "average")
}
pdf(file = paste0(out.dir,"/SampleClustering.pdf"), width = 12, height = 8)
par(mfrow=c(nSets,1))
par(mar = c(2, 6, 4, 2))
for (set in 1:nSets)
	plot(sampleTrees[[set]], main = paste("Sample clustering on all genes in", setLabels[set]), xlab="", sub="", cex = 30/exprSize$nSamples)
dev.off()

qcExpr <- t(multiExpr[[1]]$data)
cal.sample.cor <- function(qcExpr,out.prefix)
{
	qcExpr.cor <- cor(qcExpr)
	diag(qcExpr.cor) <- NA
	qcExpr.cor.summary <- t(apply(qcExpr.cor,1,function(x){ c(mean(x,na.rm = T),median(x,na.rm = T),min(x,na.rm = T),max(x,na.rm = T)) } ))
	colnames(qcExpr.cor.summary) <- c("mean","mediaan","min","max")
	qcExpr.cor.summary <- data.frame(sampleID=rownames(qcExpr.cor.summary),qcExpr.cor.summary)
	#print(qcExpr.cor.summary)
	write.table(qcExpr.cor.summary,paste0(out.prefix,".sample.cor.txt"),row.names = F,sep="\t",quote = F)
	pdf(paste0(out.prefix,".sample.cor.pdf"),width = 10,height = 8)
	plot(hclust(as.dist(1-qcExpr.cor),"average"),xlab="",sub = "",cex=30/ncol(qcExpr))
	dev.off()
}
cal.sample.cor(qcExpr,paste0(out.dir,"/RC-VST"))
rm(qcExpr)

###### network construction ######

#load(paste0(out.dir,"/network.RData"))

my.pickBeta <- function(multiExpr,nSets,out.dir,setLabels)
{
	# Choose a set of soft-thresholding powers
	powers = c(seq(1,10,by=1), seq(12,40, by=2))
	# Initialize a list to hold the results of scale-free analysis
	powerTables = vector(mode = "list", length = nSets);
	# Call the network topology analysis function for each set in turn
	for (set in 1:nSets)
	  powerTables[[set]] = list(data = pickSoftThreshold(multiExpr[[set]]$data, powerVector=powers, verbose = 2)[[2]])
	collectGarbage();
	# Plot the results:
	colors = c("black", "red")
	# Will plot these columns of the returned scale free analysis tables
	plotCols = c(2,5,6,7)
	colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity", "Max connectivity")
	# Get the minima and maxima of the plotted points
	ylim = matrix(NA, nrow = 2, ncol = 4);
	for (set in 1:nSets)
	{
	  for (col in 1:length(plotCols))
	  {
	    ylim[1, col] = min(ylim[1, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
	    ylim[2, col] = max(ylim[2, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
	  }
	}
	# Plot the quantities in the chosen columns vs. the soft thresholding power
	#sizeGrWindow(8, 6)
	pdf(file = paste0(out.dir,"/PickSoftThreshold.pdf"), width = 12, height = 12)
	par(mfcol = c(2,2));
	par(mar = c(4.2, 4.2 , 2.2, 0.5))
	cex1 = 0.7;
	for (col in 1:length(plotCols)) for (set in 1:nSets)
	{
	  if (set==1)
	  {
	    plot(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
		 xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = ylim[, col],
		 main = colNames[col]);
	    addGrid();
	  }
	  if (col==1)
	  {
	    text(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
		 labels=powers,cex=cex1,col=colors[set]);
	  } else
	    text(powerTables[[set]]$data[,1], powerTables[[set]]$data[,plotCols[col]],
		 labels=powers,cex=cex1,col=colors[set]);
	  if (col==1)
	  {
	    legend("bottomright", legend = setLabels, col = colors, pch = 20) ;
	  } else
	    legend("topright", legend = setLabels, col = colors, pch = 20) ;
	}
	dev.off()
}

if(is.null(softPower))
{
	my.pickBeta(multiExpr,nSets,out.dir,setLabels)
	save.image(paste0(out.dir,"/network.RData"))
	q()
}

net.SCS = blockwiseModules(
	multiExpr[[1]]$data, power = softPower,
	TOMType = "unsigned", minModuleSize = 20,
	reassignThreshold = 0, mergeCutHeight = 0.25,
	numericLabels = TRUE, pamRespectsDendro = FALSE,
  	maxBlockSize = 30000,
	loadTOM = loadIndividualTOM,
	saveTOMs = TRUE,
	saveTOMFileBase = paste0(out.dir,"/SCS.TOM"),
	verbose = 3)

load(paste0(out.dir,"/SCS.TOM-block.1.RData"))
TOM.SCS <- as.matrix(TOM)
rm(TOM)
collectGarbage()
save(TOM.SCS,file=paste0(out.dir,"/TOM.RData"))

#save.image(paste0(out.dir,"/network.RData"))
#### plot dendrogram  ###
my.plotDendroAndColors  <- function(out.png,tree,moduleColors)
{
	png(file = out.png, width = 1000, height = 600)
	par(cex.lab=1.5,cex.main=1.5,cex.axis=1.5)
	plotDendroAndColors(tree, moduleColors,
			    "Module colors",
			    dendroLabels = FALSE, hang = 0.03,
			    addGuide = TRUE, guideHang = 0.05,
			    cex.colorLabels = 1.2, cex.dendroLabels = 1.8, cex.rowText = 1.5,
			    main = "Gene dendrogram and module colors")
	dev.off()
}
my.plotDendroAndColors(paste0(out.dir,"/SCS.DendroAndColors.png"),net.SCS$dendrograms[[1]],labels2colors(net.SCS$colors))
###cal.sample.cor(TPM.normal,TPM.tumor,paste0(out.dir,"/TPM"))

#### ME & clinical trait correlation ###

###my.plotMECorHeatmap <- function(MEs, datTraits, out.pdf, nSamples)
###{
###	moduleTraitCor <- c()
###	moduleTraitPvalue <- c()
###
###	for(i in 1:ncol(datTraits))
###	{
###		vCor <- c()
###		vPvalue <- c()
###		if( class(datTraits[,i])=="factor" )
###		{
###			vCor <- rep(NA,ncol(MEs))
###			vPvalue <- apply(MEs,2,function(x){ summary(aov(x~datTraits[,i]))[[1]][[1,"Pr(>F)"]]  } )
###		}else
###		{
###			vCor = cor(MEs, datTraits[,i], use = "p")
###			vPvalue = corPvalueStudent(vCor, nSamples)
###		}
###		moduleTraitCor <- cbind(moduleTraitCor,vCor)
###		moduleTraitPvalue <- cbind(moduleTraitPvalue,vPvalue)
###	}
###	#moduleTraitCor <- as.data.frame(moduleTraitCor)
###	#moduleTraitPvalue <- as.data.frame(moduleTraitPvalue)
###	moduleTraitCor <- as.matrix(moduleTraitCor)
###	moduleTraitPvalue <- as.matrix(moduleTraitPvalue)
###	colnames(moduleTraitCor) <- colnames(datTraits)
###	colnames(moduleTraitPvalue) <- colnames(datTraits)
###
###
###	write.table(data.frame(sample=rownames(datTraits),datTraits,MEs),paste0(out.pdf,".dat.txt"),quote = F,sep="\t",col.names = T,row.names = F)
###	write.table(data.frame(module=rownames(moduleTraitCor),moduleTraitCor),paste0(out.pdf,".cor.txt"),quote = F,sep="\t",col.names = T,row.names = F)
###	write.table(data.frame(module=rownames(moduleTraitPvalue),moduleTraitPvalue),paste0(out.pdf,".p.txt"),quote = F,sep="\t",col.names = T,row.names = F)
###	##moduleTraitCor = cor(MEs, datTraits, use = "p")
###	##moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
###	# Will display correlations and their p-values
###	textMatrix = paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 2), ")", sep = "")
###	dim(textMatrix) = dim(moduleTraitCor)
###	pdf(file = out.pdf, width = 8, height = 8)
###	par(mar = c(8, 8.5, 3, 3));
###	# Display the correlation values within a heatmap plot
###	labeledHeatmap(Matrix = moduleTraitCor, 
###		       xLabels = names(datTraits), yLabels = names(MEs), ySymbols = names(MEs), 
###		       colorLabels = FALSE, colors = blueWhiteRed(50), 
###		       textMatrix = textMatrix, setStdMargins = FALSE, 
###		       cex.lab.y =0.7, textAdj = c(0.5, 0.5),
###		       cex.text = 0.3, zlim = c(-1,1), main = paste("Module-trait relationships"))
###	labeledHeatmap(Matrix = moduleTraitPvalue, 
###		       xLabels = names(datTraits), yLabels = names(MEs), ySymbols = names(MEs), 
###		       colorLabels = FALSE, colors = colorRampPalette(c("red", "white"))(100), 
###		       textMatrix = round(moduleTraitPvalue, 4), setStdMargins = FALSE, 
###		       cex.lab.y =0.7, textAdj = c(0.5, 0.5),
###		       cex.text = 0.3, zlim = c(0,1), main = paste("Module-trait relationships"))
###
###	dev.off()
###}
###my.plotMECorHeatmap(net.tumor$MEs,in.clinical.dat,paste0(out.dir,"/tumor.ME.traits.cor.pdf"),exprSize$nSamples)
###my.plotMECorHeatmap(net.normal$MEs,in.clinical.dat,paste0(out.dir,"/normal.ME.traits.cor.pdf"),exprSize$nSamples)


### gene & clinical trait correlation ###
###v.survival <- as.data.frame(in.clinical.dat$Survival)
###names(v.survival) <- "Survival"

my.output.geneInfo <- function(net,datExpr,tom,nSamples,out.txt)
{
	MEs <- net$MEs
	moduleLabels <- net$colors
	moduleColors <- labels2colors(net$colors)

	modNames = substring(names(MEs), 3)
	geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
	MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
	names(geneModuleMembership) = paste("MM", modNames, sep="")
	names(MMPvalue) = paste("p.MM", modNames, sep="")
	#geneTraitSignificance = as.data.frame(cor(datExpr, v, use = "p"))
	#GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
	#names(geneTraitSignificance) = paste("GS.", names(v), sep="")
	#names(GSPvalue) = paste("p.GS.", names(v), sep="")
	## identify hub genes
	moduleInfo <- as.data.frame(t(apply(geneModuleMembership,2,function(x,n){ h<-n[which.max(abs(x))]; c(h,max(abs(x)))  }, rownames(geneModuleMembership) )))
	colnames(moduleInfo) <- c("hubGeneID","hubKME")
	moduleInfo$module <- rownames(moduleInfo)
	moduleInfo$hubGeneSymbol <- entrezToXXX(moduleInfo$hubGeneID)
	moduleInfo$Size <- table(moduleLabels)[substring(moduleInfo$module,3)]
	write.table(moduleInfo,paste0(out.txt,".moduleInfo.txt"),sep="\t",row.names = F,col.names = T,quote = F)
	
	## clusterCof and connectivity
	cCoef <- clusterCoef(tom)
	k<-intramodularConnectivity(adjMat=tom,colors=moduleLabels)
	# Create the starting data frame
	geneInfo0 = data.frame(geneID = colnames(datExpr),
			       geneSymbol = entrezToXXX(colnames(datExpr)),
			       moduleLabel = moduleLabels,
			       moduleColor = moduleColors,
			       ###geneTraitSignificance,
			       ###GSPvalue,
			       highKME=apply(geneModuleMembership,1,function(x){max(abs(x))>0.9} ),
			       clusterCoef = cCoef
			       )
	geneInfo0  <- cbind(geneInfo0,k)
	## Order modules by their significance for weight
	#modOrder = order(-abs(cor(MEs, v, use = "p")))
	## Add module membership information in the chosen order
	#for (mod in 1:ncol(geneModuleMembership))
	#{
	#  oldNames = names(geneInfo0)
	#  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
	#			 MMPvalue[, modOrder[mod]])
	#  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""), paste("p.MM.", modNames[modOrder[mod]], sep=""))
	#}
	## Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
	###geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.Survival));
	###geneInfo = geneInfo0[geneOrder, ]

	write.table(geneInfo0, file = out.txt,sep = "\t",col.names = T,row.names = F,quote = F)
}

my.output.geneInfo(net.SCS,multiExpr[[1]]$data,TOM.SCS,exprSize$nSamples[1],paste0(out.dir,"/SCS.geneInfo.txt"))

### modules' GO ###

my.GO.enrichment <- function(net,datExpr,out.txt,IDType="geneID")
{
	moduleLabels <- net$colors
    if(IDType=="geneID"){
        allLLIDs <- colnames(datExpr)
    }else {
        allLLIDs <- XXXToEntrez(colnames(datExpr))
    }
	
	GOenr = GOenrichmentAnalysis(moduleLabels, allLLIDs, organism = "human", nBestP = 50)
	tab = GOenr$bestPTerms[[4]]$enrichment
	names(tab)
	write.table(tab, file = out.txt, sep = ",", quote = TRUE, row.names = FALSE)
}
my.GO.enrichment(net.SCS,multiExpr[[1]]$data,paste0(out.dir,"/SCS.module.GO.csv"),IDType=id.type)

### plot eigengene networks ###

#my.plotEigengeneNetworks <- function(net,datExpr,out.prefix,v.survival)
my.plotEigengeneNetworks <- function(net,tom,out.prefix)
{
	moduleColors <- labels2colors(net$colors)
	MEs <- net$MEs
	geneTree  <- net$dendrograms[[1]]
	out.png <- paste0(out.prefix,".TOMplot.png")
	out.pdf <- paste0(out.prefix,".plotEigengeneNetworks.pdf")

	#dissTOM = 1-tom
	#plotTOM = dissTOM^7
	#diag(plotTOM) = NA;
	#png(file = out.png, width = 1000, height = 1000)
	#TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")
	#dev.off()

	#MET = orderMEs(cbind(MEs, v.survival))
	MET = orderMEs(cbind(MEs))
	# Plot the relationships among the eigengenes and the trait
	pdf(file = out.pdf, width = 10, height = 10)
	par(cex = 1.0)
	plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.7, xLabelsAngle = 90)
	plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0), plotHeatmaps = FALSE)
	plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2), plotDendrograms = FALSE, xLabelsAngle = 90)
	dev.off()

}
my.plotEigengeneNetworks(net.SCS,TOM.SCS,paste0(out.dir,"/SCS"))

q()
### export networks ###

exportOneTOM <- function(mTOM,out.prefix,mGenes,mColors)
{
	cyt = exportNetworkToCytoscape(mTOM,
	       edgeFile = paste0(out.prefix,".CytoscapeInput.edges.txt"),
	       nodeFile = paste0(out.prefix,".CytoscapeInput.nodes.txt"),
	       weighted = TRUE,
	       threshold = 0.01,
	       nodeNames = mGenes,
	       altNodeNames = entrezToXXX(mGenes),
	       nodeAttr = mColors)
	vis = exportNetworkToVisANT(mTOM,
	      file = paste0(out.prefix,".VisANT.txt"),
	      weighted = TRUE,
	      threshold = 0.01,
	      maxNConnections=200,
	      probeToGene = data.frame(mGenes,entrezToXXX(mGenes)) )
}
my.exportModule2Cytoscape <- function(net,TOM,datExpr,modules,sampleID,out.dir)
{
	moduleColors <- labels2colors(net$colors)
	moduleLabels <- net$colors
	MEs <- net$MEs
	moduleLabels <- net$colors
	moduleColors <- labels2colors(net$colors)
	modNames = substring(names(MEs), 3)
	geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
	genes = colnames(datExpr)

	dir.create(out.dir,recursive = T,showWarnings = F)
	### all hub gene of the network
	## identify hub genes
	hub <- apply(geneModuleMembership,2,function(x){ which.max(abs(x))  } )
	hubGene <- genes[hub]
	hubTOM <- TOM[hub,hub]
	dimnames(hubTOM) <- list(hubGene,hubGene)
	exportOneTOM(hubTOM,paste0(out.dir,"/",sampleID,"hub"),hubGene,moduleColors[hub])

        highKME <- apply(geneModuleMembership,1,function(x){max(abs(x))>0.9} )
	highKMEGene <- genes[highKME]
	highKMETOM <- TOM[highKME,highKME]
	dimnames(highKMETOM) <- list(highKMEGene,highKMEGene)
	exportOneTOM(highKMETOM,paste0(out.dir,"/",sampleID,"highKMT"),highKMEGene,moduleColors[highKME])

	### per module
	for(i in 1:length(modules))
	{
		###inModule = moduleLabels %in% modules
		inModule = (moduleLabels == modules[i])
		modGenes = genes[inModule]
		modTOM = TOM[inModule, inModule]
		dimnames(modTOM) = list(modGenes, modGenes)
		exportOneTOM(modTOM,paste0(out.dir,"/",sampleID,".Modules",modules[i]),modGenes,moduleColors[inModule])
	}

}
my.exportModule2Cytoscape(net.normal,TOM.normal,multiExpr[[1]]$data,names(table(net.normal$colors)),"normal",paste0(out.dir,"/Vis"))
my.exportModule2Cytoscape(net.tumor,TOM.tumor,multiExpr[[2]]$data,names(table(net.tumor$colors)),"tumor",paste0(out.dir,"/Vis"))

## plot preservation  ###
my.modulePreservationAnalysis <- function(TOM1,TOM2,net1,net2,out.prefix,name1,name2)
{
	multiTOM <- vector(mode = "list", length = 2)
	multiTOM[[1]]=list(data=TOM1)
	multiTOM[[2]]=list(data=TOM2)
	diag(multiTOM[[1]]$data) <- 1
	diag(multiTOM[[2]]$data) <- 1

	multiColor <- vector(mode = "list", length = 2)
	multiColor[[1]]=net1$colors
	multiColor[[2]]=net2$colors

	preservationRes <- modulePreservation(multiTOM, multiColor, dataIsExpr = FALSE,
					      interpolationPlotFile = paste0(out.prefix,".modulePreservationInterpolationPlots.pdf"),
					      permutedStatisticsFile = paste0(out.prefix, ".permutedStats-Modules.RData")
					      )

	### Zsummary and medianRank
	df.x <- data.frame(module=rownames(preservationRes$preservation$observed[[1]][[2]]),preservationRes$preservation$observed[[1]][[2]])
	df.y <- data.frame(module=rownames(preservationRes$preservation$Z[[1]][[2]]),preservationRes$preservation$Z[[1]][[2]])
	df.merge <- merge(df.x,df.y,by="module")
	df.merge <- df.merge[with(df.merge, module != "0" & module !="0.1"),]
	df.merge <- df.merge[order(as.numeric(as.character(df.merge$module))),]
	write.table(df.merge,paste0(out.prefix,".modulePreservation.presObservedAndZ.txt"),row.names = F,quote = F,sep = "\t")
	pdf(paste0(out.prefix,".modulePreservation.pdf"),width = 10,height = 8)
	oldPar <- par(mfcol=c(1,2))
	plot(df.merge$moduleSize.x,df.merge$Zsummary.pres,xlab="Module Size",ylab="Zsummary",col=labels2colors(df.merge$module))
	text(df.merge$moduleSize.x,df.merge$Zsummary.pres,labels=df.merge$module,cex=0.5,pos=3,offset = 0.3)
	abline(h = 2,lty=2)
	abline(h = 10,lty=2)
	plot(df.merge$moduleSize.x,df.merge$medianRank.pres,xlab="Module Size",ylab="Median Rank",col=labels2colors(df.merge$module))
	text(df.merge$moduleSize.x,df.merge$medianRank.pres,labels=df.merge$module,cex=0.5,pos=3,offset = 0.3)
	par(oldPar)
	title(main = "Module Preservation")
	dev.off()
	### correspondence of assignment
	pdf(paste0(out.prefix,".assignment.pdf"),width=10,height=8)
	par(mar=c(8,8,4,2))
	labeledHeatmap(Matrix = -log10(preservationRes$accuracy$observedFisherPvalues[[1]][[2]]),
		xLabels = paste0("  ",labels2colors(as.numeric(names(table(net1$colors))))),
		yLabels = paste0("  ",labels2colors(as.numeric(names(table(net2$colors))))),
		colorLabels = TRUE,
		xSymbols = paste(name1, names(table(net1$colors)), "(", table(net1$colors),")", sep=""),
		ySymbols = paste(name2, names(table(net2$colors)), "(", table(net2$colors),")", sep=""),
		textMatrix = preservationRes$accuracy$observedCounts[[1]][[2]],
		colors = blueWhiteRed(100)[50:100],
		main = "Correspondence of two networks' assignment",
		cex.text = 1.0, cex.lab = 1.0, setStdMargins = FALSE)
	dev.off()

	preservationRes
}
moduleNormal.pres <- my.modulePreservationAnalysis(TOM.normal,TOM.tumor,net.normal,net.tumor,paste0(out.dir,"/normalModule"),"N","T")
moduleTumor.pres <- my.modulePreservationAnalysis(TOM.tumor,TOM.normal,net.tumor,net.normal,paste0(out.dir,"/tumorModule"),"T","N")

### specificity links
my.exportMatrix <- function(adjMat,mat1,mat2,nodeNames,altNodeNames,threshold)
{
	nRow = nrow(adjMat)
	adjDst = as.dist(adjMat)
	mat1Dst = as.dist(mat1)
	mat2Dst = as.dist(mat2)
	rowMat = matrix(c(1:nRow), nRow, nRow, byrow = TRUE)
	colMat = matrix(c(1:nRow), nRow, nRow)
	dstRows = as.dist(rowMat)
	dstCols = as.dist(colMat)
	#edges = rep(TRUE,length(adjDst))
	edges = abs(adjDst) > threshold
	nEdges = sum(edges)
	edgeData = data.frame(fromNode = nodeNames[dstRows[edges]], 
			     	toNode = nodeNames[dstCols[edges]], 
				specificity = adjDst[edges],
				weight1 = mat1Dst[edges],
				weight2 = mat2Dst[edges],
				fromAltName=altNodeNames[dstRows[edges]],
				toAltName=altNodeNames[dstCols[edges]]
			)
	print(paste0("nEdges:",nEdges))
	edgeData

}

my.specificityLink <- function(tom1,tom2,net1,net2,out.prefix)
{
	allmodules <- names(table(net1$colors))
	for(module in 2:length(allmodules))
	{
		sel <- net1$colors %in% allmodules[module]
		moduleTOM1 <- tom1[sel,sel]
		moduleTOM2 <- tom2[sel,sel]
		n <- ncol(moduleTOM1)
		gnames <- rownames(moduleTOM1)
		sel.tri <- upper.tri(moduleTOM1)
		mean1 <- mean(moduleTOM1[upper.tri(moduleTOM1)],na.rm = TRUE)
		mean2 <- mean(moduleTOM2[upper.tri(moduleTOM2)],na.rm = TRUE)
		
		out <- (moduleTOM1/mean1)/((moduleTOM1/mean1)+(moduleTOM2/mean2))
		out.df <- my.exportMatrix(out,moduleTOM1,moduleTOM2,gnames,entrezToXXX(gnames),0.8)
		
        write.table(out.df,paste0(out.prefix,".module",allmodules[module],".cytoscape.edges.txt"),sep="\t",row.names = F,quote = F)
	}
}

my.specificityLink(TOM.normal,TOM.tumor,net.normal,net.tumor,paste0(out.dir,"/specificityLink.normal"))
my.specificityLink(TOM.tumor,TOM.normal,net.tumor,net.normal,paste0(out.dir,"/specificityLink.tumor"))


######lname <- load(paste0(out.dir,"/network.RData"))
save.image(paste0(out.dir,"/network.RData"))


