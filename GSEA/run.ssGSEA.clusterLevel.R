#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-i", "--inFile", type="character", required=TRUE, help="infile")
parser$add_argument("-o", "--outPrefix", type="character", required=TRUE, help="outprefix")
parser$add_argument("-n", "--ncores", type="integer", default=4, help="number of cpu cores used [default %(default)s]")
parser$add_argument("-a", "--width", type="integer", default=8, help="pdf width [default %(default)s]")
parser$add_argument("-b", "--height", type="integer", default=8, help="pdf heights [default %(default)s]")
parser$add_argument("-c", "--nboostrap", type="integer", default=1000, help="number of boostrap [default %(default)s]")
parser$add_argument("-s", "--score", type="double", default=0.3, help="pdf heights [default %(default)s]")
parser$add_argument("-d", "--db", type="character",default="/WPSnew/zhenglt/00.database/MSigDB/msigdb_v6.1_files_to_download_locally/msigdb_v6.1_GMTs/c2.cp.v6.1.symbols.gmt", help="db file, used for external GSEA program [default %(default)s]")
parser$add_argument("-m", "--scoring_scheme", type="character", default="classic",
                    help="scoring_scheme of GSEA: classic,weighted,weighted_p2,weighted_p1.5 [default %(default)s]")
args <- parser$parse_args()
print(args)

i.file <- args$inFile
out.prefix <- args$outPrefix
args.ncores <- args$ncores
args.width <- args$width
args.height <- args$height
args.score <- args$score
args.nboostrap <- args$nboostrap
db.file <- args$db
args.scoring.scheme <- args$scoring_scheme
args.max <- 500

#i.file <- "/WPSnew/zhenglt/work/proj_zl/integrate0924/cross.patient.b20170523/functional/ssGSEA/colonC.P8.cluster.avg.txt"
#out.prefix <- "/WPSnew/zhenglt/work/proj_zl/integrate0924/cross.patient.b20170523/functional/ssGSEA/colonC.P8.cluster.avg.ssGSEA"
#args.ncores <- 8
#args.width <- 8
#args.height <- 8

### don't modify this. entrez id here; but external GSEA geneSymbol
#geneSet.file <- "/DBS/DB_temp/zhangLab/MSigDB/msigdb_v6.0_files_to_download_locally/msigdb_v6.0_GMTs/c2.cp.v6.0.entrez.gmt"
db.dir <- dirname(db.file)
geneSet.file <- sprintf("%s/%s",db.dir,gsub("symbols","entrez",basename(db.file)))
print(db.file)
print(geneSet.file)

dir.create(dirname(out.prefix),showWarnings = F,recursive = T)

suppressPackageStartupMessages(library("GSVA"))
suppressPackageStartupMessages(library("doParallel"))
suppressPackageStartupMessages(library("snow"))

if(file.exists("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.")){
    source("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")
}else{
    source("/WPSnew/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")
}

if(grepl(pattern = "\\.gmt$",geneSet.file,perl = T)){
    print("gmt format")
    geneSet.list <- readGMT(file = geneSet.file)$gSet
}else{
    geneSet.table <- read.delim(geneSet.file,header = T,check.names = F,stringsAsFactors = F,sep = "\t")
    #geneSet.table <- geneSet.table[grepl("_G|Tcm|Tem|TFH|Th17_cells|Th1_cells|Th2_cells|TReg",geneSet.table$cellType,perl=T),]
    geneSet.list <- lapply(unique(geneSet.table$cellType),function(x){
                           as.character(subset(geneSet.table,cellType==x)$geneID)
    })
    names(geneSet.list) <- unique(geneSet.table$cellType)
}

i.table <- read.table(i.file,header = T,check.names = F,stringsAsFactors = F,sep = "\t")
g.GNAME <- i.table[,2]
names(g.GNAME) <- i.table[,1]
rownames(i.table) <- i.table[,1]
i.table <- as.matrix(i.table[,c(-1,-2)])

exp.dat <- i.table[,grepl("^avg",colnames(i.table),perl = T)]

registerDoParallel(cores = args.ncores)

print(str(geneSet.list))

gsva.out <- gsva(exp.dat,geneSet.list,method="ssgsea",
                 parallel.sz=args.ncores,min.sz=5,max.sz=args.max,ssgsea.norm=T)

f.gset <- apply(gsva.out,1,function(x){ sum(x>args.score)>1 })
gsva.out.flt <- gsva.out[f.gset,]

out.df <- data.frame(gsetName=rownames(gsva.out),note="NA",stringsAsFactors = F)
out.df <- cbind(out.df,gsva.out)
write.table(out.df,file = sprintf("%s.score.txt",out.prefix),quote = F,row.names = F,col.names = T,sep = "\t")
write.table(out.df[f.gset,],file = sprintf("%s.score.flt.txt",out.prefix),quote = F,row.names = F,col.names = T,sep = "\t")

#gsva.out.2 <- gsva(exp.dat,geneSet.list,method="gsva",no.bootstraps=args.nboostrap,
#                 parallel.sz=args.ncores,min.sz=3,ssgsea.norm=T)
#if(!is.null(gsva.out.2$bootstrap$p.vals.sign)){
#    out.gsva.p.df <- data.frame(gsetName=rownames(gsva.out.2$bootstrap$p.vals.sign),note="NA",stringsAsFactors = F)
#    out.gsva.p.df <- cbind(out.gsva.p.df,gsva.out.2$bootstrap$p.vals.sign)
#    write.table(out.gsva.p.df,file = sprintf("%s.score.gsva.p.txt",out.prefix),quote = F,row.names = F,col.names = T,sep = "\t")
#    f.gset <- apply(gsva.out.2$bootstrap$p.vals.sign,1,function(x){ sum(x<0.001)>0 })
#
#    plot.dat(-log10(gsva.out.2$bootstrap$p.vals.sign[f.gset,,drop=F]),out.prefix=sprintf("%s.gsva.bootstrap.all",out.prefix),
#             mytitle="",show.number=T,pdf.width=args.width,pdf.height=20)
#
#    for(.mcls in colnames(gsva.out.2$bootstrap$p.vals.sign)){
#        plot.dat(-log10(gsva.out.2$bootstrap$p.vals.sign[gsva.out.2$bootstrap$p.vals.sign[,.mcls]<0.001,]),out.prefix=sprintf("%s.gsva.bootstrap.%s",out.prefix,.mcls),
#                 mytitle=.mcls,pdf.width = 8,pdf.height = 10)
#    }
#
#}
#save(gsva.out, gsva.out.flt,gsva.out.2,file=sprintf("%s.score.RData",out.prefix))
save(gsva.out, gsva.out.flt,file=sprintf("%s.score.RData",out.prefix))

plot.dat <- function(dat,out.prefix,mytitle="",show.number=T,pdf.width=8,pdf.height=8)
{
    suppressPackageStartupMessages(require("gplots"))
    suppressPackageStartupMessages(require("ComplexHeatmap"))
    suppressPackageStartupMessages(require("circlize"))
    suppressPackageStartupMessages(require("gridBase"))
	suppressPackageStartupMessages(require("RColorBrewer"))

    colnames(dat) <- gsub("^avg.","",colnames(dat),perl=T)
    pdf(sprintf("%s.pdf",out.prefix),width=pdf.width,height=pdf.height)
    par(mar=c(4,2,4,4))
    plot.new()
    title(main = mytitle,cex.main=2)
    #legend("topright",legend=names(colSet),fill=colSet,border=colSet,cex=1.5,inset=c(-0.03,0),xpd=T)
    ### Integrating Grid Graphics Output with Base Graphics Output
    vps <- baseViewports()
    pushViewport(vps$inner, vps$figure, vps$plot)
    tmp.var <- pretty((dat),n=8)
    z.lo <- tmp.var[1]
    z.hi <- tmp.var[length(tmp.var)]
    z.step <- tmp.var[2]-tmp.var[1]
    print(z.lo)
    print(z.hi)
    print(z.step)
    m <- ncol(dat)
    n <- nrow(dat)
    ht <- Heatmap(dat, name = "ssGSEA score", col = colorRamp2(seq(z.lo,z.hi,length=100),
                                                               colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)),
                  column_names_gp = gpar(fontsize = 12*28/max(m,32)),row_names_gp = gpar(fontsize = 10*28/max(n,32)),
                  heatmap_legend_param = list(grid_width = unit(0.6, "cm"),
                                              grid_height = unit(0.6, "cm"),
                                              at = seq(z.lo,z.hi,z.step),
                                              title_gp = gpar(fontsize = 12, fontface = "bold"),
                                              #label_gp = gpar(fontsize = 10), color_bar = "continuous"),
                                              label_gp = gpar(fontsize = 12), color_bar = "discrete"),
                  cluster_columns=F,cluster_rows=T)
    ComplexHeatmap::draw(ht, newpage= FALSE)
    dev.off()
}

plot.dat(gsva.out.flt,out.prefix,mytitle="",show.number=T,pdf.width=args.width,pdf.height=args.height)

########## use external GSEA program #############
.o.df <- data.frame(geneID=rownames(exp.dat),geneSymbol=g.GNAME[rownames(exp.dat)],stringsAsFactors = F)
.o.df <- cbind(.o.df,exp.dat)

sink(file = sprintf("%s.external.GSEA.job",out.prefix))
for(.mcls in colnames(exp.dat)){
    .dirname <- dirname(out.prefix)
    .geneSymbol.filename <- sprintf("%s.geneSymbol.%s.rnk",out.prefix,.mcls)
    .geneSymbol.filename <- sprintf("%s/%s",.dirname,gsub("[-_]",".",basename(.geneSymbol.filename)))
    write.table(.o.df[,c("geneSymbol",.mcls)],
                file = .geneSymbol.filename,
                quote = F,sep = "\t",row.names = F,col.names = F)
    .odir <- sprintf("%s.%s.%s.GSEA",out.prefix,.mcls,args.scoring.scheme)
    .odir <- sprintf("%s/%s",.dirname,gsub("[-_]",".",basename(.odir)))
    dir.create(.odir,showWarnings = F,recursive = T)
    cmd.str <- sprintf("java -cp /WPSnew/zhenglt/01.bin/GSEA/gsea2-2.2.4.jar -Xmx2048m xtools.gsea.GseaPreranked -gmx %s -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk %s -scoring_scheme %s -rpt_label %s -include_only_symbols true -make_sets true -plot_top_x 100 -rnd_seed timestamp -set_max %s -set_min 5 -zip_report false -out %s -gui false\n",
                       db.file,
                       .geneSymbol.filename,
                       args.scoring.scheme,
                       sprintf("%s",gsub("^avg.","",.mcls)),
                       args.max,.odir)
    cat(cmd.str)
    #ret <- system(cmd.str)
}
sink()

