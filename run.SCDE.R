#!/usr/env Rscript

library("scde")
# load example dataset
data(es.mef.small);
# factor determining cell types
sg <- factor(gsub("(MEF|ESC).*","\\1",colnames(es.mef.small)),levels=c("ESC","MEF"));
names(sg) <- colnames(es.mef.small); # the group factor should be named accordingly
table(sg);

# clean up the dataset
cd <- es.mef.small;
# omit genes that are never detected
cd <- cd[rowSums(cd)>0,];
# omit cells with very poor coverage
cd <- cd[,colSums(cd)>1e4]; 

# number of local process cores to use in processing
n.cores <- 4;
# calculate models
o.ifm <- scde.error.models(counts=cd,groups=sg,n.cores=n.cores,threshold.segmentation=T,save.crossfit.plots=F,save.model.plots=F,verbose=1);

head(o.ifm)

# filter out cells that don't show positive correlation with
# the expected expression magnitudes (very poor fits)
valid.cells <- o.ifm$corr.a >0;
table(valid.cells)
o.ifm <- o.ifm[valid.cells,];

#estimate gene expression prior
o.prior <- scde.expression.prior(models=o.ifm,counts=cd,length.out=400,show.plot=F)

# define two groups of cells
groups <- factor(gsub("(MEF|ESC).*","\\1",rownames(o.ifm)),levels=c("ESC","MEF")); 
names(groups) <- row.names(o.ifm);
# run differential expression tests on all genes.
ediff <- scde.expression.difference(o.ifm,cd,o.prior,groups=groups,n.randomizations=100,n.cores=n.cores,verbose=1)

# top upregulated genes (tail would show top downregulated ones)
head(ediff[order(ediff$Z,decreasing=T),])

# write out a table with all the results, showing most significantly different genes (in both directions) on top
write.table(ediff[order(abs(ediff$Z),decreasing=T),],file="results.txt",row.names=T,col.names=T,sep="\t",quote=F)

scde.test.gene.expression.difference("Tdh",models=o.ifm,counts=cd,prior=o.prior)

# view differential expression results in a browser (requires Rook and rjson packages to be installed)
scde.browse.diffexp(ediff,o.ifm,cd,o.prior,groups,geneLookupURL="http://www.informatics.jax.org/searchtool/Search.do?query={0}")

batch <- as.factor(ifelse(rbinom(nrow(o.ifm),1,0.5)==1,"batch1","batch2"))
# check the interaction between batches and cell types (shouldn't be any)
table(groups,batch)

# test the Tdh gene again
scde.test.gene.expression.difference("Tdh",models=o.ifm,counts=cd,prior=o.prior,batch=batch)

# test for all of the genes
ediff.batch <- scde.expression.difference(o.ifm,cd,o.prior,groups=groups,batch=batch,n.randomizations=100,n.cores=n.cores,return.posteriors=T,verbose=1)

# browse results
scde.browse.diffexp(ediff.batch,o.ifm,cd,o.prior,groups,batch=batch,geneLookupURL="http://www.informatics.jax.org/searchtool/Search.do?query={0}")

# calcualte joint posterior for ESCs (set return.individual.posterior.modes=T if you need p.modes)
jp <- scde.posteriors(models=o.ifm[grep("ESC",rownames(o.ifm)),],cd,o.prior,n.cores=n.cores)

# get expression magntiude estimates
o.fpm <- scde.expression.magnitude(o.ifm,counts=cd);

# get failure probabilities on the expresison range
o.fail.curves <- scde.failure.probability(o.ifm,magnitudes=log((10^o.prior$x)-1))
par(mfrow=c(1,1),mar = c(3.5,3.5,0.5,0.5), mgp = c(2.0,0.65,0), cex = 1);
plot(c(),c(),xlim=range(o.prior$x),ylim=c(0,1),xlab="expression magnitude (log10)",ylab="drop-out probability")
invisible(apply(o.fail.curves[,grep("ES",colnames(o.fail.curves))],2,function(y) lines(x=o.prior$x,y=y,col="orange")))
invisible(apply(o.fail.curves[,grep("MEF",colnames(o.fail.curves))],2,function(y) lines(x=o.prior$x,y=y,col="dodgerblue")))

# get failure probabilities on the expresison range
o.fail.curves <- scde.failure.probability(o.ifm,magnitudes=log((10^o.prior$x)-1))
# get self-fail probabilities (at a given observed count)
p.self.fail <- scde.failure.probability(models=o.ifm,counts=cd)


p.self.fail <- scde.failure.probability(models=o.ifm,counts=cd)
# simulate drop-outs
# note: using 10 sampling rounds for illustration here. ~500 or more should be used.
n.simulations <- 10; k <- 0.9;
cell.names <- colnames(cd); names(cell.names) <- cell.names;
dl <- mclapply(1:n.simulations,function(i) {
                 scd1 <- do.call(cbind,lapply(cell.names,function(nam) {
                                                      x <- cd[,nam];
                                                      # replace predicted drop outs with NAs
                                                      x[!as.logical(rbinom(length(x),1,1-p.self.fail[,nam]*k))] <- NA;
                                                      x;
                                                          }))
                   rownames(scd1) <- rownames(cd); 
                   # calculate correlation on the complete observation pairs
                   cor(log10(scd1+1),use="pairwise.complete.obs");
},mc.cores=n.cores)
# calculate average distance across sampling rounds
direct.dist <- as.dist(1-Reduce("+",dl)/length(dl))

# load boot package for the weighted correlation implementation
require(boot)
k <- 0.95;
reciprocal.dist <- as.dist(1-do.call(rbind,mclapply(cell.names,function(nam1) {
                                                      unlist(lapply(cell.names,function(nam2) {
                                                                        # reciprocal probabilities
                                                                        f1 <- scde.failure.probability(models=o.ifm[nam1,,drop=F],magnitudes=o.fpm[,nam2])
                                                                            f2 <- scde.failure.probability(models=o.ifm[nam2,,drop=F],magnitudes=o.fpm[,nam1])
                                                                            # weight factor
                                                                            pnf <- sqrt((1-f1)*(1-f2))*k +(1-k); 
                                                                                boot::corr(log10(cbind(cd[,nam1],cd[,nam2])+1),w=pnf)
                                                                                }))
                                                          },mc.cores=n.cores)),upper=F)

# reclculate posteriors with the individual posterior modes 
jp <- scde.posteriors(models=o.ifm,cd,o.prior,return.individual.posterior.modes=T,n.cores=n.cores)
# find joint posterior modes for each gene - a measure of MLE of group-average expression
jp$jp.modes <- log(as.numeric(colnames(jp$jp)))[max.col(jp$jp)]
p.mode.fail <- scde.failure.probability(models=o.ifm,magnitudes=jp$jp.modes)
# weight matrix
matw <- 1-sqrt(p.self.fail*sqrt(p.self.fail*p.mode.fail))
# magnitude matrix (using individual posterior modes here)
mat <- log10(exp(jp$modes)+1);
# weighted distance
mode.fail.dist <- as.dist(1-do.call(rbind,mclapply(cell.names,function(nam1) {
                                                     unlist(lapply(cell.names,function(nam2) {
                                                                       corr(cbind(mat[,nam1],mat[,nam2]),w=sqrt(sqrt(matw[,nam1]*matw[,nam2])))
                                                                         }))
                                                          },mc.cores=n.cores)),upper=F);





