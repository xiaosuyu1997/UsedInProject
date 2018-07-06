#!/usr/bin/env Rscript

in.file <- file("stdin")
out.file <- ""

#in.file <- "P0413.phase17.marker.exp.ann.txt"
in.table <- read.table(in.file,header = T,check.names = F,stringsAsFactors = F,sep = "\t")
out.table <- in.table

EPSILON <- 1
f.CD3 <- apply(in.table[,c("CD3D","CD3E","CD3G")],1,function(x){ mean(log2(x+EPSILON)) > log2(10+EPSILON) })
f.CD4 <- f.CD3 & in.table[,"CD4"] > 30 & apply(in.table[,c("CD8A","CD8B")],1,function(x){ mean(log2(x+EPSILON)) < log2(3+EPSILON) })
f.CD8 <- f.CD3 & in.table[,"CD4"] < 3 & apply(in.table[,c("CD8A","CD8B")],1,function(x){ mean(log2(x+EPSILON)) > log2(30+EPSILON) })
f.DP <- f.CD3 & in.table[,"CD4"] > 30 & apply(in.table[,c("CD8A","CD8B")],1,function(x){ mean(log2(x+EPSILON)) > log2(30+EPSILON) })
f.DN <- f.CD3 & in.table[,"CD4"] < 3 & apply(in.table[,c("CD8A","CD8B")],1,function(x){ mean(log2(x+EPSILON)) < log2(3+EPSILON) })
f.other <- f.CD3 & !(f.CD4 | f.CD8 | f.DP | f.DN)

##sum(f.CD4) + sum(f.CD8) + sum(f.DP) + sum(f.DN) + sum(f.other)
out.table$accessory <- "drop"
out.table[f.CD4,"accessory"] <- "CD4"
out.table[f.CD8,"accessory"] <- "CD8"
out.table[f.DP,"accessory"] <- "DP"
out.table[f.DN,"accessory"] <- "DN"
out.table[f.other,"accessory"] <- "other"

write.table(out.table,file = "",row.names = F,sep = "\t",quote = F)

