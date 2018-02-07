# figures and tables for MS
library(ggplot2)

# inference stats table
tm <- read.csv("Mean AUC and TSS trait matching.csv")
tm.f <- read.csv("Fish corrected Mean AUC and TSS trait matching.csv")
names(tm.f)[1] <- "inference"
wb <- read.csv("Webbuilder AUC, TSS, fp+fn.csv")
wb.tm <- read.csv("Wb x tm AUC, TSS, fp+fn.csv")
eu <- read.csv("PCA euclidean distance observed.csv")[,-1] # minus "X"
# put in proper order
eu <- eu[,c(1,2,3,4,6,5,7)]

# bind auc, tss, fp, fn tables
tab <- rbind(tm, tm.f, wb,wb.tm)
# add model variable
tab$model <- c("tm", "tm", "tm", "tm", "tm.f", "tm.f", "wb", "wb", "wb", "wb.tm", "wb.tm", "wb.tm")
# remove unneeded rows
tab <- tab[c(-3,-4,-5,-8,-11),]
# add euclidean distance to tab
tab$euclidean <- t(eu)
# reorder columns
tab <- tab[,c(10,1,2,3,4,5,6,7,8,9, 11)]

write.csv(tab, "figs for MS\\post poisot\\All inference summary stats.csv", row.names = F)



# plots of Adj matrices
source("C:\\Users\\Justin\\Documents\\Data\\FW modelling Petchey Github\\ttl-resources-master\\food_web\\FoodWebFunctions.r")

obs <- readRDS("observed matrices matched to inferred.RDS")
# wb.raw == initial WB inference
wb.raw <- readRDS("wb matrices matched to inferred.RDS")
# wb == neutral prune + fish correction
wb.n <- readRDS("wb for PCA.RDS")
tm <- readRDS("Initial trait matching inference.RDS")
tm.niche <- readRDS("Niche pruned trait matching inference.RDS")
tm.nn <- readRDS("tm nn fish for PCA.RDS")
names(tm.nn) <- names(obs) # need to go back to original script and fix names
# wb x trait match, no neutral
wb.tm <- readRDS("wb x tm for PCA.RDS")
# wb.tm wb * tm neutral prune fish correction
wb.tm.n <- readRDS("wb x tm x n for PCA.rds")

# dempsters, one plot ####
png(filename = "figs for MS\\post poisot\\R fig1.png",
    width = 240, height = 320, units = "mm", res =300)
layout(matrix(c(1,2,3,4,0,5,6,0,7,0,0,8), 4, 3, byrow = T))
Plot.matrix2(tm[[7]], sp.pt.ch = 18, point.cex = 1)
title(main = "Trait-matching", cex.lab = 2)
mtext(text = "A ", side = 2, cex = 1.5, las = 1, padj = 0)

Plot.matrix2(tm.niche[[7]], sp.pt.ch = 18, point.cex = 1)
title(main = "Trait-matching + Niche", cex.lab = 2)
mtext(text = "B ", side = 2, cex = 1.5, las = 1, padj = 0)

Plot.matrix2(tm.nn[[7]], sp.pt.ch = 18, point.cex = 1)
title(main = "Trait-matching + Niche + Neutral", cex.lab = 2)
mtext(text = "C ", side = 2, cex = 1.5, las = 1, padj = 0)

Plot.matrix2(wb.raw[[7]], sp.pt.ch = 18, point.cex = 1)
title(main = "WebBuilder", cex.lab = 2)
mtext(text = "D ", side = 2, cex = 1.5, las = 1, padj = 0)

Plot.matrix2(wb.n[[7]], sp.pt.ch = 18, point.cex = 1)
title(main = "WebBuilder + Neutral", cex.lab = 2)
mtext(text = "E ", side = 2, cex = 1.5, las = 1, padj = 0)

Plot.matrix2(wb.tm[[7]], sp.pt.ch = 18, point.cex = 1)
title(main = "WebBuilde:Trait-matching", cex.lab = 2)
mtext(text = "F ", side = 2, cex = 1.5, las = 1, padj = 0)

Plot.matrix2(wb.tm.n[[7]], sp.pt.ch = 18, point.cex = 1)
title(main = "WebBuilde:Trait-matching + Neutral", cex.lab = 2)
mtext(text = "G ", side = 2, cex = 1.5, las = 1, padj = 0)

Plot.matrix2(obs[[7]], sp.pt.ch = 18, point.cex = 1)
title(main = "Empirical", cex.lab = 2)
mtext(text = "H ", side = 2, cex = 1.5, las = 1, padj = 0)
dev.off()


