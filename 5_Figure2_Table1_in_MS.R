# Additional figures and tables in manuscript
library(ggplot2)

# Read in tables with AUC, TSS results
tm <- read.csv("Results/Stats/Mean_AUC_and_TSS_trait_matching.csv",
               stringsAsFactors = FALSE)
tm.f <- read.csv("Results/Stats/Fish_corrected_AUC_and_TSS_trait_match.csv",
                 stringsAsFactors = FALSE)
wb <- read.csv("Results/Stats/Webbuilder_AUC_and_TSS.csv",
               stringsAsFactors = FALSE)
wb.tm <- read.csv("Results/Stats/Wb_tm_AUC_TSS.csv")
eu <- read.csv("Results/Stats/PCA_euclidean_distance_observed.csv")
# put eu in same order as other tables
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

write.csv(tab, "Results/Figures/Table_1.csv",
          row.names = FALSE)



# plots of Adj matrices
# useful food web functions, modified from petchey
source("FoodWeb_Functions.R")

obs <- readRDS("Results/observed_matrices_matched_to_inferred.RDS")
# wb.raw == initial WB inference
wb.raw <- readRDS("Results/WebBuilder_inferred_matrices.RDS")
# wb.n == neutral prune + fish correction
wb.n <- readRDS("Results/WebBuilder_inferred_fish_correction.RDS")
# tm = initial trait matching inference
tm <- readRDS("Results/Initial_trait_matching_inference.RDS")
# tm.niche = niche pruned trait matching inference
tm.niche <- readRDS("Results/Niche_pruned_trait_matching_inference.RDS")
# tm.nn = trait matching, fish correction, niche + neutral pruned
tm.nn <- readRDS("Results/trait_match_niche_neutral_fish.RDS")
# wb x trait match composite
wb.tm <- readRDS("Results/wb_tm_initial.RDS")
# wb.tm.n = wb * tm neutral prune fish correction
wb.tm.n <- readRDS("Results/wb_tm_fish_corrected.RDS")

# dempsters, one plot ####
png(filename = "Results/Figures/fig2.png",
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


