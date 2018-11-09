# PCA analysis of food web measures
library(plyr)
library(dplyr)
library(ggplot2)
library(purrr)
library(vegan)
# useful food web functions, modified from petchey
source("Functions/FoodWeb_Functions.R")
# functions written for this manuscript
source("Functions/Inference_MS_functions.R")

obs <- readRDS("Results/observed_matrices_matched_to_inferred.RDS")
# wb.raw == initial WB inference
wb.raw <- readRDS("Results/WebBuilder_inferred_matrices.RDS")
# wb.n == neutral prune + fish correction
wb.n <- readRDS("Results/WebBuilder_inferred_fish_correction.RDS")
# tm = initial trait-matching inference
tm <- readRDS("Results/Initial_trait_matching_inference.RDS")
# tm.niche = trait-matching inference + niche pruned
tm.niche <- readRDS("Results/Niche_pruned_trait_matching_inference.RDS")
# tm.nn == trait-matching + neutral and niche pruned, after fish correction
tm.nn <- readRDS("Results/trait_match_niche_neutral_fish.RDS")
# wb.tm = initial webbuilder trait-match composite
wb.tm <- readRDS("Results/wb_tm_initial.RDS")
# wb.tm.n = wb tm composite, neutral pruned after fish correction
wb.tm.n <- readRDS("Results/wb_tm_fish_corrected.RDS")

# pca start ####
# get PCA data in order
# see get_pcdat() function in Inference_MS_functions.R for more details
pc.dat <- ldply(list(
  obs = get_pcdat(obs),
  wb.raw = get_pcdat(wb.raw),
  wb.n = get_pcdat(wb.n),
  tm = get_pcdat(tm),
  tm.niche = get_pcdat(tm.niche),
  tm.nn = get_pcdat(tm.nn),
  wb.tm = get_pcdat(wb.tm),
  wb.tm.n = get_pcdat(wb.tm.n) 
))

# add grouping variable
pc.dat$grp <- as.factor(rep(c("obs","wb.raw","wb.n", "tm", "tm.niche", "tm.nn", "wb.tm", "wb.tm.n"), each = 17))

# PCA analysis
pca.obj <- prcomp(pc.dat[,c(2:8)], center = T, scale. = T)

# biplot
# for illustration only, not used in MS
ordiplot(pca.obj)
orditorp(pca.obj, display="species",
         col="red", cex = 1.25, air = 0.0005)
ordihull(pca.obj, groups=pc.dat$grp, draw="polygon", 
         col=c("grey90", # observed
               "yellow", # webbuilder
               "gold", # webbuilder + neutral
               "gold4", #trait-matching
               "skyblue", # trait-matching + niche
               "slategray", # trait-matching + niche +neutral
               "plum", # wb tm composite
               "green"), # wb-tm + neutral
         label=F)

# can color one model at a time to compare with observed by setting color to NA
# e.g. initial trait-matching vs observed:
ordiplot(pca.obj)
orditorp(pca.obj, display="species", col="red", cex = 1.25, air = 0.0005)
ordihull(pca.obj, groups=pc.dat$grp, draw="polygon", 
         col=c("black",# observed
               "yellow", # webbuilder
               NA, # webbuilder + neutral
               NA, #trait-matching
               NA, # trait-matching + niche
               NA, # trait-matching + niche +neutral
               NA, # wb tm composite
               NA), # wb-tm + neutral
         label=F)

# distance ####
# calculate centroid for each group
obs.center = centroid(pca.obj$x[1:17,])
raw.center = centroid(pca.obj$x[18:34,])
wb.n.center = centroid(pca.obj$x[35:51,])
tm.center = centroid(pca.obj$x[52:68,])
tm.niche.center = centroid(pca.obj$x[69:85,])
tm.nn.center = centroid(pca.obj$x[86:102,])
wb.tm.center = centroid(pca.obj$x[103:119,])
wb.tm.n.center = centroid(pca.obj$x[120:136,])

# calculate distance from observed centroid to centroid of each group
obs.raw <- distance2(raw.center, obs.center)
obs.wb.n.dist <- distance2(obs.center, wb.n.center)
obs.tm.dist <- distance2(obs.center, tm.center)
obs.tm.niche.dist <- distance2(obs.center,
                               tm.niche.center)
obs.tm.nn.dist <- distance2(obs.center, tm.nn.center)
obs.wb.tm.dist <- distance2(obs.center, wb.tm.center)
obs.wb.tm.n <- distance2(obs.center, wb.tm.n.center)

# make table of distances
obs.tab <- data.frame(tm = obs.tm.dist,
                      tm.niche = obs.tm.niche.dist,
                      tm.nn = obs.tm.nn.dist,
                      wb = obs.raw, 
                      wb.tm = obs.wb.tm.dist,
                      wb.n = obs.wb.n.dist,
                      wb.tm.n = obs.wb.tm.n)
write.csv(obs.tab,
          "Results/Stats/PCA_euclidean_distance_observed.csv",
          row.names = FALSE)
