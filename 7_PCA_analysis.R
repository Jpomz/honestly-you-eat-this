# code pulled from webbuilder inference.R
# useful food web functions from petchey
# original functions can be found here: https://github.com/opetchey/ttl-resources/blob/master/food_web/FoodWebFunctions.r
source("FoodWebFunctions.r")
library(plyr)
library(tidyverse)
library(vegan)


obs <- readRDS("observed matrices matched to inferred.RDS")
# wb.raw == initial WB inference
wb.raw <- readRDS("wb matrices matched to inferred.RDS")
# wb == neutral prune + fish correction
wb.n <- readRDS("wb for PCA.RDS")
tm <- readRDS("Initial trait matching inference.RDS")
tm.niche <- readRDS("Niche pruned trait matching inference.RDS")
# tm == neutral, niche prune, fish correction
tm.nn <- readRDS("tm nn fish for PCA.RDS")
names(tm.nn) <- names(obs) # need to go back to original script and fix names
# wb x trait match, no neutral
wb.tm <- readRDS("wb x tm for PCA.RDS")
# wb.tm wb * tm neutral prune fish correction
wb.tm.n <- readRDS("wb x tm x n for PCA.rds")
names(wb.tm.n) <- names(obs) # need to go back to original script and fix names

# wb and wb.tm are basically the same
# maybe compare JUST wb*tm, without correcting fish / pruning neutral
# should be 6 total inferences

# pca start ####
# https://www.r-bloggers.com/clustering-mixed-data-types-in-r/
pc.dat <- ldply(list(
  obs = ldply(obs, function (x){
    Get.web.stats(x, which.stats = 1)
    })[,c(1, 3:7, 11:12)],
  wb.raw = ldply(wb.raw, function (x){
    Get.web.stats(x, which.stats = 1)
  })[,c(1, 3:7, 11:12)],
  wb.n = ldply(wb.n, function (x){
    Get.web.stats(x, which.stats = 1)
  })[,c(1, 3:7, 11:12)],
  tm = ldply(tm, function (x){
    Get.web.stats(x, which.stats = 1)
  })[,c(1, 3:7, 11:12)],
  tm.niche = ldply(tm.niche, function (x){
    Get.web.stats(x, which.stats = 1)
  })[,c(1, 3:7, 11:12)],
  tm.nn = ldply(tm.nn, function (x){
    Get.web.stats(x, which.stats = 1)
  })[,c(1, 3:7, 11:12)],
  wb.tm = ldply(wb.tm, function (x){
  Get.web.stats(x, which.stats = 1)
    })[,c(1, 3:7, 11:12)],
  wb.tm.n = ldply(wb.tm.n, function (x){
    Get.web.stats(x, which.stats = 1)
  })[,c(1, 3:7, 11:12)]
  ))
# add grouping variable
pc.dat$grp <- as.factor(rep(c("obs","wb.raw","wb.n", "tm", "tm.niche", "tm.nn", "wb.tm", "wb.tm.n"), each = 17))
# pc.dat$land <- rep(c("Pn", "Pn", "T","T","T","Pn","T",
#                       "T","T","T","T","Pn", "Pn", "Pn",
#                       "T","T","Pn"), 8)


pca.obj <- prcomp(pc.dat[,c(2:8)], center = T, scale. = T)
#adonis.obj <- adonis(pca.obj$x ~ grp*land, data = pc.dat, method='eu')

#biplot(pca.obj)
ordiplot(pca.obj)
orditorp(pca.obj, display="species", col="red", cex = 1.25, air = 0.0005)
ordihull(pca.obj, groups=pc.dat$grp, draw="polygon", 
         col=c("grey90","yellow", "gold", "gold4",
               "skyblue", "slategray",
               "plum", "green"),
         label=F)
# obs= grey90 
# tm = yellow; tm.niche = gold;  tm.nn= gold4
# wb.n= skyblue; wb.raw= slategray 
#wb.tm = plum; wb.tm.n = green

# pca figure of individual models
ordiplot(pca.obj)
orditorp(pca.obj, display="species", col="red", cex = 1.25, air = 0.0005)
ordihull(pca.obj, groups=pc.dat$grp, draw="polygon", 
         col=c("black","gold", NA, NA,
               NA, NA,
               NA, NA),
         label=F)
title("TM vs. Observed")

ordiplot(pca.obj)
orditorp(pca.obj, display="species", col="red", cex = 1.25, air = 0.0005)
ordihull(pca.obj, groups=pc.dat$grp, draw="polygon", 
         col=c("black", NA, "gold", NA,
               NA, NA,
               NA, NA),
         label=F)
title("TM, Niche pruned vs. Observed")

ordiplot(pca.obj)
orditorp(pca.obj, display="species", col="red", cex = 1.25, air = 0.0005)
ordihull(pca.obj, groups=pc.dat$grp, draw="polygon", 
         col=c("black", NA, NA, "gold",
               NA, NA,
               NA, NA),
         label=F)
title("TM, Niche + Neutral pruned vs. Observed")

ordiplot(pca.obj)
orditorp(pca.obj, display="species", col="red", cex = 1.25, air = 0.0005)
ordihull(pca.obj, groups=pc.dat$grp, draw="polygon", 
         col=c("black", NA, NA, NA,
               NA, "gold",
               NA, NA),
         label=F)
title("WB vs. Observed")

ordiplot(pca.obj)
orditorp(pca.obj, display="species", col="red", cex = 1.25, air = 0.0005)
ordihull(pca.obj, groups=pc.dat$grp, draw="polygon", 
         col=c("black", NA, NA, NA,
               "gold", NA,
               NA, NA),
         label=F)
title("WB, Neutral pruned vs. Observed")

ordiplot(pca.obj)
orditorp(pca.obj, display="species", col="red", cex = 1.25, air = 0.0005)
ordihull(pca.obj, groups=pc.dat$grp, draw="polygon", 
         col=c("black", NA, NA, NA, 
               NA, NA,
               "gold", NA),
         label=F)
title("WB*TM vs. Observed")

ordiplot(pca.obj)
orditorp(pca.obj, display="species", col="red", cex = 1.25, air = 0.0005)
ordihull(pca.obj, groups=pc.dat$grp, draw="polygon", 
         col=c("black", NA, NA, NA, 
               NA, NA,
               NA, "gold"),
         label=F)
title("WB*TM, Neutral pruned vs. Observed")
# distance
# surely I could come up with a function that does all of this in one go.... ####
distance2 <- function(x,y) sum((x-y)^2)
centroid <- function(x) rowMeans(x)
obs.center = centroid(pca.obj$x[1:17,])
raw.center = centroid(pca.obj$x[18:34,])
wb.n.center = centroid(pca.obj$x[35:51,])
tm.center = centroid(pca.obj$x[52:68,])
tm.niche.center = centroid(pca.obj$x[69:85,])
tm.nn.center = centroid(pca.obj$x[86:102,])
wb.tm.center = centroid(pca.obj$x[103:119,])
wb.tm.n.center = centroid(pca.obj$x[120:136,])


obs.raw <- distance2(raw.center, obs.center)
obs.wb.n.dist <- distance2(obs.center, wb.n.center)
obs.tm.dist <- distance2(obs.center, tm.center)
obs.tm.niche.dist <- distance2(obs.center,
                               tm.niche.center)
obs.tm.nn.dist <- distance2(obs.center, tm.nn.center)
obs.wb.tm.dist <- distance2(obs.center, wb.tm.center)
obs.wb.tm.n <- distance2(obs.center, wb.tm.n.center)

obs.tab <- data.frame(tm = obs.tm.dist,
                      tm.niche = obs.tm.niche.dist,
                      tm.nn = obs.tm.nn.dist,
                      wb = obs.raw, 
                      wb.tm = obs.wb.tm.dist,
                      wb.n = obs.wb.n.dist,
                      wb.tm.n = obs.wb.tm.n)
write.csv(obs.tab, "PCA euclidean distance observed.csv")
