# code pulled from webbuilder inference.R
# useful food web functions from petchey
source("C:\\Users\\Justin\\Documents\\Data\\FW modelling Petchey Github\\ttl-resources-master\\food_web\\FoodWebFunctions.r")
library(plyr)
library(tidyverse)
library(vegan)


obs <- readRDS("observed matrices matched to inferred.RDS")
# wb.raw == initial WB inference
wb.raw <- readRDS("wb matrices matched to inferred.RDS")
# wb == neutral prune + fish correction
wb.n <- readRDS("wb for PCA.RDS")
# tm == neutral, niche prune, fish correction
tm.n <- readRDS("tm nn fish for PCA.RDS")
# wb x trait match, no neutral
wb.tm <- readRDS("wb x tm for PCA.RDS")
# wb.tm wb * tm neutral prune fish correction
wb.tm.n <- readRDS("wb x tm x n for PCA.rds")

# wb and wb.tm are basically the same
# maybe compare JUST wb*tm, without correcting fish / pruning neutral
# should be 6 total inferences

# pca start ####
# https://www.r-bloggers.com/clustering-mixed-data-types-in-r/
pc.dat <- ldply(list(
  obs = ldply(obs, function (x){
    Get.web.stats(x, which.stats = 1)
    })[,c(1, 3:7, 12:13)],
  wb.raw = ldply(wb.raw, function (x){
    Get.web.stats(x, which.stats = 1)
  })[,c(1, 3:7, 12:13)],
  wb.n = ldply(wb.n, function (x){
    Get.web.stats(x, which.stats = 1)
  })[,c(1, 3:7, 12:13)],
  tm.n = ldply(tm.n, function (x){
    Get.web.stats(x, which.stats = 1)
  })[,c(1, 3:7, 12:13)],
  wb.tm = ldply(wb.tm, function (x){
  Get.web.stats(x, which.stats = 1)
    })[,c(1, 3:7, 12:13)],
  wb.tm.n = ldply(wb.tm.n, function (x){
    Get.web.stats(x, which.stats = 1)
  })[,c(1, 3:7, 12:13)]
  ))
# add grouping variable
pc.dat$grp <- as.factor(rep(c("obs","wb.raw","wb.n", "tm.n", "wb.tm", "wb.tm.n"), each = 17))
pc.dat$land <- rep(c("Pn", "Pn", "T","T","T","Pn","T",
                      "T","T","T","T","Pn", "Pn", "Pn",
                      "T","T","Pn"), 6)


pca.obj <- prcomp(pc.dat[,c(2:8)], center = T, scale. = T)
#adonis.obj <- adonis(pca.obj$x ~ grp*land, data = pc.dat, method='eu')

#biplot(pca.obj)
ordiplot(pca.obj)
orditorp(pca.obj, display="species", col="red", cex = 1.25, air = 0.0005)
ordihull(pca.obj, groups=pc.dat$grp, draw="polygon", 
         col=c("grey90","yellow", "green", "red", 
               "blue", "plum1"),
         label=F)
# obs= grey90 
# tm.n= yellow wb.n= green 
# wb.raw= red wb.tm = blue wb.tm.n = plum1

# distance
# surely I could come up with a function that does all of this in one go.... ####
distance2 <- function(x,y) sum((x-y)^2)
centroid <- function(x) rowMeans(x)
obs.center = centroid(pca.obj$x[1:17,])
raw.center = centroid(pca.obj$x[18:34,])
wb.n.center = centroid(pca.obj$x[35:51,])
tm.n.center = centroid(pca.obj$x[52:68,])
wb.tm.center = centroid(pca.obj$x[69:85,])
wb.tm.n.center = centroid(pca.obj$x[86:102,])


obs.raw <- distance2(raw.center, obs.center)
obs.wb.n.dist <- distance2(obs.center, wb.n.center)
obs.tm.n.dist <- distance2(obs.center, tm.n.center)
obs.wb.tm.dist <- distance2(obs.center, wb.tm.center)
obs.wb.tm.n <- distance2(obs.center, wb.tm.n.center)

obs.tab <- data.frame(tm.n = obs.tm.n.dist,
                      wb = obs.raw, 
                      wb.tm = obs.wb.tm.dist,
                      wb.n = obs.wb.n.dist,
                      wb.tm.n = obs.wb.tm.n)


dist.tab <- matrix(NA, 3, 4)
dimnames(dist.tab) <- list(c("Observed",
                             "Trait match",
                             "WebBuilder"),
                           c("Trait match",
                             "raw WebBuilder",
                             "WebBuilder",
                             "WebBuilder:Trait match"))
dist.tab[1,1] <- obs.tm.dist
dist.tab[1,2] <- raw.obs
dist.tab[1,3] <- obs.wb.dist
dist.tab[1,4] <- obs.wbtm.dist
dist.tab[2,2] <- raw.tm
dist.tab[2,3] <- wb.tm.dist
dist.tab[2,4] <- tm.wbtm.dist
dist.tab[3,2] <- raw.wb
dist.tab[3,4] <- wb.wbtm.dist
#write.csv(dist.tab, "PCA euclidean distance.csv")








myfun <- function(g1, g2){
  pca.pair <- prcomp(
    pc.dat[which(pc.dat$grp==g1 | pc.dat$grp == g2),2:8],
    center = T, scale. = T)
  ordiplot(pca.pair, xlim = c(-5,5))
  orditorp(pca.pair, display="species", col="red")
  ordihull(pca.pair, 
           groups = pc.dat$grp[which(pc.dat$grp==g1 | 
                                       pc.dat$grp == g2)], 
           draw="polygon", 
           col=c("blue", "red"),
           label=F)
  g1.center = pca.pair$x[1:17,]
  g2.center = pca.pair$x[18:34,]
  dist <- distance2(centroid(g1.center),
                    centroid(g2.center))
  return(dist)
}

myfun("obs", "tm")
myfun("obs", "wb.raw")
myfun("obs", "wb")
myfun("tm", "wb")


# pca by land per group
pca.dat.list <- split(pc.dat, list(pc.dat$grp))
pca.list <- llply(pca.dat.list, function (x) {
  prcomp(x[,c(2:9)], center = T, scale. = T)
})

for (g in 1:length(pca.dat.list)){
  ordiplot(pca.list[[g]])
  orditorp(pca.list[[g]],display="species",col="red",air=0.01)
  ordihull(pca.list[[g]],groups=pca.dat.list[[g]]$land,
           draw="polygon",
           col=c("blue", "green"),label=F)
  title(names(pca.list)[g])
}


