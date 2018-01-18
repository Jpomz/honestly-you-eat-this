# code pulled from webbuilder inference.R
# useful food web functions from petchey
source("C:\\Users\\Justin\\Documents\\Data\\FW modelling Petchey Github\\ttl-resources-master\\food_web\\FoodWebFunctions.r")
library(plyr)
library(tidyverse)
library(vegan)


obs <- readRDS("observed matrices matched to inferred.RDS")
# wb == neutral prune + fish correction
wb <- readRDS("wb for PCA.RDS")
# tm == neutral, niche prune, fish correction
tm <- readRDS("tm nn fish for PCA.RDS")
# wb.tm wb * tm neutral prune fish correction
wb.tm <- readRDS("wb x tm for PCA.rds")
# wb.raw == initial WB inference
wb.raw <- readRDS("wb matrices matched to inferred.RDS")

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
  wb = ldply(wb, function (x){
    Get.web.stats(x, which.stats = 1)
  })[,c(1, 3:7, 12:13)],
  tm = ldply(tm, function (x){
    Get.web.stats(x, which.stats = 1)
  })[,c(1, 3:7, 12:13)],
  wb.tm = ldply(wb.tm, function (x){
  Get.web.stats(x, which.stats = 1)
    })[,c(1, 3:7, 12:13)]
  ))
# add grouping variable
pc.dat$grp <- as.factor(rep(c("obs","wb.raw","wb", "tm", "wb.tm"), each = 17))
pc.dat$land <- rep(c("Pn", "Pn", "T","T","T","Pn","T",
                      "T","T","T","T","Pn", "Pn", "Pn",
                      "T","T","Pn"), 5)


pca.obj <- prcomp(pc.dat[,c(2:8)], center = T, scale. = T)
#adonis.obj <- adonis(pca.obj$x ~ grp*land, data = pc.dat, method='eu')

#biplot(pca.obj)
ordiplot(pca.obj, xlim = c(-5,5))
orditorp(pca.obj, display="species", col="red", air=0.01)
ordihull(pca.obj, groups=pc.dat$grp, draw="polygon", 
         col=c("grey90","red", "blue", "green", "purple"),
         label=F)
# red == TM, blue == wb, green == wb.raw, purple == wb.tm

# distance
# surely I could come up with a function that does all of this in one go.... ####
distance2 <- function(x,y) sum((x-y)^2)
centroid <- function(x) rowMeans(x)
obs.center = pca.obj$x[1:17,]
raw.center = pca.obj$x[18:34,]
wb.center = pca.obj$x[35:51,]
tm.center = pca.obj$x[52:68,]
wb.tm.center = pca.obj$x[69:85,]

obs.wb.dist <- distance2(centroid(obs.center),
                         centroid(wb.center))
obs.tm.dist <- distance2(centroid(obs.center),
                         centroid(tm.center))
obs.wbtm.dist <- distance2(centroid(obs.center),
                           centroid(wb.tm.center))
wb.wbtm.dist <- distance2(centroid(wb.tm.center),
                        centroid(wb.center))
wb.tm.dist <- distance2(centroid(tm.center),
                        centroid(wb.center))
tm.wbtm.dist <- distance2(centroid(wb.tm.center),
                          centroid(tm.center))
raw.obs <- distance2(centroid(raw.center),
                     centroid(obs.center))
raw.tm <- distance2(centroid(raw.center),
                     centroid(tm.center))
raw.wb <- distance2(centroid(raw.center),
                     centroid(wb.center))
raw.wb.tm <- distance2(centroid(raw.center),
                     centroid(wb.tm.center))



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


