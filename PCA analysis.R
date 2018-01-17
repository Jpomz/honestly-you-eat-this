# code pulled from webbuilder inference.R
# useful food web functions from petchey
source("C:\\Users\\Justin\\Documents\\Data\\FW modelling Petchey Github\\ttl-resources-master\\food_web\\FoodWebFunctions.r")
library(plyr)
library(tidyverse)
library(vegan)
require(cluster)
require(Rtsne)
obs <- readRDS("observed matrices matched to inferred.RDS")
wb <- readRDS("wb for PCA.RDS")
tm <- readRDS("tm nn fish for PCA.RDS")
wb.tm <- readRDS("wb x tm for PCA.rds")


# pca start ####
# https://www.r-bloggers.com/clustering-mixed-data-types-in-r/
pc.dat <- ldply(list(
  obs = 
    ldply(obs, function (x){
      Get.web.stats(x, which.stats = 1)
    })[,c(1, 3:7, 11:13)],
  wb = ldply(wb, function (x){
    Get.web.stats(x, which.stats = 1)
  })[,c(1, 3:7, 11:13)],
  tm = ldply(tm, function (x){
    Get.web.stats(x, which.stats = 1)
  })[,c(1, 3:7, 11:13)],
  wb.tm = 
  ldply(wb.tm, function (x){
  Get.web.stats(x, which.stats = 1)
    })[,c(1, 3:7, 11:13)]
  ))
# add grouping variable
pc.dat$grp <- rep(c("obs","wb", "tm", "wb.tm"), each = 17)
pc.dat$land <- rep(c("Pn", "Pn", "T","T","T","Pn","T",
                      "T","T","T","T","Pn", "Pn", "Pn",
                      "T","T","Pn"), 4)


pca.obj <- prcomp(pc.dat[,c(2:9)], center = T, scale. = T)
adonis.obj <- adonis(pca.obj$x ~ grp*land, data = pc.dat, method='eu')

biplot(pca.obj)
ordiplot(pca.obj, xlim = c(-5,5))
orditorp(pca.obj,display="species",col="red",air=0.01)
ordihull(pca.obj,groups=pc.dat$grp,draw="polygon",col=c("grey90", "blue", "red", "purple"),label=F)

# distance
distance2 <- function(x,y) sum((x-y)^2)
centroid <- function(x) rowMeans(x)
obs.center = pca.obj$x[1:17,]
wb.center = pca.obj$x[18:34,]
tm.center = pca.obj$x[35:51,]
wb.tm.center = pca.obj$x[52:68,]

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



myfun <- function(g){
  pca.pair <- prcomp(
    pc.dat[which(pc.dat$grp=="obs" | pc.dat$grp == g),2:9],
    center = T, scale. = T)
  ordiplot(pca.pair, xlim = c(-5,5))
  orditorp(pca.pair, display="species", col="red")
  ordihull(pca.pair, 
           groups = pc.dat$grp[which(pc.dat$grp=="obs" | 
                                       pc.dat$grp == g)], 
           draw="polygon", 
           col=c("blue", "red"),
           label=F)
  obs.center = pca.pair$x[1:17,]
  pair.center = pca.pair$x[18:34,]
  dist <- distance2(centroid(obs.center),
                    centroid(pair.center))
  return(dist)
}

myfun("tm")
myfun("wb")
myfun("wb.tm")








gower_dist <- daisy(pc.dat[,c(2:9)],
                    metric = "gower",
                    type = list(logratio = 3))
summary(gower_dist)
gower_mat <- as.matrix(gower_dist)

# most sim
pc.dat[which(gower_mat == min(
  gower_mat[gower_mat != min(gower_mat)]),
  arr.ind = TRUE)[1, ], ]
# most dissim
pc.dat[which(gower_mat == max(
  gower_mat[gower_mat != max(gower_mat)]),
  arr.ind = TRUE)[1, ], ]

pam_fit <- pam(gower_dist,
               diss = TRUE,
               k = 4)
summary(pam_fit)
pc.dat[pam_fit$medoids,]

tsne_obj <- Rtsne(gower_dist, is_distance = TRUE, perplexity = 10)

tsne_data <- tsne_obj$Y %>%
  data.frame() %>%
  setNames(c("X", "Y")) %>%
  mutate(cluster = factor(pam_fit$clustering),
         name = pc.dat$grp)

ggplot(aes(x = X, y = Y, shape = name), data = tsne_data) +
  geom_point(aes(color = cluster))








data("iris")
iris
iris_c <- scale(iris[ ,1:4])
pca <- rda(iris_c)
adonis(iris_c ~ Species, data = iris, method='eu')

plot(pca, type = 'n', display = 'sites')
cols <- c('red', 'blue', 'green')
points(pca, display='sites', col = cols[iris$Species], pch = 16)
ordihull(pca, groups=iris$Species)




pca2 <- prcomp(iris[,1:4])
biplot(pca2)
class(pca2)

adonis(pca2$x~ Species, data = iris, method='eu')




