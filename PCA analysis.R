# code pulled from webbuilder inference.R

# pca start ####
# https://www.r-bloggers.com/clustering-mixed-data-types-in-r/
pc.dat <- ldply(list(#wb = 
  #ldply(wb.matrices, function (x){
  #Get.web.stats(x)}),
  obs = 
    ldply(web.match, function (x){
      Get.web.stats(x)
    }),
  inf = ldply(inf, function (x){
    Get.web.stats(x)
  })))
pc.dat$grp <- rep(c(#"wb",
  "obs", "inf"), each = 17)

require(vegan)
require(cluster)
require(Rtsne)
pca.obj <- prcomp(pc.dat[,c(3:7, 11:12)], center = T, scale. = T)


adonis(pca.obj$x ~ grp, data = pc.dat, method='eu')

biplot(pca.obj)
ordiplot(pca.obj, xlim = c(-5,5))
orditorp(pca.obj,display="species",col="red",air=0.01)
ordihull(pca.obj,groups=pc.dat$grp,draw="polygon",col="grey90",label=F)


gower_dist <- daisy(pc.dat[,c(3:7, 11:13, 15)],
                    metric = "gower",
                    type = list(logratio = 3))
summary(gower_dist)
gower_mat <- as.matrix(gower_dist)

pc.dat[which(gower_mat == min(
  gower_mat[gower_mat != min(gower_mat)]),
  arr.ind = TRUE)[1, ], ]
pc.dat[which(gower_mat == max(
  gower_mat[gower_mat != max(gower_mat)]),
  arr.ind = TRUE)[1, ], ]

pam_fit <- pam(gower_dist,
               diss = TRUE,
               k = 3)
summary(pam_fit)
pc.dat[pam_fit$medoids,]

tsne_obj <- Rtsne(gower_dist, is_distance = TRUE, perplexity = 10)

tsne_data <- tsne_obj$Y %>%
  data.frame() %>%
  setNames(c("X", "Y")) %>%
  mutate(cluster = factor(pam_fit$clustering),
         name = pc.dat$grp)

ggplot(aes(x = X, y = Y), data = tsne_data) +
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



distance2 <- function(x,y) sum((x-y)^2)
centroid <- function(x) rowMeans(x)
x = pca.obj$x[1:17,]
y = pca.obj$x[18:34,]

distance2(centroid(x), centroid(y))
