# WebBuilder
# 12 jan 18
# justin Pomeranz
# jfpomeranz@gmail.com

# make a list of registries using all Taieri predation matrices -web trying to infer
library(plyr)
library(tidyverse)
# webbuilder functions
source("Useful WebBuilder functions.R")
# useful food web functions, modified from Owen Petchey:
#https://github.com/opetchey/ttl-resources/tree/master/food_web
source("FoodWebFunctions.r")
# useful functions 
source("Inference_MS_functions.R")

# read in all pred-prey matrices for registry construction
web.list <- readRDS("observed pred-prey.RDS")
# read in observed matrices matched to trait-matching inference. This is for reference to web names and taxa names in trait matching inference
web.match <- readRDS("observed matrices matched to inferred.RDS")

# convert predation matrix to matrix of pairs
pairs.list <- web.list %>%
  llply(function (x){
    as.data.frame(Matrix.to.list(x),
                  stringsAsFactors = FALSE)})

# make list of res-con pairs for registries
# each element in list contains all pairs EXCEPT those of
# the web they will be used to parameterize. 
# e.g. element named "blackrock" will contain all other
# pairs data [-"blackrock"], and be used to paramterize # the model in order to predict the blackrock web 
register.list <- NULL
for (i in 1:length(web.match)){
  register.list[[i]] <- ldply(pairs.list)
}


# change .id column to source.id
register.list <- register.list %>%
  llply(function (x){
    colnames(x) <- c("source.id",
                     "resource",
                     "consumer");x
  })

# read in file with taxonomy information
taxonomy <- read_csv("data/taxonomy.csv")

# add taxonomy information for resource and consumers
for (i in 1:length(register.list)){
  x <- register.list[[i]]
  x <- merge(x, taxonomy[,1:5], 
             by.x ="resource", by.y = "name", 
             all.x = T)
  colnames(x)[4:6] <- paste("res.", 
                            colnames(x[c(4:6)]), sep = "")
  x <- merge(x,taxonomy[,c(1:4,6)], by.x ="consumer",
             by.y = "name", all.x = T)         
  colnames(x)[8:10] <- paste("con.", 
                             colnames(x[c(8:10)]), sep = "")
  x$linkevidence <- "direct"
  register.list[[i]] <- x
}

# read in registry made from other sources
other.registry <- read_csv("data/fw_registry.csv")
other.registry$source.id <- as.character(other.registry$source.id)
# subset other.registry object to match with registry list
other.registry <- other.registry[,
                intersect(colnames(register.list[[1]]),
                          names(other.registry))]
# combine registry list with other registry
for (i in 1:length(register.list)){
  register.list[[i]] <- bind_rows(
    register.list[[i]], other.registry)
}

# make one complete registry to add to other registers later
complete.registry <- register.list[[1]]

# add web names
names(register.list) <- names(web.match)
# subset each element in register.list to not contain name of web trying to infer links for
for (i in names(register.list)){
  register.list[[i]] <- subset(register.list[[i]], source.id != i)
}

write.csv(complete.registry, "figs for MS/post poisot/complete webbuilder registry.csv")
# registry info ####
# write csv of registry info
# modify in excel to make table for publication
write_csv(register.list %>%
            ldply(function (x){
              data.frame(n = nrow(x))
            }),
          "C:/Users/Justin/Google Drive/Data/Predicting NZ Food Webs/figs for MS/post poisot/webbuilder registry info from R.csv")

# infer links ####
# get a list of taxa from web.list
# need to subset so only looking at webs which were also inferred using gravel method
taxa.list <- llply(web.match,
                   function (x){
                     data.frame(node = colnames(x),
                                stringsAsFactors = FALSE)
                   })

# add taxonomy information
taxa.list <- taxa.list %>% 
  llply(function (x){
    merge(x, taxonomy,by.x="node", 
          by.y="name", all.x=T)
  })

# infer links with WebBuilder function
# takes a while to run...
links <- NULL
for (i in 1:length(taxa.list)){
  links[[i]] <- WebBuilder(taxa.list[[i]],
                           register.list[[i]],
                           method = c("exact", "genus", "family", "order"))
}
names(links) <- names(web.match)

wb.matrices <- map2(taxa.list, links, pairs_to_adj)

# match wb matrices to web.match matrices
match_matr2 <- function (obj, target){
  # colnames(inf) have already been size-sorted
  index = intersect(rownames(target), rownames(obj))
  obj = obj[index, index]
  obj
}

wb.matrices <- map2(wb.matrices, web.match, match_matr2)

saveRDS(wb.matrices, file ="wb matrices matched to inferred.rds")

# neutral abundance correction
rel.ab.matr <- readRDS("relative abundance matrices.RDS")
threshold <- c(1.0e-09, 1.5e-9, 3.0e-09, 5.9e-09,
               1.0e-08, 1.5e-8, 3.0e-08, 5.9e-08,
               1.0e-07, 1.5e-7, 3.0e-07, 5.9e-07,
               1.0e-06, 1.5e-6, 3.0e-06, 5.9e-06,
               1.0e-05, 1.5e-5, 3.0e-05, 5.9e-05,
               1.0e-04, 1.5e-4, 3.0e-04, 5.9e-04,
               1.0e-03, 1.5e-3, 3.0e-03, 5.9e-03,
               1.0e-02, 1.5e-2, 3.0e-02, 5.9e-02)
               
# make a list of neutrally forbidden links at different thresholds
inf.neutral <- map(threshold, function (x){
  map(rel.ab.matr, rm_neutral, threshold = x)})
names(inf.neutral) <- threshold
# multiply inferred matrices by different neutral thresholds
inf.neutral <- map(inf.neutral, function (x){
  map2(x, wb.matrices, ~.x*.y)
})

# fish corrected abundances ####
rel.ab.fish <- readRDS("rel ab fish x 1000.RDS")

fish.neutral.list <- map(threshold, function (x){
  map(rel.ab.fish, rm_neutral, threshold = x)})
names(fish.neutral.list) <- threshold
fish.neutral.list <- map(fish.neutral.list, function (x){
  map2(x, wb.matrices, ~.x*.y)
})

# AUC ####
# initial ####
auc.init <- ldply(map2(web.match, wb.matrices, get_auc))
auc.init.mean <- mean(auc.init$V1, na.rm = TRUE)

# neutral ####
auc.neutral <- NULL
for(web in 1:length(web.match)){
  auc.web <- NULL
  for(t in 1:length(inf.neutral)){
    auc.web[[t]] <- get_auc(web.match[[web]],
                            inf.neutral[[t]][[web]])
  }
  auc.neutral[[web]] <- auc.web
}

# data frame of AUC
auc.neutral.df <- data.frame(auc = flatten_dbl(auc.neutral),
                        thresh = log10(as.numeric(threshold)),
                        site = rep(names(web.match),
                                    each = length(threshold)),
                             stringsAsFactors = FALSE)

# global max AUC Neutral
global.thresh.neutral <- auc.neutral.df %>%
  group_by(thresh) %>%
  summarize(mean.auc = mean(na.omit(auc))) %>%
  top_n(1, wt = mean.auc)

# fish corrected neutral ####
auc.f.neutral <- NULL
for(web in 1:length(web.match)){
  auc.web <- NULL
  for(t in 1:length(inf.neutral)){
    auc.web[[t]] <- get_auc(web.match[[web]],
                            fish.neutral.list[[t]][[web]])
  }
  auc.f.neutral[[web]] <- auc.web
}

# data frame of AUC Fish
auc.f.neutral.df <- data.frame(auc = 
                          flatten_dbl(auc.f.neutral),
                             thresh = 
                          log10(as.numeric(threshold)),
                             site = 
                          rep(names(web.match),
                               each = 
                          length(threshold)),
                             stringsAsFactors = FALSE)

# global max AUC Fish Neutral
global.thresh.neutral.f <- auc.f.neutral.df %>%
  group_by(thresh) %>%
  summarize(mean.auc = mean(na.omit(auc))) %>%
  top_n(1, wt = mean.auc)
# plot for suplemental
ggplot(auc.f.neutral.df, aes(x =thresh, y = auc))+
  geom_point() +
  stat_summary(aes(y = auc,group=1),
               fun.y=mean,
               colour="grey",
               geom="line",
               size = 2,
               group= 1) +
  theme_classic(base_size = 20) +
  labs(y = "AUC", x = expression(Log["10"]~Threshold))

ggsave("figs for MS\\post poisot\\wb fish corr auc.png",
       width = 420, height = 200, units = "mm")


# TSS ####
# initial
tss.init.mean <- ldply(
  map2(web.match, wb.matrices, get_tss)) %>% 
  summarize(tss = mean(V1)) %>% as.double()
tss.init.sd <- ldply(
  map2(web.match, wb.matrices, get_tss)) %>% 
  summarize(sd.tss = sd(V1)) %>% as.double()

# neutral
# threshold == 1e-9 == inf.neutral[[1]]
wb.n <- inf.neutral[[1]]
tss.n.mean <- ldply(
  map2(web.match, wb.n, get_tss)) %>% 
  summarize(tss = mean(V1)) %>% as.double()
tss.n.sd <- ldply(
  map2(web.match, wb.n, get_tss)) %>% 
  summarize(sd(V1)) %>% as.double()
# neutral fish correction * 1000
# threshold = 1.5e-5 == [[18]]
wb.f.n <- fish.neutral.list[[18]]
names(wb.f.n) <- names(web.match)
saveRDS(wb.f.n, "wb for PCA.RDS")
tss.n.f.mean <- ldply(
  map2(web.match, wb.f.n, get_tss)) %>% 
  summarize(tss = mean(V1)) %>% as.double()
tss.n.f.sd <- ldply(
  map2(web.match, wb.f.n, get_tss)) %>% 
  summarize(sd(V1)) %>% as.double()

# fp & fn ####
# initial ####
false.init <- ldply(
  map2(web.match, wb.matrices, false_prop)) %>%
  summarize(mean.fp = mean(fp), sd.fp = sd(fp),
            mean.fn = mean(fn), sd.fn = sd(fn))
# neutral ####
false.n <- ldply(
  map2(web.match, wb.n, false_prop)) %>%
  summarize(mean.fp = mean(fp), sd.fp = sd(fp),
            mean.fn = mean(fn), sd.fn = sd(fn))

# neutral fish correction ####
false.n.f <- ldply(
  map2(web.match, wb.f.n, false_prop)) %>%
  summarize(mean.fp = mean(fp), sd.fp = sd(fp),
            mean.fn = mean(fn), sd.fn = sd(fn))

false.tab <- rbind(false.init, false.n, false.n.f)
# summary table ####
wb.tab <- data.frame(inference = c("Webbuilder initial", "Neutral", "Neutral Fish correction"),
                     AUC = c(auc.init.mean,
                        as.double(global.thresh.neutral[2]),
                        as.double(global.thresh.neutral.f[2])),
                     TSS = c(tss.init.mean, tss.n.mean,
                             tss.n.f.mean),
                     sd.TSS = c(tss.init.sd, tss.n.sd,
                             tss.n.f.sd),
                     Threshold = c("NA",
                10**(as.double(global.thresh.neutral[1])),
                10**(as.double(global.thresh.neutral.f[1]))))
wb.tab <- cbind(wb.tab, false.tab)
write_csv(wb.tab, "Webbuilder AUC, TSS, fp+fn.csv")

