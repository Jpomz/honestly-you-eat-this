# webbuilder * traitmatching

library(plyr)
library(tidyverse)
# webbuilder functions
source("Useful WebBuilder functions.R")
# useful food web functions from petchey
source("C:\\Users\\Justin\\Documents\\Data\\FW modelling Petchey Github\\ttl-resources-master\\food_web\\FoodWebFunctions.r")
# useful functions 
source("Inference_MS_functions.R")
# TSS function

# data ####
web.match <- readRDS("observed matrices matched to inferred.RDS")
wb.matrices <- readRDS("wb matrices matched to inferred.rds")
tm <- readRDS("Initial trait matching inference.RDS")

wb.matrices <- map2(wb.matrices, tm, ~.x*.y)
# initial wb * tm
saveRDS(wb.matrices, "wb x tm for PCA.RDS")

# neutral abundance correction ####
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
                                     each = length(threshold)),
                               stringsAsFactors = FALSE)

# global max AUC Fish Neutral
global.thresh.neutral.f <- auc.f.neutral.df %>%
  group_by(thresh) %>%
  summarize(mean.auc = mean(na.omit(auc))) %>%
  top_n(1, wt = mean.auc)


# TSS ####
# initial
tss.init.mean <- ldply(
  map2(web.match, wb.matrices, get_tss)) %>% 
  summarize(tss = mean(V1)) %>% as.double()
tss.init.sd <- ldply(
  map2(web.match, wb.matrices, get_tss)) %>% 
  summarize(sd(V1)) %>% as.double()
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
saveRDS(wb.f.n, "wb x tm x n for PCA.RDS")
tss.n.f.mean <- ldply(
  map2(web.match, wb.f.n, get_tss)) %>% 
  summarize(tss = mean(V1)) %>% as.double()
tss.n.f.sd <- ldply(
  map2(web.match, wb.f.n, get_tss)) %>% 
  summarize(sd(V1)) %>% as.double()


# fp & fn ####
# initial
false.init <- ldply(
  map2(web.match, wb.matrices, false_prop)) %>%
  summarize(mean.fp = mean(fp), sd.fp = sd(fp),
            mean.fn = mean(fn), sd.fn = sd(fn))
# neutral
false.n <- ldply(
  map2(web.match, wb.n, false_prop)) %>%
  summarize(mean.fp = mean(fp), sd.fp = sd(fp),
            mean.fn = mean(fn), sd.fn = sd(fn))

# neutral fish correction
false.n.f <- ldply(
  map2(web.match, wb.f.n, false_prop)) %>%
  summarize(mean.fp = mean(fp), sd.fp = sd(fp),
            mean.fn = mean(fn), sd.fn = sd(fn))

false.tab <- rbind(false.init, false.n, false.n.f)
# summary table ####
wb.tab <- data.frame(inference = c("Wb * trait", "Neutral", "Neutral Fish correction"),
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
write_csv(wb.tab, "Wb x tm AUC, TSS, fp+fn.csv")




