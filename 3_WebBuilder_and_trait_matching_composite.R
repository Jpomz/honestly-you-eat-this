# webbuilder * traitmatching composite model
# Justin Pomeranz
# jfpomeranz@gmail.com

# libraries
library(plyr)
library(dplyr)
library(ggplot2)
library(purrr)
# webbuilder functions
source("Functions/Useful WebBuilder functions.R")
# useful food web functions, modified from petchey
source("Functions/FoodWeb_Functions.R")
# functions written for this manuscript
source("Functions/Inference_MS_functions.R")

# data ####
web.match <- readRDS("Results/observed_matrices_matched_to_inferred.RDS")
wb.matrices <- readRDS("Results/WebBuilder_inferred_matrices.RDS")
tm <- readRDS("Results/Initial_trait_matching_inference.RDS")

wb.tm <- map2(wb.matrices, tm, ~.x*.y)
# initial wb * tm
saveRDS(wb.tm, "Results/wb_tm_initial.RDS")

# neutral abundance correction ####
rel.ab.matr <- readRDS("Results/relative_abundance_matrices.RDS")

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
  map2(x, wb.tm, ~.x*.y)
})

# fish corrected abundances ####
rel.ab.fish <- readRDS("Results/relative_ab_fish_x_1000.RDS")

fish.neutral.list <- map(threshold, function (x){
  map(rel.ab.fish, rm_neutral, threshold = x)})
names(fish.neutral.list) <- threshold
fish.neutral.list <- map(fish.neutral.list, function (x){
  map2(x, wb.tm, ~.x*.y)
})

# initial ####
auc.init <- ldply(map2(web.match, wb.tm, get_auc))
auc.init.mean <- mean(auc.init$V1, na.rm = TRUE)

# neutral ####
auc.neutral <- get_auc_list(
  neutral.list = inf.neutral,
  obs = web.match)
# calculate numeric value of list element which has max AUC
global.thresh.neutral <- get_max_auc(x = auc.neutral,
                                     threshold = threshold,
                                     site = names(web.match))

# fish corrected neutral ####
auc.f.neutral <- get_auc_list(
  neutral.list = fish.neutral.list,
  obs = web.match)
# calculate numeric value of list element which has max AUC
fish.thresh.neutral <- get_max_auc(x = auc.f.neutral,
                                   threshold = threshold,
                                   site = names(web.match))


# TSS ####
# initial
tss.init.mean <- ldply(
  map2(web.match, wb.tm, get_tss)) %>% 
  summarize(tss = mean(V1)) %>% as.double()
tss.init.sd <- ldply(
  map2(web.match, wb.tm, get_tss)) %>% 
  summarize(sd(V1)) %>% as.double()
# neutral
# threshold == 1e-9 == inf.neutral[[1]]
wb.tm.n <- inf.neutral[[global.thresh.neutral$list.element]]
tss.n.mean <- ldply(
  map2(web.match, wb.tm.n, get_tss)) %>% 
  summarize(tss = mean(V1)) %>% as.double()
tss.n.sd <- ldply(
  map2(web.match, wb.tm.n, get_tss)) %>% 
  summarize(sd(V1)) %>% as.double()
# neutral fish correction * 1000
# threshold = 1.5e-5 == [[18]]
wb.tm.f.n <- fish.neutral.list[[fish.thresh.neutral$list.element]]

saveRDS(wb.tm.f.n, "Results/wb_tm_fish_corrected.RDS")
tss.n.f.mean <- ldply(
  map2(web.match, wb.tm.f.n, get_tss)) %>% 
  summarize(tss = mean(V1)) %>% as.double()
tss.n.f.sd <- ldply(
  map2(web.match, wb.tm.f.n, get_tss)) %>% 
  summarize(sd(V1)) %>% as.double()


# fp & fn ####
# initial
false.init <- ldply(
  map2(web.match, wb.tm, false_prop)) %>%
  summarize(mean.fp = mean(fp), sd.fp = sd(fp),
            mean.fn = mean(fn), sd.fn = sd(fn))
# neutral
false.n <- ldply(
  map2(web.match, wb.tm.n, false_prop)) %>%
  summarize(mean.fp = mean(fp), sd.fp = sd(fp),
            mean.fn = mean(fn), sd.fn = sd(fn))

# neutral fish correction
false.n.f <- ldply(
  map2(web.match, wb.tm.f.n, false_prop)) %>%
  summarize(mean.fp = mean(fp), sd.fp = sd(fp),
            mean.fn = mean(fn), sd.fn = sd(fn))

false.tab <- rbind(false.init, false.n, false.n.f)
# summary table ####
wb.tm.tab <- data.frame(
  inference = c("Wb * trait",
                "Neutral",
                "Neutral Fish correction"),
  AUC = c(auc.init.mean,
          global.thresh.neutral$auc,
          fish.thresh.neutral$auc),
  TSS = c(tss.init.mean,
          tss.n.mean,
          tss.n.f.mean),
  sd.TSS = c(tss.init.sd,
             tss.n.sd,
             tss.n.f.sd),
  Threshold = c("NA",
                global.thresh.neutral$threshold,
                fish.thresh.neutral$threshold))
wb.tm.tab <- cbind(wb.tm.tab, false.tab)
write.csv(wb.tm.tab,
          "Results/Stats/Wb_tm_AUC_TSS.csv",
          row.names = FALSE)
