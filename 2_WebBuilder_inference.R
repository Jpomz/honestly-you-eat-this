# WebBuilder inference of food webs
# justin Pomeranz
# jfpomeranz@gmail.com

# Libraries
library(plyr)
library(dplyr)
library(purrr)
library(ggplot2)

# webbuilder functions 
# source code available in the supplementary material of Gray et al. 2015. Joinng the dots: An automated method for constructing food webs from compendia of published interactions https://doi.org/10.1016/j.fooweb.2015.09.001
# if using, please cite Gray et al. 2015. 
source("Functions/Useful WebBuilder functions.R")

# useful food web functions, modified from code by Owen Petchey:
# Original source code can be found here: https://github.com/opetchey/ttl-resources/blob/master/food_web/FoodWebFunctions.r
source("Functions/Foodweb_functions.R")
# functions written for this manuscript
source("Functions/Inference_MS_functions.R")

# All observed pred-prey adjacency matrices for registry construction
web.list <- readRDS("Data/Raw_data/observed_pred_prey_A.RDS")

# Observed adjacency matrices matched to trait-matching inference. 
# This is to make comparisons between trait-mathcing and WebBuilder fair
web.match <- readRDS("Results/observed_matrices_matched_to_inferred.RDS")

# convert adjacency matrix to predator-prey pairs
pairs.list <- web.list %>%
  llply(function (x){
    as.data.frame(Matrix.to.list(x),
                  stringsAsFactors = FALSE)})

# make list of pred-prey pairs for each web to be inferred
register.list <- NULL
for (i in 1:length(web.match)){
  register.list[[i]] <- ldply(pairs.list)
  register.list[[i]]$source.id <- "Thompson and Townsend, 1999, 2004"
}

# change .id column to source.id
# this is to match data format required in WebBuilder functions 
register.list <- register.list %>%
  llply(function (x){
    colnames(x) <- c("food.web",
                     "resource",
                     "consumer",
                     "source.id");x
  })

# read in file with taxonomy information
taxonomy <- read.csv("Data/Raw_data/taxonomy.csv",
                     stringsAsFactors = FALSE)

# add taxonomy information for resource and consumers
for (i in 1:length(register.list)){
  x <- register.list[[i]]
  # add taxonomy info for resource (prey) column
  x <- merge(x, taxonomy[,1:4], 
             by.x ="resource", by.y = "name", 
             all.x = T)
  # add prefix "res." to colnames
  # this is to comply with WebBuilder functions
  colnames(x)[5:7] <- paste("res.", 
                            colnames(x[c(5:7)]), sep = "")
  # add taxonomy info for consumer (predator) column
  x <- merge(x,taxonomy[,1:4], by.x ="consumer",
             by.y = "name", all.x = T)
  # add prefix "con." to colnames
  # this is to comply with WebBuilder functions         
  colnames(x)[8:10] <- paste("con.", 
                             colnames(x[c(8:10)]), sep = "")
  x$linkevidence <- "direct"
  register.list[[i]] <- x
}

# read in registry made from other sources
other.registry <- read.csv(
  "Data/Raw_data/other_feeding_interactions.csv",
  stringsAsFactors = FALSE)

# combine registry list with other registry
for (i in 1:length(register.list)){
  register.list[[i]] <- bind_rows(
    register.list[[i]], other.registry)
}

# write csv of complete register of NZ feeding interactions
# this will be submitted to the WebBuilder registry database as described in Gray et al. 2015
write.csv(register.list[[1]][,c(-11)],
          "Data/complete_NZ_trophic_links.csv",
          row.names = FALSE)

# add web names
names(register.list) <- names(web.match)
# subset each element in register.list to not contain name of web trying to infer links for
# e.g. the list element "Blackrock" does not contain feeding pairs from Blackrock, and will be used to infer interactions at that site. 
for (i in names(register.list)){
  register.list[[i]] <- subset(register.list[[i]],
                               food.web != i)
}


# infer links ####
# get a list of taxa from each food web in web.match
# e.g. only inferring links between taxa which were also used in trait-matching method
taxa.list <- llply(web.match,
                   function (x){
                     data.frame(node = colnames(x),
                                stringsAsFactors = FALSE)
                   })

# add taxonomy information to taxa.list
# note that minimum.res/con.method column is where you could modify inferences to be more or less restrictive. 
# e.g. set to "order" to infer max number of links, and "genus" would be most specific. See Gray et al. 2015 for more details
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
# add web names to inferred links
names(links) <- names(web.match)

# convert "pairs" format to adjacency matrix
wb.matrices <- map2(taxa.list, links, pairs_to_adj)

# match webbuilder inferred matrices to observed and trait-matching so all matrices are in same order. 
wb.matrices <- map2(wb.matrices, web.match, match_matr2)

saveRDS(wb.matrices, file ="Results/WebBuilder_inferred_matrices.rds")
# neutral processes ####
# read in neutral abundance matrice created in script 1
rel.ab.matr <- readRDS("Results/relative_abundance_matrices.RDS")
# vector of thresholds, as in script 1
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


# AUC ####
# initial ####
auc.init <- ldply(map2(web.match, wb.matrices, get_auc))
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
# fish corrected abundances ####
# read in matrices with "corrected" fish abundances created in script 1b
rel.ab.fish <- readRDS("Results/relative_ab_fish_x_1000.RDS")
#make list of neutrally forbidden links with corrected fish abundances
fish.neutral.list <- map(threshold, function (x){
  map(rel.ab.fish, rm_neutral, threshold = x)})
names(fish.neutral.list) <- threshold
# multiply inferred matrices by different neutral thresholds
fish.neutral.list <- map(fish.neutral.list, function (x){
  map2(x, wb.matrices, ~.x*.y)
})

auc.f.neutral <- get_auc_list(
  neutral.list = fish.neutral.list,
  obs = web.match)
# calculate numeric value of list element which has max AUC
fish.thresh.neutral <- get_max_auc(x = auc.f.neutral,
                                     threshold = threshold,
                                     site = names(web.match))

# supplemental figure ####
auc.f.neutral.df <- data.frame(
  auc = flatten_dbl(auc.f.neutral),
  thresh = log10(as.numeric(threshold)),                        site = rep(names(web.match), each = length(threshold)),
  stringsAsFactors = FALSE)

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

ggsave("Results/Figures/S3.png",
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
wb.n <- inf.neutral[[global.thresh.neutral$list.element]]
tss.n.mean <- ldply(
  map2(web.match, wb.n, get_tss)) %>% 
  summarize(tss = mean(V1)) %>% as.double()
tss.n.sd <- ldply(
  map2(web.match, wb.n, get_tss)) %>% 
  summarize(sd(V1)) %>% as.double()

# neutral fish correction * 1000
# threshold = 1.5e-5 == [[18]]
wb.f.n <- fish.neutral.list[[fish.thresh.neutral$list.element]]
names(wb.f.n) <- names(web.match)
saveRDS(wb.f.n, "Results/WebBuilder_inferred_fish_correction.RDS")

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
wb.tab <- data.frame(
  inference = c("Webbuilder initial", "Neutral",
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
  
wb.tab <- cbind(wb.tab, false.tab)
write.csv(wb.tab,
          "Results/Stats/WebBuilder_AUC_and_TSS.csv",
          row.names = FALSE)

