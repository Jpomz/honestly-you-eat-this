# Trait-matching
# Justin Pomeranz
# jfpomeranz@gmail.com

# This script parameterizes the model presented in Gravel et al. 2013 (referred to as trait mathcing in present paper) to infer feeding interactions based on empirical observations of predator-prey bodymass. 
# Inferred feeding interactions are modified  by taking into account niche, neutral, and niche + neutral forbideen links. 
# The script is organized as follows:
# 1) Preliminary data set up
# 2) parameterizing model for each empirical food web
# 3) inferring links for each food web
# 4) modifying inferred links (e.g. niche, neutral processes)
# 5) stats on model inferences (e.g. AUC, TSS)

# libraries
library(plyr)
library(dplyr)
library(purrr)
# functions written for this manuscript
source("Functions/Inference_MS_functions.R")

# functions from Gravel et al. 2013 supporting information
# if using, please cite original publication:
# Gravel, D., Poisot, T., Albouy, C., Velez, L., & Mouillot, D. (2013). Inferring food web structure from predator-prey body size relationships. Methods in Ecology and Evolution, 4, 1083-1090. doi:10.1111/2041-210X.12103
source("Functions/gravel_functions.R")

# 1) Preliminary data setup ####
#first, download data from DataDryad  https://doi.org/10.5061/dryad.k59m37f, and put into Data/Raw_data folder. 

# all adjacency matrices should be put into "Data/Raw_data/Adjacency_matrices" folder. 

# all other data should be put into Data/Raw_data folder

# read in estimated dry weight and abundance data
dw <- read.csv("Data/Raw_data/Taxa_dry_weight_abundance.csv")
# split dw by food.web column into a list
dw <- dw %>%
  split(list(.$food.web))

# read in modified adjacency matrices to a list
# see main text in manuscript for description of how adjacency matrices have been modified from originals
# get path directories for all adjacency matrices
path.dir <- list.files("Data/Raw_data/Adjacency_matrices", 
                       full.names = TRUE)
# make vector of food web names
names.A <- gsub(pattern = ".csv",
                replacement = "",
                x = list.files(
                  "Data/Raw_data/Adjacency_matrices"))
# make empty object for files to get listed into
obs.A <- NULL
# read in each file in path directory to make list of adjacency matrices
for(i in 1:length(path.dir)){
  obs.A[[i]] <- read.csv(path.dir[[i]],
                         row.names = 1)
}
# name each element in object to match food web names
names(obs.A) <- names.A

# save obs.A as .RDS data types
# this makes loading data easier in further analyses
saveRDS(obs.A, "Data/Raw_data/observed_pred_prey_A.RDS")


# subset obs.A to only include webs with bodymass data
obs.A <- obs.A[names(obs.A) %in% names(dw)]

# pred-prey pairs ####
# empty list
A.pairs <- NULL
# for each web in obs.A
for(web in 1:length(obs.A)){
  A = obs.A[[web]]
  genus = colnames(A) # vector of taxa names
  pred = c() # empty vector for pred
  prey = c() # empty vector for prey
  for(i in 1:nrow(A)){ # i = row = prey
    for(j in 1:ncol(A)){ # j = column = pred
      if(A[i,j]==1){ # if there is a link
        pred = c(pred, genus[j]) # add col taxa to pred
        prey = c(prey, genus[i]) # add row taxa to prey
      }
    }
  }
  pairs <- matrix(c(pred,prey), ncol = 2, byrow = F)
  colnames(pairs) <- c("pred", "prey")
  A.pairs[[web]] <- pairs 
}
names(A.pairs) <- names(obs.A)

# biomass pairs of predator-prey ####
# add biomass estimate to pred and prey pairs
dw.pairs <- NULL
for(web in 1:length(dw)){
  datout <- merge(A.pairs[[web]],
                  dw[[web]][,1:2],
                  by.x = "pred",
                  by.y = "taxa",
                  all.x = T)
  datout <- merge(datout,
                  dw[[web]][,1:2],
                  by.x = "prey",
                  by.y = "taxa",
                  all.x = T)
  m.pairs <- datout[,3:4]
  colnames(m.pairs) <- c("pred", "prey")
  dw.pairs[[web]] <- m.pairs[!is.na(m.pairs[,1]) &
                               !is.na(m.pairs[,2]),]
}
names(dw.pairs) <- names(obs.A)

# 2) parameterizing model ####
#training data ####
# list, where each element contains all pred-prey biomass pairs except for the web that is being inferred
# e.g. the list element "Blackrock" does not contain feeding pairs from Blackrock, and will be used to infer interactions at that site. 
training.list <- NULL
for (web in 1:length(dw.pairs)){
  dat <- ldply(dw.pairs[-web])
  training.list[[web]] <- dat
}
names(training.list) <- names(dw.pairs)

# Gravel model ####
# get parameters based on training list
pars.list <- map(training.list, function(x){
  Bprey = log10(x$prey)
  Bpred = log10(x$pred)
  out <- reg_fn(Bprey, Bpred, quartil = c(0.01, 0.99))
})
# web params ####
# list of all body sizes at a site
Ball <- map(dw, function (x){
  log10(x$dw)
})
# calculate web paramters using get_pars_Niche() from Gravel et al 2013
web.pars <- map2(pars.list, Ball, get_pars_Niche)

# 3) inferring links ####
# calculate food web links for taieri
# predation matrix
web.links.inf  <- map(web.pars, function (x){
  L_fn(x[,"n"],
       x[,"c"],
       x[,"low"],
       x[,"high"])
})
# add taxa names to inferred matrices
web.links.inf <- map2(web.links.inf, dw,
                      function (x, y){
                        dimnames(x) <- list(y$taxa, y$taxa)
                        x
                      })
# match order of rows/cols in observed and inferred adjacency matrices
match <- map2(obs.A, web.links.inf,
              match_matr)
# observed matrices
obs <- llply(match, function (x){
  x$observed
})
# trait matching inferred matrices
inf <- llply(match, function (x){
  x$inferred
})

# save initial inferred links ####
saveRDS(inf, "Results/Initial_trait_matching_inference.RDS")
# save observed adj matrices matched to inferred
saveRDS(obs, "Results/observed_matrices_matched_to_inferred.RDS")

# 4) modifying inferences ####
# niche forbidden ####
# read in vector of niche forbidden taxa
# See supplemental table S1
taxa.forbid <- as.character(
  read.csv("Data/Raw_data/niche_forbidden_taxa.csv",
           stringsAsFactors = FALSE,
           header = FALSE)[[1]])
# save as .RDS filetype for easier loading in further analyses
saveRDS(taxa.forbid, 
        "Data/Raw_data/niche_forbidden_taxa.RDS")

# remove niche forbidden links from inferred adjacency matrices
# see rm_niche() in Inference_MS_functions.R for more info
# e.g. set named columns to 0 
inf.niche <- map(inf,
                 rm_niche,
                 taxa = taxa.forbid)

# save niche pruned trait matching ####
saveRDS(inf.niche, "Results/Niche_pruned_trait_matching_inference.RDS")

# neutral forbidden ####
# Subset dataframe with dry weight and abundance data to match taxa used in inferences
dw.sub <- NULL
for(web in 1:length(dw)){
  dw.sub[[web]] <- dw[[web]][which(dw[[web]]$taxa %in%
                                     colnames(inf[[web]])),]
}

# Make list of vector of abundances for each site
ab.vec <- llply(dw.sub, function (x){
  x$no.m2})
#vector of taxa which mathces vector of abundances
ab.taxa <- llply(dw.sub, function (x){x$taxa})

# calculate relative abundance matrices
# see get_rel_ab() in Inference_MS_functions.R for more info
rel.ab.matr <- map2(ab.vec, ab.taxa, get_rel_ab)
# match to inferred matrices
rel.ab.matr <- map2(rel.ab.matr, inf,
                    match_matr)
# rel.ab.matr$observed == relative abundance matrices
rel.ab.matr <- llply(rel.ab.matr, function (x){
  x$observed
})
names(rel.ab.matr) <- names(inf)
# save relative abundance matrices ####
saveRDS(rel.ab.matr, "Results/relative_abundance_matrices.RDS")

# threshold vector
# this vector will be used to set values in relative abundance matrices to 0 or 1. 
# e.g.values in relative abundance matrices < 1.0e-8 will be converted to 0; values > 1.0e-8 will be converted to 1
threshold <- c(1.0e-08, 1.5e-8, 3.0e-08, 5.9e-08,
               1.0e-07, 1.5e-7, 3.0e-07, 5.9e-07,
               1.0e-06, 1.5e-6, 3.0e-06, 5.9e-06,
               1.0e-05, 1.5e-5, 3.0e-05, 5.9e-05,
               1.0e-04, 1.5e-4, 3.0e-04, 5.9e-04,
               1.0e-03, 1.5e-3, 3.0e-03, 5.9e-03,
               1.0e-02, 1.5e-2, 3.0e-02, 5.9e-02)

# make a list of lists with neutrally forbidden links at different thresholds
# see rm_neutral() in Inference_MS_functions.R for more info
inf.neutral <- map(threshold, function (x){
  map(rel.ab.matr, rm_neutral, threshold = x)})
names(inf.neutral) <- threshold
# inf.neutral is a nested list  
# the first hierarchy is a list of length 28, one for each threshold. 
# each of these lists is a list of 17, one for each site. 
# all of the sites within a threshold hierarchy have had their values converted to 0 or 1, if their values was < or > than threshold, repspectively. 
# i.e. inf.neutral[["threshold"]][["site"]]
# therefore, inf.neutral[[1]][[1]] --> threshold = "1e-08", site = "AkatoreB"

# multiply inferred matrices by different neutral thresholds
inf.neutral <- map(inf.neutral, function (x){
  map2(x, inf, ~.x*.y)
})

# save Neutral forbidden trait matching ####
saveRDS(inf.neutral, "Results/Neutral_trait_matching_inference.RDS")

# prune niche forbidden links from neutral pruned matrices
inf.niche.neutral <- map(inf.neutral, function (x){
  map(x, rm_niche, taxa = taxa.forbid)})
# save niche + neutral matrices
saveRDS(inf.niche.neutral, "Results/Neutral_and_Niche_trait_matching_inference.RDS")

# 5) model stats ####
# AUC initial ####
# calculate area under the recieving operator characteristic curve (AUC) for each inference type. 
# See get_auc() in Inference_MS_functions.R for more details
auc.init <- ldply(map2(obs, inf, get_auc))
# mean AUC for initial inference
auc.init.mean <- mean(auc.init$V1, na.rm = TRUE)

# AUC niche ####
auc.niche <- ldply(map2(obs, inf.niche, get_auc))
auc.niche.mean <- mean(auc.niche$V1)

# AUC neutral ####
# calculate AUC for all webs at each threshold
auc.neutral <- get_auc_list(
  neutral.list = inf.neutral,
  obs = obs)
# numeric value representing which list element in neutral threshold resulted in highest AUC
# see TSS calculations below
neutral.threshold <- get_max_auc(
  x = auc.neutral,
  threshold = threshold,
  site = names(obs))

# AUC Niche + Neutral ####
# calculate AUC for all webs at each threshold including niche pruned links
auc.niche.neutral <- get_auc_list(inf.niche.neutral, obs)

# numeric value for niche + neutral pruned inferences that results in max AUC
nn.threshold <- get_max_auc(auc.niche.neutral,
                            threshold = threshold,
                            site = names(obs))

# TSS ####
# neutral abundance threshold which resulted in max auc = 0.00059
# e.g. inf.neutral[[20]], 
neutral <- inf.neutral[[neutral.threshold$list.element]]

# TSS initial ####
# calculate mean and SD of AUC
tss.initial <- ldply(map2(obs, inf,
                          get_tss))
tss.initial.mean <- mean(tss.initial$V1)
tss.initial.sd <- sd(tss.initial$V1)
# TSS niche ####
tss.niche <- ldply(pmap(list(obs = obs,
                             inf = inf.niche),
                        get_tss))
tss.niche.mean <- mean(tss.niche$V1)
tss.niche.sd <- sd(tss.niche$V1)
# TSS neutral ####
tss.neutral <- ldply(pmap(list(obs = obs,
                               inf = neutral),
                          get_tss))
tss.neutral.mean <- mean(tss.neutral$V1)
tss.neutral.sd <- sd(tss.neutral$V1)

# neutral and niche forbidden ####
# neutral + Niche abundance threshold 1.5e-8
neutral.niche <- inf.niche.neutral[[nn.threshold$list.element]]
tss.niche.neutral <- ldply(
  pmap(list(obs = obs,
            inf = neutral.niche),
       get_tss))
tss.nn.mean <- mean(tss.niche.neutral$V1)
tss.nn.sd <- sd(tss.niche.neutral$V1)

# false negative & false positive ####
# initial
initial.false <- ldply(map2(obs, inf, false_prop)) %>% 
  summarize(mean.fp = mean(fp), sd.fp = sd(fp),
            mean.fn = mean(fn), sd.fn = sd(fn))
# niche
niche.false <- ldply(map2(obs, inf.niche, false_prop)) %>% 
  summarize(mean.fp = mean(fp), sd.fp = sd(fp),
            mean.fn = mean(fn), sd.fn = sd(fn))
# neutral
neutral.false <- ldply(map2(obs, neutral, false_prop)) %>% 
  summarize(mean.fp = mean(fp), sd.fp = sd(fp),
            mean.fn = mean(fn), sd.fn = sd(fn))
# niche + neutral
nn.false <- ldply(map2(obs, neutral.niche, false_prop)) %>% 
  summarize(mean.fp = mean(fp), sd.fp = sd(fp),
            mean.fn = mean(fn), sd.fn = sd(fn))
# bind all fp/fn rows together
false.tab <- rbind(initial.false, niche.false, neutral.false, nn.false)
# add inference type to df
false.tab$inference <- c("Initial", "Niche", "Neutral", "Niche + Neutral")

# table of auc, tss, threshold ####
# write results
write.csv(data.frame(
  inference = c("Initial", "Niche", "Neutral",
                "Niche + Neutral"),
  AUC = c(auc.init.mean, auc.niche.mean,
          neutral.threshold$auc,
          nn.threshold$auc),
  TSS = c(tss.initial.mean, tss.niche.mean, 
          tss.neutral.mean, tss.nn.mean),
  sd.TSS = c(tss.initial.sd, tss.niche.sd,
             tss.neutral.sd, tss.nn.sd),
  Threshold = c("NA", "NA", 
                neutral.threshold$threshold,
                nn.threshold$threshold),
  mean.fp = false.tab$mean.fp,
  sd.fp = false.tab$sd.fp,
  mean.fn = false.tab$mean.fn,
  sd.fn = false.tab$sd.fn),
  "Results/Stats/Mean_AUC_and_TSS_trait_matching.csv",
  row.names = FALSE)
