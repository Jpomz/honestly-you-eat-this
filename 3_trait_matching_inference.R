# Trait-matching
# 31 Dec 2017
# Justin Pomeranz
#jfpomeranz@gmail.com

# paramterize Niche model using method presented in Gravel et al. 2013 Methods in Eco Evo. Inferring food web structure from predator-prey body size relationships

library(pracma)
library(plyr)
library(tidyverse)

# gravel functions
# from Gravel et al. 2013 supplementary information
source("gravel_functions.R")
# calculate confusion matrices
source("adj_conf_matrix function.R")
# useful functions
source("Inference_MS_functions.R")

# data ####
# invertebrate biomass and abundance
invert <- readRDS("estimated invert bodymass.RDS")
# fish biomass and abundance
# the RDS file is a list of 4 different fish sizes
# [[2]] is the mean minimum 
fish <- readRDS("estimated fish bodymass and abundance.RDS")[[2]]
# rbind invert to each fish, filter out na(dw), order by dw and then split into list of sites
dw <- bind_rows(fish, invert)
dw <- dw[!is.na(dw$dw),]
dw <- dw[order(dw$dw),]
dw <- split(dw, list(dw$site))

# observed webs
obs.A <- readRDS("observed pred-prey.RDS")
# subset obs.A to only include webs with community data
obs.A <- obs.A[names(obs.A) %in% unique(invert$site)]

# pred-prey pairs ####
# empty list
A.pairs <- NULL
# for each web in obs_A
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

# pairs biomass ####
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

#training data ####
# list, where each element contains all pred-prey biomass pairs except for the web that is being inferred
# e.g. the list element "Blackrock" does not contain feeding paris from Blackrock, and will be used to infer interactions at that site. 
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
  out <- reg_fn(Bprey, Bpred, quartil = c(0.01, 0.97))
})
# web params ####
# list of all body sizes at a site
Ball <- map(dw, function (x){
  log10(x$dw)
})
# calculate web paramters using get_pars_Niche() from Gravel et al 2013
web.pars <- map2(pars.list, Ball, get_pars_Niche)

# infer links ####
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
saveRDS(inf, "Initial trait matching inference.RDS")
# save observed adj matrices matched to inferred
saveRDS(obs, "observed matrices matched to inferred.RDS")
# sum of links per web
sum.links <- sapply(web.links.inf, sum)


# niche forbidden ####
# niche forbidden
taxa.forbid <- c("Amphipoda", "Atalophlebioides", "Austroclima", "Austrosimulium", "Coloburiscus", "Deleatidium", "Nesameletus","Oligochaetae", "Oxyethira", "Potamopyrgus", "Zephlebia")

inf.niche <- map(inf,
                 rm_niche,
                 taxa = taxa.forbid)
# # save niche pruned trait matching ####
saveRDS(inf.niche, "Niche pruned trait matching inference.RDS")


# neutral forbidden ####
# Subset dw to match inference
dw.sub <- NULL
for(web in 1:length(dw)){
  dw.sub[[web]] <- dw[[web]][which(dw[[web]]$taxa %in%
                                     colnames(inf[[web]])),]
}
# save object for rank abundance calcs, supplemental
#saveRDS(dw.sub, "ab, dw, info for sites, subset to match inference.RDS")

# vector of abundances
ab.vec <- llply(dw.sub, function (x){
  x$no.m2})
#vector of taxa
ab.taxa <- llply(dw.sub, function (x){x$taxa})

# threshold vector
threshold <- c(1.0e-08, 1.5e-8, 3.0e-08, 5.9e-08,
               1.0e-07, 1.5e-7, 3.0e-07, 5.9e-07,
               1.0e-06, 1.5e-6, 3.0e-06, 5.9e-06,
               1.0e-05, 1.5e-5, 3.0e-05, 5.9e-05,
               1.0e-04, 1.5e-4, 3.0e-04, 5.9e-04,
               1.0e-03, 1.5e-3, 3.0e-03, 5.9e-03,
               1.0e-02, 1.5e-2, 3.0e-02, 5.9e-02)

# calculate relative abundance matrices
rel.ab.matr <- map2(ab.vec, ab.taxa, get_rel_ab)
# match to inferred matrices
rel.ab.matr <- map2(rel.ab.matr, inf,
                    match_matr)
# rel.ab.matr$observed == relative abundance matrices
rel.ab.matr <- llply(rel.ab.matr, function (x){
  x$observed
})
# save relative abundance matrices ####
saveRDS(rel.ab.matr, "relative abundance matrices.RDS")

# make a list of neutrally forbidden links at different thresholds
inf.neutral <- map(threshold, function (x){
  map(rel.ab.matr, rm_neutral, threshold = x)})
names(inf.neutral) <- threshold
# multiply inferred matrices by different neutral thresholds
inf.neutral <- map(inf.neutral, function (x){
  map2(x, inf, ~.x*.y)
})

# save Neutral forbidden trait matching ####
saveRDS(inf.neutral, "Neutral trait matching inference.RDS")

# prune niche from neutral matrices
inf.niche.neutral <- map(inf.neutral, function (x){
  map(x, rm_niche, taxa = taxa.forbid)})
# save niche + neutral matrices
saveRDS(inf.niche.neutral, ("Neutral + Niche trait matching inference.RDS"))


# AUC initial ####
auc.init <- ldply(map2(obs, inf, get_auc))
auc.init.mean <- mean(auc.init$V1, na.rm = TRUE)

# AUC niche ####
auc.niche <- ldply(map2(obs, inf.niche, get_auc))
auc.niche.mean <- mean(auc.niche$V1)

# AUC neutral ####
auc.neutral <- NULL
for(web in 1:length(obs)){
  auc.web <- NULL
  for(t in 1:length(inf.neutral)){
    auc.web[[t]] <- get_auc(obs[[web]],
                            inf.neutral[[t]][[web]])
  }
  auc.neutral[[web]] <- auc.web
}

auc.neutral.df <- data.frame(auc = flatten_dbl(auc.neutral),
                 thresh = log10(as.numeric(threshold)),
                 site = rep(names(obs),
                            each = length(threshold)),
                 stringsAsFactors = FALSE)



# global max AUC Neutral
global.thresh.neutral <- auc.neutral.df %>%
  group_by(thresh) %>%
  summarize(mean.auc = mean(na.omit(auc))) %>%
  top_n(1, wt = mean.auc)


# AUC Niche + Neutral ####
auc.niche.neutral <- NULL
for(web in 1:length(obs)){
  auc.web <- NULL
  for(t in 1:length(inf.niche.neutral)){
    auc.web[[t]] <- get_auc(obs[[web]],
                            inf.niche.neutral[[t]][[web]])
  }
  auc.niche.neutral[[web]] <- auc.web
}

auc.niche.neutral.df <- data.frame(auc = 
                          flatten_dbl(auc.niche.neutral),
                        thresh =
                          log10(as.numeric(threshold)),
                        site = rep(names(obs),
                                   each =
                                     length(threshold)),
                             stringsAsFactors = FALSE)

# global max auc
global.thresh.nn <- auc.niche.neutral.df %>%
  group_by(thresh) %>%
  summarize(mean.auc = mean(na.omit(auc))) %>% 
  top_n(1, wt = mean.auc)


# TSS ####
# working with neutral abundance threshold 3.0e-04
# inf.neutral[[19]]
neutral <- inf.neutral[[19]]

# TSS initial ####
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
tss.neutral.sd <- mean(tss.neutral$V1)

# neutral and niche forbidden ####
# 1.5e-8
neutral.niche <- inf.niche.neutral[[2]]
tss.niche.neutral <- ldply(
  pmap(list(obs = obs,
            inf = neutral.niche),
       get_tss))
tss.nn.mean <- mean(tss.niche.neutral$V1)
tss.nn.sd <- sd(tss.niche.neutral$V1)



# fn & fp ####
# initial


initial.false <- ldply(map2(obs, inf, false_prop)) %>% 
  summarize(mean.fp = mean(fp), sd.fp = sd(fp),
            mean.fn = mean(fn), sd.fn = sd(fn))
niche.false <- ldply(map2(obs, inf.niche, false_prop)) %>% 
  summarize(mean.fp = mean(fp), sd.fp = sd(fp),
            mean.fn = mean(fn), sd.fn = sd(fn))
neutral.false <- ldply(map2(obs, neutral, false_prop)) %>% 
  summarize(mean.fp = mean(fp), sd.fp = sd(fp),
            mean.fn = mean(fn), sd.fn = sd(fn))
nn.false <- ldply(map2(obs, neutral.niche, false_prop)) %>% 
  summarize(mean.fp = mean(fp), sd.fp = sd(fp),
            mean.fn = mean(fn), sd.fn = sd(fn))

false.tab <- rbind(initial.false, niche.false, neutral.false, nn.false)
false.tab$inference <- c("Initial", "Niche", "Neutral", "Niche + Neutral")

false.tab <- gather(false.tab, "var", "val", 1:4) %>%
  spread(var, val)

# table of auc, tss, threshold ####
write_csv(data.frame(inference =
        c("Initial", "Niche", "Neutral", "Niche + Neutral"),
        AUC = c(auc.init.mean, auc.niche.mean,
                as.double(global.thresh.neutral[2]),
                as.double(global.thresh.nn[2])),
        TSS = c(tss.initial.mean, tss.niche.mean,
                tss.neutral.mean, tss.nn.mean),
        sd.TSS = c(tss.initial.sd, tss.niche.sd,
                   tss.neutral.sd, tss.nn.sd),
        Threshold = c("NA", "NA",
                  10^as.double(global.thresh.neutral[1]),
                  10^as.double(global.thresh.nn[1])),
        mean.fp = false.tab$mean.fp,
        sd.fp = false.tab$sd.fp,
        mean.fn = false.tab$mean.fn,
        sd.fn = false.tab$sd.fn),
        "Mean AUC and TSS trait matching.csv")

# # plots ####
# # auc ####
# # Niche ####
# # AUC ~ threshold, facet by site
# # max AUC = red
# auc.neutral.df %>% group_by(site) %>% 
#   mutate(max.auc = max(na.omit(auc)),
#          is.max = auc == max.auc) %>% 
#   ggplot(aes(x = thresh,
#              y = auc,
#              color = is.max)) +
#   facet_wrap(~site) +
#   geom_point() +
#   scale_color_manual(values = c("black", "red"))+
#   geom_hline(aes(yintercept = 0.5),
#              linetype = "dashed") +
#   theme_classic()
# # density of thresholds == max.auc
# auc.neutral.df %>% group_by(site) %>%
#   top_n(1, wt = auc) %>%
#   .[match(unique(.$site), .$site),] %>%
#   ggplot(aes(x = thresh)) +
#   geom_density() +
#   theme_classic()
# # auc ~ threshold
# # all sites combined
# # max threshold = red points
# # plot of global
# auc.neutral.df %>% 
#   ggplot(aes(x = thresh,
#              y = auc)) +
#   geom_point() +
#   stat_summary(aes(y = auc,group=1),
#                fun.y=mean,
#                colour="grey",
#                geom="line",
#                size = 2,
#                group= 1) +
#   geom_point(data = auc.neutral.df %>%
#                filter(thresh > -3.9, thresh < -3.7),
#              aes(x = thresh, y = auc), color = "red")+
#   theme_classic()
# 
# # Niche + Neutral ####
# # AUC ~ threshold, facet by site
# # max AUC = red
# # plot facet by site 
# auc.niche.neutral.df %>% group_by(site) %>% 
#   mutate(max.auc = max(na.omit(auc)),
#          is.max = auc == max.auc) %>% 
#   ggplot(aes(x = thresh,
#              y = auc,
#              color = is.max)) +
#   facet_wrap(~site) +
#   geom_point() +
#   scale_color_manual(values = c("black", "red"))+
#   geom_hline(aes(yintercept = 0.5),
#              linetype = "dashed") +
#   theme_classic()
# # density of thresholds == max.auc
# auc.niche.neutral.df %>% group_by(site) %>%
#   top_n(1, wt = auc) %>%
#   .[match(unique(.$site), .$site),] %>%
#   ggplot(aes(x = thresh)) +
#   geom_density() +
#   theme_classic()
# # auc ~ threshold
# # all sites combined
# # max threshold = red points
# # plot of global
# auc.niche.neutral.df %>% 
#   ggplot(aes(x = thresh,
#              y = auc)) +
#   geom_point() +
#   stat_summary(aes(y = auc,group=1),
#                fun.y=mean,
#                colour="grey",
#                geom="line",
#                size = 2,
#                group= 1) +
#   geom_point(data = auc.niche.neutral.df %>%
#                filter(thresh > -3.6, thresh < -3.5),
#              aes(x = thresh, y = auc), color = "red")+
#   theme_classic()
# 
# # TSS ####
# # Niche + neutral
# # TSS ~ threshold
# # colored by site
# # size = max TSS
