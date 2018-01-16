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
# calculate TSS
get_tss <- function (observed, inferred){
  # make sure adjacency matrices are same dimensions and have same colnames
  stopifnot(dim(observed) == dim(inferred), 
            identical(colnames(observed), colnames(inferred)))
  
  # subtract inferred from observed
  minus <- observed - inferred
  # multiply observed and inferred
  multiply <- observed * inferred
  # change values of 1 in multiplied matrix to 2
  multiply[multiply == 1] <- 2
  # add minus and multiply matrices
  prediction <- minus + multiply
  # prediction outcome matrix now has all 4 possibilities repreented as different integer values
  # 2 = true positive (a); links both obserevd & predicted
  # -1 = false positive (b); predicted but not observed
  # 1 = false negative (c); observed but not predicted
  # 0 = true negative (d); not predicted, not observed
    a = length(prediction[prediction==2])
    b = length(prediction[prediction==-1]) 
    c = length(prediction[prediction==1]) 
    d = length(prediction[prediction==0])
  # calculate TSS
  # TSS = (a*d - b*c)/((a+c)*(b+d))
  tss = (a*d - b*c)/((a+c)*(b+d))
  tss
}
# function to calculate relative abundance matrices
get_rel_ab <- function(vec, taxa){
  stopifnot(length(vec) == length(taxa))
  rel.ab <- vec / sum(vec)
  Nij <- matrix(0, length(vec), length(vec))
  for (i in 1:length(vec)){
    for (j in 1:length(vec)){
      Nij[i,j] <- rel.ab[i]*rel.ab[j]
    }
  }
  dimnames(Nij) <- list(taxa, taxa)
  Nij
}

# function to remove rel.abundance products < threshold
rm_neutral <- function(Nij, threshold){
  Nij[Nij > threshold] <-  1
  Nij[Nij < 1] <-  0 
  Nij
}

# function to remove links from niche forbidden taxa
rm_niche <- function(inf, taxa){
  for(name in (
    colnames(inf)[colnames(inf) %in% taxa])){
    inf[,name] <- 0
  }
  inf
}
# match observed and inferred matrices
match_matr <- function (obs, inf){
  # colnames(inf) have already been size-sorted
  index = intersect(rownames(inf), rownames(obs))
  obs = obs[index, index]
  inf = inf[index, index]
  list(observed = obs, inferred = inf)
}

# data ####
# invertebrate biomass and abundance
invert <- readRDS("estimated invert bodymass.RDS")
# fish biomass and abundance
fish <- readRDS("estimated fish bodymass and abundance.RDS")[[2]]
# rbind invert to each element in fish, filter out na(dw), order by dw and then split each list into list of sites
dw <- bind_rows(fish, invert)
dw <- dw[!is.na(dw$dw),]
dw <- dw[order(dw$dw),]
dw <- split(dw, list(dw$site))
# observed webs
obs.A <- readRDS("observed pred-prey.RDS")

# subset obs_A to only include webs with community data
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
  pred <- as.matrix(pred) # convert vec to matrix in order to add colname
  colnames(pred) <- "taxa"
  prey <- as.matrix(prey)
  colnames(prey) <- "taxa"
  pairs <- list(pred = pred, prey = prey) 
  A.pairs[[web]] <- pairs 
}
names(A.pairs) <- names(obs.A)

# pairs biomass ####
dw.pairs <- NULL
  for(web in 1:length(dw)){
    m.pred <- merge(A.pairs[[web]]$pred,
                    dw[[web]][,1:2],
                    by = "taxa",
                    all.x = T)
    m.prey <- merge(A.pairs[[web]]$prey,
                    dw[[web]][,1:2],
                    by = "taxa",
                    all.x = T)
    m.pairs <- cbind(m.pred[,2], m.prey[,2])
    colnames(m.pairs) <- c("pred", "prey")
    dw.pairs[[web]] <- m.pairs[!is.na(m.pairs[,1]) &
                        !is.na(m.pairs[,2]),]
  }
names(dw.pairs) <- names(obs.A)

#training data ####
# list of lists
training.list <- NULL
for (web in 1:length(dw.pairs)){
    dat <- ldply(dw.pairs[-web])
    training.list[[web]] <- dat
  }
names(training.list) <- names(dw.pairs)

# Gravel model ####
# get parameters based on training lists
# solution without using nested for loops
pars.list <- map(training.list, function(x){
  Bprey = log10(x$prey)
  Bpred = log10(x$pred)
  out <- reg_fn(Bprey, Bpred, quartil = c(0.01, 0.97))
})
# web params ####
# calculate web paramters using get_pars_Niche()
Ball <- map(dw, function (x){
  log10(x$dw)
})
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

match <- map2(obs.A, web.links.inf,
                  match_matr)
obs <- llply(match, function (x){
  x$observed
})
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
# vector of abundances
ab.vec <- llply(dw, function (x){x$no.m2})
#vector of taxa
ab.taxa <- llply(dw, function (x){x$taxa})
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
saveRDS(rel.ab.matr, "relative abundance matrices.RDS")

inf.neutral <- map(threshold, function (x){
  map(rel.ab.matr, rm_neutral, threshold = x)})
names(inf.neutral) <- threshold
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

# AUC logistic model ####
get_auc <- function(observed, inferred){
  require(ROCR)
  y = as.factor(as.numeric(observed))
  x = as.factor(as.numeric(inferred))
  if(length(levels(x))== 1){ 
    auc = NA
    return(auc)
  }
  mod = glm(y ~ x, family = binomial(link = "logit"))
  mod.pred = predict(mod, x, type = "response")
  prob = prediction(mod.pred, y)
  auc = performance(prob, measure = "auc")@y.values[[1]]
  return(auc)
}


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


# local max auc 
local.thresh.neutral <- auc.neutral.df %>% 
  group_by(site) %>%
  top_n(1, wt = auc) %>% 
  .[match(unique(.$site), .$site),]
# plot facet by site 
auc.neutral.df %>% group_by(site) %>% 
  mutate(max.auc = max(na.omit(auc)),
         is.max = auc == max.auc) %>% 
  ggplot(aes(x = thresh,
             y = auc,
             color = is.max)) +
  facet_wrap(~site) +
  geom_point() +
  scale_color_manual(values = c("black", "red"))+
  geom_hline(aes(yintercept = 0.5),
             linetype = "dashed") +
  theme_classic()
# density of thresholds == max.auc
auc.neutral.df %>% group_by(site) %>%
  top_n(1, wt = auc) %>%
  .[match(unique(.$site), .$site),] %>%
ggplot(aes(x = thresh)) +
  geom_density() +
  theme_classic()


global.thresh.neutral <- auc.neutral.df %>%
  group_by(thresh) %>%
  summarize(mean.auc = mean(na.omit(auc))) %>%
  top_n(1, wt = mean.auc)
# plot of global
auc.neutral.df %>% 
  ggplot(aes(x = thresh,
             y = auc)) +
  geom_point() +
  stat_summary(aes(y = auc,group=1),
               fun.y=mean,
               colour="grey",
               geom="line",
               size = 2,
               group= 1) +
  geom_point(data = auc.neutral.df %>%
               filter(thresh > -3.9, thresh < -3.7),
             aes(x = thresh, y = auc), color = "red")+
  theme_classic()

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

# local max auc 
local.thresh.nn <- auc.niche.neutral.df %>% group_by(site) %>%
  top_n(1, wt = auc) %>% 
  .[match(unique(.$site), .$site),]
# plot facet by site 
auc.niche.neutral.df %>% group_by(site) %>% 
  mutate(max.auc = max(na.omit(auc)),
         is.max = auc == max.auc) %>% 
  ggplot(aes(x = thresh,
             y = auc,
             color = is.max)) +
  facet_wrap(~site) +
  geom_point() +
  scale_color_manual(values = c("black", "red"))+
  geom_hline(aes(yintercept = 0.5),
             linetype = "dashed") +
  theme_classic()
# density of thresholds == max.auc
auc.niche.neutral.df %>% group_by(site) %>%
  top_n(1, wt = auc) %>%
  .[match(unique(.$site), .$site),] %>%
  ggplot(aes(x = thresh)) +
  geom_density() +
  theme_classic()

# global max auc
global.thresh.nn <- auc.niche.neutral.df %>%
  group_by(thresh) %>%
  summarize(mean.auc = mean(na.omit(auc))) %>% 
  top_n(1, wt = mean.auc)
# plot of global
auc.niche.neutral.df %>% 
  ggplot(aes(x = thresh,
             y = auc)) +
  geom_point() +
  stat_summary(aes(y = auc,group=1),
               fun.y=mean,
               colour="grey",
               geom="line",
               size = 2,
               group= 1) +
  geom_point(data = auc.niche.neutral.df %>%
               filter(thresh > -3.6, thresh < -3.5),
             aes(x = thresh, y = auc), color = "red")+
  theme_classic()



# TSS ####
# working with neutral abundance threshold 1.5e-04
# inf.neutral[[18]]
neutral <- inf.neutral[[18]]

# TSS initial ####
tss.initial <- ldply(map2(obs, inf,
                    get_tss))
tss.initial.mean <- mean(tss.initial$V1)
# TSS niche ####
tss.niche <- ldply(pmap(list(obs = obs,
                       inf = inf.niche),
                  get_tss))
tss.niche.meean <- mean(tss.niche$V1)

# TSS neutral ####
tss.neutral <- ldply(pmap(list(obs = obs,
                         inf = neutral),
                    get_tss))
tss.neutral.mean <- mean(tss.neutral$V1)


# neutral and niche forbidden ####
# 1e-8
neutral.niche <- inf.niche.neutral[[1]]
tss.niche.neutral <- ldply(
  pmap(list(obs = obs,
            inf = neutral.niche),
       get_tss))
tss.nn.mean <- mean(tss.niche.neutral$V1)



# TSS for niche + local neutral ####
local_tss <- function (n, inf){
  x <- sapply(inf, function (web) web[n])
  names(x) <- threshold
  out <- ldply(map(x,  get_tss, obs = obs[[n]]))
  out
}


local.tss.thresh <- NULL
for(i in 1:length(obs)){
  local.tss.thresh[[i]] <- local_tss(i, inf = inf.niche.neutral)
  names(local.tss.thresh[[i]]) <-c("thresh", "tss") 
}
names(local.tss.thresh) <- names(obs)
local.tss.thresh <- ldply(local.tss.thresh) %>%
  group_by(.id) %>%
  mutate(max.tss = max(na.omit(tss)), 
         is.max = tss == max.tss)
  
ggplot(local.tss.thresh, 
       aes(x = log10(as.numeric(thresh)),
                 y = tss,
                 color = .id,
                 size = is.max)) +
  scale_size_manual(values = c(1, 5)) +
  geom_point() +
  stat_smooth(aes(x = log10(as.numeric(thresh)),
                  y = tss, color = .id),
              alpha = 0, inherit.aes = F)+
  theme_classic()


# # total abundance
# tot.ab <- ldply(sapply(dw, function (x) sum(x$no.m2)))
# tot.ab <- left_join(tot.ab, local.thresh.nn[,c(2:3)], by = c(".id" = "site"))
# ggplot(tot.ab, aes(x = log10(V1), y = thresh)) +
#   geom_point()+
#   stat_smooth(method = "lm")
# 
# summary(lm(thresh ~ log10(V1), data = tot.ab))
# # higher abundance = smaller threshold
# # when you have more individuals, need to forbid links at smaller cross products

# table of auc, tss, threshold

write_csv(data.frame(inference =
        c("Initial", "Niche", "Neutral", "Niche + Neutral"),
        AUC = c(auc.init.mean, auc.niche.mean,
                as.double(global.thresh.neutral[2]),
                as.double(global.thresh.nn[2])),
        TSS = c(tss.initial.mean, tss.niche.mean,
                tss.neutral.mean, tss.nn.mean),
        Threshold = c("NA", "NA",
                      10^as.double(global.thresh.neutral[1]),
        10^as.double(global.thresh.nn[1]))),
        "Mean AUC and TSS trait matching.csv")
