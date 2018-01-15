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
taxa.forbid <- c("Amphipoda", "Atalophlebioides", "Austroclima", "Austrosimulium", "Blephariceridae", "Coloburiscus", "Deleatidium", "Nesameletus","Oligochaetae", "Ostracoda", "Oxyethira", "Potamopyrgus", "Zephlebia")

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
threshold <- c(#5.9e-13,
               # 1.0e-12, 1.5e-12, 3.0e-12, 5.9e-12,
               # 1.0e-11, 1.5e-11, 3.0e-11, 5.9e-11,
               # 1.0e-10, 1.5e-10, 3.0e-10, 5.9e-10,
               # 1.0e-09, 1.5e-9, 3.0e-09, 5.9e-09,
               1.0e-08, 1.5e-8, 3.0e-08, 5.9e-08,
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

inf.neutral <- map(threshold, function (x){
  map(rel.ab.matr, rm_neutral, threshold = x)})
names(inf.neutral) <- threshold

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

auc.neutral.df <- data.frame(auc =
                               flatten_dbl(auc.neutral),
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
               filter(thresh > -3.6, thresh < -3.5),
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
global.thresh.nn <- auc.niche.neutral.df %>% group_by(thresh) %>%
  summarize(mean.auc = mean(na.omit(auc))) %>% top_n(1, wt = mean.auc)
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
# working with neutral abundance threshold 3e-04
# inf.neutral[[35]]
neutral <- inf.neutral[[19]]

# TSS initial ####
tss.initial <- ldply(map2(obs, inf,
                    get_tss))
mean(tss.initial$V1)
# TSS niche ####
tss.niche <- ldply(pmap(list(obs = obs,
                       inf = inf.niche),
                  get_tss))
mean(tss.niche$V1)

# TSS neutral ####
tss.neutral <- ldply(pmap(list(obs = obs,
                         inf = neutral),
                    get_tss))
mean(tss.neutral$V1)


# neutral and niche forbidden ####
# 1.5e-4
neutral.niche <- inf.niche.neutral[[18]]
tss.niche.neutral <- ldply(
  pmap(list(obs = obs,
            inf = neutral.niche),
       get_tss))
mean(tss.niche.neutral$V1)



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




# total abundance
tot.ab <- ldply(sapply(dw, function (x) sum(x$no.m2)))
tot.ab <- left_join(tot.ab, local.thresh.nn[,c(2:3)], by = c(".id" = "site"))
ggplot(tot.ab, aes(x = log10(V1), y = thresh)) +
  geom_point()+
  stat_smooth(method = "lm")

summary(lm(thresh ~ log10(V1), data = tot.ab))
# higher abundance = smaller threshold
# when you have more individuals, need to forbid links at smaller cross products






























# fish abundance ####
# fish abundance "correction"
f_ab_corr <- function(Nij, taxa, cf){
  for(f in which(colnames(Nij) %in% taxa)){
    Nij[,f] <- Nij[,f]*cf
  }
  Nij
}
f.vec <- c("Salmo", "Galaxias", "Anguilla", "Gobiomorpus")
threshold2 <- threshold
cf <- c(10^seq(from = 1, to = 4))
auc.cf <- NULL
system.time(
  for(c in 1:length(cf)){
    auc.neutral <- NULL
    for(w in 1:length(inf)){
      auc.web <- NULL
      for(t in 1:length(threshold2)){
        N = f_ab_corr(Nij = rel.ab.matr[[w]],
                      taxa = f.vec,
                      cf = cf[c])
        Nprime = rm_neutral(N, threshold2[t])
        auc.web[[t]] = get_auc(obs[[w]], Nprime)
      }
      names(auc.web) <- as.character(threshold2)
      auc.neutral[[w]] <- auc.web
    }
    names(auc.neutral) <- names(obs)
    auc.cf[[c]] <- auc.neutral
  }
)
names(auc.cf) <- as.character(cf)

auc.cf <- llply(auc.cf, function (x){
  out = data.frame(auc = flatten_dbl(x),
                   thresh = log10(as.numeric(threshold2)),
                   site = rep(names(obs),
                              each = length(threshold2)),
                   stringsAsFactors = FALSE)})
ldply(auc.cf) %>%
ggplot(aes(x = thresh,
           y = auc, color = .id)) +
geom_point() +
facet_wrap(~site) +
stat_smooth(alpha = 0)+
theme_classic()

ldply(auc.cf) %>%
group_by(.id, thresh) %>%
summarize(mean.auc = mean(na.omit(auc))) %>%
arrange(desc(mean.auc), .id)
# fish.tss
# cf = 100, threshold = 1e-04
rel.ab.fish <- map(rel.ab.matr,
                   f_ab_corr,
                   taxa = f.vec,
                   cf = 100)
fish.neutral <- map(rel.ab.fish,
                    rm_neutral,
                    threshold = 1e-04)
tss.fish.neutral <- ldply(
  map2(obs,
       fish.neutral,
       get_tss))
tss.fish.neutral$V1 %>% mean

fish.neut.niche <- map(fish.neutral,
                       rm_niche,
                       taxa = taxa.forbid)
tss.fish.n.n <- ldply(map2(obs,
                           fish.neut.niche,
                           get_tss))
tss.fish.n.n$V1 %>% mean


# local fish neutral ####
fish.neutral.list <- map(threshold, function (x){
  map(rel.ab.fish, rm_neutral, threshold = x)})
names(fish.neutral.list) <- threshold


# local neutral ####
# # local threshold for fish abundance correction
# local_tss_f <- function (n){
#   x <- sapply(fish.neutral.list, function (web) web[n])
#   names(x) <- threshold
#   out <- ldply(map(x,  get_tss, obs = obs[[n]]))
#   out
# }
  
local.tss.f.ab <- NULL
for(i in 1:length(obs)){
  local.tss.f.ab[[i]] <- local_tss(i, inf = fish.neutral.list)
  names(local.tss.f.ab[[i]]) <-c("thresh", "tss") 
}
names(local.tss.f.ab) <- names(obs)
local.tss.f.ab <- ldply(local.tss.f.ab) %>%
  group_by(.id) %>%
  mutate(max.tss = max(na.omit(tss)), 
         is.max = tss == max.tss)

ggplot(local.tss.f.ab, 
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

global.tss.f.ab <- local.tss.f.ab %>% 
  group_by(thresh) %>%
  summarize(mean.tss = mean(tss)) %>%
  arrange(desc(mean.tss))

# local nn ####
#local fish niche neutral
fish.nn.list <- map(fish.neutral.list, function (x){
  map(x, rm_niche, taxa = taxa.forbid)})

local.tss.f.nn <- NULL
for(i in 1:length(obs)){
  local.tss.f.nn[[i]] <- local_tss(i, inf = fish.nn.list)
  names(local.tss.f.nn[[i]]) <-c("thresh", "tss") 
}
names(local.tss.f.nn) <- names(obs)
local.tss.f.nn <- ldply(local.tss.f.nn) %>%
  group_by(.id) %>%
  mutate(max.tss = max(na.omit(tss)), 
         is.max = tss == max.tss)

ggplot(local.tss.f.nn, 
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






# plot facet by site 
ggplot(local.tss.f.nn, aes(x = log10(as.numeric(thresh)),
           y = tss,
           color = is.max)) +
facet_wrap(~.id) +
geom_point() +
scale_color_manual(values = c("black", "red"))+
theme_classic()
# density of thresholds == max.auc
local.tss.f.nn %>% group_by(.id) %>%
  top_n(1, wt = tss) %>%
  .[match(unique(.$.id), .$.id),] %>%
  ggplot(aes(x = log10(as.numeric(thresh)))) +
  geom_density() +
  theme_classic()

# global max tss
global.thresh.nn <- local.tss.f.nn %>% group_by(thresh) %>%
  summarize(mean.tss = mean(na.omit(tss))) %>% top_n(1, wt = mean.tss)
# plot of global
local.tss.f.nn %>% 
  ggplot(aes(x = log10(as.numeric(thresh)),
             y = tss)) +
  geom_point() +
  stat_summary(aes(y = tss,group=1),
               fun.y=mean,
               colour="grey",
               geom="line",
               size = 2,
               group= 1) +
  geom_point(data = local.tss.f.nn %>%
               filter(thresh == "3e-05"),
             aes(x = log10(as.numeric(thresh)),
                 y = tss),
             color = "red")+
  theme_classic()





#auc fish ####
# global max auc
f.auc.neutral <- NULL
for(web in 1:length(obs)){
  auc.web <- NULL
  for(t in 1:length(fish.neutral.list)){
    auc.web[[t]] <- get_auc(obs[[web]],
                            fish.neutral.list[[t]][[web]])
  }
  f.auc.neutral[[web]] <- auc.web
}

f.auc.neutral.df <- data.frame(auc =
                               flatten_dbl(f.auc.neutral),
                             thresh = log10(as.numeric(threshold)),
                             site = rep(names(obs),
                                        each = length(threshold)),
                             stringsAsFactors = FALSE)
global.f.neutral <- f.auc.neutral.df %>%
  group_by(thresh) %>%
  summarize(mean.auc = mean(na.omit(auc))) %>%
  top_n(1, wt = mean.auc)



f.auc.nn <- NULL
for(web in 1:length(obs)){
  auc.web <- NULL
  for(t in 1:length(fish.nn.list)){
    auc.web[[t]] <- get_auc(obs[[web]],
                            fish.nn.list[[t]][[web]])
  }
  f.auc.nn[[web]] <- auc.web
}

f.auc.nn.df <- data.frame(auc =
                                 flatten_dbl(f.auc.nn),
                               thresh = log10(as.numeric(threshold)),
                               site = rep(names(obs),
                                          each = length(threshold)),
                               stringsAsFactors = FALSE)
global.f.nn <- f.auc.nn.df %>%
  group_by(thresh) %>%
  summarize(mean.auc = mean(na.omit(auc))) %>%
  top_n(1, wt = mean.auc)
