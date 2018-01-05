# Trait-matching
# 31 Dec 2017
# Justin Pomeranz
#jfpomeranz@gmail.com

# paramterize Niche model using method presented in Gravel et al. 2013 Methods in Eco Evo. Inferring food web structure from predator-prey body size relationships

library(plyr)
library(tidyverse)

# gravel functions
# from Gravel et al. 2013 supplementary information
source("gravel_functions.R")

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
# match observed and inferred matrices
match_matr <- function (obs, inf){
  # colnames(inf) have already been size-sorted
  index = intersect(rownames(inf), rownames(obs))
  obs = obs[index, index]
  inf = inf[index, index]
  list(observed = obs, inferred = inf)
}
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

# TSS ####
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

rm_niche <- function(inf, taxa){
  for(name in (
    colnames(inf)[colnames(inf) %in% taxa])){
        inf[,name] <- 0
  }
  inf
}
# # function to calc TSS from matrices matched by row/colnames
# # be careful to select appropriate obs and inf
# match_matr_tss <- function (obs, inf){
#   # colnames(inf) have already been size-sorted
#   index = intersect(rownames(inf), rownames(obs))
#   obs = obs[index, index]
#   inf = inf[index, index]
#   result = get_tss(obs, inf)
#   result
# }

# step1, biomass inference
tss.initial <- map2(obs, inf,
                  get_tss)
# niche forbidden
taxa.forbid <- c("Amphipoda", "Atalophlebioides", "Austroclima", "Austrosimulium", "Blephariceridae", "Coloburiscus", "Deleatidium", "Nesameletus", "Ostracoda", "Oxyethira", "Potamopyrgus", "Zephlebia")

inf.niche <- map(inf,
                 rm_niche,
                 taxa = taxa.forbid)
# save niche pruned trait matching ####
saveRDS(inf.niche, "Niche pruned trait matching inference.RDS")

tss.niche <- pmap(list(obs = obs,
                       inf = inf.niche),
                       get_tss)
# neutral forbidden ####
# vector of abundances
ab.vec <- llply(dw, function (x){x$dw})
#vector of taxa
ab.taxa <- llply(dw, function (x){x$taxa})
# threshold vector
threshold <- c(5.9e-13,
               1.0e-12, 1.5e-12, 3.0e-12, 5.9e-12,
               1.0e-11, 1.5e-11, 3.0e-11, 5.9e-11,
               1.0e-10, 1.5e-10, 3.0e-10, 5.9e-10,
               1.0e-09, 1.5e-9, 3.0e-09, 5.9e-09,
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

# confusion matrix
source("adj_conf_matrix function.R")
#adj_conf_matrix(observed, inferred)

conf.neutral <- map(inf.neutral, function (x){
  pmap(list(obs = obs,
            inf = x),
       adj_conf_matrix)})

test <- map(inf.neutral, function (x){
  adj_conf_matrix(observed = flatten_dbl(obs), 
                  inferred = flatten_dbl(x))
})
test <- ldply(test)
ggplot(test, aes(x = FPR, y = TPR)) +
  geom_point() +
  geom_abline(slope = 1)

trapz(test$FPR, test$TPR) 
test %>%
  summarise(AUC = -sum(diff(FPR) * (head(TPR,-1)+tail(TPR,-1)))/2) 

ROC <- llply(conf.neutral, function (x){
  ldply(x)[, c(".id", "TPR", "FPR")]})

ROC <- ROC %>% llply( function (x){
  names(x)[1] <- "site"; x
})
ROC <- ldply(ROC)
# ROC by site
ggplot(ROC, aes(x = FPR, y = TPR))+
  facet_wrap(~site)+
  geom_point(aes(color = .id)) +
  geom_line()

ggplot(ROC, aes(x = FPR, y = TPR, color = site))+
  geom_point() +
  geom_line() +
  geom_abline(slope = 1, color = "black")+
  theme_classic()
ggplot(ROC, aes(x = FPR, y = TPR))+
  facet_wrap(~.id) +
  geom_point() +
  geom_abline(slope = 1, color = "black")+
  theme_classic()





ROC.summ <- ROC %>% group_by(.id) %>%
  summarise(TPR = mean(TPR), FPR = mean(FPR))
ggplot(ROC.summ, aes(x = FPR, y = TPR))+
  geom_point(aes(color = .id)) +
  geom_line() +
  geom_abline(slope = 1)

ROC.thresh <- ROC %>%
  group_by(.id) %>%
  arrange(TPR) %>%
  summarise(AUC = sum(diff(FPR) * (head(TPR,-1)+tail(TPR,-1)))/2) 

ggplot(ROC.thresh, aes(x = log10(as.numeric(.id)),
                       y = AUC)) +
  geom_point()

require(pracma)
trapz(ROC.summ$FPR, ROC.summ$TPR)

roc.trapz <- ROC %>%
  group_by(.id) %>%
  summarise(AUC = trapz(FPR, TPR)) 
ggplot(roc.trapz, aes(x = log10(as.numeric(.id)),
                      y = -AUC)) +
  geom_point()

# 
# # ROC by threshold
# ROC %>%
#   group_by(.id) %>%
#   summarize(TPR = mean(TPR), FPR = mean(FPR)) %>%
#   ggplot(aes(x = FPR, y = TPR))+
#   geom_point(aes(color = .id), size = 2) +
#   geom_line()
# 
# 
# AUC <- llply(conf.neutral, function (x){
#   ldply(x)[, c(".id", "fn", "fp")]})
# AUC <- AUC %>% llply( function (x){
#   names(x)[1] <- "site"; x
# })
# AUC <- ldply(AUC)
# 
# 
# sutton <- NULL
# for (t in 1:length(inf.neutral)){
#   temp <- as.data.frame(as.numeric(inf.neutral[[t]][[16]]))
#   names(temp) <- "predicted"
#   temp$thresh <- as.numeric(names(inf.neutral)[t])
#   sutton <- rbind(sutton, temp)
# }
# sutton$obs <- as.numeric(obs[[16]])
# glm_obj <- glm(obs~., data = sutton, family = binomial(logit))
# pred_dat <- predict(glm_obj, newdata = sutton, type = "response")
# pred_dat <- data.frame(pred = pred_dat, thresh = sutton$thresh)
# pred_dat %>% group_by(thresh) %>% summarise(min(pred), max(pred)) %>% as.data.frame()
# roc_obj <- roc(sutton$obs, pred_dat$pred, plot = T)
# 
# nprime <- NULL
# for (web in 1:length(obs)){
#   temp <- as.data.frame(as.numeric(
#     inf.neutral[[10]][[web]]))
#   names(temp) <- "predicted"
#   temp$thresh <- as.numeric(names(inf.neutral)[10])
#   nprime <- rbind(nprime, temp)
# }
# 
# obs.vector <- llply(obs, function (x){as.numeric(x)})
# nprime$obs <- unlist(obs.vector, recursive = F)
# 
# glm_obj <- glm(~thresh, data = nprime, family = binomial(logit))
# pred_dat <- predict(glm_obj, newdata = nprime, type = "response")
# pred_dat <- data.frame(pred = pred_dat, thresh = nprime$thresh)
# pred_dat %>% group_by(thresh) %>% summarise(min(pred), max(pred)) %>% as.data.frame()
# roc_obj <- roc(sutton$obs, pred_dat$pred, plot = T)








# TSS neutral ####
tss.neutral <- map(inf.neutral, function (x){
  pmap(list(obs = obs,
            inf = x),
       get_tss)
})
tss.neutral <- ldply(flatten(tss.neutral))
tss.neutral$threshold <- rep(threshold, each = 17)
ggplot(tss.neutral, aes(x = log10(threshold), y = V1, color = .id)) +
  geom_point() +
  stat_smooth(alpha = 0.2) +
  ggtitle("Neutral forbidden") +
  theme_classic()

tss.neutral %>% group_by(.id) %>% top_n(1,wt = V1) %>% mutate(log10(threshold))

# neutral and niche forbidden
tss.niche.neutral <- map(inf.neutral, function (x){
  pmap(list(obs = obs,
            inf = map(x,rm_niche, taxa = taxa.forbid)),
       get_tss)
})

tss.niche.neutral <- ldply(flatten(tss.niche.neutral))
tss.niche.neutral$threshold <- rep(threshold, each = 17)
ggplot(tss.niche.neutral, aes(x = log10(threshold), y = V1, color = .id)) +
  geom_point() +
  stat_smooth(alpha = 0.2)+
  ggtitle("Niche and Neutral Forbidden") +
  theme_classic()
tss.niche.neutral %>% group_by(.id) %>% top_n(1,wt = V1) %>% mutate(log10(threshold))




# still need to figure this out ####
# fish abundance "correction"
f.vec <- c("Salmo", "Galaxias", "Anguilla", "Gobiomorpus")
f_ab_corr <- function(Nij, taxa, cf){
  for(f in which(colnames(Nij) %in% taxa)){
  Nij[,f] <- Nij[,f]*cf
  }
  Nij
}

rel.ab.fish <- map(rel.ab.matr,
            f_ab_corr, taxa = f.vec, cf = 1e5)
# working out how to make list of diff cf's
test1 <- map(cf, function (x) {map(Nij.list, f_ab_corr, taxa = taxa, cf = x)})

rm_neutral(f_ab_corr(Nij, taxa = "Salmo", cf = 10),
           threshold = 0.01)

test1 <- map(cf, function (x){
  pmap(list(Nij = map(Nij.list, f_ab_corr,
                      taxa = "Salmo",
                      cf = x),
            threshold = threshold),
       rm_neutral)})
test2 <- map(test1, function (x){
  pmap(list(observed = obs,
            inferred = x),
       get_tss)
})    
tss <- NULL
for(c in 1:length(cf)){
  temp.thresh <- NULL
  for(t in 1:length(threshold)){
    temp.web <- NULL
    for(web in 1:length(obs)){
      observed = obs[[web]]
      inferred = test1[[c]][[t]][[web]]
      temp.web[[web]] = get_tss(observed, inferred)
    }
    temp.thresh[[t]] <- temp.web
  }
  tss[[c]] <- temp.thresh
}
     
     
     
     
# fish neutral
fish.neutral <- map(threshold, function (x){
  map(rel.ab.fish, rm_neutral, threshold = x)})
# tss fish neutral
tss.fish.neutral <- map(fish.neutral, function (x){
  pmap(list(obs = obs.A,
            inf = x),
       match_matr_tss)
})

map(test1, map2, threshold, function(x){
  map(rm_neutral, threshold = x)
})


# invoke_map???? ####
f <- c("rm_neutral")
param <- list(
  list(Nij = Nij.list,
       threshold = threshold))
test <- invoke_map(f, param)



# example from R4ds ####
f <- c("runif", "runif", "rnorm", "rpois")
param <- list(
  list(min = -1, max = 1),
  list(min = 2, max = 10),
  list(sd = 5), 
  list(lambda = 10)
)
invoke_map(f, param, n = 2)

tss.fish.neutral <- ldply(flatten(tss.fish.neutral))
tss.fish.neutral$threshold <- rep(threshold, each = 17)
ggplot(tss.fish.neutral, aes(x = log10(threshold), y = V1, color = .id)) +
  geom_point() +
  stat_smooth(alpha = 0.2) +
  ggtitle("Neutral forbidden, Fish correction") +
  theme_classic()
tss.fish.neutral %>% group_by(.id) %>% top_n(1,wt = V1) %>% mutate(log10(threshold))

# fish + niche + neutral
tss.fish.niche.neutral <- map(fish.neutral,
                              function (x){
  pmap(list(obs = obs.A,
            inf = map(x,rm_niche, taxa = taxa.forbid)),
       match_matr_tss)
})
tss.fish.niche.neutral <- ldply(flatten(tss.fish.niche.neutral))
tss.fish.niche.neutral$threshold <- rep(threshold, each = 17)
ggplot(tss.fish.niche.neutral, aes(x = log10(threshold), y = V1, color = .id)) +
  geom_point() +
  stat_smooth(alpha = 0.2) +
  ggtitle("Neutral + Niche forbidden, Fish correction") +
  theme_classic()

tss.fish.niche.neutral %>% group_by(.id) %>% top_n(1,wt = V1) %>% mutate(log10(threshold))

tss.fish.niche.neutral %>% group_by(.id) %>% top_n(1,wt = V1) %>% .$V1 %>% mean
