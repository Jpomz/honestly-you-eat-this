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


# niche forbidden ####
# niche forbidden
taxa.forbid <- c("Amphipoda", "Atalophlebioides", "Austroclima", "Austrosimulium", "Blephariceridae", "Coloburiscus", "Deleatidium", "Nesameletus","Oligochaetae", "Ostracoda", "Oxyethira", "Potamopyrgus", "Zephlebia")

inf.niche <- map(inf,
                 rm_niche,
                 taxa = taxa.forbid)
# # save niche pruned trait matching ####
# saveRDS(inf.niche, "Niche pruned trait matching inference.RDS")


# neutral forbidden ####
# vector of abundances
ab.vec <- llply(dw, function (x){x$no.m2})
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


# AUC need to re-do with log regre ####

# AUC initial


# AUC niche


# AUC neutral


# calculate area under the curve for each threshold



# AUC Niche + Neutral


# TSS ####
# working with neutral abundance threshold 3e-05
# inf.neutral[[32]]
neutral <- inf.neutral[[33]]
# step1, biomass inference
tss.initial <- ldply(map2(obs, inf,
                    get_tss))
mean(tss.initial$V1)
# TSS niche
tss.niche <- ldply(pmap(list(obs = obs,
                       inf = inf.niche),
                  get_tss))
mean(tss.niche$V1)

# TSS neutral 
tss.neutral <- ldply(pmap(list(obs = obs,
                         inf = neutral),
                    get_tss))
mean(tss.neutral$V1)


# neutral and niche forbidden ####
tss.niche.neutral <- ldply(
  pmap(list(obs = obs,
            inf = inf.niche.neutral[[33]]),
       get_tss))
mean(tss.niche.neutral$V1)


# fish abundance ####
# fish abundance "correction"
f_ab_corr <- function(Nij, taxa, cf){
  for(f in which(colnames(Nij) %in% taxa)){
    Nij[,f] <- Nij[,f]*cf
  }
  Nij
}

f.vec <- c("Salmo", "Galaxias", "Anguilla", "Gobiomorpus")
threshold2 <- threshold[25:37]

cf <- c(10^seq(from = 1, to = 4))
cf.dat <- NULL
system.time(
for(c in 1:length(cf)){
  threshold.dat <- NULL
  for(t in 1:length(threshold2)){
    auc.dat <- NULL
    for(w in 1:length(inf)){
      N = f_ab_corr(Nij = rel.ab.matr[[w]],
                    taxa = f.vec,
                    cf = cf[c])
      Nprime = rm_neutral(N, threshold2[t])
      conf = adj_conf_matrix(obs[[w]],
                             Nprime)[,c(5,10)]
      auc.dat = rbind(auc.dat, conf)
    }
    threshold.dat[[t]] <- auc.dat %>%
      arrange(TPR, FPR) %>%
      summarize(auc = trapz(FPR, TPR))
  }
  names(threshold.dat) <- as.character(threshold2)
  cf.dat[[c]] <- threshold.dat
}
)
names(cf.dat) <- as.character(cf)

cf.dat2 <- cf.dat %>% llply(function (x){
  ldply(x)
})
cf.dat2 <- cf.dat2 %>% llply(function (x){
  names(x)[1] <- "threshold"; x
})
ldply(cf.dat2) %>% 
  ggplot(aes(x = log10(as.numeric(threshold)),
             y = auc, color = .id)) +
  geom_point() +
  stat_smooth(alpha = 0)+
  theme_classic()

ldply(cf.dat2) %>% top_n(10, wt = auc) %>% arrange(.id)

# fish.tss
# cf = 100, threshold = 5.9e-05
rel.ab.fish <- map(rel.ab.matr,
                   f_ab_corr,
                   taxa = f.vec,
                   cf = 10000)
fish.neutral <- map(rel.ab.fish,
                   rm_neutral,
                   threshold = 5.9e-05)
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

adj_conf_matrix(obs[[1]], inf[[1]])[10:11]

ldply(map2(obs, inf, adj_conf_matrix))[,10:11] %>% summarize(fnr = mean(FNR), fpr = mean(FPR))
ldply(map2(obs, inf.niche, adj_conf_matrix))[,10:11] %>% summarize(fnr = mean(FNR), fpr = mean(FPR))
ldply(map2(obs, neutral, adj_conf_matrix))[,10:11] %>% summarize(fnr = mean(FNR), fpr = mean(FPR))
ldply(map2(obs, inf.niche.neutral[[33]], adj_conf_matrix))[,10:11] %>% summarize(fnr = mean(FNR), fpr = mean(FPR))

df <- rbind(
  ldply(map2(obs, inf, adj_conf_matrix)
        )[,10:11] %>%
    summarize(fnr = mean(FNR),
              fpr = mean(FPR),
              fnr.sd = sd(FNR),
              fpr.sd = sd(FPR)),
  ldply(map2(obs, inf.niche,
             adj_conf_matrix))[,10:11] %>%
    summarize(fnr = mean(FNR),
              fpr = mean(FPR),
              fnr.sd = sd(FNR),
              fpr.sd = sd(FPR)),
  ldply(map2(obs, fish.neut.niche,
             adj_conf_matrix))[,10:11] %>%
    summarize(fnr = mean(FNR),
              fpr = mean(FPR),
              fnr.sd = sd(FNR),
              fpr.sd = sd(FPR)))
df$step <- c(1:3)
#df <- 
  gather(df, "var", "val", 1:4) %>%
  separate(var, c("mean", "sd")) 
ggplot(df, aes(x = step, y = fpr,
               color = "red")) +
  geom_point() +
  geom_errorbar(aes(ymax = fpr + fpr.sd,
                    ymin = fpr - fpr.sd)) +
  geom_line() +
  geom_point(aes(x = step, y = fnr,
             color = "black")) +
  geom_errorbar(aes(ymax = fnr + fnr.sd,
                    ymin = fnr - fnr.sd,
                    color = "black")) +
    geom_line(aes(x = step, y = fnr,
                  color = "black")) +
  theme_classic() +
  labs(y = "Proportion of links") +
  scale_colour_manual(name = 'Colour',
      values =c('black'='black','red'='red'),
      labels = c('False negatives',
                 'False positives'))




# logistic model attempt ####
# need to fix get_auc to work with all inf types!!!! ####
get_auc <- function(web, thresh){
  require(ROCR)
  y = as.factor(as.numeric(
    obs[[web]]))
  x = as.factor(as.numeric(
    inf.neutral[[thresh]][[web]]))
  if(length(levels(x))== 1){ 
    auc = NA
    return(auc)
  }
  mod = glm(y ~ x, family =
              binomial(link = "logit"))
  mod.pred = predict(mod, x,
                     type = "response")
  prob = prediction(mod.pred, y)
  auc = performance(prob, measure = "auc")@
    y.values[[1]]
  return(auc)
}
auc.all <- NULL
for(web in 1:length(obs)){
  auc.web <- NULL
  for(t in 1:length(inf.neutral)){
    auc.web[[t]] <- get_auc(web, t)
  }
auc.all[[web]] <- auc.web
}

df <- data.frame(auc = flatten_dbl(auc.all),
                 thresh = log10(as.numeric(threshold)),
                 site = rep(names(obs), each = length(threshold)),
              stringsAsFactors = FALSE)

df %>% group_by(site) %>% 
  mutate(max.auc = max(na.omit(auc)),
         is.max = auc == max.auc) %>% 
  ggplot(aes(x = thresh,
             y = auc,
             color = is.max)) +
  facet_wrap(~site) +
  geom_point() +
  scale_color_manual(values = c("black", "red"))+
  theme_classic()
 
max.auc <- df %>% group_by(site) %>%
  top_n(1, wt = auc) %>%
  .[match(unique(.$site), .$site),]
ggplot(max.auc, aes(x = thresh)) +
  geom_density()

