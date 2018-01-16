# fish abundance "correction"
# 16 Jan 2017
# Justin Pomeranz
# jfpomeranz@gmail.com

# Pruning neutrally forbidden links in trait matching inferences were resulting in low predictive power (low AUC, TSS). 
# upon examination, it was discovered that fish links were first to be removed because. 
# because they are relatively rare (e.g. <= 1 per m2), and inverts are abundant (e.g. on the scale of 10's - 1000's per m2), their cross products of relative abundances were very low. 
# however, fish have most of the observed links, and so excluding them led to low explanatory ability in feeding inferences. 
# possible explanantions: 
# 1) fish are active foragers that sample >>1m^2 in their daily lives.
# 2) fish individuals segregate into optimal feeding niches / habitats, and when sampling diets on a population scale, you are actually sampling numerous micrehabitats / niches
# 3) fish stomachs are large and can hold more prey items
# 4) fish gut contents are easier to ID, so invert stomachs are "undersampled" even though an attempt was made to allocate equivalent sample sizes. 


# libraries
library(plyr)
library(tidyverse)
# functions
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
# calc AUC using logistic model
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
# fish abundance "correction"
f_ab_corr <- function(Nij, taxa, cf){
  for(f in which(colnames(Nij) %in% taxa)){
    Nij[,f] <- Nij[,f]*cf
  }
  Nij
}

# data
# observed adjacency matrices, matched to inferences
obs <- readRDS("observed matrices matched to inferred.RDS")
# initial inferred matrices 
inf <- readRDS("Initial trait matching inference.RDS")
# relative abundance cross product matrices
rel.ab.matr <- readRDS("relative abundance matrices.RDS")
# niche forbidden
taxa.forbid <- c("Amphipoda", "Atalophlebioides", "Austroclima", "Austrosimulium", "Blephariceridae", "Coloburiscus", "Deleatidium", "Nesameletus","Oligochaetae", "Ostracoda", "Oxyethira", "Potamopyrgus", "Zephlebia")

f.vec <- c("Salmo", "Galaxias", "Anguilla", "Gobiomorpus")
threshold2 <- c(
  1.0e-08, 1.5e-8, 3.0e-08, 5.9e-08,
  1.0e-07, 1.5e-7, 3.0e-07, 5.9e-07,
  1.0e-06, 1.5e-6, 3.0e-06, 5.9e-06,
  1.0e-05, 1.5e-5, 3.0e-05, 5.9e-05,
  1.0e-04, 1.5e-4, 3.0e-04, 5.9e-04,
  1.0e-03, 1.5e-3, 3.0e-03, 5.9e-03,
  1.0e-02, 1.5e-2, 3.0e-02, 5.9e-02)
cf <- c(10^seq(from = 1, to = 4))
# correction factor ####
# # takes forever to run, commented out ####
# # examine how different correction factors influence inferences
# auc.cf <- NULL
# system.time(
#   for(c in 1:length(cf)){
#     auc.neutral <- NULL
#     for(w in 1:length(inf)){
#       auc.web <- NULL
#       for(t in 1:length(threshold2)){
#         N = f_ab_corr(Nij = rel.ab.matr[[w]],
#                       taxa = f.vec,
#                       cf = cf[c])
#         Nprime = rm_neutral(N, threshold2[t])
#         Nprime = Nprime * inf[[w]]
#         auc.web[[t]] = get_auc(obs[[w]], Nprime)
#       }
#       names(auc.web) <- as.character(threshold2)
#       auc.neutral[[w]] <- auc.web
#     }
#     names(auc.neutral) <- names(obs)
#     auc.cf[[c]] <- auc.neutral
#   }
# )
# names(auc.cf) <- as.character(cf)
# 
# auc.cf <- llply(auc.cf, function (x){
#   out = data.frame(auc = flatten_dbl(x),
#                    thresh = log10(as.numeric(threshold2)),
#                    site = rep(names(obs),
#                               each = length(threshold2)),
#                    stringsAsFactors = FALSE)})
# ldply(auc.cf) %>%
#   ggplot(aes(x = thresh,
#              y = auc, color = .id)) +
#   geom_point() +
#   facet_wrap(~site) +
#   stat_smooth(alpha = 0)+
#   theme_classic()
# 
# ldply(auc.cf) %>%
#   group_by(.id, thresh) %>%
#   summarize(mean.auc = mean(na.omit(auc))) %>%
#   arrange(desc(mean.auc), .id)

# fish corrected ####
# fish relative abundance * 1000
rel.ab.fish <- map(rel.ab.matr,
                   f_ab_corr,
                   taxa = f.vec,
                   cf = 1000)
saveRDS(rel.ab.fish, "rel ab fish x 1000.RDS")

# local fish neutral ####
fish.neutral.list <- map(threshold2, function (x){
  map(rel.ab.fish, rm_neutral, threshold = x)})
names(fish.neutral.list) <- threshold2
fish.neutral.list <- map(fish.neutral.list, function (x){
  map2(x, inf, ~.x*.y)
})

# local threshold for fish abundance correction
local_tss <- function (n, inf){
  x <- sapply(inf, function (web) web[n])
  names(x) <- threshold2
  out <- ldply(map(x,  get_tss, obs = obs[[n]]))
  out
}

# calculate tss for each threshold by site
local.tss.f.ab <- NULL
for(i in 1:length(obs)){
  local.tss.f.ab[[i]] <- local_tss(i, inf = fish.neutral.list)
  names(local.tss.f.ab[[i]]) <-c("thresh", "tss") 
}
names(local.tss.f.ab) <- names(obs)

# add a max.tss (numerical) and an is.max (logical) column
local.tss.f.ab <- ldply(local.tss.f.ab) %>%
  group_by(.id) %>%
  mutate(max.tss = max(na.omit(tss)), 
         is.max = tss == max.tss)
# calculate the threshold which gives highest mean TSS
global.thresh.ab <- local.tss.f.ab %>%
  group_by(thresh) %>%
  summarize(mean.tss = mean(na.omit(tss))) %>%
  top_n(1, wt = mean.tss)

# local fish neutral + niche ####
#local remove neutral forbidden
fish.nn.list <- map(fish.neutral.list, function (x){
  map(x, rm_niche, taxa = taxa.forbid)})

# calculate TSS for each threshold by site
local.tss.f.nn <- NULL
for(i in 1:length(obs)){
  local.tss.f.nn[[i]] <- local_tss(i, inf = fish.nn.list)
  names(local.tss.f.nn[[i]]) <-c("thresh", "tss") 
}
names(local.tss.f.nn) <- names(obs)

# add a max.tss (numerical) and an is.max (logical) column
local.tss.f.nn <- ldply(local.tss.f.nn) %>%
  group_by(.id) %>%
  mutate(max.tss = max(na.omit(tss)), 
         is.max = tss == max.tss)


# global max tss
global.thresh.nn <- local.tss.f.nn %>% group_by(thresh) %>%
  summarize(mean.tss = mean(na.omit(tss))) %>% top_n(1, wt = mean.tss)



#auc fish ####
# Neutral ####
f.auc.neutral <- NULL
for(web in 1:length(obs)){
  auc.web <- NULL
  for(t in 1:length(fish.neutral.list)){
    auc.web[[t]] <- get_auc(obs[[web]],
                            fish.neutral.list[[t]][[web]])
  }
  f.auc.neutral[[web]] <- auc.web
}

# turn results into a df
f.auc.neutral.df <- data.frame(auc =
                        flatten_dbl(f.auc.neutral),
                              thresh = 
                        log10(as.numeric(threshold2)),
                              site = 
                        rep(names(obs),
                            each = length(threshold2)),
                        stringsAsFactors = FALSE)
# global threshold == max mean auc
global.f.neutral <- f.auc.neutral.df %>%
  group_by(thresh) %>%
  summarize(mean.auc = mean(na.omit(auc))) %>%
  top_n(1, wt = mean.auc)


# niche + neutral ####
f.auc.nn <- NULL
for(web in 1:length(obs)){
  auc.web <- NULL
  for(t in 1:length(fish.nn.list)){
    auc.web[[t]] <- get_auc(obs[[web]],
                            fish.nn.list[[t]][[web]])
  }
  f.auc.nn[[web]] <- auc.web
}
# turn results into a df
f.auc.nn.df <- data.frame(auc =
                    flatten_dbl(f.auc.nn),
                          thresh =
                    log10(as.numeric(threshold2)),
                          site = 
                    rep(names(obs),
                        each = length(threshold2)),
                    stringsAsFactors = FALSE)
# global threshold = max mean.auc
global.f.nn <- f.auc.nn.df %>%
  group_by(thresh) %>%
  summarize(mean.auc = mean(na.omit(auc))) %>%
  top_n(1, wt = mean.auc)

# fn & fp ####
false_prop <- function(obs, inf){
  stopifnot(dim(obs) == dim(inf), 
            identical(colnames(obs), colnames(inf)))
  
  # subtract inferred from observed
  minus <- obs - inf
  # multiply observed and inferred
  multiply <- obs * inf
  # change values of 1 in multiplied matrix to 2
  multiply[multiply == 1] <- 2
  # add minus and multiply matrices
  prediction <- minus + multiply
  # prediction outcome matrix now has all 4 possibilities repreented as different integer values
  # 2 = true positive (a); links both obserevd & predicted
  # -1 = false positive (b); predicted but not observed
  # 1 = false negative (c); observed but not predicted
  # 0 = true negative (d); not predicted, not observed
  tss.vars <- data.frame(
    a = length(prediction[prediction==2]),
    b = length(prediction[prediction==-1]), 
    c = length(prediction[prediction==1]), 
    d = length(prediction[prediction==0]),
    S = ncol(prediction)
  )
  # calculate TSS
  # TSS = (a*d - b*c)/((a+c)*(b+d))
  # also decompose TSS to see proportion of positive and negative predictions
  # abar = proportion of links correctly predicted
  # dbar = proportion of links correctly not predicted
  # bc = proportion of links incorrectly predicted
  tss <- tss.vars %>%
    transmute(fp = b / S**2,
              fn = c / S**2)
  tss
}
f.neutral <- fish.neutral.list[[18]]
saveRDS(f.neutral, "fish correction neutral.RDS")
f.nn <- fish.nn.list[[17]]
saveRDS(f.nn, "fish correction neutral + Niche.RDS")

neutral.false <- ldply(map2(obs, f.neutral, false_prop)) %>% 
  summarize(mean.fp = mean(fp), sd.fp = sd(fp),
            mean.fn = mean(fn), sd.fn = sd(fn))
nn.false <- ldply(map2(obs, f.nn, false_prop)) %>% 
  summarize(mean.fp = mean(fp), sd.fp = sd(fp),
            mean.fn = mean(fn), sd.fn = sd(fn))

false.tab <- rbind(neutral.false, nn.false)
false.tab$inference <- c("Neutral", "Niche + Neutral")

false.tab <- gather(false.tab, "var", "val", 1:4) %>%
  spread(var, val)
# AUC TSS threshold table
write_csv(data.frame(Inference = c("Fish corrected neutral", "Fish corrected niche + Neutral"),
           AUC = c(as.double(global.f.neutral[2]),
           as.double(global.f.nn[2])),
           TSS = c(as.double(global.thresh.ab[2]),
                   as.double(global.thresh.nn[2])),
           Threshold = c(as.double(global.thresh.ab[1]),
                         as.double(global.thresh.nn[1])),
           mean.fp = false.tab$mean.fp,
           sd.fp = false.tab$sd.fp,
           mean.fn = false.tab$mean.fn,
           sd.fn = false.tab$sd.fn),
          "Fish corrected Mean AUC and TSS trait matching.csv")

# # plots ####
# # Fish correction, neutral prune ####
# # tss ~ log10(threshold)
# # each site gets a color / line
# ggplot(local.tss.f.ab, 
#        aes(x = log10(as.numeric(thresh)),
#            y = tss,
#            color = .id,
#            size = is.max)) +
#   scale_size_manual(values = c(1, 5)) +
#   geom_point() +
#   stat_smooth(aes(x = log10(as.numeric(thresh)),
#                   y = tss, color = .id),
#               alpha = 0, inherit.aes = F)+
#   theme_classic()
# 
# # fish correction, neutral + niche prune ####
# # TSS ~ log10(threshold)
# # each site gets a line / color
# ggplot(local.tss.f.nn, 
#        aes(x = log10(as.numeric(thresh)),
#            y = tss,
#            color = .id,
#            size = is.max)) +
#   scale_size_manual(values = c(1, 5)) +
#   geom_point() +
#   stat_smooth(aes(x = log10(as.numeric(thresh)),
#                   y = tss, color = .id),
#               alpha = 0, inherit.aes = F)+
#   theme_classic()
# # TSS ~ log10(threshold)
# # each facet isone site
# # max tss is colored red
# ggplot(local.tss.f.nn, aes(x = log10(as.numeric(thresh)),
#                            y = tss,
#                            color = is.max)) +
#   facet_wrap(~.id) +
#   geom_point() +
#   scale_color_manual(values = c("black", "red"))+
#   theme_classic()
# # density of thresholds which result in TSS
# local.tss.f.nn %>% group_by(.id) %>%
#   top_n(1, wt = tss) %>%
#   .[match(unique(.$.id), .$.id),] %>%
#   ggplot(aes(x = log10(as.numeric(thresh)))) +
#   geom_density() +
#   theme_classic()
# # tss~log10(threshold)
# # gray line = mean tss ~ threshold
# # red points = threshold which gives highest mean TSS
# local.tss.f.nn %>% 
#   ggplot(aes(x = log10(as.numeric(thresh)),
#              y = tss)) +
#   geom_point() +
#   stat_summary(aes(y = tss,group=1),
#                fun.y=mean,
#                colour="grey",
#                geom="line",
#                size = 2,
#                group= 1) +
#   geom_point(data = local.tss.f.nn %>%
#                filter(thresh == "0.00015"),
#              aes(x = log10(as.numeric(thresh)),
#                  y = tss),
#              color = "red")+
#   theme_classic()
# 
