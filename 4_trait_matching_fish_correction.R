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
# useful functions
source("Inference_MS_functions.R")
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

f.vec <- c("Salmo", "Galaxias", "Anguilla", "Gobiomorphus")
threshold2 <- c(
  1.0e-08, 1.5e-8, 3.0e-08, 5.9e-08,
  1.0e-07, 1.5e-7, 3.0e-07, 5.9e-07,
  1.0e-06, 1.5e-6, 3.0e-06, 5.9e-06,
  1.0e-05, 1.5e-5, 3.0e-05, 5.9e-05,
  1.0e-04, 1.5e-4, 3.0e-04, 5.9e-04,
  1.0e-03, 1.5e-3, 3.0e-03, 5.9e-03,
  1.0e-02, 1.5e-2, 3.0e-02, 5.9e-02)
cf <- c(10^seq(from = 0, to = 4))
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
# plot.corr.fact <- ldply(auc.cf) %>%
#   ggplot(aes(x = thresh,
#              y = auc, color = .id)) +
#   geom_point() +
#   scale_color_discrete(name = "cf") +
#   stat_smooth(alpha = 0, size = 1.5)+
#   theme_classic() +
#   labs(y = "AUC", x = expression(Log["10"]~Threshold))+
#   scale_colour_brewer(palette = "Set1")
# ggsave("figs for MS\\post poisot\\fish corr factor.png",
#        width = 420, height = 200, units = "mm")

# fish corrected ####
# fish relative abundance * 1000
rel.ab.fish <- map(rel.ab.matr,
                   f_ab_corr,
                   taxa = f.vec,
                   cf = 1000)
saveRDS(rel.ab.fish, "rel ab fish x 1000.RDS")

# relative abundane matrices (N) with "corrected" fish
fish.neutral.list <- map(threshold2, function (x){
  map(rel.ab.fish, rm_neutral, threshold = x)})
names(fish.neutral.list) <- threshold2
# multiply N  by inferred matrices (Ainf)
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
  summarize(mean.tss = mean(na.omit(tss)), sd = sd(tss)) %>%
  top_n(1, wt = mean.tss)

# local fish neutral + niche ####
# remove niche forbidden
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
global.thresh.nn <- local.tss.f.nn %>%
  group_by(thresh) %>%
  summarize(mean.tss = mean(na.omit(tss)), sd = sd(tss)) %>% 
  top_n(1, wt = mean.tss)



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
# neutral abundance 1.5e-4 (e.g. fish.neutral.list[[18]])
f.neutral <- fish.neutral.list[[18]]
saveRDS(f.neutral, "fish correction neutral.RDS")
# niche + neutral = same neutral threshold, 1.5e-4
f.nn <- fish.nn.list[[18]]
saveRDS(f.nn, "tm nn fish for PCA.RDS")

neutral.false <- ldply(map2(obs, f.neutral, false_prop)) %>% 
  summarize(mean.fp = mean(fp), sd.fp = sd(fp),
            mean.fn = mean(fn), sd.fn = sd(fn))
nn.false <- ldply(map2(obs, f.nn, false_prop)) %>% 
  summarize(mean.fp = mean(fp), sd.fp = sd(fp),
            mean.fn = mean(fn), sd.fn = sd(fn))

false.tab <- rbind(neutral.false, nn.false)
false.tab$inference <- c("Neutral", "Niche + Neutral")

# AUC TSS threshold table
write_csv(data.frame(inference = c("Fish corrected neutral", "Fish corrected niche + Neutral"),
           AUC = c(as.double(global.f.neutral[2]),
           as.double(global.f.nn[2])),
           TSS = c(as.double(global.thresh.ab[2]),
                   as.double(global.thresh.nn[2])),
           sd.TSS = c(as.double(global.thresh.ab[3]),
                      as.double(global.thresh.nn[3])),
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
# AUC fish corrected Neutral + Niche
(f.auc.nn.plot <- f.auc.nn.df %>%
  ggplot(aes(x = thresh, y = auc)) +
  geom_point() +
  stat_summary(aes(y = auc,group=1),
               fun.y=mean,
               colour="grey",
               geom="line",
               size = 2,
               group= 1) +
  theme_classic()) +
  labs(y = "AUC", x = expression(Log["10"]~Threshold))
ggsave("figs for MS\\post poisot\\TM fish corr nn auc.png",
       width = 420, height = 200, units = "mm")
