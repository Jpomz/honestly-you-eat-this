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


# fish abundance ####
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
# relative abundance cross product matrices
rel.ab.matr <- readRDS("relative abundance matrices.RDS")
# niche forbidden
taxa.forbid <- c("Amphipoda", "Atalophlebioides", "Austroclima", "Austrosimulium", "Blephariceridae", "Coloburiscus", "Deleatidium", "Nesameletus","Oligochaetae", "Ostracoda", "Oxyethira", "Potamopyrgus", "Zephlebia")

f.vec <- c("Salmo", "Galaxias", "Anguilla", "Gobiomorpus")
threshold2 <- c(#5.9e-13,
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

# niche neutral ####
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
