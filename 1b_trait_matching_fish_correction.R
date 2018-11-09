# fish abundance "correction"
# Justin Pomeranz
# jfpomeranz@gmail.com

# libraries
library(plyr)
library(dplyr)
library(purrr)
library(ggplot2)
# functions written for this manuscript
source("Inference_MS_functions.R")

# data
# observed adjacency matrices, matched to inferences
obs <- readRDS("Results/observed_matrices_matched_to_inferred.RDS")
# initial trait-matching inferred matrices 
inf <- readRDS("Results/initial_trait_matching_inference.RDS")
# relative abundance matrices
rel.ab.matr <- readRDS("Results/relative_abundance_matrices.RDS")
# niche forbidden taxa
taxa.forbid <- readRDS("Data/Raw_data/niche_forbidden_taxa.RDS")
# vector of fish taxa present
f.vec <- c("Salmo", "Galaxias", "Anguilla", "Gobiomorphus")

# vector of thresholds
threshold <- c(
  1.0e-08, 1.5e-8, 3.0e-08, 5.9e-08,
  1.0e-07, 1.5e-7, 3.0e-07, 5.9e-07,
  1.0e-06, 1.5e-6, 3.0e-06, 5.9e-06,
  1.0e-05, 1.5e-5, 3.0e-05, 5.9e-05,
  1.0e-04, 1.5e-4, 3.0e-04, 5.9e-04,
  1.0e-03, 1.5e-3, 3.0e-03, 5.9e-03,
  1.0e-02, 1.5e-2, 3.0e-02, 5.9e-02)
# vector of correction factors
cf <- c(10^seq(from = 0, to = 4))

# correction factor ####
# examine how correction factors influence inferences
# can take a while to run
auc.cf <- NULL
for(c in 1:length(cf)){
  auc.neutral <- NULL
  for(w in 1:length(inf)){
    auc.web <- NULL
    for(t in 1:length(threshold)){
      N = f_ab_corr(Nij = rel.ab.matr[[w]],
                    taxa = f.vec,
                    cf = cf[c])
      Nprime = rm_neutral(N, threshold[t])
      Nprime = Nprime * inf[[w]]
      auc.web[[t]] = get_auc(obs[[w]], Nprime)
    }
    names(auc.web) <- as.character(threshold)
    auc.neutral[[w]] <- auc.web
  }
  names(auc.neutral) <- names(obs)
  auc.cf[[c]] <- auc.neutral
}
names(auc.cf) <- as.character(cf)

auc.cf <- llply(auc.cf, function (x){
  out = data.frame(auc = flatten_dbl(x),
                   thresh = log10(as.numeric(threshold)),
                   site = rep(names(obs),
                              each = length(threshold)),
                   stringsAsFactors = FALSE)})
auc.cf.df <- ldply(auc.cf)
auc.cf.df$cf <- auc.cf.df$.id

# Supplemental figure 1 ####
S1 <- auc.cf.df %>%
  ggplot(aes(x = thresh,
             y = auc, color = cf)) +
  geom_point(position =
               position_jitter(width = 0.05)) +
  stat_summary(aes(y = auc),
               fun.y=mean, geom="line",
               size = 2, alpha = 0.5) +
  theme_classic(base_size = 20) +
  labs(y = "AUC", x = expression(Log["10"]~Threshold)) +
  scale_colour_brewer(palette = "Set1")
# save figure
ggsave("Results/Figures/S1.png",
       width = 420, height = 200, units = "mm")

# Correction factor = 1000
rel.ab.fish <- map(rel.ab.matr,
                   f_ab_corr,
                   taxa = f.vec,
                   cf = 1000)
saveRDS(rel.ab.fish, "Results/relative_ab_fish_x_1000.RDS")

# relative abundance matrices (N) with "corrected" fish
fish.neutral.list <- map(threshold, function (x){
  map(rel.ab.fish, rm_neutral, threshold = x)})
names(fish.neutral.list) <- threshold
# multiply N  by inferred matrices (Ainf)
fish.neutral.list <- map(fish.neutral.list, function (x){
  map2(x, inf, ~.x*.y)
})

# fish neutral + niche ####
# remove niche forbidden
fish.nn.list <- map(fish.neutral.list, function (x){
  map(x, rm_niche, taxa = taxa.forbid)})

#auc fish ####
# neutral ####
# Calculate list of AUC for each threshold:site combination
f.auc.neutral <- get_auc_list(
  neutral.list = fish.neutral.list,
  obs = obs)
# calculate numeric value of list element which has max AUC
f.neutral.threshold <- get_max_auc(x = f.auc.neutral,
            threshold = threshold,
            site = names(obs))


# niche + neutral ####
f.auc.nn <- get_auc_list(fish.nn.list, obs)

f.nn.threshold <- get_max_auc(f.auc.nn, 
                              threshold, 
                              names(obs))

# TSS ####
# TSS neutral ####
f.neutral <- fish.neutral.list[[f.neutral.threshold$list.element]]
tss.neutral <- ldply(pmap(list(obs = obs,
                               inf = f.neutral),
                          get_tss))
tss.neutral.mean <- mean(tss.neutral$V1)
tss.neutral.sd <- sd(tss.neutral$V1)

# TSS Niche + neutral ####
f.nn <- fish.nn.list[[f.nn.threshold$list.element]]
saveRDS(f.nn, "Results/trait_match_niche_neutral_fish.RDS")
tss.nn <- ldply(pmap(list(obs = obs,
                               inf = f.nn),
                          get_tss))
tss.nn.mean <- mean(tss.nn$V1)
tss.nn.sd <- sd(tss.nn$V1)

# false negative & false positive ####
neutral.false <- ldply(map2(obs, f.neutral, false_prop)) %>% 
  summarize(mean.fp = mean(fp), sd.fp = sd(fp),
            mean.fn = mean(fn), sd.fn = sd(fn))
nn.false <- ldply(map2(obs, f.nn, false_prop)) %>% 
  summarize(mean.fp = mean(fp), sd.fp = sd(fp),
            mean.fn = mean(fn), sd.fn = sd(fn))

false.tab <- rbind(neutral.false, nn.false)
false.tab$inference <- c("Neutral", "Niche + Neutral")

# AUC TSS threshold table
write.csv(data.frame(
  inference = c("Fish corrected neutral",
                "Fish corrected niche + Neutral"),
  AUC = c(f.neutral.threshold$auc, f.nn.threshold$auc),
  TSS = c(tss.neutral.mean, tss.nn.mean),
  sd.TSS = c(tss.neutral.sd, tss.nn.sd),
  Threshold = c(f.neutral.threshold$threshold,
                f.nn.threshold$threshold),
  mean.fp = false.tab$mean.fp,
  sd.fp = false.tab$sd.fp,
  mean.fn = false.tab$mean.fn,
  sd.fn = false.tab$sd.fn),
  "Results/Stats/Fish_corrected_AUC_and_TSS_trait_match.csv",
  row.names = FALSE)

# supplemental plot ####
# turn list of fish corrected niche + neutral pruned AUC's into dataframe
f.auc.nn.df <- data.frame(
  auc = flatten_dbl(f.auc.nn),
  thresh = log10(as.numeric(threshold)),                        site = rep(names(obs), each = length(threshold)),
  stringsAsFactors = FALSE)


f.auc.nn.plot <- f.auc.nn.df %>%
  ggplot(aes(x = thresh, y = auc)) +
  geom_point() +
  stat_summary(aes(y = auc,group=1),
               fun.y=mean,
               colour="grey",
               geom="line",
               size = 2,
               group= 1) +
  theme_classic(base_size = 20) +
  labs(y = "AUC", x = expression(Log["10"]~Threshold))

ggsave("Results/Figures/S2.png",
       width = 420, height = 200, units = "mm")
