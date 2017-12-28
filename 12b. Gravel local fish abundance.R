# 12 sensitivity to fish relative abundance
# 27 Jun 2017
# JPomz

# useful food web functions from Petchey
source("C:\\Users\\Justin\\Documents\\Data\\FW modelling Petchey Github\\ttl-resources-master\\food_web\\FoodWebFunctions.r")
# TSS function
source("TSS function.R")


# read in estimated Taieri dry weights 
taieri <- readRDS("estimated dw taieri webs.rds")

# note for when I want to calc relative abundances
rel.ab <- taieri %>% 
  select(taxa, no.m2) %>%
  group_by(site) %>% 
  mutate(tot.ab = sum(no.m2, na.rm = T), rel.ab = no.m2 / tot.ab) 
# make into a list
rel.ab.list <- split(rel.ab, list(rel.ab$site))


# don't have an estimate for fish abundances
# testing differnt values from 0.01 - 0.55
fish.ab <- c(0.001, 0.005, 0.01, 0.05, 0.1, 0.15, 0.2, 0.25)

# make list of lists
rel.ab.list.list <- vector("list", length(fish.ab))

# make each list (1:7) a list of rel.ab with varying fish abundance levels
for (i in 1:length(rel.ab.list.list)){
  rel.ab.list.list[[i]] <- rel.ab.list
}
# name elements of list after fish abundances
names(rel.ab.list.list) <- fish.ab

# replace NAs with fish abundance values
for (j in 1:length(rel.ab.list.list)){
  for (i in 1:length(rel.ab.list.list[[1]])){
    rel.ab.list.list[[j]][[i]]$rel.ab[is.na(rel.ab.list.list[[j]][[i]]$rel.ab)] <- fish.ab[j]
  }
}
# rel.ab.list.list[[1]] == all 10 food webs with fish ab == 0.01
# rel.ab.list.list[[2]] == all 10 food webs with fish ab == 0.05
# etc



# make list of abundance matrices
ab.matr.list <- rel.ab.list %>%
  llply(function (x){
    matrix(0, nrow(x), nrow(x))
  })
names(ab.matr.list) <- names(rel.ab.list)

# add dimnames to matrices
for (i in 1:length(ab.matr.list)){
  dimnames(ab.matr.list[[i]]) <- list(rel.ab.list[[i]]$taxa,
                                      rel.ab.list[[i]]$taxa)
}

ab.matr.list.list <- NULL
# calculate relative abundances
# this takes a while....
for (l in 1:length(rel.ab.list.list)){
  for (k in 1:length(ab.matr.list)){
    for (i in 1:nrow(rel.ab.list[[k]])){
      for (j in 1:nrow(rel.ab.list[[k]])){
      ab.matr.list[[k]][i,j] <- 
        rel.ab.list.list[[l]][[k]]$rel.ab[i] *
        rel.ab.list.list[[l]][[k]]$rel.ab[j]
      }
      ab.matr.list.list[[l]] <- ab.matr.list
    }
  }
}
names(ab.matr.list.list) <- names(rel.ab.list.list)

# # read in webbuilder inferred links, trimmed to match gravel inf
# wb.trim.list <- readRDS("webbuilder inferred list trimmed taxa.rds")
# # read in 17 observed foodweb
# obs.list <- readRDS("17 observed taieri food webs.rds")


# target vector to match all matrices
# this is also used to trim the matrices to only include taxa in this list
# e.g. matches the taxa in the gravel inferred webs
target <- readRDS("target vector to match all matrices.rds")

# make ab.matr.list match the target vector
# both trims rosw/cols and orders to target vector
for(j in 1:length(ab.matr.list.list)){
  for (i in 1:length(ab.matr.list)){
    ab.matr.list.list[[j]][[i]] <- 
      ab.matr.list.list[[j]][[i]][target[[i]], target[[i]]]
  }
}

# read in observed list, taxa trimmed and matched to target
obs.list <- readRDS("observed adj matr list trimmed taxa.rds")
# read in gravel inferred list, forbidden links pruned
grav.inf <- readRDS("gravel inferred pruned.rds")

# colnames() of obs.list, grav.inf, ab.matr.list.list [[x]], all identical

# TSS abundance####
# for loop calculating TSS for different abundance thresholds
threshold <- 5e-05

tss.grav.ab <- NULL
tss.df.list <- NULL
for (j in 1:length(ab.matr.list.list)){
  tss.df <- NULL
  for (i in 1:length(ab.matr.list.list[[1]])){
    ab <- ab.matr.list.list[[j]][[i]]
    grav <- grav.inf[[i]]
    obs <- obs.list[[i]]
    thresh.matr <- ab
    thresh.matr[thresh.matr > threshold] <- 1
    thresh.matr[thresh.matr < 1] <- 0
    grav.thresh <- grav * thresh.matr
    tss.temp <- tss(obs, grav.thresh)
    tss.temp$f.ab <- fish.ab[j] 
    tss.temp$site <- names(ab.matr.list.list[[1]])[i]
    tss.df <- rbind(tss.df, tss.temp)
  }
  tss.df.list[[j]] <- tss.df
  result <- tss.df %>%
    summarize(tss.mean = mean(tss),
              tss.sd = sd(tss),
              ad = (mean(a) + mean (d)) /mean(S**2),
              bc = mean(b +c) / mean(S**2),
              f.ab = max(f.ab))
  tss.grav.ab <- rbind(tss.grav.ab, result)
}
tss.grav.ab


# local threshold ####
local.threshold <- tss.df.list %>%
  ldply %>%
  select(tss, site, f.ab) %>%
  group_by(site) %>%
  top_n(1, tss) # i think this is wrong

tss.df.list %>% ldply() %>%
  ggplot(aes(x = f.ab, y = tss, color = site))+
  geom_point()+
  #geom_errorbar(aes(ymin = tss.mean - tss.sd, 
  #                  ymax = tss.mean + tss.sd))+
  geom_line()+
  scale_x_continuous(breaks = fish.ab) +
  theme_classic()+
  #ggtitle("Mean TSS ~ relative fish abundance") +
  labs(x = "Relative fish abundance",
       y = paste("Mean", "TSS", "\u00B1", "SD", sep = " "))
)

# plot in power point
(fish.ab.tss.plot <- ggplot(tss.grav.ab, aes(x = f.ab, y = tss.mean))+
  geom_point()+
  geom_errorbar(aes(ymin = tss.mean - tss.sd, 
                    ymax = tss.mean + tss.sd))+
  geom_line()+
  scale_x_continuous(breaks = fish.ab) +
  theme_classic()+
  #ggtitle("Mean TSS ~ relative fish abundance") +
  labs(x = "Relative fish abundance",
       y = paste("Mean", "TSS", "\u00B1", "SD", sep = " "))
)

# saveRDS(fish.ab.tss.plot, paste(getwd(), "/figs for MS/gravel tss fish ab.rds", sep = ""))



















# TSS abundance####
# for loop calculating TSS for different abundance thresholds
# threshold <- c(1e-01, 1e-02, 1e-03, 1e-04, 1e-05, 1e-06, 1e-07, 1e-08, 1e-9, 1e-10, 1e-11, 5e-02, 5e-03, 5e-04, 5e-05, 5e-06, 5e-07, 5e-08, 5e-9, 5e-10, 5e-11)
# 
# tss.grav.ab <- NULL
# tss.df.list <- NULL
# for (j in 1:length(threshold)){
#   tss.df <- NULL
#   for (i in 1:length(ab.matr.list)){
#     ab <- ab.matr.list[[i]]
#     grav <- grav.inf[[i]]
#     obs <- obs.list[[i]]
#     thresh.matr <- ab
#     thresh.matr[thresh.matr > threshold] <- 1
#     thresh.matr[thresh.matr < 1] <- 0
#     grav.thresh <- grav * thresh.matr
#     tss.temp <- tss(obs, grav.thresh)
#     tss.temp$threshold <- threshold[j]
#     tss.temp$site <- names(ab.matr.list)[i]
#     tss.df <- rbind(tss.df, tss.temp)
#   }
#   tss.df.list[[j]] <- tss.df
#   result <- tss.df %>%
#     summarize(tss.mean = mean(tss),
#               tss.sd = sd(tss),
#               ad = (mean(a) + mean (d)) /mean(S**2),
#               bc = mean(b +c) / mean(S**2),
#               threshold = max(threshold))
#   tss.grav.ab <- rbind(tss.grav.ab, result)
# }
# 
