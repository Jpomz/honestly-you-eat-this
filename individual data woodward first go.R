# Data from Woodward 2010

# useful food web functions from Petchey
source("C:\\Users\\Justin\\Documents\\Data\\FW modelling Petchey Github\\ttl-resources-master\\food_web\\FoodWebFunctions.r")

library(plyr)
library(tidyverse)
require(traitmatch)
library(gplots)
# Probability interaction new data
source("predict.niche.prob.R")
# gravel functions
source("gravel_functions.R")


# first 16 lines are description data
# skip = 16
dat <- read_csv("BroadstoneData_AERv45ch3.csv", skip = 16)

ggplot(dat, aes(x = log10(predMass), y = log10(preyMass), color = season))+
  geom_point() +
  stat_quantile(quantiles = c(0.01, .97))

ggplot(dat, aes(x = log10(predMass), y = log10(preyMass)))+
  geom_point() +
  stat_quantile(quantiles = c(0.03, .97))


##### doesn't look like it ####
# # do paramaters change seasonally?
# dat.season <- split(dat, list(dat$season))
# seas.pars <- llply(dat.season, function (x){
#   reg_fn(Bprey = log10(x$preyMass), Bpred = log10(x$predMass), quartil = c(0.03,0.97))
#   })
# result <- NULL
# for (i in 1:length(seas.pars)){
#   df <- seas.pars[[i]] %>% ldply
#   result[[i]] <- df
# }
# names(result) <- names(seas.pars)
# plot.dat <- ldply(result)
# names(plot.dat)[2] <- "Int"
# ggplot(data = plot.dat, aes(y = Int, x = .id))+
#          geom_point()
# 

# taieri ####
# read in estimated Taieri dry weights 
taieri <- readRDS("estimated dw taieri webs.rds")

# get all body sizes present at sites
# remove NAs
taieri <- taieri %>% filter(!is.na(dw))
# sort
taieri <- taieri %>% arrange(site, logdw)
# make into a list
taieri.list <- split(taieri, list(taieri$site))
# add rel abundance column
taieri.list <- llply(taieri.list, function (x){
  x %>% mutate(tot.ab = sum(no.m2, na.rm = T), 
               rel.ab = no.m2 / tot.ab)
})
# read in corrected adjaceny matrices
web.list <- readRDS("animal only adj list.rds")
web.list <- web.list[names(web.list) %in% names(taieri.list)]


# convert predation matrix to matrix of pairs
pairs.list <- web.list %>%
  llply(function (x){
    as.data.frame(Matrix.to.list(x))})

# name columns resource and consumer
pairs.list <- llply(pairs.list,
                    function (x){
                      colnames(x) <- c("resource", "consumer");x})

# make list of resource consumer pairs with dw
pairs.dw.list <- NULL
for (i in 1:length(pairs.list)){
  resource <- left_join( # add dw for resource column
    pairs.list[[i]], taieri.list[[i]][,c(2,6,8)],
    by = c("resource" = "taxa"))
  consumer <- left_join( # just consumer dw
    pairs.list[[i]], taieri.list[[i]][,c(2,6,8)],
    by = c("consumer" = "taxa"))[3:4]
  combined <- bind_cols(resource, consumer) # combine
  colnames(combined) <- c("resource", "consumer", "res.dw", "res.ab", "con.dw", "con.ab") # rename columns
  pairs.dw.list[[i]] <- combined # put in list
  rm("resource", "consumer", "combined") # remove objects
}

taieri.pairs <- ldply(pairs.dw.list)
taieri.pairs <- filter(taieri.pairs,!is.na(res.dw), !is.na(con.dw))
t.res <- taieri.pairs$res.dw
t.con <- taieri.pairs$con.dw
# ggplot(taieri.pairs, aes(x = log10(con.dw*1000),
#                          y = log10(res.dw*1000)))+
#   geom_point() +
#   stat_quantile(quantiles = c(0.03, .97))

# Trait match####
MPred <- log10(dat$predMass)
  #c(log10(dat$predMass), log10(t.con*1000))
MPrey <- log10(dat$preyMass) 
  #c(log10(dat$preyMass), log10(t.res*1000)) 

mt <- 60 #Define max.time to 60 sec to run things fast. Set to minimum 900 for a decent estimation of parameters.
pars_pre <- fit_it(integrated_model, 
                   Tlevel1 = MPrey,  
                   Tlevel2 = MPred,
                   mean_Tlevel1 = mean(MPrey),
                   sd_Tlevel1 = sd(MPrey),
                   pars = c(a0 = 0,
                            a1 = 0,
                            b0 = 0,
                            b1 = 0),
                   par_lo = c(a0 = -10,
                              a1 = 0,
                              b0 = -10,
                              b1 = -10),
                   par_hi = c(a0 = 10,
                              a1 = 10,
                              b0 = 10,
                              b1 = 10),
                   max.time = mt)


pars_pre
plot_pred(pars = pars_pre, 
          Tlevel1 = MPrey, 
          Tlevel2 = MPred, 
          xlab = "log (Predator body size)", 
          ylab = "log (prey body size)", 
          pch = "0")

pars_niche <- fit_it(niche_model, 
                     Tlevel1 = MPrey,  
                     Tlevel2 = MPred,
                     mean_Tlevel1 = mean(MPrey),
                     sd_Tlevel1 = sd(MPrey),
                     pars = c(a0 = 0,
                              a1 = 0,
                              b0 = 0,
                              b1 = 0),
                     par_lo = c(a0 = -10,
                                a1 = 0,
                                b0 = -10,
                                b1 = -10),
                     par_hi = c(a0 = 10,
                                a1 = 10,
                                b0 = 10,
                                b1 = 10),
                     max.time = mt)
pars_pre
pars_niche
plot_pred(pars = pars_niche,
          Tlevel1 = MPrey, 
          Tlevel2 = MPred,
          xlab = "log (Predator body size)", 
          ylab = "log (prey body size)",
          pch = "0")

lh_model <- -integrated_model(pars_pre,
                              MPrey,
                              MPred,
                              mean(MPrey),
                              sd(MPrey))
lh_niche <- -niche_model(pars_niche,
                         MPrey, 
                         MPred, 
                         mean(MPrey),
                         sd(MPrey))
lh_neutral <- -neutral_model(pars = NULL,
                             MPrey, 
                             MPred, 
                             mean(MPrey),
                             sd(MPrey))

barplot(c(lh_model, lh_niche, lh_neutral),
        names.arg = c("integrated",
                      "niche",
                      "neutral"))


# predict new intxn ####
# try with some of my webs
#pull out dempsters
dmp <- pairs.dw.list[[7]]
# filter out NA dw observations
dmp <- dmp %>% filter(!is.na(res.dw),
                      !is.na(con.dw))#,
                      # consumer != "Salmo",
                      # consumer != "Galaxias",
                      # consumer != "Gobiomorphus",
                      # consumer != "Anguilla")
#make vectors of 2 trophic levels
t1 <- as.vector(sort(unique(dmp$res.dw)))
t2 <- as.vector(sort(unique(dmp$con.dw)))
# combine unique values in t1 and t2
# gives you a vector of ALL body sizes
# make sure to sort the sizes
# this makes the prob.predation.matrix make sense later
all.t <- sort(unique(c(t1, t2)))
# expand grid to get all pairwise combos of body sizes
all.comb <- expand.grid(all.t, all.t)
all.comb <- log10(all.comb[with(all.comb, order( Var1)),])


# predict probability of interactions
# pars_pre
dmp.prob.pre <- predict.niche.prob(pars_pre,
                                   all.comb[[1]],
                                   all.comb[[2]],
                                   replicates = 1000)
# turn vector of probabilities into matrix
# make sure the nrows == ncols == length(all.t)
dmp.matrix.pre <- matrix(dmp.prob.pre[[1]],
                         length(all.t),
                         length(all.t))
#mybreaks <- quantile(dmp.matrix.pre, probs = seq(0,1,0.05))
mybreaks <- seq(0,1,0.1)
heatmap.2(dmp.matrix.pre,
          Rowv = NA,
          Colv = NA,
          scale = "none",
          trace = "none",
          dendrogram = "none",
          breaks = mybreaks,
          key = F,
          labRow = NA,
          labCol = NA,
          main= "Probability of interaction",
          xlab = "Consumer",
          ylab = "Resource")


# # pars_niche
# dmp.prob.niche <- predict.niche.prob(pars_niche,
#                                      all.comb[[1]],
#                                      all.comb[[2]],
#                                      replicates = 1000)
# # turn vector of probabilities into matrix
# # make sure the nrows == ncols == length(all.t)
# dmp.matrix.niche <- matrix(dmp.prob.niche[[1]],
#                            length(all.t),
#                            length(all.t))
# mybreaks <- quantile(dmp.matrix.niche, probs = seq(0,1,0.1))
# heatmap.2(dmp.matrix.niche,
#           Rowv = NA,
#           Colv = NA,
#           scale = "none",
#           trace = "none",
#           dendrogram = "none",
#           breaks = mybreaks,
#           key = F,
#           labRow = NA,
#           labCol = NA,
#           main= "Probability of interaction",
#           xlab = "Consumer",
#           ylab = "Resource")






# play with gif
library(animation)


for (i in 3:502){
  png(filename = paste(getwd(),
                       "/gif.ind.dat/dmp",
                       i,
                       ".png",
                       sep = ""))
  Plot.matrix(matrix(dmp.prob.pre[[i]], 34,34),
              point.cex = 1.5)
  dev.off()
}
gif.files <- list.files(
  paste(getwd(), "/gif.ind.dat/", sep = ""),
  full.names = F
  )
setwd(paste(getwd(),"/gif.ind.dat/",sep = ""))
ani.options(interval=.1)
im.convert(gif.files)
setwd("C:/Users/Justin/Google Drive/Data/Predicting NZ Food Webs")


