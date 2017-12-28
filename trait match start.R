# trait match package

# install.packages("devtools")
# install.packages("GenSA")
# install.packages("SDMTools")
# require(devtools)
# install_github("traitmatch", "ibartomeus")
require(traitmatch)
library(gplots)
# useful food web functions from Petchey
source("C:\\Users\\Justin\\Documents\\Data\\FW modelling Petchey Github\\ttl-resources-master\\food_web\\FoodWebFunctions.r")

# Probability interaction new data
source("predict.niche.prob.R")

# fish <- Barnes2008
# head(fish)
# 
# # Define the vectors for the predator and prey size
# MPred = log10(fish$standardised_predator_length*10)
# MPrey = log10(fish$si_prey_length*10)
# 
# mt <- 60 #Define max.time to 60 sec to run things fast. Set to minimum 900 for a decent estimation of parameters.
# pars_pre <- fit_it(integrated_model, 
#                    Tlevel1 = MPrey,  
#                    Tlevel2 = MPred,
#                    mean_Tlevel1 = mean(MPrey),
#                    sd_Tlevel1 = sd(MPrey),
#                    pars = c(a0 = 0, a1 = 0, b0 = 0, b1 = 0),
#                    par_lo = c(a0 = -10, a1 = 0, b0 = -10, b1 = -10),
#                    par_hi = c(a0 = 10, a1 = 10, b0 = 10, b1 = 10),
#                    max.time = mt)
# 
# 
# pars_pre
# plot_pred(pars = pars_pre, Tlevel1 = MPrey, 
#           Tlevel2 = MPred, xlab = "log (Predator body size)", 
#           ylab = "log (prey body size)", pch = "0")
# MPrey2 = sample(MPrey, 100)
# MPred2 = sample(MPred, 100)
# 
# integrated_model(pars = pars_pre,
#                  Tlevel1 = MPrey2, 
#                  Tlevel2 = MPred2,
#                  mean_Tlevel1 = mean(MPrey2),
#                  sd_Tlevel1 = sd(MPrey2))
# plot_pred(pars = pars_pre, Tlevel1 = MPrey2, 
#           Tlevel2 = MPred2, xlab = "log (Predator body size)", 
#           ylab = "log (prey body size)", pch = "0")
# 
# # probability of interaction given size of pred and prey
# # Equation 5 from Bartomeus et al. 2016
# # P() = exp((-(a0 = a1 * Mpred - Mprey)^2) / (2(b0 + b1 * Mpred)^2))
# # parameters from fit_it(integrated_model, ...) fxn
# a0 <- pars_pre[1]
# a1 <- pars_pre[2]
# b0 <- pars_pre[3]
# b1 <- pars_pre[4]
# 
# # probability of link given body sizes
# p.link <- exp((-(a0 + a1 * MPred - MPrey)^2) / (2*(b0 + b1 * MPred)^2))
# 
# # make df of body sizes and probability of link
# prob.df <- data.frame(Mpred = MPred,
#                       Mprey = MPrey,
#                       p.link = p.link,
#                       Pred.l = fish$standardised_predator_length,
#                       Prey.l = fish$si_prey_length)
# # mutata df, add pred:prey ratio, 
# # generate random "observed" links to work out plot below
# prob.df <- prob.df %>%
#   mutate(ratio = log10(Pred.l * 10/ Prey.l *10),
#          link.obs = as.factor(rbinom(n(), 1, .3)))
# 
# # plot link probability
# ggplot(prob.df, aes(Mpred, Mprey, color = p.link)) + 
#   geom_point() 
# 
# # plot observed links
# ggplot(subset(prob.df, link.obs ==1), aes(Mpred, Mprey)) + 
#   geom_point(color = "red")
# 
# #plot observed and probability
# ggplot(prob.df, aes(Mpred, Mprey, color = p.link)) + 
#   geom_point(shape = 15, size = 3) +
#   geom_point(data = subset(prob.df, link.obs ==1), 
#              aes(Mpred, Mprey),
#              shape = 20, size = 2, alpha = 0.5, color = "red")
# 
# # probability function based on pred:prey size
# ggplot(prob.df, aes(x = ratio, y = p.link)) +
#   stat_smooth() #+geom_point()
# 
# 
# # working out code for plotting all body size combos
# # probably not going to work as data becomes too large
# # if use expand.grid() seems to take less memory, but still...
# Mpred2 <- MPred[1:500]
# Mprey2 <- MPrey[1:500]
# # rep each value in Mpred2 length of Mprey2
# #Mpred3 <- rep(Mpred2, each = length(Mprey2))
# # rep Mprey3 vector length of Mpred2
# #Mprey3 <- rep(Mprey2, length(Mprey2))
# Mall <- expand.grid(Mpred = Mpred2, Mprey = Mprey2)
# 
# p.link2 <- exp((-(a0 + a1 * Mall$Mpred - Mall$Mprey)^2) / (2*(b0 + b1 * Mall$Mpred)^2))
# 
# prob.df2 <- data.frame(Mpred = Mall$Mpred,
#                       Mprey = Mall$Mprey,
#                       p.link = p.link2)
# # plot link probability
# ggplot(prob.df2, aes(Mpred, Mprey, color = p.link)) + 
#   geom_point() 

# useful food web functions from Petchey
source("C:\\Users\\Justin\\Documents\\Data\\FW modelling Petchey Github\\ttl-resources-master\\food_web\\FoodWebFunctions.r")
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
# should throw ~20 warnings
# this is bc output of Matrix.to.list() [line 24] is class == factor
# left_join() commands above [lines 34, 37] convert factor to character class and combine accordingly
# output is correct, can ignore warnings

# names elements in list to match site names
names(pairs.dw.list) <- names(pairs.list)

# make list into large data frame
test <- ldply(pairs.dw.list)

# filter out NA dw estimates
test <- test %>% 
  select(res.dw, con.dw, res.ab, con.ab) %>%
  filter(!is.na(res.dw), !is.na(con.dw))
# calculate relative abundance
test <- test %>%
  mutate(rel.ab = con.ab * res.ab)

# remove eels
test <- test %>% 
  filter(con.dw < 4)

# read in warburton data
# read in data
hw <- read_csv("warburton.csv")

# add pred and prey dw column
hw <- hw %>%
  mutate(pred_log = 
           pred_ln_a + 
           (pred_b * log(pred_length_mm, base = pred_base)),
         con.dw = (pred_base^pred_log)/pred_g_conversion,
         prey_log =
           prey_ln_a +
           (prey_b * log(prey_length_mm, base = prey_base)),
         res.dw =  (prey_base^prey_log)/prey_g_conversion)

# filter out pairs where no predation occured
eaten <- hw %>% filter(number_eaten_hr !=0)

# make object of just pred / prey dw
hw.pairs <- eaten %>% select(con.dw, res.dw)

# bind 2 datasets together
test <- bind_rows(test, hw.pairs)

# Define the vectors for the predator and prey size
MPred <- log10(test$con.dw)
MPrey <- log10(test$res.dw)

mt <- 60 #Define max.time to 60 sec to run things fast. Set to minimum 900 for a decent estimation of parameters.
pars_pre <- fit_it(integrated_model, 
                   Tlevel1 = MPrey,  
                   Tlevel2 = MPred,
                   mean_Tlevel1 = mean(MPrey),
                   sd_Tlevel1 = sd(MPrey),
                   pars = c(a0 = 0, a1 = 0, b0 = 0, b1 = 0),
                   par_lo = c(a0 = -10, a1 = 0, b0 = -10, b1 = -10),
                   par_hi = c(a0 = 10, a1 = 10, b0 = 10, b1 = 10),
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
                   pars = c(a0 = 0, a1 = 0, b0 = 0, b1 = 0),
                   par_lo = c(a0 = -10, a1 = 0,
                              b0 = -10, b1 = -10),
                   par_hi = c(a0 = 10, a1 = 10,
                              b0 = 10, b1 = 10),
                   max.time = mt)
pars_pre
pars_niche
plot_pred(pars = pars_niche,
          Tlevel1 = MPrey, 
          Tlevel2 = MPred,
          xlab = "log (Predator body size)", 
          ylab = "log (prey body size)",
          pch = "0")

lh_model <- -integrated_model(pars_pre, MPrey, MPred, mean(MPrey),sd(MPrey))
lh_niche <- -niche_model(pars_niche, MPrey, MPred, mean(MPrey), sd(MPrey))
lh_neutral <- -neutral_model(pars = NULL, MPrey, MPred, mean(MPrey), sd(MPrey))

barplot(c(lh_model, lh_niche, lh_neutral), names.arg = c("integrated", "niche", "neutral"))

# replicated HW ####
# replicate by number of prey eaten
n.times <- eaten$number_missing_corrected
eaten.rep <- eaten[rep(seq_len(nrow(eaten)), n.times),]
pairs.rep <- eaten.rep %>%
  select(res.dw, con.dw)

pairs.res <- log10(pairs.rep$res.dw)
pairs.con <- log10(pairs.rep$con.dw)

pars_rep_pre <- fit_it(integrated_model, 
                   Tlevel1 = pairs.res,  
                   Tlevel2 = pairs.con,
                   mean_Tlevel1 = mean(pairs.res),
                   sd_Tlevel1 = sd(pairs.res),
                   pars = c(a0 = 0, a1 = 0, b0 = 0, b1 = 0),
                   par_lo = c(a0 = -10, a1 = 0,
                              b0 = -10, b1 = -10),
                   par_hi = c(a0 = 10, a1 = 10,
                              b0 = 10, b1 = 10),
                   max.time = mt)
plot_pred(pars = pars_rep_pre, 
          Tlevel1 = pairs.res, 
          Tlevel2 = pairs.con, 
          xlab = "log (Predator body size)", 
          ylab = "log (prey body size)", 
          pch = "0")

pars_rep_niche <- fit_it(niche_model, 
                       Tlevel1 = pairs.res,  
                       Tlevel2 = pairs.con,
                       mean_Tlevel1 = mean(pairs.res),
                       sd_Tlevel1 = sd(pairs.res),
                       pars = c(a0 = 0, a1 = 0, b0 = 0, b1 = 0),
                       par_lo = c(a0 = -10, a1 = 0,
                                  b0 = -10, b1 = -10),
                       par_hi = c(a0 = 10, a1 = 10,
                                  b0 = 10, b1 = 10),
                       max.time = mt)
plot_pred(pars = pars_rep_niche, 
          Tlevel1 = pairs.res, 
          Tlevel2 = pairs.con, 
          xlab = "log (Predator body size)", 
          ylab = "log (prey body size)", 
          pch = "0")

lh_model_rep <- -integrated_model(pars_rep_pre, pairs.res, pairs.con, mean(pairs.res),sd(pairs.res))
lh_niche_rep <- -niche_model(pars_rep_niche, pairs.res, pairs.con, mean(pairs.res),sd(pairs.res))
lh_neutral_rep <- -neutral_model(pars = NULL, pairs.res, pairs.con, mean(pairs.res),sd(pairs.res))


# Probability interaction new data
source("predict.niche.prob.R")

# # sample some random body sizes
# novel.prey <- sort(sample(MPrey, 5))
# novel.pred <- sort(sample(MPred, 2))
# # expand to have all possible combinations
# all.bod <- expand.grid(novel.pred, novel.prey)
# # predict probability of interaction
# probs.int<- predict.niche.prob(pars_pre,
#                    Tlevel1 = all.bod[[2]],
#                    Tlevel2 = all.bod[[1]])
# # first item in list is probabilities of pairwise interactions
# # probs.int[[1]]
# 
# # probabilities are prey[[1]] by all pred
# # prob[[1]] == prey[[1]]*pred[[1]]
# # prob[[2]] == prey[[1]]*pred[[2]]
# # take 1st item in probability out put, 
# # make matrix with n rows and columns
# # where n = number of species
# 
# pred.a <- matrix(probs.int[[1]],5,5)
# mybreaks <- c(0.6, 0.815, 0.8499585, 0.913179, 0.9989)
# mycol <- c("red", "orange", "yellow", "white")
# 
# heatmap.2(pred.a,
#           Rowv = NA,
#           Colv = NA,
#           scale = "none",
#           trace = "none",
#           dendrogram = "none",
#           col = mycol,
#           breaks = mybreaks)
# 
# # try with bigger dataset
# novel.prey <- sort(sample(MPrey, 20))
# novel.pred <- sort(sample(MPred, 20))
# all.bod <- expand.grid(novel.pred, novel.prey)
# 
# 
# probs.int<- predict.niche.prob(pars_niche, #pars_pre
#                                Tlevel1 = all.bod[[2]],
#                                Tlevel2 = all.bod[[1]])
# a.matrix <- matrix(probs.int[[1]], 20, 20)
# # get probability distribution
# mybreaks <- quantile(a.matrix, probs = seq(0, 1, 0.1))
# 
# heatmap.2(a.matrix,
#           Rowv = NA,
#           Colv = NA,
#           scale = "none",
#           trace = "none",
#           dendrogram = "none",
#           #col = mycol,
#           breaks = mybreaks)

# try with some of my webs
#pull out dempsters
dmp <- pairs.dw.list[[7]]
# filter out NA dw observations
dmp <- dmp %>% filter(!is.na(res.dw), !is.na(con.dw))
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
mybreaks <- quantile(dmp.matrix.pre, probs = seq(0,1,0.25))
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


# pars_niche
dmp.prob.niche <- predict.niche.prob(pars_niche,
                                   all.comb[[1]],
                                   all.comb[[2]],
                                   replicates = 1000)
# turn vector of probabilities into matrix
# make sure the nrows == ncols == length(all.t)
dmp.matrix.niche <- matrix(dmp.prob.niche[[1]],
                         length(all.t),
                         length(all.t))
mybreaks <- quantile(dmp.matrix.niche, probs = seq(0,1,0.25))
heatmap.2(dmp.matrix.niche,
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




# play with gif
library(animation)

for (i in 3:502){
    png(filename = paste(getwd(),
                       "/gif/dmp",
                       i,
                       ".png",
                       sep = ""))
  Plot.matrix(matrix(dmp.prob[[i]], 34,34),
              point.cex = 1.5)
  dev.off()
}
gif.files <- list.files(paste(
  getwd(), "/gif/", sep = ""),
  full.names = F)
setwd(paste(getwd(),"/gif/",sep = ""))
ani.options(interval=.1)
im.convert(gif.files)
setwd("C:/Users/Justin/Google Drive/Data/Predicting NZ Food Webs")




dmp.niche.prob <- predict.niche.prob(pars_niche,
                   all.comb[[1]],
                   all.comb[[2]])
dmp.niche.matrix <- matrix(dmp.niche.prob[[1]],
                           length(all.t),
                           length(all.t))
mynichebreaks <- quantile(dmp.niche.matrix, probs = seq(0,1,0.03))
heatmap.2(dmp.niche.matrix,
          Rowv = NA,
          Colv = NA,
          scale = "none",
          trace = "none",
          dendrogram = "none",
          #col = mycol,
          breaks = mynichebreaks,
          key = F)

# niche model ####

probs.niche <- predict.niche.prob(pars_niche,
                    Tlevel1 = novel.prey,
                    Tlevel2 = novel.pred)




grass <- Deraison2014
head(grass)
grass <- subset(grass, Herbivory > 0)

pars_grass_bin <- fit_it(integrated_model,
                         Tlevel1 = grass$P.Leaf.dry.matter.content, 
                         Tlevel2 = grass$G.Incisive.strength,
                         mean_Tlevel1 = mean(grass$P.Leaf.dry.matter.content),
                         sd_Tlevel1 = sd(grass$P.Leaf.dry.matter.content),
                         max.time = mt)


plot_pred(pars = pars_grass_bin, 
          Tlevel1 = jitter(grass$P.Leaf.dry.matter.content), 
          Tlevel2 = jitter(grass$G.Incisive.strength), 
          xlab = "Incisive strength",
          ylab = "Leaf dry matter content")

pars_grass_bin_niche <- fit_it(niche_model,
                        Tlevel1 = grass$P.Leaf.dry.matter.content, 
                        Tlevel2 = grass$G.Incisive.strength,
                        mean_Tlevel1 = mean(grass$P.Leaf.dry.matter.content),
                        sd_Tlevel1 = sd(grass$P.Leaf.dry.matter.content),
                        max.time = mt)

plot_pred(pars = pars_grass_bin_niche, 
          Tlevel1 = jitter(grass$P.Leaf.dry.matter.content), 
          Tlevel2 = jitter(grass$G.Incisive.strength), 
          xlab = "Incisive strength", 
          ylab = "Leaf dry matter content")
          





Incisive.strength <- c()

for(i in 1:nrow(grass)){
  temp <- rep(grass$G.Incisive.strength[i], round(grass$Herbivory[i]))
  Incisive.strength <- c(Incisive.strength,temp)
}
Dry.matter <- c()
for(i in 1:nrow(grass)){
  temp <- rep(grass$P.Leaf.dry.matter.content[i], round(grass$Herbivory[i]))
  Dry.matter <- c(Dry.matter,temp)
}


pars_grass_freq <- fit_it(integrated_model,
                          Tlevel1 = Dry.matter, 
                          Tlevel2 = Incisive.strength,
                          mean_Tlevel1 = mean(grass$P.Leaf.dry.matter.content),
                          sd_Tlevel1 = sd(grass$P.Leaf.dry.matter.content),
                          max.time = mt)

plot_pred(pars = pars_grass_freq, Tlevel1 = jitter(Dry.matter,100), 
          Tlevel2 = jitter(Incisive.strength), xlab = "Incisive strength", ylab = "Leaf dry matter content")

pars_grass_freq_niche <- fit_it(niche_model, 
                      Tlevel1 = Dry.matter, 
                      Tlevel2 = Incisive.strength,
                      mean_Tlevel1 = mean(grass$P.Leaf.dry.matter.content),
                  sd_Tlevel1 = sd(grass$P.Leaf.dry.matter.content),
                                max.time = mt)

plot_pred(pars = pars_grass_freq_niche,
          Tlevel1 = jitter(Dry.matter,100),
          Tlevel2 = jitter(Incisive.strength),
          xlab = "Incisive strength",
          ylab = "Leaf dry matter content")
          