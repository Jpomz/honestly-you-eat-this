# Warburton data

# convert body lengths to dw (g)
# parameterize gravel model
# 

# read in functions
# gravel functions
source("gravel_functions.R")
# useful food web functions from Petchey
source("C:\\Users\\Justin\\Documents\\Data\\FW modelling Petchey Github\\ttl-resources-master\\food_web\\FoodWebFunctions.r")

# ===========================================

# read in data
hw <- read_csv("warburton.csv")

# add pred and prey dw column
hw <- hw %>%
  mutate(pred_log = pred_ln_a + (pred_b * log(pred_length_mm, base = pred_base)),
         pred_dw_g = (pred_base^pred_log)/pred_g_conversion,
         pred_log_dw = log10(pred_dw_g),
         prey_log = prey_ln_a + (prey_b * log(prey_length_mm, base = prey_base)),
         prey_dw_g =  (prey_base^prey_log)/prey_g_conversion,
         prey_log_dw = log10(prey_dw_g))

# filter out pairs where no predation occured
eaten <- hw %>% filter(number_eaten_hr !=0)

# make object of just pred / prey dw
hw.pairs <- eaten %>% select(pred_dw_g, prey_dw_g)
#saveRDS(hw.pairs, file = "hw.pred.prey.dw.pairs.rds")

# parameterize model ####
Bprey <- eaten$prey_log_dw
Bpred <- eaten$pred_log_dw
# plot of log prey ~ log pred with quantiles
ggplot(eaten, aes(x = pred_log_dw, y = prey_log_dw)) +
  geom_point() +
  stat_quantile(quantiles = c(0.000000001, .99))

# Calculate the parameters
pars_reg <-  reg_fn(Bprey,Bpred,quartil = c(0.00000001,.99))

# # make up some body sizes
# Ball <- seq(-5, 12.2, by = 0.05)
# # get niche parameters for made up body sizes
# pars_niche <-  get_pars_Niche(pars_reg,Ball)
# # Treshold  <-  which(pars_niche[,"n"]<=1.367915)
# # pars_niche[,"low"][Treshold] <- 0
# # pars_niche[,"high"][Treshold] <- 0
# 
# # make predation matrix using niche paramters 
# L <-  L_fn(pars_niche[,1],pars_niche[,2],pars_niche[,3],pars_niche[,4])
# # plot predation matrix
# Plot.matrix(L)

# =========================================

# apply to Taieri data ####
# read in estimated Taieri dry weights 
taieri <- readRDS("estimated dw taieri webs.rds")
# remove NAs
taieri <- taieri %>% filter(!is.na(dw))
# sort
taieri <- taieri %>% arrange(site, logdw)
taieri.list <- split(taieri, list(taieri$site))


# calculate niche parameters
taieri.list.pars <- llply(taieri.list, function (x){
  pars = as.data.frame(get_pars_Niche(pars_reg, Ball = x$logdw))
  taxa = x$taxa
  cbind(taxa, pars)
})

# make a threshold for small bugs
treshold.list <- llply(taieri.list.pars, function (x){
  which(x[,"n"]<=-4)
})

# set small bugs hi/low niche to 0
for (i in 1:length(taieri.list.pars)){
  taieri.list.pars[[i]][,"low"][treshold.list[[i]]] <- 0
  taieri.list.pars[[i]][,"high"][treshold.list[[i]]] <- 0
}

# make "wider" niche range for fish  
max.list <- llply(taieri.list.pars, function (x){
  which(x[,"n"]>-2.5)
})
# set fish niche range (currently arbitrary number)
for (i in 1:length(taieri.list.pars)){
  taieri.list.pars[[i]][,"low"][max.list[[i]]] <- -5
  taieri.list.pars[[i]][,"high"][max.list[[i]]] <- -1
}

# calculate food web links for taieri
taieri.list.links <- llply(taieri.list.pars, function (x){
  L_fn(x[,"n"],x[,"c"],x[,"low"],x[,"high"])
})
# add taxa names to row/col of predation matrix
for (i in 1:length(taieri.list.links)){
  dimnames(taieri.list.links[[i]]) <- list(
    taieri.list.pars[[i]][,"taxa"],
    taieri.list.pars[[i]][,"taxa"]
  )
}

# plot predation matrices
for (i in 1:length(taieri.list.links)){
  Plot.matrix(taieri.list.links[[i]])
}

# taieri.list.pars %>% llply(function (x){
#   ggplot(x, aes(x = n, y = c))+
#     geom_line(lty = 2) +
#     geom_line(aes(x = n, y = low), lty = 4) +
#     geom_line(aes(x = n, y = high), lty = 3) +
#     geom_line(aes(x = n, y = n), lty = 1)
# })


taieri.list.pars %>% bind_rows() %>% 
  ggplot( aes(x = n, y = n))+
  geom_point() +
  geom_line(aes(x = n, y = low), lty = 4) +
  geom_line(aes(x = n, y = high), lty = 3) +
  geom_line(aes(x = n, y = c), lty = 1)
