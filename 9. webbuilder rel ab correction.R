# useful food web functions from Petchey
source("C:\\Users\\Justin\\Documents\\Data\\FW modelling Petchey Github\\ttl-resources-master\\food_web\\FoodWebFunctions.r")
# TSS function
source("TSS function.R")


# different ab.matr lists
# determined that Npred * Nprey gives best tss score, and is simplest model to explain
# determined 23 june 2017
# script 8 has code to read in all saved plots
ab.matr.list <- readRDS("17 relative abundance matrices list.rds")
#ab.matr.list <- readRDS("17 relative abundance matrices prey2.rds")
#ab.matr.list <- readRDS("17 relative abundance matrices prey only.rds")

# webbuilder inferred ####
# calculated using all observations in links
# this is likely what I will use in manuscript
# read in webbuilder inferred links, trimmed to match gravel inf
wb.trim.list <- readRDS("webbuilder inferred list trimmed taxa.rds")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# filtered links ####
# webbuilder inferred from filtered links > 1 
# just checking to see if rel ab correction helps this
# likely not going to use this in manuscript....
#wb.trim.list <- readRDS(file ="webbuilder filtered links inferred matr.rds")

# checked this on 24 june 2017
# TSS of filtered links does NOT improve with relative abundance correction
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# read in 17 observed foodweb
obs.list <- readRDS("17 observed taieri food webs.rds")


# target vector to match all matrices
# this is also used to trim the matrices to only include taxa in this list
# e.g. matches the taxa in the gravel inferred webs
target <- readRDS("target vector to match all matrices.rds")

# make obs.list match order of wb.trim.list
# this also trims 
for (i in 1:length(obs.list)){
  obs.list[[i]] <- obs.list[[i]][
    match(target[[i]], rownames(obs.list[[i]])),
    match(target[[i]], colnames(obs.list[[i]]))]
}


# TSS abundance####
# for loop calculating TSS for different abundance thresholds

threshold <- c(1e-01, 1e-02, 1e-03, 1e-04, 1e-05, 1e-06, 1e-07,  5.9e-02, 5.9e-03, 5.9e-04, 5.9e-05, 5.9e-06, 5.9e-07, 3e-02, 3e-03, 3e-04, 3e-05, 3e-06, 3e-07)

tss.wb.ab <- NULL
tss.df.list <- NULL
for (j in 1:length(threshold)){
  tss.df <- NULL
  for (i in 1:length(ab.matr.list)){
    ab <- ab.matr.list[[i]]
    grav <- wb.trim.list[[i]]
    obs <- obs.list[[i]]
    thresh.matr <- ab
    thresh.matr[thresh.matr > threshold[j]] <- 1
    thresh.matr[thresh.matr < 1] <- 0
    grav.thresh <- grav * thresh.matr
    tss.temp <- tss(obs, grav.thresh)
    tss.temp$threshold <- threshold[j]
    tss.temp$site <- names(ab.matr.list)[i]
    tss.df <- rbind(tss.df, tss.temp)
  }
  tss.df.list[[j]] <- tss.df
  result <- tss.df %>%
    summarize(tss.mean = mean(tss),
              tss.sd = sd(tss),
              ad = (mean(a) + mean (d)) /mean(S**2),
              bc = mean(b +c) / mean(S**2),
              threshold = max(threshold))
  tss.wb.ab <- rbind(tss.wb.ab, result)
}
tss.wb.ab %>% arrange(tss.mean)
# 1.0e-5 increases tss to 0.690414
# 0.690 - 0.667 = 0.023

# tss.df.list %>% ldply %>%
#   ggplot(aes(x = log10(threshold), y = tss, color = site)) +
#   geom_point() + 
#   stat_smooth(alpha = 0.1) +
#   theme_classic() +
#   ggtitle("WebBuilder inferred food webs")



# log10(threshold)
# this is the plot in power point ####
(wb.relab.tss.plot <- ggplot(tss.wb.ab, aes(x = log10(threshold), y = tss.mean))+
  geom_point()+
  geom_line() +
  geom_errorbar(aes(ymin = tss.mean - tss.sd, 
                  ymax = tss.mean + tss.sd)) +
  #ggtitle("NPred * Nprey") +
  labs(x = expression(Log[10]~Threshold), parse = T,
       y = paste("Mean", "TSS", "\u00B1", "SD", sep = " ")) +
  theme_classic()
)
saveRDS(wb.relab.tss.plot, paste(getwd(), "/figs for MS/wb rel ab tss plot.rds", sep = ""))


# # this is not working, not sure if it ever did, or more importantly, that it is necessary
  # ggplot(tss.wb.ab, aes(x = log10(threshold), y = tss))+
  # geom_point() +
  # stat_smooth(method = "lm", formula = y~ poly(x, 4))
  # 
  # ggplot(tss.wb.ab, aes(x = log10(threshold), y = tss))+
  # geom_point()+
  # stat_smooth()


