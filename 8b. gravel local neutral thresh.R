# gravel inferred matrices abundance correction

# useful food web functions from Petchey
source("C:\\Users\\Justin\\Documents\\Data\\FW modelling Petchey Github\\ttl-resources-master\\food_web\\FoodWebFunctions.r")
# TSS function
source("TSS function.R")

# read in 17 observed foodweb
obs.list <- readRDS("17 observed taieri food webs.rds")
# read in gravel inferred list, forbidden links pruned
grav.inf <- readRDS("gravel inferred pruned.rds")

# target vector to match all matrices
# this is also used to trim the matrices to only include taxa in this list
# e.g. matches the taxa in the gravel inferred webs
target <- readRDS("target vector to match all matrices.rds")

# make obs.list match order of target
# this also trims
for (i in 1:length(obs.list)){
  obs.list[[i]] <- obs.list[[i]][
    match(target[[i]], rownames(obs.list[[i]])),
    match(target[[i]], colnames(obs.list[[i]]))]
}


# make grav.inf match order of target
# this also trims 
for (i in 1:length(grav.inf)){
  grav.inf[[i]] <- grav.inf[[i]][
    match(target[[i]], rownames(grav.inf[[i]])),
    match(target[[i]], colnames(grav.inf[[i]]))]
}

# read in matrices with relative abundances
# made in Relative abundance...R script
ab.matr.list <- readRDS("17 relative abundance matrices list.rds")
#ab.matr.list <- readRDS("10 relative abundance matrices prey2.rds")
#ab.matr.list <- readRDS("10 relative abundance matrices prey only.rds")

# TSS abundance####
# for loop calculating TSS for different abundance thresholds
threshold <- c(1e-01, 1e-02, 1e-03, 1e-04, 1e-05, 1e-06, 1e-07,  5.9e-02, 5.9e-03, 5.9e-04, 5.9e-05, 5.9e-06, 5.9e-07, 3e-02, 3e-03, 3e-04, 3e-05, 3e-06, 3e-07)



tss.grav.ab <- NULL
tss.df.list <- NULL
for (j in 1:length(threshold)){
  tss.df <- NULL
  for (i in 1:length(ab.matr.list)){
    ab <- ab.matr.list[[i]]
    grav <- grav.inf[[i]]
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
  tss.grav.ab <- rbind(tss.grav.ab, result)
}
# arrange by mean TSS
tss.grav.ab %>% arrange(tss.mean)
# highest tss = threshold 5.9e-05


# local threshold ####
local.threshold <- tss.df.list %>%
  ldply %>%
  select(tss, site, threshold) %>%
  group_by(site) %>%
  top_n(1, tss)
# density plot
ggplot(local.threshold,
       aes(x = log10(threshold))) +
  geom_density()+
  coord_cartesian(xlim = c(-7,-1))
# histogram
ggplot(local.threshold,
       aes(x = as.factor(threshold))) +
  geom_histogram(stat = "count")

# Tss ~ threshold + site 
tss.df.list %>% ldply %>%
  ggplot(aes(x = log10(threshold), y = tss, color = site)) +
  geom_point() +
  stat_smooth(alpha = 0.1) +
  theme_classic() +
  ggtitle("Gravel inferred food webs")

# this is the plot in power point ####
(gravel.relab.plot <- ggplot(tss.grav.ab, aes(x = log10(threshold), y = tss.mean))+
  geom_point()+
  geom_line() +
  geom_errorbar(aes(ymin = tss.mean - tss.sd, 
                    ymax = tss.mean + tss.sd)) +
  #ggtitle("NPred * Nprey") +
  labs(y = paste("Mean", "TSS", "\u00B1", "SD", sep = " "),
       x = expression(Log[10]~Threshold), parse = T) +
  theme_classic()
)

saveRDS(gravel.relab.plot, paste(getwd(), "/figs for MS/gravel rel ab tss plot.rds", sep = ""))


# read in tss for 2 steps
tss.2.step <- readRDS(paste(getwd(),"/figs for MS/gravel tss 2 steps.rds", sep = ""))
# select threshold that gives highest tss
# 5.9e-05
thresh.tss <- tss.df.list[[11]]
# rename variables, mutate so matches 2 step table above
thresh.tss <- thresh.tss %>% 
  rename(web = site, step = threshold) %>% 
  mutate(step = "abundance")

# join two tables
tss.3.step <- bind_rows(tss.2.step, thresh.tss)
saveRDS(tss.3.step, file = paste(getwd(), "/figs for MS/gravel tss 3 steps.rds", sep =""))

# tss.3.step <- readRDS(file = paste(getwd(), "/figs for MS/gravel tss 3 steps.rds", sep =""))

# get means / sd
tss.3.step.summ <- tss.3.step %>%
  group_by(step) %>%
  summarize(tss.mean = mean(tss), 
            tss.sd = sd(tss), 
            a.mean = mean(abar),
            a.sd = sd(abar),
            d.mean = mean(dbar), 
            d.sd = sd(dbar),
            b.mean = mean(bbar),
            b.sd = sd(bbar),
            c.mean = mean(cbar), 
            c.sd = sd(cbar))

# make plots 
# change order of factor levels for plotting
tss.3.step.summ$step <- factor(tss.3.step.summ$step,
                               levels = c("model", "prune",
                                          "abundance"))
# # mean tss, a, d, bc
# ggplot(tss.3.step.summ, aes(x = step, y = tss.mean, color = "tss"))+
#   geom_point(shape = 0, size = 4) +
#   geom_point(aes(x = step, y = a.mean, color = "a"), size = 4) +
#   geom_point(aes(x = step, y = b.mean, color = "b"), 
#              shape = 10, size = 4) +
#   geom_point(aes(x = step, y = c.mean, color = "c"), 
#              shape = 15, size = 4) +
#   geom_point(aes(x = step, y = d.mean, color = "d"), size = 4, shape = 17) +
#   labs(y = "Proportion")+
#   theme_classic() +
#   theme(panel.border = element_rect(colour = "black",
#                                     fill=NA, size=1))
# # all tss, a, bc, d
# # change order of factor levels for plotting
# tss.3.step$step <- factor(tss.3.step$step,
#                   levels = c("model", "prune","abundance"))
# 
# ggplot(tss.3.step, aes(x = step, y = tss, color = "tss"))+
#   geom_jitter(shape = 15, size = 4, width = 0.2) +
#   geom_jitter(aes(x = step, y = abar, color = "a"), size = 4, width = 0.2, shape = 1) +
#   geom_jitter(aes(x = step, y = bbar, color = "b"), 
#              shape = 10, size = 4, width = 0.2) +
#   geom_jitter(aes(x = step, y = cbar, color = "c"), 
#               shape = 10, size = 4, width = 0.2) +
#   geom_jitter(aes(x = step, y = dbar, color = "d"), size = 4, shape = 2, width = 0.2) +
#   labs(y = "Proportion")+
#   theme_classic() +
#   theme(panel.border = element_rect(colour = "black",
#                                     fill=NA, size=1))

# tpr and fp
ggplot(tss.3.step.summ, aes(x = step, y = a.mean, color = "a"))+
  geom_point(shape = 0, size = 4) +
  geom_point(aes(x = step, y = b.mean, color = "b"), 
             shape = 10, size = 4) +
  labs(y = "Proportion")+
  theme_classic() +
  theme(panel.border = element_rect(colour = "black",
                                    fill=NA, size=1))




# calculate food web stats for selected threshold
selected.threshold <- threshold[12]
# grav.inf[[1]]
# ab.matr.list[[1]]

grav.ab.inf <- NULL
for(i in 1:length(ab.matr.list)){
  ab <- ab.matr.list[[i]]
  grav <- grav.inf[[i]]
  thresh.matr <- ab
  thresh.matr[thresh.matr > selected.threshold] <- 1
  thresh.matr[thresh.matr < 1] <- 0
  grav.ab.inf[[i]] <- grav * thresh.matr
}
names(grav.ab.inf) <- names(ab.matr.list)

saveRDS(grav.ab.inf, file = "gravel inferred pruned rel ab.rds")

grav.primary <- readRDS("gravel inferred initial.rds")
grav.pruned <- grav.inf
rm(grav.inf)
grav.ab.inf <- readRDS("gravel inferred pruned rel ab.rds")


primary.stats <- ldply(grav.primary, function (x){
  Get.web.stats(x)}
)
primary.stats$step <- 1

pruned.stats <- ldply(grav.pruned, function (x){
    Get.web.stats(x)}
)
pruned.stats$step <- 2

ab.stats <- ldply(grav.ab.inf, function (x){
  Get.web.stats(x)}
)
ab.stats$step <- 3

gravel.fw.stats <- bind_rows(primary.stats, pruned.stats, ab.stats)

# add model column so it can match with other stats below
gravel.fw.stats$model <- "gravel"
saveRDS(gravel.fw.stats, file = "food web stats gravel 3 steps.rds")








obs.wb.fw.stats <- readRDS(file = "food web stats observed and webbuilder.rds")


gravel.fw.stats <- readRDS("food web stats gravel 3 steps.rds")


all.stats <- bind_rows(gravel.fw.stats, obs.wb.fw.stats)

saveRDS(all.stats, file = "food web stats all.rds")
