# 21 Gravel neutral model threshold accuracy
# JPomz

#**************************************************
#**************************************************
#**************************************************
#**************************************************
#### I think this is not a good analysis
# Accuracy increases to ~>90% when nearly all links are removed e.g. threshold == 0.1
# I think this is due to the fact that most possible links go unrealised e.g. v low connectance in empirical webs
# therefore max accuracy occurs when there are essentially no links
# going to scrap this and just go with TSS
# 03 dec 2017
#**************************************************
#**************************************************
#**************************************************
#**************************************************


# re-running neutral model threshold sensitivity analysis using accuracy (ACC) as the response variable. 
# T Poisot's suggestion, confusion matrix derivation from
# https://en.wikipedia.org/wiki/Confusion_matrix
# accessed 3 dec 2017

# just ACC function
source("get_ACC function.R")
# useful food web functions from Petchey
source("C:\\Users\\Justin\\Documents\\Data\\FW modelling Petchey Github\\ttl-resources-master\\food_web\\FoodWebFunctions.r")

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
threshold <- c(1e-01, 1e-02, 1e-03, 1e-04, 1e-05, 1e-06, 1e-07, 5.9e-02, 5.9e-03, 5.9e-04, 5.9e-05, 5.9e-06, 5.9e-07, 3e-02, 3e-03, 3e-04, 3e-05, 3e-06, 3e-07)




acc.grav.ab <- NULL
acc.df.list <- NULL
for (j in 1:length(threshold)){
  acc.df <- NULL
  for (i in 1:length(ab.matr.list)){
    ab <- ab.matr.list[[i]]
    grav <- grav.inf[[i]]
    obs <- obs.list[[i]]
    thresh.matr <- ab
    thresh.matr[thresh.matr > threshold[j]] <- 1
    thresh.matr[thresh.matr < 1] <- 0
    grav.thresh <- grav * thresh.matr
    acc.temp <- data.frame(acc = get_ACC(obs,
                                         grav.thresh))
    acc.temp$threshold <- threshold[j]
    acc.temp$site <- names(ab.matr.list)[i]
    acc.df <- rbind(acc.df, acc.temp)
  }
  acc.df.list[[j]] <- acc.df
  result <- acc.df %>%
    summarize(acc.mean = mean(acc),
              acc.sd = sd(acc),
              threshold = max(threshold))
  acc.grav.ab <- rbind(acc.grav.ab, result)
}
# arrange by mean TSS
acc.grav.ab %>% arrange(acc.mean)
# highest acc = threshold 1.0e-01

ggplot(acc.grav.ab, aes(x = log10(threshold),
                        y = acc.mean))+
    geom_point()+
    geom_line() +
    geom_errorbar(aes(ymin = acc.mean - acc.sd, 
                      ymax = acc.mean + acc.sd)) +
    #ggtitle("NPred * Nprey") +
    labs(y = paste("Mean", "acc", "\u00B1", "SD", sep = " "),
         x = expression(Log[10]~Threshold), parse = T) +
    theme_classic()



selected.threshold <- threshold[1]
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
