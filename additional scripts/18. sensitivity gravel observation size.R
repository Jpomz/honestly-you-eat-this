# 18. sensitivity analysis Gravel database size
# JPomz 29 Nov 2017

# this script investigates the sensitivity of
# results obtained when removing varying proportions
# of the pairwise interactions in the global 
# data registry

# run analyses on largest web - Dempsters
# remove varying proportions of links 0-99%
# run each simulation 50 times

# libraries / functions ####
# TSS function
source("TSS function.R")
# gravel functions
source("gravel_functions.R")
# WebBuilder function
source("Useful WebBuilder functions.R")
#fw fxnx Petchey
source("C:\\Users\\Justin\\Documents\\Data\\FW modelling Petchey Github\\ttl-resources-master\\food_web\\FoodWebFunctions.r")

# Empirical web ####
# Dempsters empirical web
obs.demp <- as.matrix(readRDS("animal only adj list.rds")[["Dempsters"]])
# list of taxa in Dempsters, ordered by body size
target <- readRDS("target vector to match all matrices.RDS")[[7]] # list not named, but double checked that [[7]] == Dempsters
# order rows/cols by increasing body size
# of empirical web
obs.demp <- obs.demp[match(target,
                           rownames(obs.demp)),
                     match(target,
                           colnames(obs.demp))]

# training data ####
# complete data training list for gravel method
training.list <- readRDS("complete gravel training list.rds")
# vector of resource logdw from dempsters 
# training data
# the Dempsters list has already had the empirical
# links from the dempsters web removed
# e.g. script 4 
Bprey <- training.list[["Dempsters"]]$logres
# vector of consumer logdw from dempsters training data
Bpred <- training.list[["Dempsters"]]$logcon

# body size data ####
# body size information for taieri food webs
# read in estimated Taieri dry weights 
taieri <- as.data.frame(readRDS("estimated dw taieri webs.rds"))
# subset Dempsters data
demp.B <- taieri[taieri$site=="Dempsters",]
# subset only taxa in target vector
demp.B <- demp.B[demp.B$taxa %in% target, ]
# vector of all logdw from Dempsters
Ball <- sort(demp.B$logdw)

# forbidden taxa ####
# read in data with FFG information
ffg <- read_csv("taxa.ffg.csv")
# join taxa target with FFG
target.ffg <- left_join(as.data.frame(target),
                        ffg,
                        c("target" = "name"))
# make a list of "forbidden taxa"
forbid.taxa <- target.ffg %>%
  filter(FFG == "SC" | FFG == "G" | FFG == "FF") %>% select(target)
# make index of forbidden columns
forbid.col <- which(target %in% forbid.taxa$target)

# neutral processes ####
# read in matrices with relative abundances
ab.matr <- readRDS("17 relative abundance matrices list.rds")[["Dempsters"]]
# "turn on" pairwise products > 5.9e-05
# threshold determined in script 8
ab.matr[ab.matr>5.9e-05] <- 1
# turn off all other pairs
ab.matr[ab.matr<1] <- 0

# sensitivity analysis ####
#vector of proportions
prop <- seq(.01,1,.01)
# number of iterations
niter <- 50

# empty object for results
result <- NULL

# for loop
system.time(
for(p in seq_along(prop)){ # each proportion value
  tss.prop <- numeric(niter)
  for(i in 1:niter){ # number of iterations
    # make an index for sampling data
    row.prop <- sample(length(Bprey),
                       prop[p]*length(Bprey))
    # quantile regression
    reg <- reg_fn(Bprey[row.prop],
                  Bpred[row.prop],
                  quartil = c(0.01,0.97))
    # niche paramters from quant reg
    niche <- get_pars_Niche(reg, Ball)
    # infer links
    links <- L_fn(niche[,"n"],
                  niche[,"c"],
                  niche[,"low"],
                  niche[,"high"])
    # name dims of links (for tss() below)
    dimnames(links) <- list(target, target)
    # niche forbidden taxa
    links[,forbid.col] <- 0
    # neutral forbidden
    links <- links * ab.matr
    # calc tss for each iteration
    tss.prop[i] <- tss(obs.demp, links)[[6]]
  }
  # mean tss for each proportion
  tss.mean <- mean(tss.prop)
  # sd tss for each proportion
  tss.sd <- sd(tss.prop)
  # combine mean and sd, and put in list
  result[[p]] <- c(tss.mean, tss.sd)
})

# name result list elements
names(result) <- as.character(prop)
# collapse result into df
result.df <- ldply(result)
# save results
saveRDS(result.df, "Gravel TSS db size df.rds")
# readRDS("Gravel TSS db size df.rds")

# make plot of tss ~ proportion of links in db
(tss.db.size <- ggplot(result.df,
       aes(x = as.numeric(.id),
           y = V1))+
  geom_line(size = 1) +
  geom_ribbon(aes(x = as.numeric(.id),
                  ymin = V1 - V2,
                  ymax = V1 + V2),
              alpha = 0.5) +
  labs(x = "Proportion of links in database",
       y = "TSS") +
  theme_classic()
)
# save plot
saveRDS(tss.db.size, "Gravel TSS db size plot.rds")
# readRDS("Gravel TSS db size plot.rds")
