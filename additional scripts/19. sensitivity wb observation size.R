# 19. sensitivity WB observation size
# JPomz 30 Nov 2017

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
# pairs to adjacency function
source("pairs_to_adj function.R")

# Empirical web ####
# Dempsters empirical web
obs.demp <- as.matrix(readRDS("animal only adj list.rds")[["Dempsters"]])
# list of taxa in Dempsters, ordered by body size
target <- readRDS("target vector to match all matrices.RDS")[[7]] # list not named, but double checked that [[7]] == Dempsters
# order rows/cols by increasing body size
# of empirical adjacency matrix
obs.demp <- obs.demp[match(target,
                           rownames(obs.demp)),
                     match(target,
                           colnames(obs.demp))]

# prep for WebBuilder ####
# make df of "nodes"
nodes <- data.frame(node = target,
                    stringsAsFactors = FALSE)
# csv with taxonomy information
taxonomy <- read_csv("taxonomy.csv")
# merge nodes with taxonmoy
nodes <- merge(nodes,
               taxonomy,
               by.x = "node",
               by.y = "name",
               sort = FALSE)
# paste "minimum" to res./con.method
# (requirement of WebBuilder)
colnames(nodes)[5:6] <- paste("minimum.",
                              colnames(nodes)[5:6],
                              sep = "")
# read in complete registry
complete.registry <- read_csv("complete webbuilder registry.csv")
# filter out Dempsters from registry
registry <- as.data.frame(complete.registry[
  complete.registry$source.id!="Dempsters" &
  complete.registry$source.id!="DempstersSp" &
  complete.registry$source.id!="DempstersAu",])


# sensitivity analysis ####
#vector of proportions
prop <- seq(.01,1,.01)
# number of iterations
niter <- 50
# empty object for results
result <- NULL
# for loop
# each proportion (.01-1) * niter (50)
system.time(
  for(p in seq_along(prop)){ # each proportion value
    tss.prop <- numeric(niter)
    for(i in 1:niter){ # number of iterations
      # sample p fraction of registry
      regist <- registry[sample(nrow(registry),
                         prop[p]*nrow(registry)),]
      # infer links with fractioned regist
      link <- WebBuilder(nodes, regist,
                        methods = c("exact",
                                    "genus",
                                    "family", 
                                    "order"))
      # convert links pairs to adjacency matrix
      # necessary for tss()
      adj <- pairs_to_adj(nodes, link)
      # calc tss for each iteration
      tss.prop[i] <- tss(obs.demp, adj)[[6]]
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
# calc number of obs col
db.size <- ceiling(as.numeric(result.df$.id) * 1252)
# add number of obs to df
result.df$n.obs <- db.size
# save results
saveRDS(result.df, "WebBuilder TSS db size df.rds")
# readRDS("WebBuilder TSS db size df.rds")

# plot ####
# make plot of tss ~ proportion of links in db
(tss.db.size <- ggplot(result.df,
                       aes(x = db.size,#as.numeric(.id),
                           y = V1))+
    geom_line(size = 1) +
    geom_ribbon(aes(x = db.size, #as.numeric(.id),
                    ymin = V1 - V2,
                    ymax = V1 + V2),
                alpha = 0.5) +
    labs(x = "Number of links in database",
         y = "TSS") +
    theme_classic()
)
# save plot
saveRDS(tss.db.size, "WebBuilder TSS db size plot.rds")
# readRDS("Gravel TSS db size plot.rds")




