# Trait-matching
# 31 Dec 2017
# Justin Pomeranz
#jfpomeranz@gmail.com

# paramterize Niche model using method presented in Gravel et al. 2013 Methods in Eco Evo. Inferring food web structure from predator-prey body size relationships

library(plyr)
#library(tidyverse)

# gravel functions
# from Gravel et al. 2013 supplementary information
source("gravel_functions.R")

# data ####
# invertebrate biomass and abundance
invert <- readRDS("estimated invert bodymass.RDS")
# fish biomass and abundance
fish <- readRDS("estimated fish bodymass and abundance.RDS")
# rbind invert to each element in fish and then split each list into list of sites
dw <- llply(fish, function (x){
  y = rbind(x, invert)
  #y = y[!is.na(y$dw),]
  y = split(y, list(y$site))
  y
})

# observed webs
obs.A <- readRDS("observed pred-prey.RDS")

# subset obs_A to only include webs with community data
obs.A <- obs.A[names(obs.A) %in% unique(invert$site)]

# pred-prey pairs ####
# empty list
A.pairs <- NULL
# for each web in obs_A
for(web in 1:length(obs.A)){
  A = obs.A[[web]]
  genus = colnames(A) # vector of taxa names
  pred = c() # empty vector for pred
  prey = c() # empty vector for prey
  for(i in 1:nrow(A)){ # i = row = prey
    for(j in 1:ncol(A)){ # j = column = pred
      if(A[i,j]==1){ # if there is a link
        pred = c(pred, genus[j]) # add col taxa to pred
        prey = c(prey, genus[i]) # add row taxa to prey
      }
    }
    
  }
  pred <- as.matrix(pred)
  colnames(pred) <- "taxa"
  prey <- as.matrix(prey)
  colnames(prey) <- "taxa"
  pairs <- list(pred = pred, prey = prey) 
  A.pairs[[web]] <- pairs 
}
names(A.pairs) <- names(obs.A)

# pairs biomass ####
# list of 4 fish sizes
  # each fish size has 17 elements
    # each element is paired biomass for pred-prey
dw.pairs <- NULL
for(f in 1:length(dw)){
  temp <- NULL
  for(web in 1:length(dw[[f]])){
    m.pred <- merge(A.pairs[[web]]$pred,
                    dw[[f]][[web]][,1:2],
                    by = "taxa",
                    all.x = T)
    m.prey <- merge(A.pairs[[web]]$prey,
                    dw[[f]][[web]][,1:2],
                    by = "taxa",
                    all.x = T)
    m.pairs <- cbind(m.pred[,2], m.prey[,2])
    colnames(m.pairs) <- c("pred", "prey")
    temp[[web]] <- m.pairs[!is.na(m.pairs[,1]) &
                        !is.na(m.pairs[,2]),]
  }
  names(temp) <- names(A.pairs)
  dw.pairs[[f]] <- temp
}
names(dw.pairs) <- names(dw)




