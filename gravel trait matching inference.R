# Trait-matching
# 31 Dec 2017
# Justin Pomeranz
#jfpomeranz@gmail.com

# paramterize Niche model using method presented in Gravel et al. 2013 Methods in Eco Evo. Inferring food web structure from predator-prey body size relationships

library(plyr)
library(tidyverse)

# gravel functions
# from Gravel et al. 2013 supplementary information
source("gravel_functions.R")

tss <- function (observed, inferred){
  # make sure adjacency matrices are same dimensions and have same colnames
  stopifnot(dim(observed) == dim(inferred), 
            identical(colnames(observed), colnames(inferred)))
  
  # subtract inferred from observed
  minus <- observed - inferred
  # multiply observed and inferred
  multiply <- observed * inferred
  # change values of 1 in multiplied matrix to 2
  multiply[multiply == 1] <- 2
  # add minus and multiply matrices
  prediction <- minus + multiply
  # prediction outcome matrix now has all 4 possibilities repreented as different integer values
  # 2 = true positive (a); links both obserevd & predicted
  # -1 = false positive (b); predicted but not observed
  # 1 = false negative (c); observed but not predicted
  # 0 = true negative (d); not predicted, not observed
    a = length(prediction[prediction==2])
    b = length(prediction[prediction==-1]) 
    c = length(prediction[prediction==1]) 
    d = length(prediction[prediction==0])
  # calculate TSS
  # TSS = (a*d - b*c)/((a+c)*(b+d))
  tss = (a*d - b*c)/((a+c)*(b+d))
  tss
}

# data ####
# invertebrate biomass and abundance
invert <- readRDS("estimated invert bodymass.RDS")
# fish biomass and abundance
fish <- readRDS("estimated fish bodymass and abundance.RDS")
# rbind invert to each element in fish, filter out na(dw), order by dw and then split each list into list of sites
dw <- llply(fish, function (x){
  y = bind_rows(x, invert)
  y = y[!is.na(y$dw),]
  y = y[order(y$dw),]
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
  pred <- as.matrix(pred) # convert vec to matrix in order to add colname
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

#training data ####
# list of lists
training.list <- NULL
for (f in 1:length(dw.pairs)){
  temp <- NULL
  for (web in 1:length(dw.pairs[[f]])){
    dat <- ldply(dw.pairs[[f]][-web])
    temp[[web]] <- dat
  }
  training.list[[f]] <- temp
  names(training.list[[f]]) <- names(dw.pairs[[f]])
}
names(training.list) <- names(dw.pairs)

# Gravel model ####
# get parameters based on training lists
# solution without using nested for loops
pars.list <- map(training.list, map, function(x){
  Bprey = log10(x$prey)
  Bpred = log10(x$pred)
  out <- reg_fn(Bprey, Bpred, quartil = c(0.03, 0.97))
})
# coef plot ####
# plot parameters for 4 different fish sizes
# funciton to extract parameters
pull_params <- function(dat){
  data.frame(B0hi = dat[[1]][1],
             B0center = dat[[2]][1],
             B0lo = dat[[3]][1],
             B1hi = dat[[1]][2],
             B1center = dat[[2]][2],
             B1lo = dat[[3]][2])
  }

# list of params for 4 fish sizes
param.coef.list <- map(pars.list, ldply,
                       function (x){
                         pull_params(x)
})

# mean/SD of paramters for 4 fish sizes
param.summary.list <- llply(param.coef.list,
                            function (x){
  data.frame(mean = apply(x[,-1],2, mean),
             sd = apply(x[,-1], 2, sd),
             coef = colnames(x)[-1])
})
# add names for plotting below
names(param.summary.list) <- names(pars.list)

# plot mean coef +- SD
ggplot(ldply(param.summary.list), aes(x = coef, y = mean, color = .id)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean -sd,
                    ymax = mean +sd))


# web params ####
# calculate web paramters using get_pars_Niche()
Ball <- map(dw, map, function (x){
  log10(x$dw)
})
web.pars <- map2(pars.list, Ball,
                 map2, get_pars_Niche)
# infer links ####

# calculate food web links for taieri
# predation matrix
web.links.inf  <- map(web.pars, map, function (x){
  L_fn(x[,"n"],
       x[,"c"],
       x[,"low"],
       x[,"high"])
})

# sum of links per web
# maybe I can use this to "pick" fish size??
sum_fn <- function (x){
  sapply(x, sum)
}
sum.links <- ldply(web.links.inf, function (x){
  sum_fn(x)})

# add taxa names to inferred matrices
web.links.inf <- map2(web.links.inf, dw, map2,
             function (x, y){
  dimnames(x) <- list(y$taxa, y$taxa)
  x
})

# taxa index ####
# make index of taxa names to match matrices
get_dim_index <- function(object, target){
  # subset object to names in target
  index <- intersect(rownames(target),
                     rownames(object))
  index
}

# taxa names index
dim.index <- map2(rep(list(obs.A),
                      length(web.links.inf)),
                  web.links.inf,
                  map2,
                  get_dim_index)
names(dim.index) <- names(web.links.inf)

# TSS ####
# step1, biomass inference
tss.step1 <- NULL
for(f in 1:length(web.links.inf)){
  temp <- NULL
  for(web in 1:length(web.links.inf[[f]])){
    index <- dim.index[[f]][[web]]
    obs <- obs.A[[web]][index, index]
    inf <- web.links.inf[[f]][[web]][index, index]
    temp[[web]] <- tss(obs, inf)
  }
  tss.step1[[f]] <-temp
}


# step2, prune niche forbidden links
# e.g. taxa that cannot eat prey due to mouthparts, scrapers, filter feeders, etc
taxa.forbid <- c("Amphipoda", "Atalophlebioides", "Austroclima", "Austrosimulium", "Blephariceridae", "Coloburiscus", "Deleatidium", "Nesameletus", "Ostracoda", "Oxyethira", "Potamopyrgus", "Zephlebia")

tss.step2 <- NULL
for(f in 1:length(web.links.inf)){
  temp <- NULL
  for(web in 1:length(web.links.inf[[f]])){
    index <- dim.index[[f]][[web]]
    obs <- obs.A[[web]][index, index]
    inf <- web.links.inf[[f]][[web]][index, index]
    for(name in (colnames(inf)[colnames(inf) %in%
                               taxa.forbid])){
      inf[,name] <- 0
    }
    temp[[web]] <- tss(obs, inf)
  }
  tss.step2[[f]] <-temp
}

unlist(lapply(tss.step1, mean))
unlist(lapply(tss.step2, mean))
unlist(lapply(tss.step1, sd))
unlist(lapply(tss.step2, sd))

get_rel_ab <- function(vec){
  rel.ab <- vec / sum(vec)
  Nij <- matrix(0, length(vec), length(vec))
  for (i in 1:length(vec)){
    for (j in 1:length(vec)){
      Nij[i,j] <- rel.ab[i]*rel.ab[j]
    }
  }
  Nij
}

rm_neutral <- function(Nij, threshold){
    Nij[Nij > threshold] <-  1
    Nij[Nij < 1] <-  0 
    Nij
}
out <- NULL
for(t in 1:length(threshold)){
  out[[t]] <- rm_neutral(nij.test, threshold[t])
}

# apply function to 1 nested list
vec <- rnorm(5, 25, 10)
vecs <- rep(list(vec), 5)
map(vecs, get_rel_ab)
vecs2 <- rep(list(vecs),2)
vecs2.out <- map(vecs2, map, get_rel_ab)

# playing around with applying fxn to 2 nested lists
Nij <- matrix(rnorm(25, 1e-4, 1e-5), ncol = 5 )
Nij2 <- matrix(rnorm(25, 1e-4, 1e-5), ncol = 5 )
test2 <- list(rep(list(Nij), 5), rep(list(Nij2), 5)) 
threshold <- rnorm(5, 9e-5, 8e06)
thresh2 <- rep(list(threshold), 2)
out <- map2(test2, thresh2, map2, rm_neutral)
