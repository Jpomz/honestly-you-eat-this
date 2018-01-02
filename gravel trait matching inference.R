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

get_tss <- function (observed, inferred){
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
fish <- readRDS("estimated fish bodymass and abundance.RDS")[[2]]
# rbind invert to each element in fish, filter out na(dw), order by dw and then split each list into list of sites
dw <- bind_rows(fish, invert)
dw <- dw[!is.na(dw$dw),]
dw <- dw[order(dw$dw),]
dw <- split(dw, list(dw$site))
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
dw.pairs <- NULL
  for(web in 1:length(dw)){
    m.pred <- merge(A.pairs[[web]]$pred,
                    dw[[web]][,1:2],
                    by = "taxa",
                    all.x = T)
    m.prey <- merge(A.pairs[[web]]$prey,
                    dw[[web]][,1:2],
                    by = "taxa",
                    all.x = T)
    m.pairs <- cbind(m.pred[,2], m.prey[,2])
    colnames(m.pairs) <- c("pred", "prey")
    dw.pairs[[web]] <- m.pairs[!is.na(m.pairs[,1]) &
                        !is.na(m.pairs[,2]),]
  }
names(dw.pairs) <- names(obs.A)

#training data ####
# list of lists
training.list <- NULL
for (web in 1:length(dw.pairs)){
    dat <- ldply(dw.pairs[-web])
    training.list[[web]] <- dat
  }
names(training.list) <- names(dw.pairs)

# Gravel model ####
# get parameters based on training lists
# solution without using nested for loops
pars.list <- map(training.list, function(x){
  Bprey = log10(x$prey)
  Bpred = log10(x$pred)
  out <- reg_fn(Bprey, Bpred, quartil = c(0.01, 0.97))
})
# web params ####
# calculate web paramters using get_pars_Niche()
Ball <- map(dw, function (x){
  log10(x$dw)
})
web.pars <- map2(pars.list, Ball, get_pars_Niche)
# infer links ####

# calculate food web links for taieri
# predation matrix
web.links.inf  <- map(web.pars, function (x){
  L_fn(x[,"n"],
       x[,"c"],
       x[,"low"],
       x[,"high"])
})
# add taxa names to inferred matrices
web.links.inf <- map2(web.links.inf, dw,
          function (x, y){
            dimnames(x) <- list(y$taxa, y$taxa)
            x
                      })
# sum of links per web
sum.links <- sapply(web.links.inf, sum)

# TSS ####
# function to calc TSS from matrices matched by row/colnames
# be careful to select appropriate obs and inf
match_matr_tss <- function (obs, inf, forbidden.taxa = NULL, ab.vec = NULL, ab.taxa = NULL, ab.threshold = NULL){
  # colnames(inf) have already been size-sorted
  index = intersect(rownames(inf), rownames(obs))
  obs = obs[index, index]
  inf = inf[index, index]
    if(is.character(forbidden.taxa)){
      for(name in (colnames(inf)[colnames(inf) %in% forbidden.taxa])){
      inf[,name] <- 0
      }
      # else(warning("\n***************\nNo forbidden taxa supplied\nTSS calculated for unmodified inf object\n***************"))
    }
  if(is.numeric(ab.vec)){
    Nij = get_rel_ab(ab.vec, ab.taxa)
    Nij = rm_neutral(Nij, ab.threshold)
    Nij = Nij[index, index]
    inf = inf * Nij
  }
  result = get_tss(obs, inf)
  result
}

# step1, biomass inference
tss.step1 <- map2(obs.A, web.links.inf,
                  match_matr_tss)

# step2, prune niche forbidden links
# e.g. taxa that cannot eat prey due to mouthparts, scrapers, filter feeders, etc
taxa.forbid <- c("Amphipoda", "Atalophlebioides", "Austroclima", "Austrosimulium", "Blephariceridae", "Coloburiscus", "Deleatidium", "Nesameletus", "Ostracoda", "Oxyethira", "Potamopyrgus", "Zephlebia")

tss.step2 <- pmap(list(obs = obs.A,
                  inf = web.links.inf),
                  match_matr_tss,
                  forbidden.taxa = taxa.forbid)


# calculate relative abundance matrices
get_rel_ab <- function(vec, taxa){
  stopifnot(length(vec) == length(taxa))
  rel.ab <- vec / sum(vec)
  Nij <- matrix(0, length(vec), length(vec))
  for (i in 1:length(vec)){
    for (j in 1:length(vec)){
      Nij[i,j] <- rel.ab[i]*rel.ab[j]
    }
  }
  dimnames(Nij) <- list(taxa, taxa)
  Nij
}
# function to remove interactions < threshold
rm_neutral <- function(Nij, threshold){
    Nij[Nij > threshold] <-  1
    Nij[Nij < 1] <-  0 
    Nij
}

threshold <- c(1e-03, 1e-04, 1e-05, 1e-06, 1e-07, 1e-08, 1e-09, 5.9e-03, 5.9e-04, 5.9e-05, 5.9e-06, 5.9e-07, 5.9e-08, 5.9e-9, 3e-03, 3e-04, 3e-05, 3e-06, 3e-07, 3e-8, 3e-9, 1e-10, 3e-10, 5.9e-10, 1e-11, 3e-11, 5.9e-11, 1e-12, 3e-12, 5.9e-12, 1e-13, 3e-13, 5.9e-13,1e-14, 3e-14, 5.9e-14,1e-15, 3e-15, 5.9e-15, 1e-2, 3e-2, 5.9e-2)

ab.vec <- llply(dw, function (x){x$dw})
ab.taxa <- llply(dw, function (x){x$taxa})
test <- map(threshold, function (x){
                 pmap(list(obs = obs.A,
                      inf = web.links.inf,
                      ab.vec = ab.vec,
                      ab.taxa = ab.taxa),
                 match_matr_tss,
                 forbidden.taxa = taxa.forbid,
                 ab.threshold = x)
  }
)

#***************************************************************
# clean up "test" object, select local and global thresholds
# try and correct for fish abundance

#***************************************************************
test <- ldply(flatten(test))
test$threshold <- rep(threshold, each = 17)
ggplot(test, aes(x = log10(threshold), y = V1, color = .id)) +
  geom_point() +
  stat_smooth(alpha = 0.2)
test %>% group_by(.id) %>% top_n(1,wt = V1) %>% mutate(log10(threshold))

# fish abundance "correction"
fish.corr <- NULL
for (f in 1:length(web.links.inf)){
  t.temp <- NULL
  for (t in 1:length(threshold)){
    web.temp <- NULL
    for(web in 1:length(web.links.inf[[f]])){
      index <- dim.index[[f]][[web]]
      Nij <- get_rel_ab(vec = dw[[f]][[web]]$dw,
                        taxa = dw[[f]][[web]]$taxa)
      f.vec <- c("Salmo", "Galaxias", "Anguilla", "Gobiomorpus")
      for(fv in (colnames(Nij)[colnames(Nij) %in%                   f.vec])){
        Nij[,fv] <- Nij[,fv]*1e10
      }
      Nij <- rm_neutral(Nij, threshold[t])
      Nij <- Nij[index, index]
      inf <- web.links.inf[[f]][[web]][index, index]
      obs <- obs.A[[web]][index, index]
      for(name in (colnames(inf)[colnames(inf) %in%                   taxa.forbid])){
        inf[,name] <- 0
      }
      inf <- inf * Nij
      web.temp[[web]] <- tss(obs, inf)
    }
    t.temp[[t]] <- web.temp
  }
  names(t.temp) <- threshold
  fish.corr[[f]] <- t.temp
}
names(fish.corr) <- names(web.links.inf)

fish.corr.summary <- map(fish.corr, ldply, c(mean, sd))

fish.corr.summary <- ldply(fish.corr.summary)
fish.corr.summary$size <- c(rep("min.min", length(threshold)), rep("mean.min", length(threshold)),rep("mean.max", length(threshold)),rep("max.max", length(threshold)))

ggplot(fish.corr.summary, aes(x = log10(as.numeric(.id)), y = V1, color = size)) +
   geom_point() #+
  # geom_errorbar(aes(ymin = V1 - V2,
  #               ymax = V1 +V2))


unlist(lapply(tss.step1, mean))
unlist(lapply(tss.step2, mean))

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
