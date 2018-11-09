# Functions for: Pomeranz, Thompson, Poisot, Harding. (in press). Inferring predator-prey interactions in food webs. Methods Ecology and Evolution
# Justin Pomeranz
# jfpomeranz@gmail.com


# calculate TSS
get_tss <- function (observed, inferred){
  # observed and inferred adjacency matrices
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
  tss = (a*d - b*c)/((a+c)*(b+d))
  tss
}

# function to calculate relative abundance matrices
# vec = vector of abundances
# taxa = vector of taxa names in same order as abundance vector
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

# function to remove rel.abundance products < threshold
# Nij = relative abundance matrix
# threshold = numeric value
rm_neutral <- function(Nij, threshold){
  Nij[Nij > threshold] <-  1
  Nij[Nij < 1] <-  0 
  Nij
}

# function to remove inferred links in named columns that match vector of niche forbidden taxa
# inf = inferred adjacency matrix with named columns
# taxa = vector of niche forbidden taxa
# must make sure names in 'inf' match names in 'taxa' !!!
rm_niche <- function(inf, taxa){
  for(name in (
    colnames(inf)[colnames(inf) %in% taxa])){
    inf[,name] <- 0
  }
  inf
}

# match order of named columns and rows between two matrices 
# obs and inf = adjacency matrices
# This function mathces obs to inf, make sure inf has already been sorted as needed (e.g. by body size)
match_matr <- function (obs, inf){
  index = intersect(rownames(inf), rownames(obs))
  obs = obs[index, index]
  inf = inf[index, index]
  list(observed = obs, inferred = inf)
}

# modified function to match webbuilder inferred matrices to trait-matching inferred matrices
match_matr2 <- function (obj, target){
  # columns in inf have already been size-sorted
  index = intersect(rownames(target), rownames(obj))
  obj = obj[index, index]
  obj
}

# AUC logistic model
# function to calculate AUC 
# requires the ROCR package
# see help for 'ROCR' package for more details
get_auc <- function(observed, inferred){
  # observed and inferred are adjacency matrices
  # need to be identical dimensions and dimnames
  stopifnot(length(observed) == length(inferred))
  stopifnot(colnames(observed) == colnames(inferred))
  require(ROCR)
  # convert observed adjacency matrix to vector of factors
  y = as.factor(as.numeric(observed))
  # convert inferred adjacency matrix to numeric vector
  x = as.numeric(inferred)
  # if all values in inferred matrix are 0 or 1, auc == NA
  if(length(levels(x))== 1){ 
    auc = NA
    return(auc)
  }
  # logistic regression
  mod = glm(y ~ x, family = binomial(link = "logit"))
  # model prediction
  mod.pred = predict(mod, x = x, type = "response")
  # create prediction object (necessary for ROCR package)
  prob = prediction(mod.pred, y)
  # calculate AUC from prediction object
  auc = performance(prob, measure = "auc")@y.values[[1]]
  return(auc)
}

# function to get list of AUC for site:threshold combination
get_auc_list <- function(neutral.list, obs){
  #neutral.list must be a nested list with 2 levels. 
  # level 1 = threshold to be tested
  # level 2 = matrix for each site
  out <- NULL
  for(web in 1:length(obs)){
    auc.web <- NULL
    for(t in 1:length(neutral.list)){
      auc.web[[t]] <- get_auc(obs[[web]], 
                              neutral.list[[t]][[web]])
    }
    out[[web]] <- auc.web
  }
  return(out)
}

get_max_auc <- function(x, threshold, site){
  # this is a helper function to clean up analysis
  # it takes a list of results from get_auc() function, calculates the mean AUC by threshold, and selects threshold which results in greatest AUC. It then provides a numeric value to select list element which results in greatest AUC
  # dependencies = dplyr, purrr
  # x = list, output from get_auc_list()
  # threshold = vector of thresholds tested in get_auc_list()
  # site = vector of site names
  auc <- flatten_dbl(x)
  
  stopifnot(length(auc) == length(threshold) * length(site))
  
  if(class(threshold) == "character")
    threshold <- as.numeric(threshold)
  
  df <- data.frame(auc = auc,
                   thresh = threshold,
                   site = rep(site, each = length(threshold)),
    stringsAsFactors = FALSE)
  top.auc <- df %>%
    group_by(thresh) %>%
    summarize(mean.auc = mean(na.omit(auc))) %>%
    top_n(1, wt = mean.auc)
  auc.value <- as.numeric(top.auc[2])
  thresh.value <- as.numeric(top.auc[1])
  list.element <- which(threshold == thresh.value)
  return(list(auc = auc.value, threshold = thresh.value, list.element = list.element))
}
  
# calculate proportion of false positive and false negative predictions
false_prop <- function(obs, inf){
  # obs and inf are adjacency matrices
  # check that they are equal dimensions and have same colnames
  stopifnot(dim(obs) == dim(inf), 
            identical(colnames(obs), colnames(inf)))
  
  # subtract inferred from observed
  minus <- obs - inf
  # multiply observed and inferred
  multiply <- obs * inf
  # change values of 1 in multiplied matrix to 2
  multiply[multiply == 1] <- 2
  # add minus and multiply matrices
  prediction <- minus + multiply
  # prediction outcome matrix now has all 4 possibilities repreented as different integer values
  # 2 = true positive (a); links both obserevd & predicted
  # -1 = false positive (b); predicted but not observed
  # 1 = false negative (c); observed but not predicted
  # 0 = true negative (d); not predicted, not observed
  vars <- data.frame(
    a = length(prediction[prediction==2]),
    b = length(prediction[prediction==-1]), 
    c = length(prediction[prediction==1]), 
    d = length(prediction[prediction==0]),
    S = ncol(prediction)
  )
  
  wrong <- vars %>%
    transmute(fp = b / S**2,
              fn = c / S**2)
  wrong
}

# function to make adjacency matrix from dataframe of pairwise species interactions (e.g predator and prey cols).  
pairs_to_adj <- function(taxa, pairs){
  # taxa is a data.frame which must have a "node" column and contains a single row for each taxa present
  # if duplicated nodes, could use unique(taxa$node)
  # pairs is a data frame that has all biotic interactions in paired "resource" and "consumer" columns
  S = nrow(taxa)
  A = matrix(0, nrow = S, ncol = S)
  dimnames(A) = list(taxa$node, taxa$node)
  for (i in 1:nrow(pairs))
    A[pairs$resource[i], pairs$consumer[i]] = 1
  A[match(taxa, rownames(A)),
    match(taxa, colnames(A))]
  return(A)
}

# fish abundance "correction"
f_ab_corr <- function(Nij, taxa, cf){
  for(f in which(colnames(Nij) %in% taxa)){
    Nij[,f] <- Nij[,f]*cf
  }
  Nij
}

# get_pcdat() 
# this function calculates the food web measures for the PCA analysis. 
# x = list of adjacency matrices 
# e.g. all observed adjacency matrices
# Get.web.stats calculates a suite of food web measures (see Foodweb_functions.R for details)
# This function calculates all stats for each adjacency matrix within a collection, as a list
# it then collapses the list to a dataframe (e.g. ldply())
# It then subsets to the variables we are interested in
# i.e. [,c(1, 3:7, 11:12)]
get_pcdat <- function(x){
  ldply(x, function(z){
    Get.web.stats(z, which.stats = 1)
  })[,c(1, 3:7, 11:12)]
}

# pca distance functions
# calculate centroid of cloud of points
centroid <- function(x) rowMeans(x)
# calculate distance between 2 centroids
distance2 <- function(x,y) sum((x-y)^2)

