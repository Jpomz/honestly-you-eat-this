# WebBuilder
# 12 jan 18
# justin Pomeranz
# jfpomeranz@gmail.com

# make a list of registries using all Taieri predation matrices -web trying to infer
library(plyr)
library(tidyverse)
# webbuilder functions
source("Useful WebBuilder functions.R")
# useful food web functions from petchey
source("C:\\Users\\Justin\\Documents\\Data\\FW modelling Petchey Github\\ttl-resources-master\\food_web\\FoodWebFunctions.r")
# function to convert list of pairs to adjacency matrix
source("pairs_to_adj function.R")
# TSS function
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

# read in corrected adjaceny matrices
web.list <- readRDS("observed pred-prey.RDS")
web.match <- readRDS("observed matrices matched to inferred.RDS")

# convert predation matrix to matrix of pairs
pairs.list <- web.list %>%
  llply(function (x){
    as.data.frame(Matrix.to.list(x), stringsAsFactors = FALSE)})



# make list of res-con pairs for registries
# each element in list contains all pairs EXCEPT those of
# the web they will be used to parameterize. 
# e.g. element named "blackrock" will contain all other
# pairs data [-"blackrock"], and be used to paramterize # the model in order to predict the blackrock web 
register.list <- NULL
for (i in 1:length(web.match)){
  register.list[[i]] <- ldply(pairs.list)
}
# make one complete registry to add to other registers later
complete.registry <- register.list[[1]]
colnames(complete.registry) <- c("source.id", "resource", "consumer")

# add web names
names(register.list) <- names(web.match)
# subset each element in register.list to not contain name of web trying to infer links for
for (i in names(register.list)){
  register.list[[i]] <- subset(register.list[[i]], .id != i)
}

# change .id column to source.id
register.list <- register.list %>%
  llply(function (x){
    colnames(x) <- c("source.id", "resource", "consumer");x
  })

# read in registry made from other sources
other.registry <- read_csv("fw_registry.csv")
other.registry$source.id <- as.character(other.registry$source.id)

# read in file with taxonomy information
taxonomy <- read_csv("taxonomy.csv")

# add taxonomy information for resource and consumers
for (i in 1:length(register.list)){
  x <- register.list[[i]]
  x <- merge(x, taxonomy[,1:5], 
             by.x ="resource", by.y = "name", 
             all.x = T)
  colnames(x)[4:6] <- paste("res.", 
                            colnames(x[c(4:6)]), sep = "")
  x <- merge(x,taxonomy[,c(1:4,6)], by.x ="consumer",
             by.y = "name", all.x = T)         
  colnames(x)[8:10] <- paste("con.", 
                             colnames(x[c(8:10)]), sep = "")
  x$linkevidence <- "direct"
  register.list[[i]] <- x
}






# intersect of register.list and other register names
col.names <- intersect(colnames(register.list[[1]]), names(other.registry))

# subset other.registry object to match with registry list
other.registry <- other.registry[,col.names]

# combine registry list with other registry
for (i in 1:length(register.list)){
  register.list[[i]] <- bind_rows(
    register.list[[i]], other.registry)
}

# complete registry add taxonomy
complete.registry <- merge(complete.registry,
                           taxonomy[,1:5],
                           by.x ="resource",
                           by.y = "name",
                           all.x = T)
colnames(complete.registry)[4:6] <- paste("res.",
                                          colnames(complete.registry[c(4:6)]),
                                          sep = "")
complete.registry <- merge(complete.registry,
                           taxonomy[,c(1:4,6)],
                           by.x ="consumer",
                           by.y = "name",
                           all.x = T)         
colnames(complete.registry)[8:10] <- paste("con.", 
                                           colnames(complete.registry[c(8:10)]),
                                           sep = "")
complete.registry$linkevidence <- "direct"

#merge complete with others
complete.registry <- bind_rows(complete.registry, other.registry)

write.csv(complete.registry, "complete webbuilder registry.csv")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# make table of registry information for manuscript / SI
# want site, n pairwise intxn
# na pairwise intxns
n.pairs <- register.list %>%
  ldply(function (x){
    data.frame(n = nrow(x))
  })

# write csv of registry size
# modify in excel to make table for publication
write_csv(n.pairs, "C:/Users/Justin/Google Drive/Data/Predicting NZ Food Webs/figs for MS/webbuilder registry info from R.csv")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# get a list of taxa from web.list
# need to subset so only looking at webs which were also inferred using gravel method
taxa.list <- llply(web.match,
                   function (x){
                     data.frame(node = colnames(x),
                                stringsAsFactors = FALSE)
                   })

# add taxonomy information
taxa.list <- taxa.list %>% 
  llply(function (x){
    merge(x, taxonomy,by.x="node", 
          by.y="name", all.x=T)
  })



# infer links with WebBuilder function
# takes a while to run...
links <- NULL
for (i in 1:length(taxa.list)){
  links[[i]] <- WebBuilder(taxa.list[[i]],
                           register.list[[i]],
                           method = c("exact", "genus", "family", "order"))
}
names(links) <- names(web.match)

wb.matrices <- map2(taxa.list, links, pairs_to_adj)

# match wb matrices to web.match matrices
match_matr2 <- function (obj, target){
  # colnames(inf) have already been size-sorted
  index = intersect(rownames(target), rownames(obj))
  obj = obj[index, index]
  obj
}

wb.matrices <- map2(wb.matrices, web.match, match_matr2)

saveRDS(wb.matrices, file ="wb matrices matched to inferred.rds")

wb.tss <- map2(web.match, wb.matrices, get_tss)
