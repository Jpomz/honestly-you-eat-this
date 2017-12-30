# Cleaning data
# 29 Dec 2017
# Justin Pomeranz
# jfpomeranz@gmail.com

# # This script takes the raw Taieri adjacency matrices (obtained from Ross Thompson November 2016 at New Zealand Ecological Society meeting, Hamilton NZ) and translates the names to the generic level (family / subfamily for coleoptera, Diptera, etc). It converts basal resources to categorical names (e.g. diatoms, algae, etc)
# also cleans up names in Taieri_community_data.csv and calculates dry weight (dw in grams) for genera

# load plyr first to not mask tidyverrse functions
library(plyr)
library(tidyverse)

# recoder function 
# found this online somewhere, I think the reference information is in one of my oiginal webbuilder scripts
# takes a character vector from "data", and replaces "oldvalue" with "newvalue"
# I have a translation "dictionary" that I made by hand
recoderFunc <- function(data, oldvalue, newvalue) {
  # convert any factors to characters
  if (is.factor(data)) data <- as.character(data)
  if (is.factor(oldvalue)) oldvalue <- as.character
    (oldvalue)
  if (is.factor(newvalue)) newvalue <- as.character
    (newvalue)
  # create the return vector
  newvec <- data
  # put recoded values into the correct position in the return vector
  for (i in unique(oldvalue)) newvec[data == i] <-
    newvalue[oldvalue == i]
  newvec
  }


# read in original adj_matr
# make list of adj_matr
adj.list <- NULL
for(f in list.files(path="All Thompson 2004", full.names=T)){
  A = as.matrix(read.csv(f, h=T, row.names=1))
  row.names(A) <- colnames(A)
  adj.list[f] <- list(A)
}

# clean up adjacency names
# remove all text up to and inclusing "/"
names(adj.list) <- gsub(".*/", "", names(adj.list))
# remove ".csv" from names
names(adj.list) <- gsub(".csv", "", names(adj.list))

# make list of all taxa names
adj.list.names <- llply(adj.list,
                        function(x){
                          colnames(x)
                        })

# translation "dictionary" for recoder function 
translate <- read_csv("translation.csv") 
# species to genus / category data
sp.gen <- read_csv("species genus category ffg.csv")

# translate list of names
names.cor <- llply(adj.list.names, function (x){
  recoderFunc(x, translate$Wrong, translate$Corrected)})
#make list of generic names with recoderFunc
generic.names <- llply(names.cor, function (x) {
    recoderFunc(x, sp.gen$Species, sp.gen$Genus)
  })

# basal categories ####
# truncate sp.gen to basal resources
basal.gen <- sp.gen %>% 
  filter(ffg == "producer")
# truncate sp.gen 
# convert basal genera to category
new.names <- generic.names %>%
  llply(function (x){ 
    recoderFunc(x, basal.gen$Genus, basal.gen$Category)
  })
# should give warning that number of items to replace is not multiple of replacement length 
# This is because only replacing the basal.gen names


# change dimnames of adj.list to new corrected names
for (i in 1:length(adj.list)){
  dimnames(adj.list[[i]]) <- list(new.names[[i]],
                                  new.names[[i]])
}

# combine observations in duplicate names
adj.list <- adj.list %>% 
  llply(function (x){
    # convert to table, data.frame; get frequencies
    y <- as.data.frame(as.table(x)) %>%
      xtabs(Freq ~ Var1 + Var2, .)
    # convert all values > 1 in table to 1
    y[y>1] <- 1
    # convert back to matrix
    as.matrix(as.data.frame.matrix(y))
    })

# remove basal categories to make pred-prey interactions only 
# vector of basal categories, or taxa that are unknown 
basal.cat <- c("Detritus", "Terrestrial.invertebrates", "Macrophyte", "Diatom", "Algae", "Moss", "Meiofauna", "Pelicypod", "Amphora", "Plant.material")
# # integer vector of col/row numbers to remove for basal categories
# rm.basal.col.num <- llply(adj.list, function (x){
#   which(colnames(x) %in% basal.cat)
# })

# remove basal resources from adj_matr
for (i in 1:length(adj.list)){
  basal.col <- which(colnames(x) %in% basal.cat)
  adj.list[[i]] <-  adj.list[[i]][-basal.col,
                                  -basal.col]
}

# removing probable data errors
# eg, non-predatory taxa shown to be consuming prey potamopyrgus eating austrosim, deleatidum etc deleatidium consuming animal prey the mouthparts and habits of these animals make it extremely unlikely that they are consuming animals
errors <- c("Amphipoda", "Atalophlebioides", "Austroclima", "Austrosimulium", "Coloburiscus", "Deleatidium", "Nesameletus", "Oxyethira", "Potamopyrgus", "Zephlebia")


# make non=predatory taxa columns == 0
for (i in 1:length(adj.list)){
  x <- adj.list[[i]]
  forbid <- which(colnames(x) %in% errors)
  x[, forbid] <- 0
  adj.list[[i]] <- x
}

saveRDS(adj.list, "observed pred-prey.RDS")
