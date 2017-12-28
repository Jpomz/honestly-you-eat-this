# Cleaning data

# # This script takes the raw Taieri adjacency matrices (obtained from Ross Thompson November 2016 at New Zealand Ecological Society meeting, Hamilton NZ) and translates the names to the generic level (family / subfamily for coleoptera, Diptera, etc). It converts basal resources to categorical names (e.g. diatoms, algae, etc)

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
# look at new names to make sure replacement is correct


# change dimnames of adj.list to new corrected names
for (i in 1:length(adj.list)){
  dimnames(adj.list[[i]]) <- list(new.names[[i]],
                                  new.names[[i]])
}

# compile duplicate names
test <- adj.list %>% 
  llply(function (x){
    # convert to data.frame to get frequencies
    y <- as.data.frame(as.table(x)) %>%
      xtabs(Freq ~ Var1 + Var2, .)
    # convert all values > 1 in table 
    y[y>1] <- 1
    # convert back to matrix
    as.matrix(as.data.frame.matrix(y))
    })


