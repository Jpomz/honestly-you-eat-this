# Relative abundances
# modifying webbuilder inferred links with relative abundances correction
# JPomz
# june 2017


# useful food web functions from Petchey
source("C:\\Users\\Justin\\Documents\\Data\\FW modelling Petchey Github\\ttl-resources-master\\food_web\\FoodWebFunctions.r")
# TSS function
source("TSS function.R")


# read in estimated Taieri dry weights 
taieri <- readRDS("estimated dw taieri webs.rds")

# note for when I want to calc relative abundances
rel.ab <- taieri %>% 
  select(taxa, no.m2) %>%
  group_by(site) %>% 
  mutate(tot.ab = sum(no.m2, na.rm = T), rel.ab = no.m2 / tot.ab) 



#rel.ab.in <- readRDS("rel.ab with only invert.RDS")

# don't have an estimate for fish abundances
# species specific correction factor ####
# need to do some research to inform correction factor value
rel.ab[rel.ab$taxa == "Salmo",]$rel.ab <- 
  rel.ab[rel.ab$taxa == "Salmo",]$rel.ab * 10^5

rel.ab[rel.ab$taxa == "Gobiomorphus",]$rel.ab <- 
  rel.ab[rel.ab$taxa == "Gobiomorphus",]$rel.ab * 10^5

rel.ab[rel.ab$taxa == "Anguilla",]$rel.ab <- 
  rel.ab[rel.ab$taxa == "Anguilla",]$rel.ab * 10^5

rel.ab[rel.ab$taxa == "Galaxias",]$rel.ab <- 
  rel.ab[rel.ab$taxa == "Galaxias",]$rel.ab * 10^5

# Arbitrary ####
# making an arbitrary value for rel.ab
# fish are highly mobile, and in nearly all 
# webs fish consume ~ all prey
# therefore making a relatively high rel.ab
# fish abundance sensitivity indicates that any abundance > 0.20 equivalent TSS score
#rel.ab$rel.ab[is.na(rel.ab$rel.ab)] <- 0.2


# make into a list
rel.ab.list <- split(rel.ab, list(rel.ab$site))


# make list of abundance matrices
ab.matr.list <- rel.ab.list %>%
  llply(function (x){
    matrix(0, nrow(x), nrow(x))
  })

# add dimnames to matrices
for (i in 1:length(ab.matr.list)){
  dimnames(ab.matr.list[[i]]) <- list(rel.ab.list[[i]]$taxa,
                                      rel.ab.list[[i]]$taxa)
}

# calculate relative abundances
for (k in 1:length(ab.matr.list)){
  for (i in 1:nrow(rel.ab.list[[k]])){
      for (j in 1:nrow(rel.ab.list[[k]])){
        ab.matr.list[[k]][i,j] <- 
          rel.ab.list[[k]]$rel.ab[i] * rel.ab.list[[k]]$rel.ab[j]
      }
  }
}
  
# # read in webbuilder inferred links, trimmed to match gravel inf
# wb.trim.list <- readRDS("webbuilder inferred list trimmed taxa.rds")
# # read in 10 observed foodweb
# obs.list <- readRDS("10 observed taieri food webs.rds")


# target vector to match all matrices
# this is also used to trim the matrices to only include taxa in this list
# e.g. matches the taxa in the gravel inferred webs
target <- readRDS("target vector to match all matrices.RDS")

# make ab.matr.list match the target vector
# both trims rosw/cols and orders to target vector
for (i in 1:length(ab.matr.list)){
  ab.matr.list[[i]] <- ab.matr.list[[i]][target[[i]],
                                         target[[i]]]
}


# save abundance matrices to use use for gravel inferred
# Npred * Nprey
saveRDS(ab.matr.list, "17 relative abundance matrices list.rds")


# test out different rel abundance models
# Npred * Nprey**2

# make list of abundance matrices
ab.matr.list <- rel.ab.list %>%
  llply(function (x){
    matrix(0, nrow(x), nrow(x))
  })

# add dimnames to matrices
for (i in 1:length(ab.matr.list)){
  dimnames(ab.matr.list[[i]]) <- list(rel.ab.list[[i]]$taxa,
                                      rel.ab.list[[i]]$taxa)
}

# calculate relative abundances
for (k in 1:length(ab.matr.list)){
  for (i in 1:nrow(rel.ab.list[[k]])){
    for (j in 1:nrow(rel.ab.list[[k]])){
      ab.matr.list[[k]][i,j] <- 
        rel.ab.list[[k]]$rel.ab[i]^2 * rel.ab.list[[k]]$rel.ab[j]
    }
  }
}

# # read in webbuilder inferred links, trimmed to match gravel inf
# wb.trim.list <- readRDS("webbuilder inferred list trimmed taxa.rds")
# # read in 17 observed foodweb
# obs.list <- readRDS("17 observed taieri food webs.rds")


# target vector to match all matrices
# this is also used to trim the matrices to only include taxa in this list
# e.g. matches the taxa in the gravel inferred webs
#target <- readRDS("target vector to match all matrices.rds")

# make ab.matr.list match the target vector
# both trims rosw/cols and orders to target vector
for (i in 1:length(ab.matr.list)){
  ab.matr.list[[i]] <- ab.matr.list[[i]][target[[i]],
                                         target[[i]]]
}


# save abundance matrices to use use for gravel inferred
saveRDS(ab.matr.list, "17 relative abundance matrices prey2.rds")




# test out different rel abundance models
# Nprey

# make list of abundance matrices
ab.matr.list <- rel.ab.list %>%
  llply(function (x){
    matrix(0, nrow(x), nrow(x))
  })

# add dimnames to matrices
for (i in 1:length(ab.matr.list)){
  dimnames(ab.matr.list[[i]]) <- list(rel.ab.list[[i]]$taxa,
                                      rel.ab.list[[i]]$taxa)
}

# calculate relative abundances
for (k in 1:length(ab.matr.list)){
  for (i in 1:nrow(rel.ab.list[[k]])){
    for (j in 1:nrow(rel.ab.list[[k]])){
      ab.matr.list[[k]][i,j] <- 
        rel.ab.list[[k]]$rel.ab[i] * 1
    }
  }
}

# # read in webbuilder inferred links, trimmed to match gravel inf
# wb.trim.list <- readRDS("webbuilder inferred list trimmed taxa.rds")
# # read in 17 observed foodweb
# obs.list <- readRDS("17 observed taieri food webs.rds")


# target vector to match all matrices
# this is also used to trim the matrices to only include taxa in this list
# e.g. matches the taxa in the gravel inferred webs
#target <- readRDS("target vector to match all matrices.rds")

# make ab.matr.list match the target vector
# both trims rosw/cols and orders to target vector
for (i in 1:length(ab.matr.list)){
  ab.matr.list[[i]] <- ab.matr.list[[i]][target[[i]],
                                         target[[i]]]
}


# save abundance matrices to use use for gravel inferred
saveRDS(ab.matr.list, "17 relative abundance matrices prey only.rds")
