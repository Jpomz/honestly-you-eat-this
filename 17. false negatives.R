# script 17, false negatives
# JPOMZ october 2017

# script examining the false negatives from food webs inferred using the Gravel et al. 2013 trait-matching and Gray et al. 2015 WebBuilder methods

source("false.neg function.R")

grav <- readRDS("gravel inferred pruned rel ab.rds")
wb <- readRDS("webbuilder inferred list trimmed taxa.rds")
obs <- readRDS("observed adj matr list trimmed taxa.rds")

grav.false <- NULL
for (i in 1:length(grav)){
  grav.false[[i]] <- false.neg(obs[[i]], grav[[i]])
}
names(grav.false) <- names(grav)

wb.false <- NULL
for (i in 1:length(wb)){
  wb.false[[i]] <- false.neg(obs[[i]], wb[[i]])
}
names(wb.false) <- names(wb)

for(i in 1:length(grav.false)){
  print(list(grav.false[[i]], wb.false[[i]]))
}
