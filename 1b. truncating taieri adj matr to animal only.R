# truncating taieri adj_matr to just pred:prey
# J POMZ
# June 2017

# cut out basal categories, turn adj_matr to just animals
# start with list of adj matrices corrected in "fixing taxa...R" script
# remove basal categories
# save resulting adjacency matrices as animal only 
# file = "animal only adj list.rds"

# useful food web functions from Petchey
source("C:\\Users\\Justin\\Documents\\Data\\FW modelling Petchey Github\\ttl-resources-master\\food_web\\FoodWebFunctions.r")
# gravel functions
source("gravel_functions.R")

# read in adj.list corrected in "fixing taxa...R" script
animal.list <- readRDS("list adj matr corr names basal cat.rds")

# vector of basal categories
# or taxa that are unknown 
basal.cat <- c("Detritus", "Terrestrial.invertebrates", "Macrophyte", "Diatom", "Algae", "Moss", "Meiofauna", "Pelicypod", "Amphora", "Plant.material")
# integer vector of col/row numbers to remove for basal categories
rm.basal.col.num <- llply(animal.list, function (x){
  which(colnames(x) %in% basal.cat)
})

# remove basal resources from adj_matr
for (i in 1:length(animal.list)){
  animal.list[[i]] <-  animal.list[[i]][-rm.basal.col.num[[i]], -rm.basal.col.num[[i]]]
}
# removing probable data errors
# eg, non-predatory taxa shown to be consuming prey
# potamopyrgus eating austrosim, deleatidum etc
# deleatidium consuming animal prey
# the mouthparts and habits of these animals make it extremely unlikely that they are consuming animals
errors <- c("Amphipoda", "Atalophlebioides", "Austroclima", "Austrosimulium", "Coloburiscus", "Deleatidium", "Nesameletus", "Oxyethira", "Potamopyrgus", "Zephlebia")

# make a list of columns which match forbidden taxa
forbid.col.list <- animal.list %>% 
  llply(function (x){
    which(colnames(x) %in% errors)
  })

# make forbidden taxa in columns = 0
for (i in 1:length(animal.list)){
  x <- animal.list[[i]]
  x[, forbid.col.list[[i]]] <- 0
  animal.list[[i]] <- x
}



# save object
saveRDS(animal.list, file = "animal only adj list.rds")
