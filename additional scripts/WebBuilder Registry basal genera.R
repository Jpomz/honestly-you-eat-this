# Making webbuilder registry
# JPOMZ 26 Nov 2017

# Making "one registry to rule them all!"
# Registry with complete species

source("C:\\Users\\Justin\\Documents\\Data\\FW modelling Petchey Github\\ttl-resources-master\\food_web\\FoodWebFunctions.r")

# recoder function 
# found this online somewhere, I think the reference information is in one of my oiginal webbuilder scripts
recoderFunc <- function(data, oldvalue, newvalue) {
  # convert any factors to characters
  if (is.factor(data))     data     <- as.character(data)
  if (is.factor(oldvalue)) oldvalue <- as.character(oldvalue)
  if (is.factor(newvalue)) newvalue <- as.character(newvalue)
  # create the return vector
  newvec <- data
  # put recoded values into the correct position in the return vector
  for (i in unique(oldvalue)) newvec[data == i] <- newvalue[oldvalue == i]
  newvec}

# read in original adj_matr
# set upper.dir
upper.dir <- "C:\\Users\\Justin\\Google Drive\\Data\\Predicting NZ Food Webs\\All Thompson 2004"
# full path names of adj_matr
full.name <- list.files(upper.dir, full.names = T)
# file names
a.names <- list.files(upper.dir)

# make list of adj_matr
adj.list <- NULL
for (i in 1:length(a.names)){
  j <- read.csv(full.name[i], header = T)
  adj.list[i] <- list(j)
  rm(j)
}

names(adj.list) <- c(a.names[1:length(a.names)])
# remove ".csv" from names
names(adj.list) <- gsub(".csv", "", names(adj.list))

##remove "X" column
adj.list <- llply(adj.list, 
                  function(x) { x["X"] <- NULL; x })
#summary(adj.list)

##add rownames to elements within list
adj.list <- llply(adj.list,function (x){
  row.names(x) <- colnames(x);x})

adj.list.names <- llply(adj.list,
                        function(x){
                          colnames(x)
                        })

# fix names, typos, misnomers, etc. 
translate <- read.csv("translation.csv") 
# translate list of names
adj.list.names.cor <- llply(adj.list.names,
                            function (x){
                              recoderFunc(x, translate$Wrong, translate$Corrected)})
# rename old webs with corrected names
for (i in 1:length(adj.list)){
  colnames(adj.list[[i]]) <- adj.list.names.cor[[i]]
}
## coerce elements to matrices to allow for duplicate rownames
adj.list <- llply(adj.list,
                  function(x) {
                    as.matrix(x)})
##add rownames to elements within list
adj.list <- llply(adj.list,
                  function (x){
                    row.names(x) <- colnames(x);x})
# genera ####
# species to genus translation
sp.gen <- read_csv("species genus category ffg.csv")


#make list of "old" names
old.names <- adj.list %>% 
  llply( function (x){colnames(x)}) 

#make list of generic names with recoderFunc
generic.names <- old.names %>% 
  llply(function (x) {
    recoderFunc(x, sp.gen$Species, sp.gen$Genus)
  })

#rename colnames in tt1 list of x tabs
for (i in 1:length(adj.list)){
  colnames(adj.list[[i]]) <- generic.names[[i]]
}
# make webs matrix (for duplicate row names)
adj.list <- llply(adj.list,
                  function (x){
                    as.matrix(x)
                  })

#make rownames match colnames
adj.list <- llply(adj.list,
                  function (x){
                    rownames(x) <- colnames(x);x
                  })




#need to re-compile with duplicated generic names
adj.list <- adj.list %>% 
  llply(function (x){
    as.data.frame(as.table(x)) %>%
      xtabs(Freq ~ Var1 + Var2, .)})

adj.list <- adj.list %>% 
  llply(function (x){x[x>1] <- 1; x})

adj.list <- adj.list %>% 
  llply(function (x){as.data.frame.matrix(x)})


# removing probable data errors
# eg, non-predatory taxa shown to be consuming prey
# potamopyrgus eating austrosim, deleatidum etc
# deleatidium consuming animal prey
# the mouthparts and habits of these animals make it extremely unlikely that they are consuming animals
errors <- c("Amphipoda", "Atalophlebioides", "Austroclima", "Austrosimulium","Green.algae", "Coloburiscus", "Deleatidium", "Nesameletus", "Oxyethira", "Potamopyrgus", "Zephlebia")

# make a list of columns which match forbidden taxa
forbid.col.list <- adj.list %>% 
  llply(function (x){
    which(colnames(x) %in% errors)
  })

# make forbidden taxa in columns = 0
for (i in 1:length(adj.list)){
  x <- adj.list[[i]]
  x[, forbid.col.list[[i]]] <- 0
  adj.list[[i]] <- x
}


# convert predation matrix to matrix of pairs
pairs.list <- adj.list %>%
  llply(function (x){
    as.data.frame(Matrix.to.list(x))})

register <- ldply(pairs.list)
colnames(register) <- c("source.id", "resource", "consumer")

# read in registry made from other sources
other.registry <- read_csv("fw_registry.csv")
other.registry$source.id <- as.character(other.registry$source.id)

# read in file with taxonomy information
taxonomy <- read_csv("C:\\Users\\Justin\\Google Drive\\Data\\Rose's Tui FW inference\\taxonomy.csv")
#taxonomy <- read_csv("C:\\Users\\Justin\\Documents\\R\\fish diatom and invert order family genus.csv")


# add taxonomy information for resource and consumers
register <- merge(register, taxonomy[,1:5],
                  by.x = "resource", by.y = "name",
                  all.x = T)
colnames(register)[4:6] <- paste("res.",
          colnames(register[c(4:6)]), sep = "")
register <- merge(register, taxonomy[,c(1:4,6)],
                  by.x = "consumer", by.y = "name",
                  all.x = T)
colnames(register)[8:10] <- paste("con.", 
          colnames(register[c(8:10)]), sep = "")

register <- merge(register, sp.gen[,2:3], by.x = "consumer", by.y = "Genus", all.x = T)
colnames(register)[12] <- "con.category"

register <- merge(register, sp.gen[,2:3], by.x = "resource", by.y = "Genus", all.x = T)
colnames(register)[13] <- "res.category"
register$linkevidence <- "Direct"

other.registry <- other.registry[intersect(
                              names(other.registry),
                  names(register))]
register <- bind_rows(register, other.registry)

write_csv(register, "WebBuilder registery basal genera.csv")
