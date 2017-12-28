# fixing Taieri adjacency matrices
# Jun 2017
# JPomz

# This script takes the raw Taieri adjacency matrices (obtained from Ross Thompson November 2016 at New Zealand Ecological Society meeting, Hamilton NZ) and translates the names to the generic level (family / subfamily for coleoptera, Diptera, etc). It converts basal resources to categorical names (e.g. diatoms, algae, etc)

# corrected adjacency matrices are saved as a list. These corrected matrices are then used to parameterize the gravel model, as well as make individual registries for the WebBuilder function (2 separate scripts)



# food web functions from Petchey
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



#rename colnames in tt1 list of x tabs
for (i in 1:length(adj.list)){
  colnames(adj.list[[i]]) <- new.names[[i]]
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

saveRDS(adj.list, file = "list adj matr corr names basal cat.rds")
