# change method to family for all bugs

# from script 6
# Taieri food webs WebBuilder
# JPomz
# June 2017

# make a list of registries using all Taieri predation matrices -web trying to infer

source("Useful WebBuilder functions.R")
source("C:\\Users\\Justin\\Documents\\Data\\FW modelling Petchey Github\\ttl-resources-master\\food_web\\FoodWebFunctions.r")
# TSS function
source("TSS function.R")

# GENUS ####
# read in corrected adjaceny matrices
web.list <- readRDS("animal only adj list.rds")
# read in inferred list to subset by names
# only using names of these webs
inf.list <- readRDS("inferred adj matr gravel method.rds")

#subset web.list to only those also in inf.list
web.subset <- web.list[names(web.list) %in% names(inf.list)]

# convert predation matrix to matrix of pairs
pairs.list <- web.list %>%
  llply(function (x){
    as.data.frame(Matrix.to.list(x))})



# make list of res-con pairs for registries
# each element in list contains all pairs EXCEPT those of
# the web they will be used to parameterize. 
# e.g. element named "blackrock" will contain all other
# pairs data [-"blackrock"], and be used to paramterize # the model in order to predict the blackrock web 
register.list <- NULL
for (i in 1:length(inf.list)){
  register.list[[i]] <- ldply(pairs.list)
}

# add web names
names(register.list) <- names(inf.list)
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
taxonomy <- read_csv("C:\\Users\\Justin\\Documents\\R\\fish diatom and invert order family genus.csv")

# res.meth <- taxonomy$res.method
# con.meth <- taxonomy$con.method
# res.meth[res.meth == "genus"] <- "family"
# con.meth[con.meth == "genus"] <- "family"
# taxonomy$res.method <- res.meth
# taxonomy$con.method <- con.meth

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
taxa.list <- llply(web.subset,
                   function (x){
                     data.frame(node = colnames(x))
                   })

# add taxonomy information
taxa.list <- taxa.list %>% 
  llply(function (x){
    merge(x, taxonomy,by.x="node", 
          by.y="name", all.x=T)
  })

# add "minimum" to method columns
taxa.list <- taxa.list  %>%
  llply(function (x) {
    colnames(x)[5:6] <- paste("minimum.", colnames(x[c(5:6)]), sep = "");x
  })

# make node as.character
# thought this was necessary, but maybe not...
taxa.list <- taxa.list  %>%
  llply(function (x) {
    x$node <- as.character(x$node);x
  })

# infer links with WebBuilder function
# takes a while to run...
links.genus <- NULL
for (i in 1:length(taxa.list)){
  links.genus[[i]] <- WebBuilder(taxa.list[[i]],
                               register.list[[i]],
                               method = c("genus","family","order"))
}
names(links.genus) <- names(web.subset)

# look at N.records
# ldply(links, function (x){x %>% filter(N.records >10) %>% select(resource, consumer, N.records)})


# make empty matrices
wb.matrices <- taxa.list %>% llply(function (x){
  matrix(0, nrow(x), nrow(x))
})
#add dimnames to matrices
for (i in seq_along(wb.matrices)){
  dimnames(wb.matrices[[i]]) <- list(
    taxa.list[[i]]$node,
    taxa.list[[i]]$node)
}

# put 1's in at resource / consumer index
#this part is wrong...
for(i in seq_along(wb.matrices)) #for each matrix
{ 
  for (j in seq_along(links.genus[[i]]$resource)) #for each row in links.genus list
  {
    wb.matrices[[i]][as.character(links.genus[[i]]$resource[j]),
                     as.character(links.genus[[i]]$consumer[j])] <-  1
  }
}

wb.trim.list <- readRDS("webbuilder inferred list trimmed taxa.rds")

target <- readRDS("target vector to match all matrices.rds")

# match order of matrices to target vector
# this orders matrices in order of increasing biomass
# observed adj
for (i in 1:length(wb.matrices)){
  wb.matrices[[i]] <- wb.matrices[[i]][
    match(target[[i]], rownames(wb.matrices[[i]])),
    match(target[[i]], colnames(wb.matrices[[i]]))
    ]}

obs <- readRDS("observed adj matr list trimmed taxa.rds")

# plot of dempsters
Plot.matrix(wb.matrices[[7]])
title("Dempsters, Genus", cex.main = 1.5)
# plot of little
Plot.matrix(wb.matrices[[11]])
title("Little Kye Burn, Genus", cex.main = 1.5)


# # Plots of predicted and observed food webs
# for (i in 1:length(obs)){
#   Plot.matrix(wb.matrices[[i]], pt.col = "tan", point.cex = 1.25)
#   par(new = T)
#   Plot.matrix(obs[[i]], point.cex = 1)
# }

tss.genus <- NULL
for(i in 1:length(obs)){
  temp <- tss(obs[[i]], wb.matrices[[i]])
  tss.genus <- rbind(tss.genus, temp)
}
tss.genus




# FAMILY ####
# read in corrected adjaceny matrices
web.list <- readRDS("animal only adj list.rds")
# read in inferred list to subset by names
# only using names of these webs
inf.list <- readRDS("inferred adj matr gravel method.rds")

#subset web.list to only those also in inf.list
web.subset <- web.list[names(web.list) %in% names(inf.list)]

# convert predation matrix to matrix of pairs
pairs.list <- web.list %>%
  llply(function (x){
    as.data.frame(Matrix.to.list(x))})



# make list of res-con pairs for registries
# each element in list contains all pairs EXCEPT those of
# the web they will be used to parameterize. 
# e.g. element named "blackrock" will contain all other
# pairs data [-"blackrock"], and be used to paramterize # the model in order to predict the blackrock web 
register.list <- NULL
for (i in 1:length(inf.list)){
  register.list[[i]] <- ldply(pairs.list)
}

# add web names
names(register.list) <- names(inf.list)
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
taxonomy <- read_csv("C:\\Users\\Justin\\Documents\\R\\fish diatom and invert order family genus.csv")

res.meth <- taxonomy$res.method
con.meth <- taxonomy$con.method
res.meth[res.meth == "genus"] <- "family"
con.meth[con.meth == "genus"] <- "family"
taxonomy$res.method <- res.meth
taxonomy$con.method <- con.meth

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
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# # make table of registry information for manuscript / SI
# # want site, n pairwise intxn 
# # na pairwise intxns
# n.pairs <- register.list %>%
#   ldply(function (x){
#     data.frame(n = nrow(x))
#   })
# 
# # write csv of registry size
# # modify in excel to make table for publication
# write_csv(n.pairs, "C:/Users/Justin/Google Drive/Data/Predicting NZ Food Webs/figs for MS/webbuilder registry info from R.csv")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# get a list of taxa from web.list
# need to subset so only looking at webs which were also inferred using gravel method
taxa.list <- llply(web.subset,
                   function (x){
                     data.frame(node = colnames(x))
                   })

# add taxonomy information
taxa.list <- taxa.list %>% 
  llply(function (x){
    merge(x, taxonomy,by.x="node", 
          by.y="name", all.x=T)
  })

# add "minimum" to method columns
taxa.list <- taxa.list  %>%
  llply(function (x) {
    colnames(x)[5:6] <- paste("minimum.", colnames(x[c(5:6)]), sep = "");x
  })

# make node as.character
# thought this was necessary, but maybe not...
taxa.list <- taxa.list  %>%
  llply(function (x) {
    x$node <- as.character(x$node);x
  })

# infer links with WebBuilder function
# takes a while to run...
links.fam <- NULL
for (i in 1:length(taxa.list)){
  links.fam[[i]] <- WebBuilder(taxa.list[[i]],
                           register.list[[i]],
                           method = c("family","order"))
}
names(links.fam) <- names(web.subset)

# look at N.records
# ldply(links, function (x){x %>% filter(N.records >10) %>% select(resource, consumer, N.records)})


# make empty matrices
wb.matrices <- taxa.list %>% llply(function (x){
  matrix(0, nrow(x), nrow(x))
})
#add dimnames to matrices
for (i in seq_along(wb.matrices)){
  dimnames(wb.matrices[[i]]) <- list(
    taxa.list[[i]]$node,
    taxa.list[[i]]$node)
}

# put 1's in at resource / consumer index
#this part is wrong...
for(i in seq_along(wb.matrices)) #for each matrix
{ 
  for (j in seq_along(links.fam[[i]]$resource)) #for each row in links.fam list
  {
    wb.matrices[[i]][as.character(links.fam[[i]]$resource[j]),
                     as.character(links.fam[[i]]$consumer[j])] <-  1
  }
}

wb.trim.list <- readRDS("webbuilder inferred list trimmed taxa.rds")

target <- readRDS("target vector to match all matrices.rds")

# match order of matrices to target vector
# this orders matrices in order of increasing biomass
# observed adj
for (i in 1:length(wb.matrices)){
  wb.matrices[[i]] <- wb.matrices[[i]][
    match(target[[i]], rownames(wb.matrices[[i]])),
    match(target[[i]], colnames(wb.matrices[[i]]))
    ]}

obs <- readRDS("observed adj matr list trimmed taxa.rds")

# plot of dempsters
Plot.matrix(wb.matrices[[7]])
title("Dempsters, Family", cex.main = 1.5)
# plot of little
Plot.matrix(wb.matrices[[11]])
title("Little Kye Burn, Family", cex.main = 1.5)


# # Plots of predicted and observed food webs
# for (i in 1:length(obs)){
#   Plot.matrix(wb.matrices[[i]], pt.col = "tan", point.cex = 1.25)
#   par(new = T)
#   Plot.matrix(obs[[i]], point.cex = 1)
# }

tss.fam <- NULL
for(i in 1:length(obs)){
  temp <- tss(obs[[i]], wb.matrices[[i]])
  tss.fam <- rbind(tss.fam, temp)
}
tss.fam





# ORDER ####
register.list <- NULL
for (i in 1:length(inf.list)){
  register.list[[i]] <- ldply(pairs.list)
}

# add web names
names(register.list) <- names(inf.list)
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
taxonomy <- read_csv("C:\\Users\\Justin\\Documents\\R\\fish diatom and invert order family genus.csv")

res.meth <- taxonomy$res.method
con.meth <- taxonomy$con.method
res.meth[res.meth == "genus"| res.meth =="family"] <- "order"
con.meth[con.meth == "genus"| con.meth =="family"] <- "order"
taxonomy$res.method <- res.meth
taxonomy$con.method <- con.meth

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

# get a list of taxa from web.list
# need to subset so only looking at webs which were also inferred using gravel method
taxa.list <- llply(web.subset,
                   function (x){
                     data.frame(node = colnames(x))
                   })

# add taxonomy information
taxa.list <- taxa.list %>% 
  llply(function (x){
    merge(x, taxonomy,by.x="node", 
          by.y="name", all.x=T)
  })

# add "minimum" to method columns
taxa.list <- taxa.list  %>%
  llply(function (x) {
    colnames(x)[5:6] <- paste("minimum.", colnames(x[c(5:6)]), sep = "");x
  })

# make node as.character
# thought this was necessary, but maybe not...
taxa.list <- taxa.list  %>%
  llply(function (x) {
    x$node <- as.character(x$node);x
  })

links.ord <- NULL
for (i in 1:length(taxa.list)){
  links.ord[[i]] <- WebBuilder(taxa.list[[i]],
                               register.list[[i]],
                               method = c("order"))
}
names(links.ord) <- names(web.subset)

# look at N.records
# ldply(links, function (x){x %>% filter(N.records >10) %>% select(resource, consumer, N.records)})


# make empty matrices
wb.matrices <- taxa.list %>% llply(function (x){
  matrix(0, nrow(x), nrow(x))
})
#add dimnames to matrices
for (i in seq_along(wb.matrices)){
  dimnames(wb.matrices[[i]]) <- list(
    taxa.list[[i]]$node,
    taxa.list[[i]]$node)
}

# put 1's in at resource / consumer index
#this part is wrong...
for(i in seq_along(wb.matrices)) #for each matrix
{ 
  for (j in seq_along(links.ord[[i]]$resource)) #for each row in links.ord list
  {
    wb.matrices[[i]][as.character(links.ord[[i]]$resource[j]),
                     as.character(links.ord[[i]]$consumer[j])] <-  1
  }
}

wb.trim.list <- readRDS("webbuilder inferred list trimmed taxa.rds")

target <- readRDS("target vector to match all matrices.rds")

# match order of matrices to target vector
# this orders matrices in order of increasing biomass
# observed adj
for (i in 1:length(wb.matrices)){
  wb.matrices[[i]] <- wb.matrices[[i]][
    match(target[[i]], rownames(wb.matrices[[i]])),
    match(target[[i]], colnames(wb.matrices[[i]]))
    ]}

obs <- readRDS("observed adj matr list trimmed taxa.rds")

# plot of dempsters
Plot.matrix(wb.matrices[[7]])
title("Dempsters, Order", cex.main = 1.5)
# plot of little
Plot.matrix(wb.matrices[[11]])
title("Little Kye Burn, Order", cex.main = 1.5)


# Plots of predicted and observed food webs
for (i in 1:length(obs)){
  Plot.matrix(wb.matrices[[i]], pt.col = "tan", point.cex = 1.25)
  par(new = T)
  Plot.matrix(obs[[i]], point.cex = 1)
}

tss.ord <- NULL
for(i in 1:length(obs)){
  temp <- tss(obs[[i]], wb.matrices[[i]])
  tss.ord <- rbind(tss.ord, temp)
}
tss.ord

tss.exact <- readRDS("webbuilder tss model inferred.rds")

tss.exact %>% summarize(mean(tss))
tss.genus%>% summarize(mean(tss))
tss.fam %>% summarize(mean(tss))
tss.ord %>% summarize(mean(tss))

# for outline power point
# tss for dempters 
tss.genus[7,6]
tss.fam[7,6]
tss.ord[7,6]
# tss for little
tss.genus[11,6]
tss.fam[11,6]
tss.ord[11,6]
