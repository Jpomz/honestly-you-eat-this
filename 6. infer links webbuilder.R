# Taieri food webs WebBuilder
# JPomz
# June 2017

# make a list of registries using all Taieri predation matrices -web trying to infer

source("Useful WebBuilder functions.R")
source("C:\\Users\\Justin\\Documents\\Data\\FW modelling Petchey Github\\ttl-resources-master\\food_web\\FoodWebFunctions.r")
# TSS function
source("TSS function.R")

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
# make one complete registry to add to other registers later
complete.registry <- register.list[[1]]
colnames(complete.registry) <- c("source.id", "resource", "consumer")

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
# should throw ~32 warnings for coercing factor to character vector

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
links <- NULL
for (i in 1:length(taxa.list)){
  links[[i]] <- WebBuilder(taxa.list[[i]],
        register.list[[i]],
        method = c("exact", "genus", "family", "order"))
}
names(links) <- names(web.subset)

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
  for (j in seq_along(links[[i]]$resource)) #for each row in links list
  {
    wb.matrices[[i]][as.character(links[[i]]$resource[j]),
          as.character(links[[i]]$consumer[j])] <-  1
  }
}

# need to sort wb.matrices by order in obs matr
target <- NULL
for (i in 1:length(web.subset)){
  target[[i]] <- colnames(web.subset[[i]])
} # the %in% command subsets the target only to taxa that occur in both
# if you don't do this the following step throws an error (e.g. undefined columns selected)

# match order of matrices to target vector
# this orders matrices in order of increasing biomass
# observed adj
for (i in 1:length(wb.matrices)){
  wb.matrices[[i]] <- wb.matrices[[i]][
    match(target[[i]], rownames(wb.matrices[[i]])),
    match(target[[i]], colnames(wb.matrices[[i]]))
    ]}


# # Plots of predicted and observed food webs
# for (i in 1:length(web.subset)){
#   Plot.matrix(wb.matrices[[i]], pt.col = "tan", point.cex = 1.25)
#   par(new = T)
#   Plot.matrix(web.subset[[i]], point.cex = 1)
# }


# # TSS ####
# # this calculates true tss for observed v inferred
# # need to trim taxa so that it is comparable to gravel method
# tss.wb.df <- NULL
# for (i in 1:length(web.subset)){
#   j <- tss(web.subset[[i]], wb.matrices[[i]])
#   j$web <- names(web.subset)[[i]]
#   j$step <- "model"
#   tss.wb.df <- rbind(tss.wb.df, j)
#   rm(j)
# }

# tss.mean <- tss.wb.df %>%
#   summarize(ad = (mean(a) + mean (d)) /mean(S**2),
#             bc = mean(b +c) / mean(S**2))

# ad = 0.818346
# bc = 0.181654

# # plot TSS
# ggplot(tss.wb.df, aes(x = web, y = tss)) +
#   geom_point(shape = 0, size = 2) 
# # plot a, d, and b+c
# ggplot(tss.wb.df, aes(x = web, y = abar)) +
#   geom_point(aes(color = "a"), size = 2) +
#   geom_point(aes(x = web, y = dbar, color = "d"), shape = 9, size = 2)+
#   geom_point(aes(x = web, y = bc, color = "bc"), shape = 1, size = 2) +
#   scale_shape(breaks = c("abar", "dbar", "bc"))


# subset taxa ####
# in order to compare to the lnks inferred with Gravel method, need to trim matrices to only include taxa that occured in gravel matrices
# taxa were missing for a variaety of reasons, mostly that there weren't biomass estimates for all taxa

trimmed.taxa <- readRDS("trimmed taxa in gravel method.rds")

# subset inf and obs adj matrices by trimmed.taxa
web.trim <- NULL
for (i in 1:length(web.subset)){
  web.trim[[i]] <- web.subset[[i]][trimmed.taxa[[i]], trimmed.taxa[[i]]]}
names(web.trim) <- names(web.subset)

wb.matr.trim <- NULL
for (i in 1:length(wb.matrices)){
  wb.matr.trim[[i]] <- wb.matrices[[i]][trimmed.taxa[[i]], trimmed.taxa[[i]]]}
names(wb.matr.trim) <- names(wb.matrices)

# read in target list
# list made in script 5 compare obs v...
target <- readRDS("target vector to match all matrices.RDS")
# match order of matrices to target vector
# this orders matrices in order of increasing biomass
# matches order of gravel inferred matrices

# trimmed observed 
for (i in 1:length(web.trim)){
  web.trim[[i]] <- web.trim[[i]][
    match(target[[i]], rownames(web.trim[[i]])),
    match(target[[i]], colnames(web.trim[[i]]))
    ]}
# webbuilder inferred trimmed 
for (i in 1:length(wb.matr.trim)){
  wb.matr.trim[[i]] <- wb.matr.trim[[i]][
    match(target[[i]], rownames(wb.matr.trim[[i]])),
    match(target[[i]], colnames(wb.matr.trim[[i]]))
    ]}
saveRDS(wb.matr.trim, file ="webbuilder inferred list trimmed taxa.rds")
saveRDS(web.trim, file ="observed adj matr list trimmed taxa.rds")
# web.trim <- readRDS("observed adj matr list trimmed taxa.rds")

# TSS trimmed matrices ####
tss.trim.wb.df <- NULL
for (i in 1:length(web.trim)){
  j <- tss(web.trim[[i]], wb.matr.trim[[i]])
  j$web <- names(web.trim)[[i]]
  j$step <- "model"
  tss.trim.wb.df <- rbind(tss.trim.wb.df, j)
  rm(j)
}

saveRDS(tss.trim.wb.df, file = "webbuilder tss model inferred.rds")

# tss.trim.wb.df <- readRDS(file = "webbuilder tss model inferred.rds")

tss.mean <- tss.trim.wb.df %>%
  summarize(tss.mean = mean(tss),
            tss.sd = sd(tss),
            ad = (mean(a) + mean (d)) /mean(S**2),
            bc = mean(b +c) / mean(S**2))
# tss.mean = 0.66747
# tss.sd = 0.0897
# ad = 0.8141299
# bc = 0.1858701

# plot trim TSS
ggplot(tss.trim.wb.df, aes(x = web, y = tss)) +
  geom_point(shape = 0, size = 2) 
# plot a, d, and b+c
ggplot(tss.trim.wb.df, aes(x = web, y = abar)) +
  geom_point(aes(color = "a"), size = 2) +
  geom_point(aes(x = web, y = dbar, color = "d"), shape = 9, size = 2)+
  geom_point(aes(x = web, y = bc, color = "bc"), shape = 1, size = 2) +
  scale_shape(breaks = c("abar", "dbar", "bc"))



# for (i in 1:length(web.trim)){
#   Plot.matrix(wb.matr.trim[[i]], pt.col = "tan", point.cex = 1.25)
#   par(new = T)
#   Plot.matrix(web.trim[[i]], point.cex = 1)
# }

# filter links > 1 ####
# filter links by N.records
filtered.links <- links %>%
  llply(function (x){
    x %>% filter(N.records >1)
  })

# make empty matrices
wb.filt.matr <- taxa.list %>% llply(function (x){
  matrix(0, nrow(x), nrow(x))
})
#add dimnames to matrices
for (i in seq_along(wb.filt.matr)){
  dimnames(wb.filt.matr[[i]]) <- list(
    taxa.list[[i]]$node,
    taxa.list[[i]]$node)
}

# put 1's in at resource / consumer index
#this part is wrong...
for(i in seq_along(wb.filt.matr)) #for each matrix
{ 
  for (j in seq_along(filtered.links[[i]]$resource)) #for each row in filtered.links list
  {
    wb.filt.matr[[i]][as.character(filtered.links[[i]]$resource[j]),
                     as.character(filtered.links[[i]]$consumer[j])] <-  1
  }
}

# match order of matrices to target vector
# this orders matrices in order of increasing biomass
# observed adj
for (i in 1:length(wb.filt.matr)){
  wb.filt.matr[[i]] <- wb.filt.matr[[i]][
    match(target[[i]], rownames(wb.filt.matr[[i]])),
    match(target[[i]], colnames(wb.filt.matr[[i]]))
    ]}
# # subset inf and obs adj matrices by trimmed.taxa
# web.trim <- NULL
# for (i in 1:length(web.subset)){
#   web.trim[[i]] <- web.subset[[i]][trimmed.taxa[[i]], trimmed.taxa[[i]]]}
# names(web.trim) <- names(web.subset)

wb.filt.trim <- NULL
for (i in 1:length(wb.filt.matr)){
  wb.filt.trim[[i]] <- wb.filt.matr[[i]][trimmed.taxa[[i]], trimmed.taxa[[i]]]}
names(wb.filt.trim) <- names(wb.filt.matr)

# match order of matrices to target vector
# this orders matrices in order of increasing biomass
# matches order of gravel inferred matrices

# # trimmed observed 
# for (i in 1:length(web.trim)){
#   web.trim[[i]] <- web.trim[[i]][
#     match(target[[i]], rownames(web.trim[[i]])),
#     match(target[[i]], colnames(web.trim[[i]]))
#     ]}
# filtered webbuilder inferred trimmed 
for (i in 1:length(wb.filt.trim)){
  wb.filt.trim[[i]] <- wb.filt.trim[[i]][
    match(target[[i]], rownames(wb.filt.trim[[i]])),
    match(target[[i]], colnames(wb.filt.trim[[i]]))
    ]}
 saveRDS(wb.filt.trim, file ="webbuilder filtered links inferred matr.rds")
# wb.filt.trim <- readRDS("webbuilder filtered links inferred matr.rds")

# TSS filtered matrices ####
tss.filt.wb.df <- NULL
for (i in 1:length(web.trim)){
  j <- tss(web.trim[[i]], wb.filt.trim[[i]])
  j$web <- names(web.trim)[[i]]
  j$step <- "filter"
  tss.filt.wb.df <- rbind(tss.filt.wb.df, j)
  rm(j)
}

wb.tss.2step <- bind_rows(tss.filt.wb.df, tss.trim.wb.df)

wb.tss.2step %>% group_by(step) %>% summarise(mean(tss))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# calculate food web stats ####

# wb.matr.trim <- readRDS("webbuilder inferred list trimmed taxa.rds")
# web.subset <- readRDS("observed adj matr list trimmed taxa.rds")

saveRDS(wb.matr.trim, file ="webbuilder inferred list trimmed taxa.rds")
saveRDS(web.trim, file ="observed adj matr list trimmed taxa.rds")


# make web.list / web.trim into matrices
web.subset <- llply(web.subset,
                  function(x) {
                    as.matrix(x)})
# web.trim <- llply(web.trim,
#                   function(x) {
#                     as.matrix(x)})
wb.matr.trim <- llply(wb.matr.trim,
                  function(x) {
                    as.matrix(x)})


obs.stats <- ldply(web.subset, function (x){
  Get.web.stats(x)}
)
obs.stats$model <- "observed"
obs.stats$step <- NA

wb.trim.stats <- ldply(wb.matr.trim, function (x){
  Get.web.stats(x)
})
wb.trim.stats$model <- "webbuilder"
wb.trim.stats$step <- NA

obs.wb.fw.stats <- bind_rows(obs.stats, wb.trim.stats)
saveRDS(obs.wb.fw.stats, file = "food web stats observed and webbuilder.rds")


# this object isn't made until script 8. Need to fix the order of this
# gravel.fw.stats <- readRDS("food web stats gravel 3 steps.rds")
# 
# all.stats <- bind_rows(gravel.fw.stats, obs.wb.fw.stats)
# saveRDS(all.stats, file = "food web stats all.rds")

