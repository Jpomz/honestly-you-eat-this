# 16. webbuilder including basal sp
# J Pomz aug 2017




# This script takes the raw Taieri adjacency matrices (obtained from Ross Thompson November 2016 at New Zealand Ecological Society meeting, Hamilton NZ) and translates the names to the generic level (family / subfamily for coleoptera, Diptera, etc). 

# corrected adjacency matrices are saved as a list. 
# this will then be used to infer food web strucutre using the webbuilder method, using ALL taxa present. 
# these will be compared with observed networks (names resolved to genera) 



# food web functions from Petchey
source("C:\\Users\\Justin\\Documents\\Data\\FW modelling Petchey Github\\ttl-resources-master\\food_web\\FoodWebFunctions.r")
# webbuilder fxns
source("Useful WebBuilder functions.R")
# TSS function
source("TSS function.R")
# pairs to adj function
source("pairs_to_adj function.R")

library(vegan)
library(betalink)
library(igraph)

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

adj.list <- adj.list %>%
  llply( function (x){
    if("Green.algae" %in% colnames(x) ==T)
      x[x[,"Green.algae"]>0] <- 0; x
  })

# adj.list %>%
#   llply( function (x){
#     if("Green.algae" %in% colnames(x) ==T)
#       x[x[,"Green.algae"]>0]
#   })

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

saveRDS(adj.list, file = "list adj matr corr names basal genera.rds")
# adj.list <- readRDS("list adj matr corr names basal genera.rds")

# infer links ####

web.list <- llply(adj.list, function (x){
  as.matrix(x)
  }
)

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
for (i in 1:length(pairs.list)){
  register.list[[i]] <- ldply(pairs.list)
}

# add web names
names(register.list) <- names(web.list)
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
# should throw ~42 warnings for coercing factor to character vector

# get a list of taxa from web.list
# need to subset so only looking at webs which were also inferred using gravel method
taxa.list <- llply(web.list,
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
                           method = c("exact",
                                      "genus",
                                      "family",
                                      "order"))
}
names(links) <- names(web.list)
# saveRDS(links, "wb full community links.rds")
# links <- readRDS("wb full community links.rds")

# make matrices with new function
wb.matrices <- NULL
for (i in 1: length(links)){
  wb.matrices[[i]] <- pairs_to_adj(taxa.list[[i]],
                                 links[[i]])
}
names(wb.matrices) <- names(web.list)


target.df <- read.csv("C:/Users/Justin/Google Drive/Data/FW modelling Petchey Github/food web compilation/FW compilation/order.for.matrix.tl.csv")
target.df <- target.df$Genus %>% unique

target <- NULL
for (i in 1:length(wb.matrices)){
  target[[i]] <- target.df[
    target.df %in% colnames(wb.matrices[[i]])]
}

# mathc wb.matrices to target
for (i in 1:length(wb.matrices)){
  wb.matrices[[i]] <- wb.matrices[[i]][
    match(target[[i]], rownames(wb.matrices[[i]])),
    match(target[[i]], colnames(wb.matrices[[i]]))
    ]}

# # plot all wb + obs webs
# for (i in 1:length(wb.matrices)){
#   Plot.matrix(wb.matrices[[i]],
#               pt.col = "red", point.cex = 2)
#   par(new = T)
#   Plot.matrix(web.list[[i]], point.cex = 1)
# }

# match web.list to target
for (i in 1:length(web.list)){
  web.list[[i]] <- web.list[[i]][
    match(target[[i]], rownames(web.list[[i]])),
    match(target[[i]], colnames(web.list[[i]]))
    ]}

# save wb and observed matrices
saveRDS(wb.matrices, file = "wb full comm matrices.rds")
saveRDS(web.list, file = "obs matr full comm.rds")

# TSS ####
obs.df <- llply(web.list, function (x){
  as.data.frame(x)
  }
)

wb.df <- llply(wb.matrices, function (x){
  as.data.frame(x)
})

tss.wb.full.df <- NULL
for (i in 1:length(wb.df)){
  j <- tss(obs.df[[i]], wb.df[[i]])
  j$web <- names(web.list)[[i]]
  tss.wb.full.df <- rbind(tss.wb.full.df, j)
  rm(j)
}


tss.wb.full.df %>% summarize(tss = mean(tss), a = mean(abar), d = mean(dbar))
# mean TSS = 0.711
# mean a = 0.064
# d = 0.779

tss.wb.full.df %>%
  summarize(ad = (mean(a) + mean (d)) /mean(S**2),
            bc = mean(b +c) / mean(S**2))


# betalink ####
# obs.mat <- obs %>% llply(function (x){as.matrix(x)})
obs.graph <- web.list %>%
  llply( function (x){
    graph_from_adjacency_matrix(x)
  })
obs.beta <- beta_os_prime(obs.graph)
plot(density(obs.beta), main = "B'os observed")
abline(v = mean(obs.beta))

wb.full.graph <- wb.matrices %>%
  llply( function (x){
    graph_from_adjacency_matrix(x)
  })
wb.full.beta <- beta_os_prime(wb.full.graph)
plot(density(wb.full.beta), main = "B'os webbuilder")
abline(v = mean(wb.full.beta))

# food web stats ####
# webbuilder
wb.full.comm.stats <- wb.matrices %>% 
  ldply(function (x){
    Get.web.stats(x)
  }
)
# observed
obs.full.comm.stats <- web.list %>% 
  ldply(function (x){
    Get.web.stats(x)
  }
)

# make df with landuse cat
site.cat <- data.frame(.id = wb.full.comm.stats$.id)
site.cat$cat <- c("Pn", "Pn", "Pn", "P", "P", "P", "Pn", "T", NA, NA, "T", "T", "T", "T", "Pn", "N", "N", "T", "T", NA, "Pn")

# left_join landuse cat to food web stats
wb.full.comm.stats <- left_join(wb.full.comm.stats, site.cat, by = ".id")
# add "model" variable
wb.full.comm.stats$model <- "webbuilder"

# left_join landuse cat to food web stats
obs.full.comm.stats <- left_join(obs.full.comm.stats, site.cat, by = ".id")
# add "model" variable
obs.full.comm.stats$model <- "observed"

# remove NA landuse (seasonal webs from 1999)
wb.full.comm.stats <- wb.full.comm.stats %>%
  filter(!is.na(cat))
obs.full.comm.stats <- obs.full.comm.stats %>%
  filter(!is.na(cat))

# combine both stats into one df
both.stats <- rbind(wb.full.comm.stats, obs.full.comm.stats)
# saverds
saveRDS(both.stats, file = "both full comm fw stats.rds")

# add L/S
both.stats <- both.stats %>%
  mutate(LS = L / S)


# clean up for better ggplotting
# make the model variable as a factor
both.stats$model <- as.factor(both.stats$model)

# saverds
saveRDS(both.stats, file = "both full comm fw stats.rds")

# ANOVA ####
# one way anova testing different methods
# for differences among land use types
# using aov() in order to perform 
# posthoc Tukey test

# function for getting one way aov
one_way_aov <- function(d, x, y){
  aov(d[[y]] ~ d[[x]])
}


l <- dlply(both.stats,
           .(model),
           one_way_aov,
           x = "cat", y = "L")

ls <-  dlply(both.stats,
             .(model),
             one_way_aov,
             x = "cat", y = "LS")
c <-  dlply(both.stats,
            .(model),
            one_way_aov,
            x = "cat", y = "C")
mcl <-  dlply(both.stats,
              .(model),
              one_way_aov,
              x = "cat", y = "mean.TL")
b <-  dlply(both.stats,
            .(model),
            one_way_aov,
            x = "cat", y = "B")
int <-  dlply(both.stats,
              .(model),
              one_way_aov,
              x = "cat", y = "I")
top <-  dlply(both.stats,
              .(model),
              one_way_aov,
              x = "cat", y = "T")
g <-  dlply(both.stats,
            .(model),
            one_way_aov,
            x = "cat", y = "Gensd")
v <- dlply(both.stats,
           .(model),
           one_way_aov,
           x = "cat", y = "Vulsd")
msim <-  dlply(both.stats,
               .(model),
               one_way_aov,
               x = "cat", y = "Maxsim")
maxtl <-  dlply(both.stats,
                .(model),
                one_way_aov,
                x = "cat", y = "max.TL")

aov.results <- list(l, ls, c, mcl, b, int,
                    top, g, v, msim, maxtl)


# length(aov.res...) == 11
# length each aov.res element == 2
# variables == 11
# each variable has three methods
aov.table <- matrix(0, 11, 2)
for (i in 1:length(aov.results)){
  for (j in 1:length(aov.results[[1]])){
    aov.table[i,j] <- anova(aov.results[[i]][[j]])[1,5]
  }
}
dimnames(aov.table) <- list(
  c("L", "LS", "C", "meanTL", "B", "I", "T", "Gensd", "Vulsd", "Maxsim", "max.TL"),
  c("observed", "webbuilder"))
aov.table <- signif(aov.table, 5)
# significant results:
# Obs:L, LS, C, I, T, Vul, Maxsim
# wb: L, LS,       T, Vul, Maxsim

# Tukey L ####
# NS
TukeyHSD(l$observed)
# T - Pn, P > Pn (p=0.0577)
TukeyHSD(l$webbuilder)
# T > Pn, P > PN

# Tukey L/S ####
TukeyHSD(ls$observed)
# T > Pn, P > PN, T > N (p=0.0690)
TukeyHSD(ls$webbuilder)
# T > Pn, Pn > P, T > N, P > N (p=0.0727)

# Tukey C ####
TukeyHSD(c$observed)
# T > Pn, P > Pn, T > N
TukeyHSD(c$webbuilder)
# NA

# Tukey meanTL ####
# No significance
TukeyHSD(mcl$observed)
TukeyHSD(mcl$webbuilder)

# Tukey B ####
# NS
TukeyHSD(b$observed)
TukeyHSD(b$webbuilder)

# Tukey I ####
TukeyHSD(int$observed)
# T > Pn 
TukeyHSD(int$webbuilder)
# NA

# Tukey T ####
TukeyHSD(top$observed)
# Pn > T, Pn > P; N > T (p=0.0592), 
# webbuilder
TukeyHSD(top$webbuilder)
# Pn > T; Pn > P, N > T 

# Tukey Gensd ####
TukeyHSD(g$observed)
TukeyHSD(g$webbuilder)

# Tukey Vulsd ####
TukeyHSD(v$observed)
# Pn > T, Pn > P, N > T 
TukeyHSD(v$webbuilder)
# N > T, 

# Tukey maxsim ####
TukeyHSD(msim$observed)
# T > Pn, T > N (p=0.0658)
TukeyHSD(msim$webbuilder)
# T > Pn, T > N, P > N (p=0.0612) 

# Tukey maxtl ####
# NS
TukeyHSD(maxtl$observed)
TukeyHSD(maxtl$webbuilder)


# plot food web metrics ####
# functions for plots
plot_bar <- geom_bar(position = position_dodge(),
                     stat = "identity")
plot_error <- geom_errorbar(aes(ymin = avg - sd,
                                ymax = avg +sd),
                            position = position_dodge())
# links
both.stats %>%
  group_by(model, cat) %>%
  summarize(avg = mean(L), sd = sd(L)) %>%
  ggplot(aes(model, avg, fill = cat)) +
  plot_bar +
  plot_error

# L/S
both.stats %>%
  group_by(model, cat) %>%
  summarize(avg = mean(L/S), sd = sd(L/S)) %>%
  ggplot(aes(model, avg, fill = cat)) +
  plot_bar +
  plot_error

# Connectance
both.stats %>%
  group_by(model, cat) %>%
  summarize(avg = mean(C), sd = sd(C)) %>%
  ggplot(aes(model, avg, fill = cat)) +
  plot_bar +
  plot_error

# proportion Top species
both.stats %>%
  group_by(model, cat) %>%
  summarize(avg = mean(T), sd = sd(T)) %>%
  ggplot(aes(model, avg, fill = cat)) +
  plot_bar +
  plot_error

# mean.TL
both.stats %>%
  group_by(model, cat) %>%
  summarize(avg = mean(mean.TL), sd = sd(mean.TL)) %>%
  ggplot(aes(model, avg, fill = cat)) +
  plot_bar +
  plot_error


# ordinations ####

# make cat variable a factor
both.stats$cat <- as.factor(both.stats$cat)

# split stats into list
data.list <- dlply(both.stats, .(model))

# RDA ####
# make ordinations based on method
ord <- dlply(both.stats, .(model),
             function (x){
               rda(x[c(2:7,11:16,19)],
                   scale = T)
             })

# variables for ordination plots
scl <- 3 
colvec <- c("red2", "green4", "mediumblue", "black")

# observed ordination
plot(ord$observed, type = "n", scaling = scl)
with(data.list$observed,
     points(ord$observed,
            display = "sites",
            col = colvec[cat],
            scaling = scl,
            pch = 21,
            bg = colvec[cat]))
text(ord$observed, display = "species", scaling = scl, cex = 0.8, col = "darkcyan")
with(data.list$observed,
     legend("topright",
            legend = levels(cat),
            bty = "n",
            col = colvec,
            pch = 21,
            pt.bg = colvec))
title("Observed networks RDA")


# webbuilder ordination
plot(ord$webbuilder, type = "n", scaling = scl)
with(data.list$webbuilder,
     points(ord$webbuilder,
            display = "sites",
            col = colvec[cat],
            scaling = scl,
            pch = 21,
            bg = colvec[cat]))
text(ord$webbuilder, display = "species", scaling = scl, cex = 0.8, col = "darkcyan")
with(data.list$webbuilder,
     legend("topright",
            legend = levels(cat),
            bty = "n",
            col = colvec,
            pch = 21,
            pt.bg = colvec))
title("webbuilder networks RDA")



# MDS ####
# make data into distance matrix
mat.list <- llply(data.list, function (x){
  as.matrix(dist(x[c(2:7, 11:16, 19)]))
})
# make ordinations based on method
mds <- llply(mat.list,
             function (x){
               wcmdscale(x, k = 2,
                         eig = T,
                         w = rep(1, nrow(x)),
                         x.ret = TRUE)
             })
# observed ordination
plot(mds$observed, type = "n")
with(data.list$observed,
     points(mds$observed$points,
            col = colvec[cat],
            pch = 21,
            bg = colvec[cat]))
with(data.list$observed,
     legend("bottomright",
            legend = levels(cat),
            bty = "n",
            col = colvec,
            pch = 21,
            pt.bg = colvec))
title("Observed networks MDS")

# webbuilder ordination
plot(mds$webbuilder, type = "n")
with(data.list$webbuilder,
     points(mds$webbuilder$points,
            col = colvec[cat],
            pch = 21,
            bg = colvec[cat]))
with(data.list$webbuilder,
     legend("bottomleft",
            legend = levels(cat),
            bty = "n",
            col = colvec,
            pch = 21,
            pt.bg = colvec))
title("webbuilder networks MDS")
