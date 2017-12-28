# comparing inferred to observed
# JPomz
# june 2017

# comparing adjacency matrices inferred using Gravel method (Parameterizing Gravel ...R script ) to animal only matrices observed in taieri food webs (truncating taieri...R script)

# useful food web functions from Petchey
source("C:\\Users\\Justin\\Documents\\Data\\FW modelling Petchey Github\\ttl-resources-master\\food_web\\FoodWebFunctions.r")
# TSS function
source("TSS function.R")

# read in animal only links from observed adjacency matrices
obs.list <- readRDS("animal only adj list.rds")
 
# read in Gravel inferred link matrices
inf.list <- readRDS("inferred adj matr gravel method.rds")
# subset web.list to match only names in taieri.list
obs.list <- obs.list[names(obs.list) %in% names(inf.list)]

# subsetting / ordering adj ####

# some names in taieri adjacency matrices don't match names in dataset with dryweight estimates (e.g. Hydora in adj, elmidae in dw data)
# make separate lists of vectors with taxa names for observed and inferred matrices 
  # *** inf.taxa, obs.taxa
# make a list that contains taxa names which occur in both 
  # *** taxa.intersect <- intersect(inf.taxa, obs.taxa)
# subset both inf and obs adj matr to taxa that intersect
# make target vector of taxa names to match adj matrices to
  # first subset 

# make list of names in inferred matrices 
inf.taxa <- llply(inf.list, 
                        function (x){
                          colnames(x)
                        })
# make list of names in observed matrices
obs.taxa <- llply(obs.list, 
                  function (x){
                    colnames(x)
                  })
# make list of names that occur in both inf and obs matr
taxa.intersect <- NULL
for (i in 1:length(obs.taxa)){
  taxa.intersect[[i]] <- intersect(obs.taxa[[i]], inf.taxa[[i]])
}
# save this object for use in webbuilder analysis
saveRDS(taxa.intersect, "trimmed taxa in gravel method.rds")

# subset inf and obs adj matrices by taxa.intersect
for (i in 1:length(obs.list)){
  obs.list[[i]] <- obs.list[[i]][taxa.intersect[[i]], taxa.intersect[[i]]]}

for (i in 1:length(inf.list)){
  inf.list[[i]] <- inf.list[[i]][taxa.intersect[[i]], taxa.intersect[[i]]]}

# make target vector to match order
# this will allow obs and inf adj to both be sorted by increasing biomass
target <- NULL
for (i in 1:length(inf.taxa)){
  target[[i]] <- inf.taxa[[i]][inf.taxa[[i]] %in% taxa.intersect[[i]]]
} 

saveRDS(target, "target vector to match all matrices.RDS")

# the %in% command subsets the target only to taxa that occur in both
# if you don't do this the following step throws an error (e.g. undefined columns selected)

# match order of matrices to target vector
# this orders matrices in order of increasing biomass
# observed adj
for (i in 1:length(obs.list)){
  obs.list[[i]] <- obs.list[[i]][
    match(target[[i]], rownames(obs.list[[i]])),
    match(target[[i]], colnames(obs.list[[i]]))
    ]}
#saveRDS(obs.list, "observed webs names and taxa matched.rds")
# inf
for (i in 1:length(inf.list)){
  inf.list[[i]] <- inf.list[[i]][
    match(target[[i]], rownames(inf.list[[i]])),
    match(target[[i]], colnames(inf.list[[i]]))
    ]}

saveRDS(inf.list, file = "gravel inferred initial.rds")
# # Plots of predicted and observed food webs
# for (i in 1:length(inf.list)){
#   Plot.matrix(inf.list[[i]], pt.col = "tan", point.cex = 1.25)
#   par(new = T)
#   Plot.matrix(obs.list[[i]], point.cex = 1)
# }

# plot inferred v obs web 


# #TSS ####
# tss() function written June 19 2017

# TSS = (ad - bc)/[(a+c)(b+d)]
# a = TPR = number of links predicted and observed = matrix value == 2
# b = predicted not observed; matrix value == -1
# c = observed but not predicted; matrix value == 1
# d = not predicted, not observed; matrix value ==0


tss.df <- NULL
for (i in 1:length(obs.list)){
  j <- tss(obs.list[[i]], inf.list[[i]])
  j$web <- names(obs.list)[[i]]
  j$step <- "model"
  tss.df <- rbind(tss.df, j)
  rm(j)
}

# tss.df <- tss.df %>%
#   mutate(mean.tss = mean(tss),
#          mean.abar = mean(abar),
#          mean.bbar = mean(b),
#          mean.cbar = mean(cbar),
#          mean.dbar = mean(dbar))

# tss.mean <- tss.df %>%
#   summarize(ad = (mean(a) + mean (d)) /mean(S**2),
#             bc = mean(b +c) / mean(S**2))
# ad = 0.4375731
# bc = 0.5624269

# plot TSS
ggplot(tss.df, aes(x = web, y = tss)) +
  geom_point(shape = 0, size = 2) 
# plot a, d, and b+c
ggplot(tss.df, aes(x = web, y = abar)) +
  geom_point(aes(color = "a"), size = 2) +
  geom_point(aes(x = web, y = dbar, color = "d"), shape = 9, size = 2)+
  geom_point(aes(x = web, y = bc, color = "bc"), shape = 1, size = 2)



# pruning forbidden links ####
# get list of all taxa in inferred webs
inf.taxa.unique <- inf.taxa %>%
  unlist %>%
  unique() %>%
  sort
inf.taxa.df <- inf.taxa.unique %>% as.data.frame()
colnames(inf.taxa.df) <- "taxa"

# read in data with FFG information
ffg <- read_csv("taxa.ffg.csv")

# join inferred taxa with FFG
inf.taxa.df <- left_join(inf.taxa.df, ffg, c("taxa" = "name"))

# make a list of "forbidden taxa"
# e.g. taxa that cannot eat prey due to mouthparts,
# e.g. scrapers or filter feeders
forbid.taxa <- inf.taxa.df %>%
  filter(FFG == "SC" | FFG == "G" | FFG == "FF") %>% select(taxa)
# write_csv(forbid.taxa, path = paste(getwd(), "/figs for MS/pruned taxa.csv", sep =""))
# convert to character vector
forbid.taxa <- forbid.taxa %>% as_vector()

# make new list of webs to prune
inf.forbid.list <- inf.list
# make a list of columns which match forbidden taxa
forbid.col.list <- inf.forbid.list %>% 
  llply(function (x){
    which(colnames(x) %in% forbid.taxa)
  })

# make forbidden taxa in columns = 0
for (i in 1:length(inf.forbid.list)){
  x <- inf.forbid.list[[i]]
  x[, forbid.col.list[[i]]] <- 0
  inf.forbid.list[[i]] <- x
}
saveRDS(inf.forbid.list, "gravel inferred pruned.rds")

# # Plots of pruned and observed food webs
# for (i in 1:length(inf.forbid.list)){
#   Plot.matrix(inf.forbid.list[[i]], pt.col = "tan", point.cex = 1.25)
#   par(new = T)
#   Plot.matrix(obs.list[[i]], point.cex = 1)
# }

#TSS ####
tss.prune.df <- NULL
for (i in 1:length(obs.list)){
  j <- tss(obs.list[[i]], inf.forbid.list[[i]])
  j$web <- names(obs.list)[[i]]
  j$step <- "prune"
  tss.prune.df <- rbind(tss.prune.df, j)
  rm(j)
}

gravel.tss.2step <- bind_rows(tss.df, tss.prune.df)
saveRDS(gravel.tss.2step, paste(getwd(), "/figs for MS/gravel tss 2 steps.rds", sep =""))


tss.prune.mean <- tss.prune.df %>%
  summarize(ad = (mean(a) + mean (d)) /mean(S**2),
            bc = mean(b +c) / mean(S**2))


# plot TSS
ggplot(tss.prune.df, aes(x = web, y = tss)) +
  geom_point(shape = 0, size = 2) 
# plot a, d, and b+c
ggplot(tss.prune.df, aes(x = web, y = abar)) +
  geom_point(aes(color = "a"), size = 2) +
  geom_point(aes(x = web, y = dbar, color = "d"), shape = 9, size = 2)+
  geom_point(aes(x = web, y = bc, color = "bc"), shape = 1, size = 2)

