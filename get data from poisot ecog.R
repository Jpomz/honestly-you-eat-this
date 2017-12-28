#! /usr/bin/env Rscript

library(rglobi)
library(taxize)
library(spocc)
library(plyr)
library(igraph)
library(RColorBrewer)

genus_name = function(n) {
  n = strsplit(n, " ")
  n = laply(n, function(x) x[1])
  return(n)
}

interactions_from_globi = function(n) {
   globi_prey = c()
   globi_pred = c()
   #cat(paste("Getting interactions for", n,"\n"))
   gl_int = get_interactions_by_taxa(sourcetaxon=n)
   gl_int = subset(gl_int, interaction_type %in% c("preysOn", "preyedUponBy"))
   if(NROW(gl_int) > 0){
      for(i in c(1:nrow(gl_int))){
         r = gl_int[i,]
         if(r$interaction_type == "preysOn"){
            globi_pred = c(globi_pred, r$source_taxon_name)
            globi_prey = c(globi_prey, r$target_taxon_name)
         } else {
            globi_prey = c(globi_prey, r$source_taxon_name)
            globi_pred = c(globi_pred, r$target_taxon_name)
         }
      }
      globi_interactions = unique(
      data.frame(i = genus_name(globi_pred), j = genus_name(globi_prey))
      )
      globi_interactions$i = as.vector(globi_interactions$i)
      globi_interactions$j = as.vector(globi_interactions$j)
      return(globi_interactions)
   }
   return(NA)
}

multi_int_from_taxa = function(genus) {
   globi_int = data.frame()
   for(g in genus) {
      g_int = interactions_from_globi(g)
      if(!is.na(g_int)) globi_int = rbind(globi_int, g_int)
   }
   globi_int$i = as.vector(globi_int$i)
   globi_int$j = as.vector(globi_int$j)
   return(globi_int)
}

cat("COLLATING ORIGINAL DATA\n")
cat(date(), "\n")

# Get original data
prey = c()
pred = c()
for(f in list.files(path = "All Thompson 2004", full.names=T))
  #path="webs", ## original argument in poisot script
{
   A = read.csv(f, h=T, sep=',', row.names=1)
   A = A[,-ncol(A)]
   u_sp = rownames(A)
   genus = genus_name(u_sp)
   for(i in c(1:nrow(A))){
      for(j in c(1:ncol(A))){
         if(A[i,j]==1){
            pred = c(pred, genus[i])
            prey = c(prey, genus[j])
         }
      }
   }
}

interactions = unique(data.frame(i = pred, j = prey))

INVALID_NAMES = c("Unidentified", "Unidentifiable", "Unknown", "cf.", "Terrestrial")

for(INV_NAME in INVALID_NAMES){
   interactions = subset(interactions, (i != INV_NAME)&(j != INV_NAME))
}

interactions$i = as.vector(interactions$i)
interactions$j = as.vector(interactions$j)

genus = sort(unique(c(interactions$i, interactions$j)))

cat("CURRENT STEP: ")
cat(length(genus))
cat(" - ")
cat(nrow(interactions))
cat("\n")

# Match genus names using gnr_resolve
# No matched names = NA
# If not, returns best match
cat("MATCHING AND CLEANING GENERA NAMES\n")
cat(date(), "\n")
snames = c()
cnames = c()
score = c()
for(g in genus) {
   #cat("Matching genus ")
   #cat(g)
   #cat("\n")
   nr = gnr_resolve(g)
   # if (nrow(nr$results) > 1) { # original
   if (nrow(nr) > 1) {
      snames = c(snames, g)
      # original
      # cnames = c(cnames, genus_name(nr$results$matched_name[1]))
      # score = c(score, nr$results$score[1])
      # modified
      cnames = c(cnames, genus_name(nr$matched_name[1]))
      score = c(score, nr$score[1])
   } else {
      snames = c(snames, g)
      cnames = c(cnames, "NONE")
      score = c(score, 0)
   }
}

names(cnames) = snames
interactions$i = cnames[interactions$i]
interactions$j = cnames[interactions$j]
interactions = subset(interactions, (i!='NONE')&(j!='NONE'))
interactions = unique(interactions)

genus = sort(unique(c(interactions$i, interactions$j)))

cat("CURRENT STEP: ")
cat(length(genus))
cat(" - ")
cat(nrow(interactions))
cat("\n")

# Add data from Globi
cat("GETTING INTERACTION DATA\n")
cat(date(), "\n")
globi_int = multi_int_from_taxa(genus)
all_int = rbind(interactions, globi_int)
all_gen = unique(c(all_int$i, all_int$j))
globi_gen = sort(unique(c(globi_int$i, globi_int$j)))
new_globi_gen = globi_gen[!(globi_gen %in% genus)]
new_globi_int = multi_int_from_taxa(new_globi_gen)
new_globi_int = subset(new_globi_int, (i %in% all_gen) & (j %in% all_gen))
interactions = unique(rbind(all_int, new_globi_int))

# Get interactions in igraph
g = graph.data.frame(interactions)

cat("CURRENT STEP: ")
cat(length(V(g)))
cat(" - ")
cat(length(E(g)))
cat("\n")

# Check taxonomy
cat("GET TAXONOMY DATA\n")
cat(date(), "\n")
invalid_taxa = c()
nv = length(V(g)$name)
for(i in c(1:nv)) {
   tax = V(g)$name[i]
   taxid = NA
   if(is.na(taxid)) taxid = get_gbifid(tax, rows=1)
   if(is.na(taxid)) {
      invalid_taxa = c(invalid_taxa, i)
   } else {
      cl = classification(taxid, db='gbif')[[1]]
      tdepth = nrow(cl)
      if(cl$rank[tdepth] != 'genus') invalid_taxa = c(invalid_taxa, i)
   }
}

g = delete.vertices(g, invalid_taxa)
g = delete.vertices(g, which(degree(g)==0))

cat("CURRENT STEP: ")
cat(length(V(g)))
cat(" - ")
cat(length(E(g)))
cat("\n")

# Get occurence data from gbif using spocc
sources = c('gbif', 'bison')
dat = occ(V(g)$name, from=sources, limit=2000, gbifopts=list(hasCoordinate=TRUE))

# Assemble occurences in a data frame
ppd = data.frame()
for(source in dat){
   for(sp in names(source$data)){
      X = source$data[[sp]]
      if(nrow(X)>0){
         if('latitude' %in% colnames(X)) ppd = rbind(ppd, X[,c('name', 'latitude', 'longitude')])
      }
   }
}

ppd$genus = genus_name(ppd$name)
ppd = unique(ppd[,c('genus', 'latitude', 'longitude')])

ppd = subset(ppd, genus %in% V(g)$name)

save(ppd, file="ppd.Rdata")
save(g, file="g.Rdata")
