# all fish eat all prey all the time??
library(plyr)
library(tidyverse)
inf.niche <- readRDS("Niche pruned trait matching inference.RDS")
inf.neutral <- readRDS("Neutral trait matching inference.RDS")
inf.niche.neutral <- readRDS("Neutral + Niche trait matching inference.RDS")
obs <- readRDS("observed matrices matched to inferred.RDS")
# wb.raw == initial WB inference
wb.raw <- readRDS("wb matrices matched to inferred.RDS")
# wb == neutral prune + fish correction
wb.n <- readRDS("wb for PCA.RDS")
tm <- readRDS("Initial trait matching inference.RDS")
tm.niche <- readRDS("Niche pruned trait matching inference.RDS")
# tm == neutral, niche prune, fish correction
tm.nn <- readRDS("tm nn fish for PCA.RDS")
# wb x trait match, no neutral
wb.tm <- readRDS("wb x tm for PCA.RDS")
# wb.tm wb * tm neutral prune fish correction
wb.tm.n <- readRDS("wb x tm x n for PCA.rds")

f.taxa = c("Salmo", "Galaxias", "Gobiomorphus", "Anguilla")
fish.prey <- function(web){
  prop = c()
  for(f in which(colnames(web) %in% f.taxa)){
    i <- sum(web[,f]) / length(web[,f])
    prop = c(prop, i)
  }
  prop
}


sum(web.match[[17]][,25])/length(web.match[[17]][,25])


which(colnames(web.match[[17]]) %in% f.taxa)
fish.prey(web.match[[10]])

mean(flatten_dbl(map(obs, fish.prey)))
mean(flatten_dbl(map(wb.raw, fish.prey)))
mean(flatten_dbl(map(tm.nn, fish.prey)))
mean(flatten_dbl(map(tm.niche, fish.prey)))
mean(flatten_dbl(map(wb.tm, fish.prey)))
mean(flatten_dbl(map(wb.tm.n, fish.prey)))

mean(flatten_dbl(map(wb.n, fish.prey)))


data.frame(obs = map(obs, fish.prey) %>% sapply(mean),
           tm.nn = map(tm.nn, fish.prey) %>% sapply(mean),
           wb.raw = map(wb.raw, fish.prey) %>% sapply(mean)) %>%
  plot(wb.raw~obs, data = .)



# saveRDS(inf.niche, "Niche pruned trait matching inference.RDS")
# saveRDS(inf.neutral, "Neutral trait matching inference.RDS")
# saveRDS(inf.niche.neutral, ("Neutral + Niche trait matching inference.RDS"))


mean(flatten_dbl(map(inf.niche, fish.prey)))
mean(flatten_dbl(map(inf.neutral, fish.prey)))
mean(flatten_dbl(map(inf.niche.neutral, fish.prey)))
