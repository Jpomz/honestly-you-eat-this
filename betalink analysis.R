# betalink 

# install.packages("betalink") # one time only
library(plyr)
library(tidyverse)
library(betalink)
library(igraph)

obs <- readRDS("observed matrices matched to inferred.RDS")
wb <- readRDS("wb for PCA.RDS")
tm <- readRDS("tm nn fish for PCA.RDS")
wb.tm <- readRDS("wb x tm for PCA.rds")

# betalink(graph1, graph2)

obs.mat <- obs %>% llply(function (x){as.matrix(x)})
obs.graph <- obs.mat %>%
  llply( function (x){
    graph_from_adjacency_matrix(x)
  })
obs.beta <- beta_os_prime(obs.graph)
plot(density(obs.beta), main = "B'os observed")
abline(v = mean(obs.beta))

tm.mat <- tm %>%
  llply(function (x){as.matrix(x)})
tm.graph <- tm.mat %>%
  llply( function (x){
    graph_from_adjacency_matrix(x)
  })
tm.beta <- beta_os_prime(tm.graph)
plot(density(tm.beta), main = "B'os Trait matching")
abline(v = mean(tm.beta))

wb.mat <- wb %>%
  llply(function (x){as.matrix(x)})
wb.graph <- wb.mat %>%
  llply( function (x){
    graph_from_adjacency_matrix(x)
  })
wb.beta <- beta_os_prime(wb.graph)
plot(density(wb.beta), main = "B'os WebBuilder")
abline(v = mean(wb.beta))

wb.tm.mat <- wb.tm %>%
  llply(function (x){as.matrix(x)})
wb.tm.graph <- wb.tm.mat %>%
  llply( function (x){
    graph_from_adjacency_matrix(x)
  })
wb.tm.beta <- beta_os_prime(wb.tm.graph)
plot(density(wb.tm.beta), main = "B'os WebBuilder")
abline(v = mean(wb.tm.beta))


# one plot
plot(density(obs.beta), lty = 1, main = NA,
     xlim = c(0,0.8), ylim = c(0, 6),
     xlab = NA)
abline(v = mean(obs.beta), lty = 1)
par(new = T)
plot(density(tm.beta), lty = 2, main = NA,
     xlim = c(0,0.8), ylim = c(0, 6),
     xlab = NA)
abline(v = mean(tm.beta), lty = 2)
par(new = T)
plot(density(wb.beta), main = NA, lty = 3,
     xlim = c(0,0.8), ylim = c(0, 6),
     xlab = NA)
abline(v = mean(wb.beta), lty = 3)
par(new = T)
plot(density(wb.tm.beta), main = NA, lty = 3,
     xlim = c(0,0.8), ylim = c(0, 6),
     xlab = NA)
abline(v = mean(wb.tm.beta), lty = 3)







beta_os <- as.data.frame(rbind(obs.beta, grav.beta, wb.beta))
beta_os$model <- rownames(beta_os)
rownames(beta_os) <- NULL
  
beta_os <- beta_os %>% gather("site", "dis", 1:10)

ggplot(beta_os, aes(x = site, y = dis, color = model))+
  geom_point()

beta_os %>% group_by(model) %>% summarise(mean(dis))

# #betalink(obs.graph[[1]], grav.graph[[1]]) %>% unlist 
# 
# # dissim ####
# # not sure what the point of these dissim calcs was...
# dissim.grav <- NULL
# for (i in 1:length(obs.graph)){
#   dis <- betalink(obs.graph[[i]], grav.graph[[i]]) %>%
#     unlist
#   dissim.grav <- rbind(dissim.grav, dis)
#   rownames(dissim.grav)[i] <- names(obs.graph)[i]
# }
# 
# 
# dissim.wb <- NULL
# for (i in 1:length(obs.graph)){
#   dis <- betalink(obs.graph[[i]], wb.graph[[i]]) %>%
#     unlist
#   dissim.wb <- rbind(dissim.wb, dis)
#   rownames(dissim.wb)[i] <- names(obs.graph)[i]
# }
# 
# dissim <- cbind(dissim.grav, dissim.wb)[,c(2,6)]
# colnames(dissim) <- c("gravel", "wb")
# dissim %>% as.data.frame() %>% mutate(dif = gravel - wb)
# 
# dissim.inf <- NULL
# for (i in 1:length(wb.graph)){
#   dis <- betalink(grav.graph[[i]], wb.graph[[i]]) %>%
#     unlist
#   dissim.inf <- rbind(dissim.inf, dis)
#   rownames(dissim.inf)[i] <- names(wb.graph)[i]
# }
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# 
# # remove rare taxa ####
# # remove interactions of relatively rare taxa from observed matrices
# # re-run beta_os_prime and see if there is a difference
# # does the interaction variability of the observed webs disappear when rare taxa links are removed?
# # Is there evidence of interaction variability between common species?
# 
# ab.matr.list <- readRDS("18 relative abundance matrices list.rds")
# names(ab.matr.list)[12] <- "Little"
# 
# obs.mat <- obs %>% llply(function (x){as.matrix(x)})
# obs.graph <- obs.mat %>%
#   llply( function (x){
#     graph_from_adjacency_matrix(x)
#   })
# obs.beta <- beta_os_prime(obs.graph)
# 
# threshold <- c(1e-01, 1e-02, 1e-03, 1e-04, 1e-05, 1e-06, 1e-07, 1e-08, 5e-02, 5e-03, 5e-04, 5e-05, 5e-06, 5e-07, 5e-08)
# threshold <- sort(threshold)
# 
# # this takes a long time...
# # ~ 20 - 30 sec
# ptm <- proc.time()
# obs.ab.bos.list <- NULL
# thresh.list <- NULL
# for (j in 1:length(threshold)){
#   thresh <- NULL
#   for (i in 1:length(ab.matr.list)){
#     ab <- ab.matr.list[[i]]
#     obs.mat <- obs[[i]]
#     thresh.matr <- ab
#     thresh.matr[thresh.matr > threshold[j]] <- 1
#     thresh.matr[thresh.matr < 1] <- 0
#     obs.ab <- obs.mat * thresh.matr
#     obs.ab.mat <- as.matrix(obs.ab)
#     obs.ab.graph <- graph_from_adjacency_matrix(obs.ab.mat)
#     thresh[[i]] <- obs.ab.graph
#     names(thresh)[i] <- names(obs)[i]
#   }
#   obs.ab.bos.list[[j]] <- beta_os_prime(thresh)
# }
# proc.time() - ptm
# names(obs.ab.bos.list) <- as.character(threshold)
# 
# for (i in 1:length(obs.ab.bos.list)){
#   plot(density(na.omit(obs.ab.bos.list[[i]])),
#        main = paste("Observed Beta, removing links at threshold =",
#                     names(obs.ab.bos.list)[i],
#                     sep = " ")
#   )
# }
# 
# obs.ab.bos.df <- ldply(obs.ab.bos.list)
# obs.ab.bos.df <- obs.ab.bos.df %>%
#   gather("site", "bos", 2:11)
# obs.ab.bos.df$.id <- as.numeric(obs.ab.bos.df$.id)
# ggplot(obs.ab.bos.df, aes(x = log10(.id),
#                           y=bos,
#                           color = site )) +
#   geom_point(size = 4)
# 
# # distribution of Bos
# # need to look at poisot again and 
# # read interpretation of Z score, CVD 
# z.score <- obs.ab.bos.df %>% group_by(.id) %>%
#   mutate(mean = mean(bos),
#             sd = sd(bos),
#             dist = dnorm(bos),
#          z = (mean(dist) - mean)/sd)
# 
# ggplot(z.score, 
#        aes(x = log10(.id),
#            y = z,
#            color = site)) +
#   geom_point()
# 
# ggplot(z.score, 
#        aes(x = bos,
#            y = dist,
#        color = as.factor(.id))) +
#   geom_line()
# 
# 
# 
# 
# 
# cvd <- (1+1/(4 * 10)) * (sd(obs.ab.bos.df %>% filter(.id == 1e-04)%>% select(bos) %>% .[[1]]) / mean(obs.ab.bos.df %>% filter(.id == 1e-04)%>% select(bos) %>% .[[1]]))
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # # TSS of webbuilder V gravel
# # wb.grav <- NULL
# # for (i in 1:length(wb)){
# #   temp <- tss(wb[[i]], grav3[[i]])
# #   wb.grav <- rbind(wb.grav, temp)
# # }
# # grav.wb <- NULL
# # for (i in 1:length(wb)){
# #   temp <- tss(grav3[[i]], wb[[i]])
# #   grav.wb <- rbind(grav.wb, temp)
# # }
