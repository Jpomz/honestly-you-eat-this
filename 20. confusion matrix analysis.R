# 20 confusion matrix analysis
# JPomz 2 Dec 2017

# running confusion matrices analyses per T poisot's
# recommendations. 

# adjacency confusion matrix function
source("adj_conf_matrix function.R")
# just ACC function
source("get_ACC function.R")
# initial trait-match
initial <- readRDS("gravel inferred initial.rds")
# niche pruned trait-match
niche.prune <- readRDS("gravel inferred pruned.rds")
# neutral pruned trait-match
neutral.prune <- readRDS("gravel inferred pruned rel ab.rds")
# webbuilder inferred
wb <- readRDS("webbuilder inferred list trimmed taxa.rds")
# observed
obs <- readRDS("observed webs names and taxa matched.rds")

conf.initial <- ldply(map2(obs, initial, adj_conf_matrix))
conf.neutral <-ldply(map2(obs, neutral.prune, adj_conf_matrix))
conf.niche <- ldply(map2(obs, niche.prune, adj_conf_matrix))

initial.mean <- lapply(na.omit(conf.initial[-1]), mean)
initial.sd <- lapply(na.omit(conf.initial[-1]), sd)
niche.mean <- lapply(na.omit(conf.niche[-1]), mean)
niche.sd <- lapply(na.omit(conf.niche[-1]), sd)
neutral.mean <- lapply(na.omit(conf.neutral[-1]), mean)
neutral.sd <- lapply(na.omit(conf.neutral[-1]), sd)

trait.conf.der <- cbind(initial.mean, initial.sd, niche.mean, niche.sd,  neutral.mean, neutral.sd)

conf.wb <- ldply(map2(obs, wb, adj_conf_matrix))
wb.mean <- lapply(conf.wb[-1], mean)
wb.sd <- lapply(conf.wb[-1], sd)
wb.conf.der <- cbind(wb.mean, wb.sd)

data.frame(mean = lapply(conf.wb[-1], mean),
           sd = lapply(conf.wb[-1], sd))



acc <- data.frame(initial = unsplit(map2(obs,
                            initial,
                            get_ACC),
                          names(obs)),
                  niche = unsplit(map2(obs,
                            niche.prune,
                            get_ACC),
                          names(obs)),
                  neutral = unsplit(map2(obs,
                           neutral.prune,
                           get_ACC),
                          names(obs)),
                  wb = unsplit(map2(obs,
                            wb,
                            get_ACC),
                          names(obs))
)
acc.sum <- data.frame(mean = sapply(acc, mean),
                      sd = sapply(acc, sd))
