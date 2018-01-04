# learning AUC
# 4 Jan 2018


# https://www.r-bloggers.com/calculating-auc-the-area-under-a-roc-curve/

library(plyr)
library(pROC)
library(tidyverse)

# example from my data
inf.list <- readRDS("Initial trait matching inference.RDS")
# observed webs
obs.A <- readRDS("observed matrices matched to inferred.RDS")

obs <- flatten_dbl(obs.A)
inf <- flatten_dbl(inf.list)
initial.AUC <- roc(inf, obs, plot = T)





# formula from stock
df$dif <- df$obs - df$inf
auc <- function(dat, group){
  pos = dat[dat>0,]
  neg = dat[dat<0,]
  tot.pos = sum(pos)
  tot.neg = sum(neg)
  m.pos = mean(pos)
  m.neg = mean(neg)
  auc = m.pos * m.neg / 
    (abs(tot.pos) * abs(tot.neg)) 
  auc
}

df %>% group_by(threshold) %>% select(dif) %>% summarize(auc(.))

df2 <- df %>% filter(threshold == 0.5)
auc(df2$dif)

glm_test <- glm(obs ~ threshold, data = df, family = binomial(logit))
pred <- predict(glm_test, df, type = "response")
roc(df$obs, pred, plot = T)

# trying to calculate probabilities based on 
# https://ms.mcmaster.ca/~bolker/emdbook/book.pdf
source("adj_conf_matrix function.R")
c.matrix <- adj_conf_matrix(test$observed, test$inferred)

pred.prob <- data.frame(plus = sum(test$inferred) / 
                          (length(test$inferred)*2),
                        neg = (length(test$inferred) -
                                 sum(test$inferred)) / 
                          (length(test$inferred) * 2),
                        o = sum(test$observed) / 
                          (length(test$inferred) * 2),
                        u = (length(test$observed) - 
                               sum(test$observed)) / 
                          (length(test$inferred) *2))

tp.prob <- (pred.prob$plus * pred.prob$plus *
              pred.prob$o) / (pred.prob$plus *
                                pred.prob$o)
tn.prob <- (pred.prob$neg * pred.prob$neg *
              pred.prob$u) / (pred.prob$neg *
                                pred.prob$u)
fp.prob <- (pred.prob$neg * pred.prob$neg *
              pred.prob$o) / (pred.prob$neg *
                                pred.prob$o)
fn.prob <- (pred.prob$plus * pred.prob$plus *
              pred.prob$u) / (pred.prob$plus *
                                pred.prob$u)
