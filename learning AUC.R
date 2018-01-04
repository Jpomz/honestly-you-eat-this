# learning AUC
# 4 Jan 2018


# https://www.r-bloggers.com/calculating-auc-the-area-under-a-roc-curve/

library(pROC)
category <- c(1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0)
prediction <- rev(seq_along(category))
prediction[9:10] <- mean(prediction[9:10])

roc_obj <- roc(category, prediction)
auc(roc_obj)

roc_df <- data.frame(
  TPR=rev(roc_obj$sensitivities), 
  FPR=rev(1 - roc_obj$specificities), 
  labels=roc_obj$response, 
  scores=roc_obj$predictor)

# help page example
data("aSAH")
roc(aSAH$outcome, aSAH$s100b,
    levels=c("Good", "Poor"))



# example from my data
inf.list <- readRDS("Initial trait matching inference.RDS")
# observed webs
obs.A <- readRDS("observed pred-prey.RDS")
# subset obs_A to only include webs with community data
obs.A <- obs.A[names(obs.A) %in% names(inf.list)]

match_matr <- function (obs, inf){
  # colnames(inf) have already been size-sorted
  index = intersect(rownames(inf), rownames(obs))
  obs = obs[index, index]
  inf = inf[index, index]
  list(observed = obs, inferred = inf)
}

test <- match_matr(obs.A[[11]], inf.list[[11]])

inf.prob <- rbinom(length(test$observed), 1, 0.1)
roc(as.numeric(test$observed), inf.prob)


simp_roc <- roc(as.numeric(test$observed), as.numeric(test$inferred))
with(simp_roc, lines(1 - FPR, TPR, col="blue", lty=2))





obs <- as.numeric(obs.A[[11]])
inf1 <- rbinom(length(obs), 1, 0)
inf2 <- rbinom(length(obs), 1, 0.25)
inf3 <- rbinom(length(obs), 1, 0.5)
inf4 <- rbinom(length(obs), 1, 0.75)
inf5 <- rbinom(length(obs), 1, 1)

df <- data.frame(obs = rep(obs, 5), 
           inf = c(inf1, inf2, inf3, inf4, inf5),
           threshold = rep(c(0, .25, .5, .75, 1),
                           each = length(obs)))
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
