# calculating AUC
# 5 Jan 2018
# pomeranz

# read in data
# data has already been sorted by TPR & FPR within each threshold
dat <- read.csv("Table for AUC.csv")
# split into list
dat.list <- split(dat, list(dat$Threshold))
# function to calculate AUC using trapezoidal method
get_auc <- function(dat, x, y){
  sum(diff(dat[[x]]) * (head(dat[[y]], -1)+ tail(dat[[y]], -1)))/ 2
}
# apply fxn to list
auc.base <- lapply(dat.list, get_auc, x = "FPR", y = "TPR")


# using pracma::trapz
library(pracma)
auc.trapz <- lapply(dat.list, function (.data){trapz(
  x = .data$FPR, y = .data$TPR)})
  

# compare
dat2 <- data.frame(auc.base = unlist(auc.base),
                   auc.trapz = unlist(auc.trapz))
difference <- dat2$auc.base - dat2$auc.trapz
# mean difference <<< small
mean(difference)


# plot
dat2$threshold <- log10(as.numeric(rownames(dat2)))
plot(auc.base ~ threshold, data = dat2)

# max auc at threshold of 1e-05
dat2[order(dat2$auc.base),]
