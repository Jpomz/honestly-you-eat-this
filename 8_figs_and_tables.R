# figures and tables for MS


# inference stats table
tm <- read.csv("Mean AUC and TSS trait matching.csv")
tm.f <- read.csv("Fish corrected Mean AUC and TSS trait matching.csv")
names(tm.f)[1] <- "inference"
wb <- read.csv("Webbuilder AUC, TSS, fp+fn.csv")
names(wb)[4] <- "Threshold"
wb.tm <- read.csv("Wb x tm AUC, TSS, fp+fn.csv")
names(wb.tm)[4] <- "Threshold"
eu <- read.csv("PCA euclidean distance observed.csv")[,-1] # minus "X"
# put in proper order
eu <- eu[,c(1,2,3,4,6,5,7)]

# bind auc, tss, fp, fn tables
tab <- rbind(tm, tm.f, wb,wb.tm)
# add model variable
tab$model <- c("tm", "tm", "tm", "tm", "tm.f", "tm.f", "wb", "wb", "wb", "wb.tm", "wb.tm", "wb.tm")
# remove unneeded rows
tab <- tab[c(-3,-4,-5,-8,-11),]
# add euclidean distance to tab
tab$euclidean <- t(eu)
# reorder columns
tab <- tab[,c(9,1,2,3,4,5,6,7,8,10)]

write.csv(tab, "figs for MS\\post poisot\\All inference summary stats.csv", row.names = F)
