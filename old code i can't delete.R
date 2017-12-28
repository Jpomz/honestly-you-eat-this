# old / random code

# working out rel ab threshold, WebBuilder

# threshold ####
#********************************************************
# make a threshold value
ab.threshold <- 1e-05
#********************************************************

# make new matrices with abundances converted to binary
# turn values above/below threshold in ab.matr to 1/0 respectively
thresh.matr <- ab.matr.list
for (i in 1:length(thresh.matr)){
  thresh.matr[[i]][thresh.matr[[i]] < ab.threshold] <- 0
  thresh.matr[[i]][thresh.matr[[i]] > ab.threshold] <- 1
}

wb.thresh <- thresh.matr
for (i in 1:length(wb.thresh)){
  wb.thresh[[i]] <- thresh.matr[[i]] * wb.trim.list[[i]]
}



# # plot matrices
# for (i in 1:length(wb.thresh)){
#   Plot.matrix(wb.thresh[[i]], pt.col = "tan", point.cex = 1.7)
#   par(new = T)
#   Plot.matrix(obs.list[[i]], point.cex = 1)
# }

# #TSS ####
# tss() function written June 19 2017

# TSS = (ad - bc)/[(a+c)(b+d)]
# a = TPR = number of links predicted and observed = matrix value == 2
# b = predicted not observed; matrix value == -1
# c = observed but not predicted; matrix value == 1
# d = not predicted, not observed; matrix value ==0

# calculate TSS for
# 1) webbuilder inferred webs [obs.list, wb.trim.list]
# 2) abundance corrected [obs.list, wb.thresh]


tss.wb <- NULL
for (i in 1:length(obs.list)){
  j <- tss(obs.list[[i]], wb.trim.list[[i]])
  j$web <- names(obs.list)[[i]]
  tss.wb <- rbind(tss.wb, j)
  rm(j)
}

tss.wb.mean <- tss.wb %>%
  summarize(tss = mean(tss),
            ad = (mean(a) + mean (d)) /mean(S**2),
            bc = mean(b +c) / mean(S**2))
# webbuilder inferred, not corrected for abundance
# tss = 0.6860638
# ad = 0.8351997
# bc = 0.1648003

tss.wb.ab <- NULL
for (i in 1:length(obs.list)){
  j <- tss(obs.list[[i]], wb.thresh[[i]])
  j$web <- names(obs.list)[[i]]
  tss.wb.ab <- rbind(tss.wb.ab, j)
  rm(j)
}

tss.wb.ab.mean <- tss.wb.ab %>%
  summarize(tss = mean(tss),
            ad = (mean(a) + mean (d)) /mean(S**2),
            bc = mean(b +c) / mean(S**2))
# for threshold == 1e-05
# tss = 0.702392
# ad = 0.8599365
# bc = 0.1400635


# old ####
# code template for saving plots
# saveRDS(ggplot(tss.wb.ab, aes(x = log10(threshold), y = tss))+
#   geom_point()+
#   geom_line() +
#   ggtitle("Npred * Nprey **2 "),
#   file = "webbuilder tss ~ ab threshold Nprey.rds")




# rel ab gravel ####
# old code working out gravel relative abundance correction 


# tss.wb.ab <- NULL
# tss.df.list <- NULL
# for (j in 1:length(threshold)){
#   tss.df <- NULL
#   for (i in 1:length(ab.matr.list)){
#     ab <- ab.matr.list[[i]]
#     grav <- wb.trim.list[[i]]
#     obs <- obs.list[[i]]
#     thresh.matr <- ab
#     thresh.matr[thresh.matr > threshold[j]] <- 1
#     thresh.matr[thresh.matr < 1] <- 0
#     grav.thresh <- grav * thresh.matr
#     tss.temp <- tss(obs, grav.thresh)
#     tss.temp$threshold <- threshold[j]
#     tss.temp$site <- names(ab.matr.list)[i]
#     tss.df <- rbind(tss.df, tss.temp)
#   }
#   tss.df.list[[j]] <- tss.df
#   result <- tss.df %>%
#     summarize(tss = mean(tss),
#               ad = (mean(a) + mean (d)) /mean(S**2),
#               bc = mean(b +c) / mean(S**2),
#               threshold = max(threshold))
#   tss.wb.ab <- rbind(tss.wb.ab, result)
# }



# log10(threshold)
ggplot(tss.grav.ab, aes(x = log10(threshold), y = tss))+
  geom_point()+
  geom_line() +
  ggtitle("NPred Nprey")

ggplot(tss.grav.ab, aes(x = log10(threshold), y = tss))+
  geom_point() +
  stat_smooth(method = "lm", formula = y~ poly(x, 4))
ggplot(tss.grav.ab, aes(x = log10(threshold), y = tss))+
  geom_point()+
  stat_smooth()

# raw threshold
ggplot(tss.grav.ab, aes(x = threshold, y = tss))+
  geom_point() +
  coord_cartesian(xlim = c(0, 0.005))
geom_line()

# threshold ####
#********************************************************
# make a threshold value
ab.threshold <- 1e-04
#********************************************************

# make new matrices with abundances converted to binary
# turn values above/below threshold in ab.matr to 1/0 respectively
thresh.matr <- ab.matr.list
for (i in 1:length(thresh.matr)){
  thresh.matr[[i]][thresh.matr[[i]] < ab.threshold] <- 0
  thresh.matr[[i]][thresh.matr[[i]] > ab.threshold] <- 1
}

# multiply gravel inferred matrices by abundance threshold matrices
grav.thresh <- thresh.matr
for (i in 1:length(grav.thresh)){
  grav.thresh[[i]] <- thresh.matr[[i]] * grav.inf[[i]]
}


# #TSS ####
# tss() function written June 19 2017

# TSS = (ad - bc)/[(a+c)(b+d)]
# a = TPR = number of links predicted and observed = matrix value == 2
# b = predicted not observed; matrix value == -1
# c = observed but not predicted; matrix value == 1
# d = not predicted, not observed; matrix value ==0

# calculate TSS for
# 1) webbuilder inferred webs [obs.list, wb.trim.list]
# 2) abundance corrected [obs.list, wb.thresh]


tss.grav <- NULL
for (i in 1:length(obs.list)){
  j <- tss(obs.list[[i]], grav.inf[[i]])
  j$web <- names(obs.list)[[i]]
  tss.grav <- rbind(tss.grav, j)
  rm(j)
}

tss.grav.mean <- tss.grav %>%
  summarize(tss = mean(tss),
            ad = (mean(a) + mean (d)) /mean(S**2),
            bc = mean(b +c) / mean(S**2))
# webbuilder inferred, not corrected for abundance
# tss = 0.4760068
# ad = 0.5594184
# bc = 0.4405816

tss.grav.ab.thresh <- NULL
for (i in 1:length(obs.list)){
  j <- tss(obs.list[[i]], grav.thresh[[i]])
  j$web <- names(obs.list)[[i]]
  tss.grav.ab.thresh <- rbind(tss.grav.ab.thresh, j)
  rm(j)
}

tss.grav.ab.thresh.mean <- tss.grav.ab.thresh %>%
  summarize(tss = mean(tss),
            ad = (mean(a) + mean (d)) /mean(S**2),
            bc = mean(b +c) / mean(S**2))
# for threshold == 1e-05
# tss = 0.5944057
# ad = 0.6884506
# bc = 0.3115494
tss.grav.ab.mean

# plot matrices
# uncorrected
for (i in 1:length(obs.list)){
  Plot.matrix(grav.thresh[[i]], pt.col = "tan", point.cex = 1.7)
  par(new = T)
  Plot.matrix(obs.list[[i]], point.cex = 1, title = "Abundance corrected")
}


# old code used to determine which relative abundance formula gave best  tss score
# code template used to make plots
# saveRDS(ggplot(tss.grav.ab, aes(x = log10(threshold), y = tss))+
#           geom_point()+
#           geom_line() +
#           ggtitle("NPred Nprey**2"),
#         file = "gravel tss ~ ab threshold Npred Nprey2.rds")


# all plots saved for comparison 
# wb.pred.prey <- readRDS("webbuilder tss ~ ab threshold Npred Nprey.rds")
# wb.pred.prey2 <- readRDS("webbuilder tss ~ ab threshold Npred Nprey2.rds")
# wb.prey <- readRDS("webbuilder tss ~ ab threshold Nprey.rds")
# grav.pred.prey <- readRDS("gravel tss ~ ab threshold Npred Nprey.rds")
# grav.pred.prey2 <- readRDS("gravel tss ~ ab threshold Npred Nprey2.rds")
# grav.prey <- readRDS("gravel tss ~ ab threshold Nprey.rds")
# 
# wb.pred.prey
# wb.pred.prey2
# wb.prey
# grav.pred.prey
# grav.pred.prey2
# grav.prey

