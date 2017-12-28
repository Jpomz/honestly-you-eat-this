# figs for manuscript

source("C:\\Users\\Justin\\Documents\\Data\\FW modelling Petchey Github\\ttl-resources-master\\food_web\\FoodWebFunctions.r")
# Gravel ####
# read in data for plotting TSS, a, b, c, d
tss.3.step <- readRDS(paste(getwd(), "/figs for MS/gravel tss 3 steps.rds", sep =""))

# change "step to numeric"
step.no <- data.frame(step.no = c(1, 2, 3),
                      step = c("model", "prune", "abundance"))

tss.3.step <- left_join(tss.3.step, step.no)

tss.3.step %>% filter(step.no == 3) %>%
ggplot(aes(x = fct_reorder(web, S), y = abar, color = "True Pos")) +
  geom_point() +
  geom_point(aes(y = dbar, color = "True neg"))+
  geom_point(aes(y = bbar, color = "False pos"))+
  geom_point(aes(y = cbar, color = "False neg"))

# calc FPR FNR via wiki page
# https://en.wikipedia.org/wiki/Confusion_matrix
# 2 Dec 2017
rates <- tss.3.step %>% group_by(step.no)  %>% mutate(FPR = b / (b + d), FNR = c / (c+a)) %>% summarize(FPR.m = mean(FPR), FPR.sd = sd(FPR), FNR.m = mean(FNR), FNR.sd = sd(FNR)) 
ggplot(rates, aes(x = step.no, y = FPR.m))+
  geom_point(shape = 15, size = 4) +
  geom_errorbar(aes(x = step.no,
                    ymax = FPR.m + FPR.sd,
                    ymin = FPR.m - FPR.sd),
                width = 0.1) +
  geom_point(aes(x = step.no, y = FNR.m), shape = 10, size = 4) +
  geom_errorbar(aes(x = step.no,
                    ymax = FNR.m + FNR.sd,
                    ymin = FNR.m - FNR.sd),
                width = 0.1) 

# get means / sd
tss.3.step.summ <- tss.3.step %>%
  group_by(step.no) %>%
  summarize(tss.mean = mean(tss), 
            tss.sd = sd(tss), 
            a.mean = mean(abar),
            a.sd = sd(abar),
            d.mean = mean(dbar), 
            d.sd = sd(dbar),
            b.mean = mean(bbar),
            b.sd = sd(bbar),
            c.mean = mean(cbar), 
            c.sd = sd(cbar))

# make table of trait-matching false/true positive
tss.3.step.summ %>% select(a.mean, a.sd, b.mean, b.sd) %>% write.csv("trait.matching true false pos.csv")

# make plots 
#TSS 
gravel.tss.mean.plot <- ggplot(tss.3.step.summ,
                               aes(x = step.no,
                            y = tss.mean)) +
  geom_point(shape = 15, size = 4) +
  geom_errorbar(aes(x = step.no,
                    ymax = tss.mean + tss.sd,
                    ymin = tss.mean - tss.sd),
                width = 0.1) +
  scale_x_continuous(breaks = c(1, 2, 3)) +
  geom_line() +
  theme_classic()+
  theme(panel.border = element_rect(colour = "black",
                                    fill=NA, size=1))+
  labs(y = "Mean TSS", x = "Step") +
  ylim(c(0,1))
  
 saveRDS(gravel.tss.mean.plot, file = paste(getwd(), "/figs for MS/gravel mean tss plot.rds", sep = ""))

# tss <- readRDS(file = paste(getwd(), "/figs for MS/gravel mean tss plot.rds", sep = ""))

# # all decomp TSS
# (gravel.mean.tpr.fp.plot <- ggplot(tss.3.step.summ, 
#        aes(x = step.no, y = a.mean, color = "TPR")) +
#   geom_point(size = 4, shape = 15) +
#   geom_line() +
#   geom_errorbar(aes(x = step.no, ymax = a.mean + a.sd,
#                     ymin = a.mean - a.sd),
#                 width = 0.1) +
#   geom_point(aes(x = step.no, y = b.mean, color = "FP"), 
#              size = 4, shape = 0) +
#   geom_line(aes(step.no, b.mean, color = "FP")) +
#   geom_errorbar(aes(x = step.no, ymax = b.mean + b.sd,
#                     ymin = b.mean - b.sd, color = "FP"),
#                 width = 0.1) +
#   geom_point(aes(x = step.no, y = c.mean, color = "FN"), 
#              size = 4, shape = 0) +
#   geom_line(aes(step.no, c.mean, color = "FN")) +
#   geom_errorbar(aes(x = step.no, ymax = c.mean + c.sd,
#                     ymin = c.mean - c.sd, color = "FN"),
#                 width = 0.1) +
#   geom_point(aes(x = step.no, y = d.mean, color = "TN"), 
#              size = 4, shape = 17) +
#   geom_line(aes(step.no, d.mean, color = "TN")) +
#   geom_errorbar(aes(x = step.no, ymax = d.mean + d.sd,
#                     ymin = d.mean - d.sd, color = "TN"),
#                 width = 0.1) +
#   scale_x_continuous(breaks = c(1, 2, 3)) +
#   theme_classic() +
#   theme(panel.border = element_rect(colour = "black",
#                        fill=NA, size=1))+
#   labs(y = "Proportion of links", x = "Step") +
#   ylim(c(0,1)))
# saveRDS(gravel.mean.tpr.fp.plot, file = paste(getwd(), "/figs for MS/gravel mean decomposed tss plot.rds", sep = ""))

# gravel 3 step
# False Positive / False Negative
ggplot(tss.3.step.summ, 
       aes(x = step.no, y = c.mean, color = "False negative")) +
  geom_line()  +
  geom_errorbar(aes(x = step.no, ymax = c.mean + c.sd,
                    ymin = c.mean - c.sd),
                width = 0.1) +
  geom_line(aes(step.no, b.mean, color = "False positive")) +
  geom_errorbar(aes(x = step.no, ymax = b.mean + b.sd,
                    ymin = b.mean - b.sd, color = "False positive"),
                width = 0.1) +
  scale_x_continuous(breaks = c(1, 2, 3)) +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black",
                                    fill=NA, size=1))+
  labs(y = "Proportion of links", x = "Step") +
  ylim(c(-0.0005,1))


# decomp.tss <- readRDS(file = paste(getwd(), "/figs for MS/gravel mean decomposed tss plot.rds", sep = ""))

# webbuilder TSS

wb.tss <- readRDS( file = "webbuilder tss model inferred.rds")

# webbuilder decomp TSS ~ web ####
ggplot(wb.tss, aes(x = fct_reorder(web, S), y = abar, color = "True Pos")) +
  geom_point() +
  geom_point(aes(y = dbar, color = "True neg"))+
  geom_point(aes(y = bbar, color = "False pos"))+
  geom_point(aes(y = cbar, color = "False neg"))

# compare decomp TSS grav / wb by site
tss.3.step[47,12] <- "Little"
wb.grav.tss <- tss.3.step %>% filter(step.no == 3) %>% bind_rows(wb.tss)
wb.grav.tss$step[wb.grav.tss$step == "abundance"] <- "gravel step 3"
wb.grav.tss$step[wb.grav.tss$step == "model"] <- "WebBuilder"

ggplot(wb.grav.tss, aes(x = fct_reorder(web, S),
                        y = abar, 
                        shape = step,
                        color = step)) +
  geom_point(size = 4) +
  ggtitle("True Positive") +
  theme_classic() +
  theme(legend.title=element_blank())

ggplot(wb.grav.tss, aes(x = fct_reorder(web, S),
                        y = bbar,
                        shape = step, 
                        color = step)) +
  geom_point(size = 4) +
  ggtitle("False Positive") +
  theme_classic() +
  theme(legend.title=element_blank())

ggplot(wb.grav.tss, aes(x = fct_reorder(web, S),
                        y = cbar, 
                        shape = step, 
                        color = step)) +
  geom_point(size = 4) +
  ggtitle("False Negative")+
  theme_classic() +
  theme(legend.title=element_blank())

ggplot(wb.grav.tss, aes(x = fct_reorder(web, S),
                        y = dbar, 
                        shape = step,
                        color = step)) +
  geom_point(size = 4) +
  ggtitle("True Negative")+
  theme_classic() +
  theme(legend.title=element_blank())

ggplot(wb.grav.tss, aes(x = fct_reorder(web, S),
                        y = tss, 
                        shape = step,
                        color = step)) +
  geom_point(size = 4) +
  ggtitle("TSS")+
  theme_classic() +
  theme(legend.title=element_blank())


# 
#   geom_point(aes(y = dbar, color = "True neg"))+
#   geom_point(aes(y = bbar, color = "False pos"))+
#   geom_point(aes(y = cbar, color = "False neg"))


# wb TSS summary ####
wb.tss.summ <- wb.tss %>%
  summarize(tss.mean = mean(tss), 
            tss.sd = sd(tss), 
            a.mean = mean(abar),
            a.sd = sd(abar),
            d.mean = mean(dbar), 
            d.sd = sd(dbar),
            b.mean = mean(bbar),
            b.sd = sd(bbar),
            c.mean = mean(cbar), 
            c.sd = sd(cbar))


# # plot of tss for gravel 3 steps, 
# # with tss for WB overlaid 
# ggplot(tss.3.step.summ, aes(x = step.no,
#                             y = tss.mean)) +
#   geom_point(shape = 15, size = 4) +
#   geom_errorbar(aes(x = step.no, ymax = tss.mean + tss.sd,
#                     ymin = tss.mean - tss.sd),
#                 width = 0.1) +
#   scale_x_continuous(breaks = c(1, 2, 3)) +
#   geom_line() +
#   theme_classic()+
#   theme(panel.border = element_rect(colour = "black",
#                                     fill=NA, size=1))+
#   labs(y = "Mean TSS", x = "Step") +
#   ylim(c(0,1)) +
#   geom_point(data = wb.tss.summ, aes(x = 3, y = tss.mean)) +
#   geom_errorbar(data = wb.tss.summ, 
#     aes(x = 3,
#         ymax = tss.mean + tss.sd,
#         ymin = tss.mean - tss.sd),
#                 width = 0.1, color = "red")

# # Decomposed TSS for webbuilder
# ggplot(wb.tss.summ, 
#        aes(x = 1, y = a.mean, color = "TPR")) +
#   geom_point(size = 4, shape = 15) +
#   geom_errorbar(aes(x = 1, ymax = a.mean + a.sd,
#                     ymin = a.mean - a.sd),
#                 width = 0.01) +
#   geom_point(aes(x = 1, y = b.mean, color = "FP"), 
#              size = 4, shape = 0) +
#   geom_errorbar(aes(x = 1, ymax = b.mean + b.sd,
#                     ymin = b.mean - b.sd, color = "FP"),
#                 width = 0.01) +
#   geom_point(aes(x = 1, y = c.mean, color = "FN"), 
#              size = 4, shape = 15) +
#   geom_errorbar(aes(x = 1, ymax = c.mean + c.sd,
#                     ymin = c.mean - c.sd, color = "FN"),
#                 width = 0.01) +
#   geom_point(aes(x = 1, y = d.mean, color = "TN"), 
#              size = 4, shape = 2) +
#   geom_errorbar(aes(x = 1, ymax = d.mean + d.sd,
#                     ymin = d.mean - d.sd, color = "TN"),
#                 width = 0.01) +
#   scale_x_continuous(breaks = c(1)) +
#   theme_classic() +
#   theme(panel.border = element_rect(colour = "black",
#                                     fill=NA, size=1))+
#   labs(y = "Proportion of links", x = "Step") +
#   ylim(c(0,1))


# webbuilder bar plot of all TSS components
ggplot(wb.tss.summ, aes(x = 1, y = a.mean, fill = "TPR")) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(x = 1, ymax = a.mean + a.sd,
                    ymin = a.mean - a.sd), width = 0.1) +
  geom_bar(aes(x = 2, y = b.mean, fill = "FP"), 
           stat = "identity") +
  geom_errorbar(aes(x = 2, ymax = b.mean + b.sd,
                    ymin = b.mean - b.sd), width = 0.1) +
  geom_bar(aes(x = 3, y = c.mean, fill = "FN"), stat = "identity") +
  geom_errorbar(aes(x = 3, ymax = c.mean + c.sd,
                    ymin = c.mean - c.sd), width = 0.1) +
  geom_bar(aes(x = 4, y = d.mean, fill = "TN"), stat = "identity") +
  geom_errorbar(aes(x = 4, ymax = d.mean + d.sd,
                    ymin = d.mean - d.sd),  width = 0.1) +
  scale_x_continuous(breaks = c(1, 2, 3, 4)) +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black",
                                    fill=NA, size=1))+
  labs(y = "Proportion of links") +
  ylim(c(0,1))

# # webbuilder bar plot of all False Pos False Neg components
# ggplot(wb.tss.summ, aes(x = 1, y = b.mean, fill = "False positive")) +
#   geom_bar(stat = "identity") +
#   geom_errorbar(aes(x = 1, ymax = b.mean + b.sd,
#                     ymin = b.mean - b.sd), width = 0.1) +
#   geom_bar(aes(x = 2, y = c.mean, fill = "False negative"), stat = "identity") +
#   geom_errorbar(aes(x = 2, ymax = c.mean + c.sd,
#                     ymin = c.mean - c.sd), width = 0.1) +
#   scale_x_continuous(breaks = c(1, 2)) +
#   theme_classic() +
#   theme(panel.border = element_rect(colour = "black",
#                                     fill=NA, size=1))+
#   labs(y = "Proportion of links") +
#   ylim(c(0,1))



# combine tss.sum for both models
tss.combined <- bind_rows(tss.3.step.summ, wb.tss.summ)


# write_csv(tss.combined, paste(getwd(), "/figs for MS/mean tss both models RAW.csv", sep = ""))



# FIX DATA FOR EASIER PLOTTING
names(tss.combined)
tss.combined$model <- c(rep("Gravel",3), "WebBuilder")








# Plot of gravel TSS ~ relative abundance threshold values
# from script 8
# this is the plot in power point ####
# tss.grav.ab <- readRDS("gravel tss real ab correction.rds")
# ggplot(tss.grav.ab, aes(x = log10(threshold), y = tss.mean))+
#   geom_point()+
#   geom_line() +
#   geom_errorbar(aes(ymin = tss.mean - tss.sd, 
#                     ymax = tss.mean + tss.sd)) +
#   ggtitle("NPred * Nprey") +
#   labs(y = "Mean TSS (+- SD)") +
#   theme_classic()





# # food web stats ####

# modelled vs observed stats ####
# compare empirical and modeled stats
all.stats <- readRDS(file = "food web stats all.rds")
all.stats <- all.stats %>% mutate(LS = L / S)
# remove gravel steps 2 and 3 
model.stats <- all.stats[35:85,]
model.stats <- model.stats %>% arrange(.id, model)
write.csv(model.stats, "model v obs fw stats.csv")
# split into a list
stats.list <- model.stats %>% split(., .$model)
# remove row numbers for combining later
stats.list <- llply(stats.list, function (x){
  rownames(x) <- NULL; x
})
# make object of observed stats
obs <- stats.list$observed
# paste "obs" at beginning of all colnames
# allows for easier plotting below
colnames(obs) <- paste("obs", colnames(obs), sep = ".") 
# remove observed item from list
stats.list <- stats.list[-2]

# bind rows for gravel and wb models ([[1]] and [[2]], respectively)
# then bind columns of obseved stats
exp.obs <- rbind(stats.list[[1]], stats.list[[2]]) %>% cbind(obs)

# all of the following plots have the same outlay
# x = observed metric
# y = modelled metric
# red line = gravel
# blue line = webbuilder
# black line = 1:1
ggplot(exp.obs, aes(x = obs.C, y = C, shape = model, color = model)) +
  geom_point() +
  #stat_smooth(method = "lm", alpha = 0.1) +
  geom_abline(aes(intercept = 0, slope = 1)) +
  coord_cartesian(ylim = c(0, 0.4)) +
  theme_classic()

ggplot(exp.obs, aes(x = obs.L, y = L, shape = model, color = model)) +
  geom_point() +
  #stat_smooth(method = "lm", alpha = 0.1) +
  geom_abline(aes(intercept = 0, slope = 1)) +
  coord_cartesian(ylim = c(0, 400)) +
  theme_classic()
# removing highest link web
ggplot(exp.obs[c(-7,-24),], aes(x = obs.L, y = L, shape = model, color = model)) +
  geom_point() +
#  stat_smooth(method = "lm", alpha = 0.1) +
  geom_abline(aes(intercept = 0, slope = 1)) +
  coord_cartesian(ylim = c(0, 300)) +
  theme_classic()

ggplot(exp.obs, aes(x = obs.LS, y = LS, shape = model, color = model)) +
  geom_point() +
#  stat_smooth(method = "lm", alpha = 0.1) +
  geom_abline(aes(intercept = 0, slope = 1)) +
  coord_cartesian(ylim = c(0, 11)) +
  theme_classic()

b.plot <- ggplot(exp.obs, aes(x = obs.B, y = B, shape = model, color = model)) +
  geom_point() +
#  stat_smooth(method = "lm", alpha = 0.1) +
  geom_abline(aes(intercept = 0, slope = 1)) +
  coord_cartesian(ylim = c(0, 0.4)) +
  theme_classic()

i.plot <- ggplot(exp.obs, aes(x = obs.I, y = I, shape = model, color = model)) +
  geom_point() +
#  stat_smooth(method = "lm", alpha = 0.1) +
  geom_abline(aes(intercept = 0, slope = 1)) +
  coord_cartesian(ylim = c(0, 0.9)) +
  theme_classic()

t.plot <- ggplot(exp.obs, aes(x = obs.T, y = T, shape = model, color = model)) +
  geom_point() +
  stat_smooth(method = "lm", alpha = 0.1) +
  geom_abline(aes(intercept = 0, slope = 1)) +
  coord_cartesian(ylim = c(0, 0.3)) +
  theme_classic()

library(gridExtra)

grid.arrange(b.plot, i.plot, t.plot)

ggplot(exp.obs, aes(x = obs.Gensd, y = Gensd, shape = model, color = model)) +
  geom_point() +
#  stat_smooth(method = "lm", alpha = 0.1) +
  geom_abline(aes(intercept = 0, slope = 1)) +
  coord_cartesian(xlim = c(1, 4.1),
                  ylim = c(0.5, 2)) +
  theme_classic()

ggplot(exp.obs, aes(x = obs.Vulsd, y = Vulsd, shape = model, color = model)) +
  geom_point() +
#  stat_smooth(method = "lm", alpha = 0.1) +
  geom_abline(aes(intercept = 0, slope = 1)) +
  coord_cartesian(xlim = c(0.5, 2.6)) +
  theme_classic()

ggplot(exp.obs, aes(x = obs.Maxsim, y = Maxsim, shape = model, color = model)) +
  geom_point() +
#  stat_smooth(method = "lm", alpha = 0.1) +
  geom_abline(aes(intercept = 0, slope = 1)) +
  coord_cartesian(ylim = c(0.25, 0.9)) +
  theme_classic()

(mean.tl.plot <- ggplot(exp.obs, aes(x = obs.mean.TL, y = mean.TL, shape = model, color = model)) +
  geom_point() +
  stat_smooth(method = "lm", alpha = 0.1) +
  geom_abline(aes(intercept = 0, slope = 1)) +
  coord_cartesian(ylim = c(1, 7.5)) +
  theme_classic()
)

(max.tl.plot <- ggplot(exp.obs, aes(x = obs.max.TL, y = max.TL, shape = model, color = model)) +
  geom_point() +
  stat_smooth(method = "lm", alpha = 0.1) +
  geom_abline(aes(intercept = 0, slope = 1)) +
  coord_cartesian(ylim = c(2, 11)) +
  theme_classic()
)
(sd.tl.plot <- ggplot(exp.obs, aes(x = obs.sd.TL, y = sd.TL, shape = model, color = model)) +
  geom_point() +
#  stat_smooth(method = "lm", alpha = 0.1) +
  geom_abline(aes(intercept = 0, slope = 1)) +
  coord_cartesian(ylim = c(0, 4)) +
  theme_classic()
)
grid.arrange(max.tl.plot, mean.tl.plot, sd.tl.plot)

# stats relationship ####
# relationship between obs / mod stats
# c.mod <- lm(C~obs.C*model, data= exp.obs)
# l.mod <- lm(L ~ obs.L*model, data = exp.obs)
# i.mod <- lm(I~obs.I*model, data = exp.obs)
# t.mod <- lm(T~obs.T*model, data = exp.obs)

# for loop calculating statistical relationships between observed and inferred statistics

# make numeric vector of column variables want to test
col.vars <- c(3:16, 19)
# for loop for webbuilder inferred
wb.df.result <- NULL
for (i in col.vars){
  j = i + 19
  dat = exp.obs[exp.obs$model == "webbuilder",]
  y = dat[,i]
  x = dat[,j]
  mod = lm(y~x)
  x.var = colnames(dat)[i]
  p = summary(mod)$coefficients[8]
  result <- c(x.var,p)
  wb.df.result <- rbind(wb.df.result, result)
}
# save results
write.csv(wb.df.result, file = paste(getwd(), "/figs for MS/wb fw stats obs pred.csv", sep =""))
# for loop for gravel inferred
grav.df.result <- NULL
for (i in col.vars){
  j = i + 19
  dat = exp.obs[exp.obs$model == "gravel",]
  y = dat[,i]
  x = dat[,j]
  mod = lm(y~x)
  x.var = colnames(dat)[i]
  p = summary(mod)$coefficients[8]
  result <- c(x.var,p)
  grav.df.result <- rbind(grav.df.result, result)
}
# save results
write.csv(grav.df.result, file = paste(getwd(), "/figs for MS/grav fw stats obs pred.csv", sep =""))


# predation matrices ####
source("C:\\Users\\Justin\\Documents\\Data\\FW modelling Petchey Github\\ttl-resources-master\\food_web\\FoodWebFunctions.r")
grav.primary <- readRDS("gravel inferred initial.rds")
grav.pruned <- readRDS("gravel inferred pruned.rds")
grav.ab.inf <- readRDS("gravel inferred pruned rel ab.rds")
wb.matr <- readRDS("webbuilder inferred list trimmed taxa.rds")
obs <- readRDS("observed adj matr list trimmed taxa.rds")

# plot largest, Dempsters; [[7]]
Plot.matrix2(grav.primary[[7]], sp.pt.ch = 18, point.cex = 1.5, title = "Biomass", axes.labels = T)

Plot.matrix2(grav.pruned[[7]], sp.pt.ch = 18, point.cex = 1.5, title = "Pruned", axes.labels = T)

Plot.matrix2(grav.ab.inf[[7]], sp.pt.ch = 18, point.cex = 1.5, title = "Relative abundance", axes.labels = T)

Plot.matrix2(wb.matr[[7]], sp.pt.ch = 18, point.cex = 1.5, title = "WebBuilder", axes.labels = T)

Plot.matrix2(obs[[7]], sp.pt.ch = 18, point.cex = 1.5, title = "Observed", axes.labels = T)

# dempsters, one plot ####
png(filename = "C:\\Users\\Justin\\Documents\\UC PhD Thesis\\Manuscripts\\inferring food webs\\R fig1.png", width = 240,
    height = 180, units = "mm", res =300)
#layout parameters
# 3 x 3 matrix
# 1,2,3 == first three plots in firts row
# 0,0,4 == one plot in row 2, col 3
# 0,0,5 == one plot in row 3, col3
layout(matrix(c(1,2,3,0,0,4,0,0,5), 3, 3, byrow = T))
#biomass inferred
Plot.matrix2(grav.primary[[7]], sp.pt.ch = 18, point.cex = 1)
title(main = "Biomass inference", cex.lab = 2)
mtext(text = "A ", side = 2, cex = 1.5, las = 1, padj = 0)
# niche pruned
Plot.matrix2(grav.pruned[[7]], sp.pt.ch = 18, point.cex = 1)
title(main = "Niche forbidden", cex.lab = 2)
mtext(text = "B ", side = 2, cex = 1.5, las = 1, padj = 0)
# neutrally pruned
Plot.matrix2(grav.ab.inf[[7]], sp.pt.ch = 18, point.cex = 1)
title(main = "Neutrally forbidden", cex.lab = 2)
mtext(text = "C ", side = 2, cex = 1.5, las = 1, padj = 0)
# webbuilder
Plot.matrix2(wb.matr[[7]], sp.pt.ch = 18, point.cex = 1)
title(main = "WebBuilder", cex.lab = 2)
mtext(text = "D ", side = 2, cex = 1.5, las = 1, padj = 0)
#empirical
Plot.matrix2(obs[[7]], sp.pt.ch = 18, point.cex = 1)
title(main = "Empirical", cex.lab = 2)
mtext(text = "E ", side = 2, cex = 1.5, las = 1, padj = 0)
dev.off() # turns off png() function

# plot smallest; Little; [[11]]
Plot.matrix2(grav.primary[[11]], sp.pt.ch = 18, point.cex = 1.5, title = "Biomass", axes.labels = T)

Plot.matrix2(grav.pruned[[11]], sp.pt.ch = 18, point.cex = 1.5, title = "Pruned", axes.labels = T)

Plot.matrix2(grav.ab.inf[[11]], sp.pt.ch = 18, point.cex = 1.5, title = "Relative abundance", axes.labels = T)

Plot.matrix2(wb.matr[[11]], sp.pt.ch = 18, point.cex = 1.5, title = "WebBuilder", axes.labels = T)

Plot.matrix2(obs[[11]], sp.pt.ch = 18, point.cex = 1.5, title = "Observed", axes.labels = T)

# plot grav step3 webbuilder Observed together
for(i in 1:length(grav.ab.inf)){
  Plot.matrix(grav.ab.inf[[i]], pt.col = "tan", point.cex =2.5)
  par(new = T)
  Plot.matrix(wb.matr[[i]], pt.col = "red", point.cex = 1.5)
  par(new = T)
  Plot.matrix(obs[[i]], point.cex = 1)
}


# supplementary
fish.size <- readRDS(paste(getwd(), "/figs for MS/mean tss 4 fish sizes.rds", sep = ""))

fish.params <- readRDS(paste(getwd(), "/figs for MS/gravel params 4 fish sizes.rds", sep = ""))

fish.ab <- readRDS(paste(getwd(), "/figs for MS/gravel tss fish ab.rds", sep = ""))

# gravel relative abundance TSS
gravel.relab.plot <- readRDS(paste(getwd(), "/figs for MS/gravel rel ab tss plot.rds", sep = ""))

# webbuilder relative abundance TSS
wb.relab.tss.plot <- readRDS(paste(getwd(), "/figs for MS/wb rel ab tss plot.rds", sep = ""))

# Mean TSS ~ relative fish abundance
# data from script 12
# plot in power point
ggplot(tss.grav.ab, aes(x = f.ab, y = tss.mean))+
  geom_point()+
  geom_errorbar(aes(ymin = tss.mean - tss.sd, 
                    ymax = tss.mean + tss.sd))+
  geom_line()+
  scale_x_continuous(breaks = fish.ab) +
  theme_classic()+
  ggtitle("Mean TSS ~ relative fish abundance") +
  labs(x = "Relative fish abundance",
       y = "Mean TSS")

