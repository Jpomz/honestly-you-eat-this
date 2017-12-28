# method sensitivity to land use effects
# J Pomz August 2017

library(vegan)


# read in stats
all.stats <- readRDS(file = "food web stats all.rds")
# add linkage density (L/S)
all.stats <- all.stats %>% mutate(LS = L / S)
# remove gravel steps 2 and 3 
model.stats <- all.stats[35:85,]


# make data frame with sites + landuse category
sites <- data.frame(.id = unique(model.stats$.id))
# add category (cat)
# make sure that these are in correct order 
sites$cat <- c("Pn", "Pn", "P", "P", "P", "Pn","T",
               "T", "T", "T", "T", "Pn", "N", "N",
               "T", "T", "Pn")
# join the landuse cat to the stats file
model.stats <- left_join(model.stats, sites, by = ".id")

# clean up for better ggplotting
# make the model variable as a factor
model.stats$model <- as.factor(model.stats$model)
# reorder factor for plots
model.stats$model <- factor(
  model.stats$model,
  levels(model.stats$model)[c(2,1,3)])

# ANOVA ####
# one way anova testing different methods
# for differences among land use types
# using aov() in order to perform 
# posthoc Tukey test

# function for getting one way aov
one_way_aov <- function(d, x, y){
  aov(d[[y]] ~ d[[x]])
}


l <- dlply(model.stats,
            .(model),
            one_way_aov,
            x = "cat", y = "L")

ls <-  dlply(model.stats,
             .(model),
             one_way_aov,
             x = "cat", y = "LS")
c <-  dlply(model.stats,
            .(model),
            one_way_aov,
            x = "cat", y = "C")
mcl <-  dlply(model.stats,
              .(model),
              one_way_aov,
              x = "cat", y = "mean.TL")
b <-  dlply(model.stats,
            .(model),
            one_way_aov,
            x = "cat", y = "B")
int <-  dlply(model.stats,
              .(model),
              one_way_aov,
              x = "cat", y = "I")
top <-  dlply(model.stats,
              .(model),
              one_way_aov,
              x = "cat", y = "T")
g <-  dlply(model.stats,
            .(model),
            one_way_aov,
            x = "cat", y = "Gensd")
v <- dlply(model.stats,
           .(model),
           one_way_aov,
           x = "cat", y = "Vulsd")
msim <-  dlply(model.stats,
               .(model),
               one_way_aov,
               x = "cat", y = "Maxsim")
maxtl <-  dlply(model.stats,
                .(model),
                one_way_aov,
                x = "cat", y = "max.TL")

aov.results <- list(l, ls, c, mcl, b, int,
                    top, g, v, msim, maxtl)

# length(aov.res...) == 11
# length each aov.res element == 3
# variables == 11
# each variable has three methods
aov.table <- matrix(0, 11, 3)
for (i in 1:length(aov.results)){
  for (j in 1:length(aov.results[[1]])){
    aov.table[i,j] <- anova(aov.results[[i]][[j]])[1,5]
  }
}
dimnames(aov.table) <- list(
  c("L", "LS", "C", "meanTL", "B", "I", "T", "Gensd", "Vulsd", "Maxsim", "max.TL"),
  c("observed", "gravel", "webbuilder"))
aov.table <- signif(aov.table, 4)
# significant results:
# Obs:    LS, C, B, I, T, Vul, Maxsim
# gravel: LS, C, B,       Vul,        Gen
# wb:                  T, Vul

# TukeyHSD ####

# Tukey L ####
# NS
TukeyHSD(l$observed)
TukeyHSD(l$gravel)
TukeyHSD(l$webbuilder)

# Tukey L/S ####
TukeyHSD(ls$observed)
# T > Pn
TukeyHSD(ls$gravel)
# T > N; T > Pn (p=0.103)
TukeyHSD(ls$webbuilder)

# Tukey C ####
TukeyHSD(c$observed)
# T > Pn
TukeyHSD(c$gravel)
# T > N, Pn > N (p=0.0633)
TukeyHSD(c$webbuilder)

# Tukey meanTL ####
# No significance
TukeyHSD(mcl$observed)
TukeyHSD(mcl$gravel)
TukeyHSD(mcl$webbuilder)

# Tukey B ####
TukeyHSD(b$observed)
# T > N; T > Pn
TukeyHSD(b$gravel)
# N > Pn, T > Pn, P > Pn; P > T; 
TukeyHSD(b$webbuilder)

# Tukey I ####
TukeyHSD(int$observed)
# T > Pn (p = 0.0577);  N > Pn (p=0.0772)
TukeyHSD(int$gravel)
# Pn > P (p=0.0983)
TukeyHSD(int$webbuilder)

# Tukey T ####
TukeyHSD(top$observed)
# Pn > P; Pn > T 
TukeyHSD(top$gravel)
# webbuilder
TukeyHSD(top$webbuilder)
# Pn > T; Pn > N (p=0.0856); Pn > P (p=0.0698)

# Tukey Gensd ####
TukeyHSD(g$observed)
TukeyHSD(g$gravel)
# N > P; N > Pn ; N > T; 
TukeyHSD(g$webbuilder)

# Tukey Vulsd ####
TukeyHSD(v$observed)
# Pn > N
TukeyHSD(v$gravel)
# N > P; N > T; Pn > T (p=0.0782)
TukeyHSD(v$webbuilder)
# N > T (p=0.0826); Pn > T (0.0679)

# Tukey maxsim ####
TukeyHSD(msim$observed)
# T > Pn
TukeyHSD(msim$gravel)
TukeyHSD(msim$webbuilder)

# Tukey maxtl ####
# NS
TukeyHSD(maxtl$observed)
TukeyHSD(maxtl$gravel)
TukeyHSD(maxtl$webbuilder)

# plot food web metrics ####
# functions for plots
plot_bar <- geom_bar(position = position_dodge(),
                     stat = "identity")
plot_error <- geom_errorbar(aes(ymin = avg - sd,
                                ymax = avg +sd),
                            position = position_dodge())
# links
model.stats %>%
  group_by(model, cat) %>%
  summarize(avg = mean(L), sd = sd(L)) %>%
  ggplot(aes(model, avg, fill = cat)) +
  plot_bar +
  plot_error
  
# L/S
model.stats %>%
  group_by(model, cat) %>%
  summarize(avg = mean(L/S), sd = sd(L/S)) %>%
  ggplot(aes(model, avg, fill = cat)) +
  plot_bar +
  plot_error
  
# Connectance
model.stats %>%
  group_by(model, cat) %>%
  summarize(avg = mean(C), sd = sd(C)) %>%
  ggplot(aes(model, avg, fill = cat)) +
  plot_bar +
  plot_error
  
# mean.TL
model.stats %>%
  group_by(model, cat) %>%
  summarize(avg = mean(mean.TL), sd = sd(mean.TL)) %>%
  ggplot(aes(model, avg, fill = cat)) +
  plot_bar +
  plot_error

# ordinations ####
# make cat variable a factor
model.stats$cat <- as.factor(model.stats$cat)

# split stats into list
data.list <- dlply(model.stats, .(model))

# RDA ####
# make ordinations based on method
ord <- dlply(model.stats, .(model),
             function (x){
               rda(x[c(2:7,11:16,19)],
                      scale = T)
})

# variables for ordination plots
scl <- 3 
colvec <- c("red2","mediumblue", "green4", "black")

# observed ordination
plot(ord$observed, type = "n", scaling = scl)
with(data.list$observed,
     points(ord$observed,
            display = "sites",
            col = colvec[cat],
            scaling = scl,
            pch = 21,
            bg = colvec[cat]))
text(ord$observed, display = "species", scaling = scl, cex = 0.8, col = "darkcyan")
with(data.list$observed,
     legend("topright",
            legend = levels(cat),
            bty = "n",
            col = colvec,
            pch = 21,
            pt.bg = colvec))
title("Observed networks RDA")

# gravel step 3 ordination
plot(ord$gravel, type = "n", scaling = scl)
with(data.list$gravel,
     points(ord$gravel,
            display = "sites",
            col = colvec[cat],
            scaling = scl,
            pch = 21,
            bg = colvec[cat]))
text(ord$gravel, display = "species", scaling = scl, cex = 0.8, col = "darkcyan")
with(data.list$gravel,
     legend("topright",
            legend = levels(cat),
            bty = "n",
            col = colvec,
            pch = 21,
            pt.bg = colvec))
title("gravel networks RDA")

# webbuilder ordination
plot(ord$webbuilder, type = "n", scaling = scl)
with(data.list$webbuilder,
     points(ord$webbuilder,
            display = "sites",
            col = colvec[cat],
            scaling = scl,
            pch = 21,
            bg = colvec[cat]))
text(ord$webbuilder, display = "species", scaling = scl, cex = 0.8, col = "darkcyan")
with(data.list$webbuilder,
     legend("topright",
            legend = levels(cat),
            bty = "n",
            col = colvec,
            pch = 21,
            pt.bg = colvec))
title("webbuilder networks RDA")



# MDS ####
# make data into distance matrix
mat.list <- llply(data.list, function (x){
  as.matrix(dist(x[c(2:7, 11:16, 19)]))
})
# make ordinations based on method
mds <- llply(mat.list,
             function (x){
               wcmdscale(x, k = 2,
                         eig = T,
                         w = rep(1, nrow(x)),
                         x.ret = TRUE)
             })
# observed ordination
plot(mds$observed, type = "n")
with(data.list$observed,
     points(mds$observed$points,
            col = colvec[cat],
            pch = 21,
            bg = colvec[cat]))
with(data.list$observed,
     legend("bottomright",
            legend = levels(cat),
            bty = "n",
            col = colvec,
            pch = 21,
            pt.bg = colvec))
title("Observed networks MDS")

# gravel ordination
plot(mds$gravel, type = "n")
with(data.list$gravel,
     points(mds$gravel$points,
            col = colvec[cat],
            pch = 21,
            bg = colvec[cat]))
with(data.list$gravel,
     legend("bottomleft",
            legend = levels(cat),
            bty = "n",
            col = colvec,
            pch = 21,
            pt.bg = colvec))
title("gravel networks MDS")

# webbuilder ordination
plot(mds$webbuilder, type = "n")
with(data.list$webbuilder,
     points(mds$webbuilder$points,
            col = colvec[cat],
            pch = 21,
            bg = colvec[cat]))
with(data.list$webbuilder,
     legend("bottomleft",
            legend = levels(cat),
            bty = "n",
            col = colvec,
            pch = 21,
            pt.bg = colvec))
title("webbuilder networks MDS")


# Gravel 3 steps ####
# looking at other gravel steps
gravel.stats <- all.stats[-52:-85,]
gravel.stats <- left_join(gravel.stats, sites, by = ".id")

# # clean up for better ggplotting
# # make the model variable as a factor
# gravel.stats$model <- as.factor(gravel.stats$model)
# # reorder factor for plots
# gravel.stats$model <- factor(
#   gravel.stats$model,
#   levels(gravel.stats$model)[c(2,1,3)])

# ANOVA ####
# one way anova testing different methods
# for differences among land use types
# using aov() in order to perform 
# posthoc Tukey test

# function for getting one way aov
one_way_aov <- function(d, x, y){
  aov(d[[y]] ~ d[[x]])
}


l <- dlply(gravel.stats,
           .(step),
           one_way_aov,
           x = "cat", y = "L")

ls <-  dlply(gravel.stats,
             .(step),
             one_way_aov,
             x = "cat", y = "LS")
c <-  dlply(gravel.stats,
            .(step),
            one_way_aov,
            x = "cat", y = "C")
# mcl <-  dlply(gravel.stats,
#               .(step),
#               one_way_aov,
#               x = "cat", y = "mean.TL")
b <-  dlply(gravel.stats,
            .(step),
            one_way_aov,
            x = "cat", y = "B")
int <-  dlply(gravel.stats,
              .(step),
              one_way_aov,
              x = "cat", y = "I")
top <-  dlply(gravel.stats,
              .(step),
              one_way_aov,
              x = "cat", y = "T")
g <-  dlply(gravel.stats,
            .(step),
            one_way_aov,
            x = "cat", y = "Gensd")
v <- dlply(gravel.stats,
           .(step),
           one_way_aov,
           x = "cat", y = "Vulsd")
msim <-  dlply(gravel.stats,
               .(step),
               one_way_aov,
               x = "cat", y = "Maxsim")
# maxtl <-  dlply(gravel.stats,
#                 .(step),
#                 one_way_aov,
#                 x = "cat", y = "max.TL")

aov.results2 <- list(l, ls, c, b, int,
                    top, g, v, msim)

# length(aov.res...) == 9
# length each aov.res element == 3
# variables == 9
# each variable has three steps
aov.table2 <- matrix(0, 9, 3)
for (i in 1:length(aov.results2)){
  for (j in 1:length(aov.results2[[1]])){
    aov.table2[i,j] <- anova(aov.results2[[i]][[j]])[1,5]
  }
}
dimnames(aov.table2) <- list(
  c("L", "LS", "C", "B", "I", "T", "Gensd", "Vulsd", "Maxsim"),
  c("1", "2", "3"))
aov.table2 <- signif(aov.table2, 4)
# significant results:
# 1:                  Vul
# 2:        B, I, Gen, Vul,
# 3: LS, C, B,    Gen, Vul
