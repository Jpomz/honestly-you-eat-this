# 22 WebBuilder && Trait-matching
# JPomz
# 3 dec 2017

# multiplying webbuilder inferred and trait-matching inferred adjacency matrices together 

source("C:\\Users\\Justin\\Documents\\Data\\FW modelling Petchey Github\\ttl-resources-master\\food_web\\FoodWebFunctions.r")
# TSS function
source("TSS function.R")
source("adj_conf_matrix function.R")


obs <- readRDS("observed webs names and taxa matched.rds")
grav.in <- readRDS("gravel inferred initial.rds")
grav.niche <- readRDS("gravel inferred pruned.rds")
grav.neutral <- readRDS("gravel inferred pruned rel ab.rds")
wb <- readRDS("webbuilder inferred list trimmed taxa.rds")

# grav initial * wb
grav.in.wb <- map2(grav.in, wb, ~.x*.y)
# grav niche * wb
grav.niche.wb <- map2(grav.niche, wb, ~.x*.y)
# grav neutral * wb
grav.neutral.wb <- map2(grav.neutral, wb, ~.x*.y)




grav.wb.tss <- map2(obs,grav.in.wb, tss)
ldply(grav.wb.tss) %>% summarize(tss.mean = mean(tss),
                                 tss.sd = sd(tss))

# gravel initial
g1 <- map2(obs,grav.in, tss) %>%
  ldply() %>%
  summarize(tss.mean = mean(tss),
            tss.sd = sd(tss))
# graven niche pruned
g2 <- map2(obs, grav.niche, tss) %>%
  ldply() %>%
  summarize(tss.mean = mean(tss),
            tss.sd = sd(tss))
# gravel neutral
g3 <- map2(obs, grav.neutral, tss) %>%
  ldply() %>%
  summarize(tss.mean = mean(tss),
            tss.sd = sd(tss))
# wb
wb1 <- map2(obs, wb, tss) %>%
  ldply() %>%
  summarize(tss.mean = mean(tss),
            tss.sd = sd(tss))
# grav.in * wb
g1wb <- map2(obs, grav.in.wb, tss) %>%
  ldply() %>%
  summarize(tss.mean = mean(tss),
            tss.sd = sd(tss))
# grav.niche * wb
g2wb <- map2(obs, grav.niche.wb, tss) %>%
  ldply() %>%
  summarize(tss.mean = mean(tss),
            tss.sd = sd(tss))
# grav.neutral * wb
g3wb <- map2(obs, grav.neutral.wb, tss) %>%
  ldply() %>%
  summarize(tss.mean = mean(tss),
            tss.sd = sd(tss))
rbind(g1, g2, g3, wb1, g1wb, g2wb, g3wb)


# confusion matrix
# gravel initial
g1.acc <- map2(obs,grav.in, adj_conf_matrix) %>%
  ldply() %>%
  summarize(acc.mean = mean(ACC),
            tpr.mean = mean(TPR),
            fpr.mean = mean(FPR),
            fnr.mean = mean(FNR))
# graven niche pruned
g2.acc <- map2(obs,grav.niche, adj_conf_matrix) %>%
  ldply() %>%
  summarize(acc.mean = mean(ACC),
            tpr.mean = mean(TPR),
            fpr.mean = mean(FPR),
            fnr.mean = mean(FNR))
# gravel neutral
g3.acc <- map2(obs,grav.neutral, adj_conf_matrix) %>%
  ldply() %>%
  summarize(acc.mean = mean(ACC),
            tpr.mean = mean(TPR),
            fpr.mean = mean(FPR),
            fnr.mean = mean(FNR))
# wb
wb.acc <- map2(obs,wb, adj_conf_matrix) %>%
  ldply() %>%
  summarize(acc.mean = mean(ACC),
            tpr.mean = mean(TPR),
            fpr.mean = mean(FPR),
            fnr.mean = mean(FNR))
# grav.in * wb
g1.wb.acc <- map2(obs,grav.in.wb,
                  adj_conf_matrix) %>%
  ldply() %>%
  summarize(acc.mean = mean(ACC),
            tpr.mean = mean(TPR),
            fpr.mean = mean(FPR),
            fnr.mean = mean(FNR))
# grav.niche * wb
g2.wb.acc <- map2(obs,grav.niche.wb,
                  adj_conf_matrix) %>%
  ldply() %>%
  summarize(acc.mean = mean(ACC),
            tpr.mean = mean(TPR),
            fpr.mean = mean(FPR),
            fnr.mean = mean(FNR))
# grav.neutral * wb
g3.wb.acc <- map2(obs,grav.neutral.wb,
                  adj_conf_matrix) %>%
  ldply() %>%
  summarize(acc.mean = mean(ACC),
            tpr.mean = mean(TPR),
            fpr.mean = mean(FPR),
            fnr.mean = mean(FNR))
rbind(g1.acc, g2.acc, g3.acc, wb.acc, g1.wb.acc, g2.wb.acc, g3.wb.acc)
