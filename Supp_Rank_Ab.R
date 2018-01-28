library(plyr)
library(tidyverse)
# invertebrate biomass and abundance
#invert <- readRDS("estimated invert bodymass.RDS")

invert <- readRDS("ab, dw, info for sites, subset to match inference.RDS")
invert <- ldply(invert)

invert %>% group_by(site) %>%
  mutate(tot = sum(no.m2),
         a = no.m2 / tot) %>%
  group_by(taxa) %>%
  summarize(ma = mean(a), sd(a)) %>%
  arrange(desc(ma)) %>%
  View
  
  
invert %>% group_by(site) %>%
  mutate(tot = sum(no.m2),
         a = no.m2 / tot) %>%
  ggplot(aes(x = taxa, y = a))+
  geom_boxplot()


invert %>% group_by(site) %>%
  mutate(tot = sum(no.m2),
         a = no.m2 / tot) %>%
  group_by(taxa) %>%
  summarize(ma = max(a), sd = sd(a)) %>%
  ggplot(aes(x = ma, y = (sd))) +
  geom_point() +
  stat_smooth()


invert %>% group_by(site) %>%
  mutate(tot = sum(no.m2),
         a = no.m2 / tot) %>%
  group_by(taxa) %>%
  summarize(ma = max(a), sd = sd(a)) %>%
  mutate(m.s = ma / sd) %>%
  ggplot(aes(ma, m.s)) +
  geom_point()


invert %>% group_by(site) %>%
  mutate(tot = sum(no.m2),
         a = no.m2 / tot) %>%
  mutate(rank = rank(-no.m2, ties.method = "first")) %>% #View
  group_by(taxa) %>%
  summarize(n(), mean = mean(rank), min = min(rank), max = max(rank)) %>%
  arrange(mean) %>% View


invert[order(desc(rank(invert$no.m2, ties = "first"))),]
invert[312,]


invert %>% 
  mutate(tot = sum(no.m2),
          rel.ab = no.m2 / tot) %>%
  group_by(taxa) %>%
  summarize(mean = mean(rel.ab), sd = sd(rel.ab), n = n()) %>% 
  arrange(desc(mean))





invert %>% group_by(site) %>%
  mutate(tot = sum(no.m2),
         a = no.m2 / tot) %>%
  mutate(rank = rank(-no.m2, ties.method = "average")) %>%
  filter(rank <20) %>%
  ggplot(aes(rank))+
  facet_wrap(~taxa) +
  geom_histogram()

invert %>% group_by(site) %>%
  mutate(tot = sum(no.m2),
         a = no.m2 / tot) %>%
  mutate(rank = rank(-no.m2, ties.method = "average")) %>%
  filter(rank <8) %>%
  ggplot(aes(rank, fill = taxa))+
  geom_histogram(position = "fill", bins = 8)

invert %>% group_by(site) %>%
  mutate(tot = sum(no.m2),
         a = no.m2 / tot) %>%
  mutate(rank = rank(-no.m2, ties.method = "average")) %>%
  group_by(taxa) %>%
  count(rank) %>%
  arrange(rank, taxa) %>% View



# make a graph of M ~ mean(M)

invert %>% group_by(taxa) %>%
  summarize(mean = mean(dw), sd = sd(dw)) %>%
  na.omit %>%
  ggplot(aes(log10(mean), log10(mean)))+
  geom_point()+
  geom_errorbar(aes(ymin = log10(mean) - log10(sd),
                    ymax = log10(mean) + log10(sd)))


# local site variation from mean biomass
invert %>%
  group_by(taxa) %>%
  mutate(mean.dw = mean(dw)) %>%
  #group_by(site, taxa) %>%
  mutate(x = (dw - mean.dw)^2) %>%
  group_by(taxa) %>%
  summarize(max.x = log10(max(x))) %>%
  arrange(desc(max.x))

  
invert %>%
    group_by(taxa) %>%
    mutate(mean.no.m2 = mean(log10(no.m2))) %>%
    #group_by(site, taxa) %>%
    mutate(x = (no.m2 - mean.no.m2)^2) %>%
    group_by(taxa) %>%
  summarize(max.x = log10(max(x))) %>%
  arrange(desc(max.x))



# boxplots ####
# filtering out fish (same values at all sites)
# reducing to only taxa that were observed >= 7 sites
# dw
dplot <- invert %>% filter(dw<.10) %>%
  group_by(taxa) %>%
  mutate(n = n()) %>% 
  filter(n > 7) %>%
  ggplot(aes(fct_reorder(taxa, no.m2), log10(dw))) +
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1))
# N
nplot <- invert %>% filter(dw<.10) %>%
  group_by(taxa) %>%
  mutate(n = n()) %>% 
  filter(n > 7) %>%
  ggplot(aes(fct_reorder(taxa, no.m2), log10(no.m2))) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90,
                                 hjust = 1))

library(gridExtra)
grid.arrange(dplot, nplot)

# geom_density
invert %>% filter(dw<.10) %>%
  group_by(taxa) %>%
  mutate(n = n()) %>% 
  filter(n > 5) %>%
  ggplot(aes(log10(dw))) +
  geom_density()+
  facet_wrap(~taxa)

invert %>% filter(dw<.10) %>%
  group_by(taxa) %>%
  mutate(n = n()) %>% 
  filter(n > 5) %>%
  ggplot(aes(log10(no.m2))) +
  geom_density()+
  facet_wrap(~taxa)

invert %>% filter(dw<.10) %>%
  group_by(taxa) %>%
  #mutate(n = n()) %>% 
  filter(log10(dw) < 5) %>%
  ggplot(aes(log10(dw))) +
  geom_density()+
  facet_wrap(~taxa)

invert %>% filter(dw<.10) %>%
  group_by(site) %>%
  mutate(tot = sum(no.m2),
         a = no.m2 / tot) %>%
  mutate(rank = rank(-no.m2, ties.method = "first")) %>%
  group_by(taxa) %>%
  mutate(m.rank = mean(rank)) %>%
  filter(m.rank < 10) %>%
  ggplot(aes(log10(no.m2))) +
  geom_density()+
  facet_wrap(~taxa)

