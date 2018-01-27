library(plyr)
library(dplyr)
library(ggplot2)
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
