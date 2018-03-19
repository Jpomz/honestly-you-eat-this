# consumption rate
# 25 Jan 2018
#jfpomeranz.@gmail.com

# consumption rate scales with bodymass as c ~ xM**0.78
# where x is a constant and M is consumer body mass
# from Pawar et al 2012 Nature doi:10.1038/nature11131


library(plyr)
library(dplyr)
library(ggplot2)
invert <- readRDS("ab, dw, info for sites, subset to match inference.RDS")
invert <- ldply(invert)


invert %>% transmute(log10(dw**0.78)) %>%
  range()
#~ 4 orders of magnitude