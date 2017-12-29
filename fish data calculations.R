# Fish length, biomass, abundance estimates
# 29 Dec 2017
# Justin Pomeranz
# jfpomeranz@gmail.com

# fish data in originals only includes presence / absence
#downloaded all fish data from Taireri Catchment
# downloaded January 2017 from https://www.niwa.co.nz/our-services/online-services/freshwater-fish-database
# "Taieri NZ Fish datbase.csv" ~3200 rows, 22 columns

# convert fish lengths to biomass estimates 
# Formulas from Jellyman et al. NZJMFR "Does one size fit all? An evaluation of length–weight relationships for New Zealand's freshwater fish species"

# estimate numerical abundance using metabolic theory based on body size
# formula in Tang et al. 2014 Ecology Letters "Correlation between interaction strengths ..." Supporting Information, section 2

# libraries 
library(tidyverse)
# functions
#log10(w) = log(a) + b * log(l)
get_fish_dw <- function(l, a, b){
  x = log(a) + (b * log(l))
  dw = exp(1)^x
  dw
}


#downloaded fish data
fish <- read_csv("Taieri NZ Fish datbase.csv")

# vector of spcodes to keep
# species codes are in "NZFFD species codes.csv"
# No observations in origingal data:
# "angaus", "galarg", "galeld", "galfas", "galmac", "galpul", "gobcot", "gobiom",
fishcode <- c("angdie", "galano",  "galaxi", "galbre", "galdep", "gobbre", "saltru")

# filter out "0" observations and other spcodes
# only keep spcode and min/max length observations
fish <- fish %>%
  filter(minl >0 & maxl >0, spcode %in% fishcode) %>%
  select(spcode, minl, maxl) 
# only keep first three characters of spcode to get genus average
fish$spcode <- substr(fish$spcode, 1, 3)

fish.summary <-  fish %>%
  group_by(spcode) %>%
  summarize(min.min = min(minl),
            mean.min = mean(minl),
            mean.max = mean(maxl),
            max.max = max(maxl))

# formulas used for length weight regression
# Jellyman et al. 2013 table 2, dry weight equations
# Anguilla dieffenbachii
# Galaxias brevipinnis
# Gobiomorphus breviceps
# Salmo trutta
# vector of coefficients for length - weight regression
a <- c(1.39e-8, 1.4e-8, 2.155e-7, 2.294e-6)
b <- c(3.641, 4.082, 3.616, 3.05)
  
# make a list from fish summary columns
fish.list <- as.list(fish.summary)

# make an empty list for dry weight calculations
dw.list <- NULL
for(i in 2:length(fish.list)){
  dw.list[[i]] = data.frame(taxa = fish.list[[1]],
    dw = get_fish_dw(fish.list[[i]], a, b))

}
names(dw.list) <- names(fish.list)
dw.list[[1]] <- NULL

