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
# calculate dry weight of fish based on length measure
# Jellyman et al. 2013
get_fish_dw <- function(l, a, b){
  x = log(a) + (b * log(l))
  dw = exp(1)^x
  dw
}

#downloaded fish data
fish <- read_csv("data/Taieri NZ Fish datbase.csv")

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
# data frame of sites with fish
fish.sites <- data.frame(site =c("Berwick",
          "Blackrock", "Broad", "Canton", "Dempsters",
          "German", "Kyeburn", "Little", "NorthCol",
          "Powder", "Venlaw", "Sutton","Catlins",
          "Dempsters", "German", "Healy", "Little",
          "Narrowdale", "Stony", "Venlaw", "Dempsters",
          "Little", "Dempsters"),
             taxa = c(rep("Salmo",12),
                    rep("Galaxias", 8),
                    rep("Anguilla", 2),
                    rep("Gobiomorphus", 1)),
             stringsAsFactors = F)

# make an empty list for dry weight calculations
# estimate fish abundance using metabolic scaling
# N ~ M^b
# b generally in the range [-0.675 - -1.2] when using individual data
# estimated exponent b (linear slope of log10 transformed data) as ~-0.33 (see script "supplemental fish calculations.R")
dw.list <- NULL
for(i in 2:length(fish.list)){
  dw.list[[i]] = data.frame(taxa = c("Anguilla", "Galaxias", "Gobiomorphus", "Salmo"),
# calculate body mass 
          dw = get_fish_dw(fish.list[[i]], a, b),
      stringsAsFactors = F)
# estimate abundances
  dw.list[[i]] <- dw.list[[i]] %>%
    mutate(no.m2 = dw^-1) %>%
    right_join(fish.sites, by = "taxa")
}
names(dw.list) <- names(fish.list)
dw.list[[1]] <- NULL

saveRDS(dw.list, "estimated fish bodymass and abundance.RDS")

# scratch calcs for seeing if exponent changes estimated abundances 

mean.min <- fish.list$mean.min
get_fish_dw(mean.min, a, b)
fish.ab <- data.frame(taxa = c("Anguilla", "Galaxias",
                               "Gobiomorphus", "Salmo"),
                      dw = get_fish_dw(mean.min, a, b))

fish.ab <- fish.ab %>%
  mutate(no.05 = round(dw^-0.5, 4),
         no.06 = round(dw^-0.6, 4) ,
         no.07 = round(dw^-0.7, 4),
         no.08 = round(dw^-0.8, 4),
         no.09 = round(dw^-0.9, 4),
         no.1 = round(dw^-1, 4),
         no.11 = round(dw^-1.1, 4),
         no.12 = round(dw^-1.2, 4)) 

