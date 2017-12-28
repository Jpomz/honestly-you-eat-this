#Fish data
#J Pomz
# 31 Jan 2017

#formula for converting length to DW
#from Jellyman et al 2013 NZJMFR
#see citation for coefficients
#log10(w) = log(a) + b * log(l)

#downloaded all fish data from Taireri Catchment
#downloaded from https://www.niwa.co.nz/our-services/online-services/freshwater-fish-database
#~3200 rows, 22 columns

fish <- read.csv("Taieri NZ Fish datbase.csv")

fish %>%  
  filter(minl >0 & maxl >0) %>%
  group_by(locality, spcode) %>% 
  summarize(mean.min = mean(minl),
            mean.max = mean(maxl)) %>%
  select(locality, spcode, mean.min,mean.max) %>% arrange(locality) #%>% write_csv("taieri mean min max fish.csv")



fish %>%  
  filter(minl >0 & maxl >0) %>%
  group_by( spcode) %>% 
  summarize(mean.min = mean(minl),
            mean.max = mean(maxl)) %>%
  select( spcode, mean.min,mean.max) %>%
  as.data.frame



# salmo trutta ####
# a = 2.294e-6
# ln(a) = -12.9852
#mean max saltru = 171.6028
#mean min saltru = 73.85357
-12.985 + (3.05 * log(73.85357, base = exp(1)))
#log(dw) = 0.1363573
exp(1)^0.1363573
# dw = 1.146091
log10(1.146091)
# 0.0592191




-12.985 + 3.05 * log(130)
#log(dw) = 1.86098
exp(1)^1.86098
# dw = 6.430035
log10(6.430035)
#0.8082133

# galax
# mean.min = 53.28232
-14.77496979 + (3.417 * log(53.28232,	base = exp(1)))
# log dw = -1.190329
exp(1)^-1.190329
# dw = 0.3041212
log10(0.3041212)
# -0.5169533
                
# gobio 
# mean length 50.333
# G. brevicepps, jellyman 2013
# a = 2.155e-07
# b = 3.616
log(2.155e-07)
# ln(a) = -15.3503
-15.3503 + 3.616 * log(50.333, base = exp(1))
# log(dw) = -1.180422
exp(1)^-1.180422
# dw = 0.3071491
log10(0.3071491)
# log10dw = -0.5126508

# Anguilla
# mean length 417.375
# A dieffenbachii, jellyman 2013
# a = 1.390e-8
# b = 3.641
log(1.390e-8)
# ln(a) = -18.09138
-18.09138 + 3.641 * log(417.375)
# log(dw) = 3.87836
exp(1)^3.87836
# dw = 48.34486
log10(48.34486)
# log10dw = 1.68435

# estimating abundance using metabolic scaling
# Tang et al. 2014
# xstar = 10^(x0 + 3*gamma + eps)*mi^(1+gamma)
# mi needs to be in kg
# gamma = -0.675
# x0 = -1.16
# ignoring eps for now

fish.dw <- data.frame(species = c("Salmo",
                      "Galaxias","Anguilla",
                      "Gobiomorphus"),
                      # dw calc above
                      dw = c(1.145851, 0.3041212,
                             48.34485, 0.3071564))

fish.dw <- fish.dw %>% mutate(dw.kg = dw / 10^3)

# biomass at equilibrium
get_xstar <- function(mi, x0 = -1.16, gamma = -0.675){
  xstar = 10^(x0 + 3 * gamma)*mi^(1 + gamma)
  return(xstar)
}

fish.dw <- fish.dw %>%
  mutate(xstar = get_xstar(mi = .$dw.kg),
         # numerical abundance = eq biomass / body mass
         no.m2 = xstar / dw.kg) 
