# Translate Taieri_community_data.csv

# Taieri_community_data.csv has Ross's old taxa codes in it
# converting to genera names
# add formulas for biomass estimation
# estimate sp average biomass

# ***** if you want individual dry weights, need to convert files in
# C:\Users\Justin\Documents\Data\Taieri food web data\Body lengths of taxa by food web
# to .csv's with individual lengths replace "Taieri_community_data.csv"
# with list of files, then repeat rest of process. 
# see "biomass taxa names.R" in gravel project for more details


#recoderFunc from S Johnston
recoderFunc <- function(data, oldvalue, newvalue) {
  # convert any factors to characters
  if (is.factor(data))     data     <- as.character(data)
  if (is.factor(oldvalue)) oldvalue <- as.character(oldvalue)
  if (is.factor(newvalue)) newvalue <- as.character(newvalue)
  # create the return vector
  newvec <- data
  # put recoded values into the correct position in the return vector
  for (i in unique(oldvalue)) newvec[data == i] <- newvalue[oldvalue == i]
  newvec
}

#translation file
# fix names, typos, misnomers, etc. 
translate <- read.csv("translation.csv")

# read in Taieri csv
taieri.comm <- read_csv("Taieri_community_data.csv")
taieri.comm$avg.mm.BL <- taieri.comm$avg.mm.BL %>%
  as.numeric()

#replacing " " with a "." to match old tranlsation file which read spaces in as a period
taieri.comm$taxa <- taieri.comm %$%
  taxa %>%
  gsub(" ", "\\.", .)

# double check that all names are in translate file
# setdiff(taieri.comm$taxa, translate$Wrong) %>% sort()

# re code taxa to Corrected
taieri.comm$taxa <- recoderFunc(taieri.comm$taxa, translate$Wrong, translate$Corrected)

# species to genus translation
sp.gen <- read_csv("species genus category ffg.csv")

# re code taxa to genus
taieri.comm$taxa <- recoderFunc(taieri.comm$taxa, sp.gen$Species, sp.gen$Genus)

# combine genus names that are now duplicated
taieri.comm <- taieri.comm %>% 
  group_by(site, taxa) %>% 
  summarize(no.m2 = sum(no.m2),
            avg.mm.bl = mean(avg.mm.BL, na.rm = T))

# biomass formula ####
# read in formula
# file from Helen, modified containing all variable values
formula <- read_csv("C:\\Users\\Justin\\Documents\\Data\\Length DW conversion\\FWsurvey\\invert_meas_CJP\\length_weight_formulas.csv") 
formula <- formula %>%.[,c(4,5,6,7)] %>% distinct

# merge tables
taieri.length <- left_join(taieri.comm, formula, by = c("taxa" = "Name"))


# estimate sp average dw
# dw is in grams
taieri.dw <- taieri.length %>% 
  mutate(log = ln_a + (b * log(avg.mm.bl, base = base)), 
         dw = (base^log)/1000,
         logdw = log10(dw))

taieri <- taieri.dw %>% 
  select(site, taxa, no.m2, avg.mm.bl, logdw, dw)


# add fish biomass and abundance
# calculated mean min fish size from Taieri database  
# in 2. fish data.R
# dw is in grams
# abundance estimated using metabolic scaling theory
# c.f. Tang et al. 2014 ecology letters
salmo <- data.frame(site =c("Berwick", "Blackrock", "Broad", "Canton", "Dempsters", "German", "Kyeburn", "Little", "NorthCol", "Powder", "Venlaw", "Sutton"), 
                    taxa = "Salmo", 
                    no.m2 = 0.06311, 
                    avg.mm.bl = 73.85357, 
                    logdw = 0.0592191, 
                    dw = 1.145851)

galax <- data.frame(site = c("Catlins", "Dempsters", "German", "Healy", "Little", "Narrowdale", "Stony", "Venlaw"), 
                    taxa = "Galaxias", 
                    no.m2 = 0.15451, 
                    avg.mm.bl = 53.28232, 
                    logdw = -0.5169533, 
                    dw = 0.3041212)

eel <- data.frame(site = c("Dempsters", "Little"), 
                  taxa = "Anguilla", 
                  no.m2 = 0.00505, 
                  avg.mm.bl = 417.375, 
                  logdw = 1.68435, 
                  dw = 48.34485)

gobio <- data.frame(site = "Dempsters", 
                    taxa = "Gobiomorphus", 
                    no.m2 = 0.15347, 
                    avg.mm.bl = 50.333, 
                    logdw = -0.5126508, 
                    dw = 0.3071564)

# bind fish data rows to taieri
taieri <- bind_rows(taieri, salmo, galax, gobio, eel)

# save RDS
saveRDS(taieri, file = "estimated dw taieri webs.rds")


