# does fish biomass affect gravel parameters


# 1) select range of fish sizes
# 2) estimate dry weight for range of sizes
# 3) incorporate fish size ranges into Taieri biomass estimates
# 4) parameterize gravel model with different fish estimates
# 5) copmare parameters to see if they change

# gravel functions
source("gravel_functions.R")
# useful food web functions from Petchey
source("C:\\Users\\Justin\\Documents\\Data\\FW modelling Petchey Github\\ttl-resources-master\\food_web\\FoodWebFunctions.r")
# TSS function
source("TSS function.R")



#downloaded all fish data from Taireri Catchment
#downloaded from https://www.niwa.co.nz/our-services/online-services/freshwater-fish-database
#~3200 rows, 22 columns

fish <- read.csv("Taieri NZ Fish datbase.csv")

# fish %>%  
#   filter(minl >0 & maxl >0) %>%
#   group_by(locality, spcode) %>% 
#   summarize(mean.min = mean(minl),
#             mean.max = mean(maxl)) %>%
#   select(locality, spcode, mean.min,mean.max) %>% arrange(locality) #%>% write_csv("taieri mean min max fish.csv")
# 
# 
# 
# fish %>%  
#   filter(minl >0 & maxl >0) %>%
#   group_by( spcode) %>% 
#   summarize(mean.min = mean(minl),
#             mean.max = mean(maxl),
#             min.min = min(minl),
#             max.max = max(maxl)) %>%
#   select( spcode, mean.min,mean.max, min.min, max.max) %>%
#   as.data.frame

sp.code <- fish$spcode %>% as.character() %>% unique %>% sort()
genus.number <- c(2:4, 6:15, 17:20, 33:34)

sp.subset <- sp.code[genus.number]
genus <- c("Ang", "Ang","Ang", "galax", "galax", "galax", "galax", "galax", "galax", "galax", "galax", "galax", "galax", "gobio", "gobio", "gobio", "gobio", "salmo", "salmo")


sp.genus <- data.frame(spcode = sp.subset, genus)

fish <- left_join(fish, sp.genus, by = "spcode")

fish.l <- fish %>%  
  filter(minl >0 & maxl >0) %>%
  group_by(genus) %>% 
  summarize(mean.min = mean(minl),
            mean.max = mean(maxl),
            min.min = min(minl),
            max.max = max(maxl)) %>%
  select(genus, mean.min,mean.max, min.min, max.max) %>%
  as.data.frame()


# get dw funtions ####
#formula for converting length to DW
#from Jellyman et al 2013 NZJMFR
#see citation for coefficients
#log10(w) = log(a) + b * log(l)
get.salmo.dw <- function(x, lna = -12.98521, b = 3.05, base = exp(1)){
  logdw = lna + (b * log(x, base = base))
  dw = base^logdw
  return(dw)
}

get.galax.dw <- function(x, lna = -14.77496979, b = 3.417, base = exp(1)){
  logdw = lna + (b * log(x, base = base))
  dw = base^logdw
  return(dw)
}

get.gobio.dw <- function(x, lna = -15.3503, b = 3.616, base = exp(1)){
  logdw = lna + (b * log(x, base = base))
  dw = base^logdw
  return(dw)
}

get.angu.dw <- function(x, lna = -18.09138, b = 3.641, base = exp(1)){
  logdw = lna + (b * log(x, base = base))
  dw = base^logdw
  return(dw)
}

# make vectors of estimated dw's
angu.dw <- get.angu.dw(fish.l[1,-1])
galax.dw <- get.galax.dw(fish.l[2,-1])
gobio.dw <- get.gobio.dw(fish.l[3,-1])
salmo.dw <- get.salmo.dw(fish.l[4,-1])


# make lists with different fish dw estimates
angu.list <- NULL
for (i in 1:length(angu.dw)){
  angu.list[[i]] <- data.frame(site = c("Dempsters", "Little"), 
                        taxa = "Anguilla", 
                        no.m2 = NA, 
                        avg.mm.bl = fish.l[1,i+1], 
                        logdw = NA, 
                        dw = angu.dw[i])
}
galax.list <- NULL
for (i in 1:length(galax.dw)){
  galax.list[[i]] <- data.frame(site = c("Catlins", "Dempsters", "German", "Healy", "Little", "Narrowdale", "Stony", "Venlaw"),
                                taxa = "Galaxias", 
                               no.m2 = NA, 
                               avg.mm.bl = fish.l[2, i+1], 
                               logdw = NA, 
                               dw = galax.dw[i])
}
gobio.list <- NULL
for (i in 1:length(gobio.dw)){
  gobio.list[[i]] <- data.frame(site = "Dempsters",
                                taxa = "Gobiomorphus", 
                                no.m2 = NA, 
                                avg.mm.bl = fish.l[3, i+1], 
                                logdw = NA, 
                                dw = gobio.dw[i])
}
salmo.list <- NULL
for (i in 1:length(salmo.dw)){
  salmo.list[[i]] <- data.frame(site = c("Berwick", "Blackrock", "Broad", "Canton", "Dempsters", "German", "Kyeburn", "Little", "NorthCol", "Powder", "Venlaw", "Sutton"), 
                                taxa = "Salmo",
                                no.m2 = NA, 
                                avg.mm.bl = fish.l[4, i+1], 
                                logdw = NA, 
                                dw = salmo.dw[i])
}


# code from 3. estimating dw ...R
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
  summarize(no.m2 = sum(no.m2), avg.mm.bl = mean(avg.mm.BL, na.rm = T))

# biomass formula ####
# read in formula
# file from Helen, modified containing all variable values
formula <- read_csv("C:\\Users\\Justin\\Documents\\Data\\Length DW conversion\\FWsurvey\\invert_meas_CJP\\length_weight_formulas.csv") 
formula <- formula %>%.[,c(4,5,6,7)] %>% distinct

# merge tables
taieri.length <- left_join(taieri.comm, formula, by = c("taxa" = "Name"))


# estimate sp average dw
taieri.dw <- taieri.length %>% 
  mutate(log = ln_a + (b * log(avg.mm.bl, base = base)), 
         dw = (base^log)/1000,
         logdw = log10(dw))

taieri <- taieri.dw %>% 
  select(site, taxa, no.m2, avg.mm.bl, logdw, dw)

# make list (length = 4) with all taieri comm data and 4 diff fish estimates
taieri.fish.list <- NULL
for (i in 1:length(salmo.list)){
  taieri.fish.list[[i]] <- bind_rows(taieri,
                                     salmo.list[[i]],
                                     galax.list[[i]],
                                     gobio.list[[i]],
                                     angu.list[[i]])
}

# get all body sizes present at sites
# remove NAs
taieri.fish.list <- taieri.fish.list %>%
  llply(function (x){
    x %>% filter(!is.na(dw))
  }
)
# sort
taieri.fish.list <- taieri.fish.list %>%
  llply(function (x){
    x %>% arrange(site, logdw)
  }
)

# make into a list of lists
taieri.list.list <- NULL
for (i in 1:length(taieri.fish.list)){
  taieri.list.list[[i]] <- split(taieri.fish.list[[i]],
                                 list(taieri.fish.list[[i]]$site))
}

# read in corrected adjaceny matrices
web.list <- readRDS("17 observed taieri food webs.rds")


# convert predation matrix to matrix of pairs
pairs.list <- web.list %>%
  llply(function (x){
    as.data.frame(Matrix.to.list(x))})

# name columns resource and consumer
pairs.list <- llply(pairs.list,
                    function (x){
                      colnames(x) <- c("resource", "consumer");x})

# make list of lists of consumer - resource pairs with estimated dw
pairs.dw.list.list <- NULL
pairs.dw.list <- NULL
for( j in 1:length(taieri.list.list)){
  for (i in 1:length(pairs.list)){
  resource <- left_join( # add dw for resource column
    pairs.list[[i]], taieri.list.list[[j]][[i]][,c(2,6)],
    by = c("resource" = "taxa"))
  consumer <- left_join( # just consumer dw
    pairs.list[[i]], taieri.list.list[[j]][[i]][,c(2,6)],
    by = c("consumer" = "taxa"))[3]
  combined <- bind_cols(resource, consumer) # combine
  colnames(combined) <- c("resource", 
                          "consumer", "res.dw", "con.dw") # rename columns
  pairs.dw.list[[i]] <- combined # put in list
  rm("resource", "consumer", "combined") # remove objects
  }
  pairs.dw.list.list[[j]] <- pairs.dw.list
}


# make list of res-con biomass pairs for model parameterizations
# each element in list contains all pairs EXCEPT those of
# the web they will be used to parameterize. 
# e.g. element named "blackrock" will contain all other
# pairs data [-"blackrock"], and be used to paramterize 
# the model in order to predict the blackrock web 

training.list.list <- NULL
training.list <- NULL
for( j in 1:length(taieri.list.list)){
  for (i in 1:length(pairs.dw.list)){
  a <- ldply(pairs.dw.list.list[[j]][-i])
  training.list[[i]] <- a
  rm(a)
  }
  training.list.list[[j]] <- training.list
}

# select estimated dry weights
# filter out dryweight NA's
# log10 transform dry weights, test if res > con
# filter out res > con == TRUE
for (j in 1:length(training.list.list)){
  training.list.list[[j]] <- training.list.list[[j]] %>%
    llply(function (x){
      x %>% 
        select(res.dw, con.dw) %>%
        filter(!is.na(res.dw), !is.na(con.dw)) %>%
        mutate(logres = log10(res.dw), 
               logcon = log10(con.dw),
               res.greater = res.dw > con.dw) %>%
        filter(res.greater == FALSE)
    })  
}

names(training.list.list) <- colnames(fish.l)[-1]

# parameterize gravel model ####
# parametrize model for each food web
# make a list of regression coefficients
pars.list.list <- NULL
pars.list <- NULL
for (j in 1:length(training.list.list)){
  for (i in 1:length(training.list.list[[j]])){
      Bprey <- training.list.list[[j]][[i]]$logres
      Bpred <- training.list.list[[j]][[i]]$logcon
      pars.list[[i]] <- reg_fn(Bprey,Bpred,quartil = c(0.01,0.97))
  }
  pars.list.list[[j]] <- pars.list
}

param.coef.list <- NULL
for (i in 1:length(pars.list.list)){
  param.coef.list[[i]] <- pars.list.list[[i]] %>% 
  ldply(function (x){
    data.frame(B0hi = x[[1]][1],
               B0center = x[[2]][1],
               B0lo = x[[3]][1],
               B1hi = x[[1]][2],
               B1center = x[[2]][2],
               B1lo = x[[3]][2])
  })
}

names(param.coef.list) <- colnames(fish.l)[-1]

param.coef.df <- param.coef.list %>% ldply

# compare food webs using different fish values
# dempsters
# get parameters for dempsters for all fish val
dempster.pars <- NULL
for (i in 1:length(pars.list.list)){
  dempster.pars[[i]] <- get_pars_Niche(
    pars.list.list[[i]][[8]],
    Ball = log10(taieri.list.list[[i]][[8]]$dw))
}

dempster.links.inf <- llply(dempster.pars, function (x){
  L_fn(x[,"n"],x[,"c"],x[,"low"],x[,"high"])
})
names(dempster.links.inf) <- names(training.list.list)

# for (i in 1:length(dempster.links.inf)){
#   Plot.matrix(dempster.links.inf[[i]])
#   title(names(dempster.links.inf)[i])
# }

# black rock
# get parameters for dempsters for all fish val
blackrock.pars <- NULL
for (i in 1:length(pars.list.list)){
  blackrock.pars[[i]] <- get_pars_Niche(
    pars.list.list[[i]][[4]],
    Ball = log10(taieri.list.list[[i]][[4]]$dw))
}

blackrock.links.inf <- llply(blackrock.pars, function (x){
  L_fn(x[,"n"],x[,"c"],x[,"low"],x[,"high"])
})
names(blackrock.links.inf) <- names(training.list.list)

# for (i in 1:length(blackrock.links.inf)){
#   Plot.matrix(blackrock.links.inf[[i]])
#   title(names(blackrock.links.inf)[i])
# }


# min.min.fish ####
# get paramters using min.min fish
min.pars<- NULL
for (i in 1:length(pars.list.list[[3]])){
  min.pars[[i]] <- get_pars_Niche(
    pars.list.list[[3]][[i]],
    Ball = log10(taieri.list.list[[3]][[i]]$dw))
}
# add taxa names to min.pars
for (i in 1:length(min.pars)){
  min.pars[[i]] <- as.data.frame(min.pars[[i]])
  min.pars[[i]]$taxa <- taieri.list.list[[3]][[i]]$taxa
}

# infer links min.min
min.links.inf <- llply(min.pars, function (x){
  L_fn(x[,"n"],x[,"c"],x[,"low"],x[,"high"])
})

# add taxa names to inferred matrices
for (i in 1:length(min.links.inf)){
  dimnames(min.links.inf[[i]]) <- list(min.pars[[i]]$taxa,
                                       min.pars[[i]]$taxa)
}

# read in observed
obs.list <- web.list
# read in target vectors
target <- readRDS("target vector to match all matrices.RDS")
# subset / match observed to target
for (i in 1:length(obs.list)){
  obs.list[[i]] <- obs.list[[i]][
    match(target[[i]], rownames(obs.list[[i]])),
    match(target[[i]], colnames(obs.list[[i]]))
    ]}
# subset / match minlinks.inf to target
for (i in 1:length(min.links.inf)){
  min.links.inf[[i]] <- min.links.inf[[i]][
    match(target[[i]], rownames(min.links.inf[[i]])),
    match(target[[i]], colnames(min.links.inf[[i]]))
    ]}
# TSS min.min model ####


tss.df <- NULL
for (i in 1:length(obs.list)){
  j <- tss(obs.list[[i]], min.links.inf[[i]])
  j$web <- names(obs.list)[[i]]
  j$step <- "model"
  tss.df <- rbind(tss.df, j)
  rm(j)
}

tss.3.step <- readRDS(paste(getwd(), "/figs for MS/gravel tss 3 steps.rds", sep = ""))

data.frame(mean.max = tss.3.step %>% 
             filter(step == "model") %>%
             select(tss),
           min.min = tss.df$tss) %>%
  summarize_all(mean)


# make list of taxa names in inferred matrices 
inf.taxa <- llply(min.links.inf, 
                  function (x){
                    colnames(x)
                  })
# pruning forbidden links ####
# get list of all taxa in inferred webs
inf.taxa.unique <- inf.taxa %>%
  unlist %>%
  unique() %>%
  sort
inf.taxa.df <- inf.taxa.unique %>% as.data.frame()
colnames(inf.taxa.df) <- "taxa"

# read in data with FFG information
ffg <- read_csv("taxa.ffg.csv")

# join inferred taxa with FFG
inf.taxa.df <- left_join(inf.taxa.df, ffg, c("taxa" = "name"))

# make a list of "forbidden taxa"
# e.g. taxa that cannot eat prey due to mouthparts,
# e.g. scrapers or filter feeders
forbid.taxa <- inf.taxa.df %>%
  filter(FFG == "SC" | FFG == "G" | FFG == "FF") %>% select(taxa)
# convert to character vector
forbid.taxa <- forbid.taxa %>% as_vector()

# make new list of webs to prune
inf.forbid.list <- min.links.inf
# make a list of columns which match forbidden taxa
forbid.col.list <- inf.forbid.list %>% 
  llply(function (x){
    which(colnames(x) %in% forbid.taxa)
  })

# make forbidden taxa in columns = 0
for (i in 1:length(inf.forbid.list)){
  x <- inf.forbid.list[[i]]
  x[, forbid.col.list[[i]]] <- 0
  inf.forbid.list[[i]] <- x
}

#TSS min.min prune ####
tss.prune.df <- NULL
for (i in 1:length(obs.list)){
  j <- tss(obs.list[[i]], inf.forbid.list[[i]])
  j$web <- names(obs.list)[[i]]
  j$step <- "prune"
  tss.prune.df <- rbind(tss.prune.df, j)
  rm(j)
}

min.min.tss.2step <- bind_rows(tss.df, tss.prune.df)

data.frame(mean.max = tss.3.step %>% 
             filter(step == "prune") %>%
             select(tss),
           min.min = tss.prune.df$tss) %>%
  summarize_all(mean)





# mean.min fish####
# get paramters using mean.min fish
min.mean.pars<- NULL
for (i in 1:length(pars.list.list[[1]])){
  min.mean.pars[[i]] <- get_pars_Niche(
    pars.list.list[[1]][[i]],
    Ball = log10(taieri.list.list[[1]][[i]]$dw))
}
# add taxa names to web.pars
for (i in 1:length(min.mean.pars)){
  min.mean.pars[[i]] <- as.data.frame(min.mean.pars[[i]])
  min.mean.pars[[i]]$taxa <- taieri.list.list[[1]][[i]]$taxa
}


min.mean.links.inf <- llply(min.mean.pars, function (x){
  L_fn(x[,"n"],x[,"c"],x[,"low"],x[,"high"])
})

# add taxa names to predation matrices
for (i in 1:length(min.mean.links.inf)){
  dimnames(min.mean.links.inf[[i]]) <- list(
    min.mean.pars[[i]]$taxa,
    min.mean.pars[[i]]$taxa)
}

# subset / match minlinks.inf to target
for (i in 1:length(min.mean.links.inf)){
  min.mean.links.inf[[i]] <- min.mean.links.inf[[i]][
    match(target[[i]], rownames(min.mean.links.inf[[i]])),
    match(target[[i]], colnames(min.mean.links.inf[[i]]))
    ]}

tss.mean.min.df <- NULL
for (i in 1:length(obs.list)){
  j <- tss(obs.list[[i]], min.mean.links.inf[[i]])
  j$web <- names(obs.list)[[i]]
  j$step <- "model"
  tss.mean.min.df <- rbind(tss.mean.min.df, j)
  rm(j)
}

tss.3.step <- readRDS(paste(getwd(), "/figs for MS/gravel tss 3 steps.rds", sep = ""))
mean.max.tss.2step <- readRDS(paste(getwd(), "/figs for MS/gravel tss 2 steps.rds", sep = ""))

data.frame(mean.max = tss.3.step %>% 
             filter(step == "model") %>%
             select(tss),
           mean.min = tss.mean.min.df$tss) %>%
  summarize_all(mean)


# make new list of mean.min webs to prune
inf.forbid.mean.min.list <- min.mean.links.inf
# make a list of columns which match forbidden taxa
forbid.col.list <- inf.forbid.mean.min.list %>% 
  llply(function (x){
    which(colnames(x) %in% forbid.taxa)
  })

# make forbidden taxa in columns = 0
for (i in 1:length(inf.forbid.mean.min.list)){
  x <- inf.forbid.mean.min.list[[i]]
  x[, forbid.col.list[[i]]] <- 0
  inf.forbid.mean.min.list[[i]] <- x
}

#TSS mean.min prune ####
tss.mean.min.prune.df <- NULL
for (i in 1:length(obs.list)){
  j <- tss(obs.list[[i]], inf.forbid.mean.min.list[[i]])
  j$web <- names(obs.list)[[i]]
  j$step <- "prune"
  tss.mean.min.prune.df <- rbind(tss.mean.min.prune.df, j)
  rm(j)
}

mean.min.tss.2step <- bind_rows(tss.mean.min.df, tss.mean.min.prune.df)

data.frame(mean.max = tss.3.step %>% 
             filter(step == "prune") %>%
             select(tss),
           mean.min = tss.mean.min.prune.df$tss) %>%
  summarize_all(mean)


# max.max fish####
# get paramters using max.max fish
max.max.pars<- NULL
for (i in 1:length(pars.list.list[[4]])){
  max.max.pars[[i]] <- get_pars_Niche(
    pars.list.list[[4]][[i]],
    Ball = log10(taieri.list.list[[4]][[i]]$dw))
}
# add taxa names to web.pars
for (i in 1:length(max.max.pars)){
  max.max.pars[[i]] <- as.data.frame(max.max.pars[[i]])
  max.max.pars[[i]]$taxa <- taieri.list.list[[4]][[i]]$taxa
}


max.max.links.inf <- llply(max.max.pars, function (x){
  L_fn(x[,"n"],x[,"c"],x[,"low"],x[,"high"])
})

# add taxa names to predation matrices
for (i in 1:length(max.max.links.inf)){
  dimnames(max.max.links.inf[[i]]) <- list(
    max.max.pars[[i]]$taxa,
    max.max.pars[[i]]$taxa)
}

# subset / match minlinks.inf to target
for (i in 1:length(max.max.links.inf)){
  max.max.links.inf[[i]] <- max.max.links.inf[[i]][
    match(target[[i]], rownames(max.max.links.inf[[i]])),
    match(target[[i]], colnames(max.max.links.inf[[i]]))
    ]}

tss.max.max.df <- NULL
for (i in 1:length(obs.list)){
  j <- tss(obs.list[[i]], max.max.links.inf[[i]])
  j$web <- names(obs.list)[[i]]
  j$step <- "model"
  tss.max.max.df <- rbind(tss.max.max.df, j)
  rm(j)
}

data.frame(max.max = tss.3.step %>% 
             filter(step == "model") %>%
             select(tss),
           max.max = tss.max.max.df$tss) %>%
  summarize_all(mean)


# make new list of max.max webs to prune
inf.forbid.max.max.list <- max.max.links.inf
# make a list of columns which match forbidden taxa
forbid.col.list <- inf.forbid.max.max.list %>% 
  llply(function (x){
    which(colnames(x) %in% forbid.taxa)
  })

# make forbidden taxa in columns = 0
for (i in 1:length(inf.forbid.max.max.list)){
  x <- inf.forbid.max.max.list[[i]]
  x[, forbid.col.list[[i]]] <- 0
  inf.forbid.max.max.list[[i]] <- x
}

#TSS max.max prune ####
tss.max.max.prune.df <- NULL
for (i in 1:length(obs.list)){
  j <- tss(obs.list[[i]], inf.forbid.max.max.list[[i]])
  j$web <- names(obs.list)[[i]]
  j$step <- "prune"
  tss.max.max.prune.df <- rbind(tss.max.max.prune.df, j)
  rm(j)
}

max.max.tss.2step <- bind_rows(tss.max.max.df, tss.max.max.prune.df)

data.frame(mean.max = tss.3.step %>% 
             filter(step == "prune") %>%
             select(tss) %>% rename(mean.max = tss),
           max.max = tss.max.max.prune.df$tss) %>%
  summarize_all(funs(mean, sd))

tss.4.fish.size <- data.frame(max.max = 
                   max.max.tss.2step$tss,
                   mean.max = tss.3.step %>% 
                     filter(step == "model") %>%
                     select(tss) %>%
                     rename(mean.max = tss),
                   mean.min = tss.mean.min.df$tss,
                   min.min = tss.df$tss) 


tss.4.fish.size <- data.frame(step = max.max.tss.2step$step,
           max.max = max.max.tss.2step$tss,
           mean.max = mean.max.tss.2step$tss,
           mean.min = mean.min.tss.2step$tss,
           min.min = min.min.tss.2step$tss) %>%
  group_by(step) %>%
  summarize_at(c(2:5), funs(mean, sd))


tss.4.fish.size <- tss.4.fish.size %>% gather("fish.size", "tss", 2:9)
tss.4.fish.size <- tss.4.fish.size %>%
  separate(fish.size, c("fish.size", "var"), sep = "_")
tss.4.fish.size <- tss.4.fish.size %>% spread(var, tss)
# plot 4 fish sizes ####
# plot of mean TSS for 4 fish sizes
mean.tss.4.fish.plot <- ggplot(tss.4.fish.size,
                               aes(x = step, y = mean,
                                   color = fish.size)) +
  geom_jitter(width = 0.15) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean +sd)) +
  theme_classic() +
  labs(y = "Mean (+- SD) TSS")

saveRDS(mean.tss.4.fish.plot, paste(getwd(), "/figs for MS/mean tss 4 fish sizes.rds", sep = ""))

# plot fish size ~ tss
(tss.4.fish.plot <- ggplot(tss.4.fish.size %>%
                            filter(step == "model"),
                               aes(x = fct_reorder(fish.size, mean),
                                   y = mean)) +
  geom_point(size = 3)+
  geom_errorbar(aes(ymin = mean - sd, ymax = mean +sd), size = 1) +
  theme_classic() +
  labs(y = paste("Mean", "TSS", "\u00B1", "SD", sep = " "),
       x = "Fish length value used") 
) 

saveRDS(tss.4.fish.plot, paste(getwd(), "/figs for MS/4 fish size TSS.rds", sep = ""))

# mean.tss.4.fish.plot <- readRDS( paste(getwd(), "/figs for MS/mean tss 4 fish sizes.rds", sep = ""))

# thoughts 25 june 2017 ####
# coefficients for minimum min fish length generally higher
# compare TSS for mean.max (what I've been working with so far) to TSS for minimum min fish length

# average TSS for min.min is lower than mean.max
# After pruning:
# TSS is higher for some webs, lower for others
# mean min.min TSS is still lower than mean.max 

# looking at mean.min V mean.max
# TSS is higher and lower for some webs
# mean tss for mean.min is higher than mean.max
# 0.3455853 V 0.3558779
# after pruning, mean.min > mean.max
# 0.4760068 V 0.4847927

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# plot parameter coefficients
# calculate coeff mean +-sd
param.coef.summ <- param.coef.df %>% group_by(.id) %>%
  summarise_all(funs(mean, sd))
# get data into order for easier plotting
param.coef.summ <- param.coef.summ %>% 
  gather("coef", "value", 2:13)
param.coef.summ <- param.coef.summ %>%
  separate(coef, c("parameter", "var"))
param.coef.summ <- param.coef.summ %>% spread(var, value)
param.coef.summ <- param.coef.summ %>% rename(fish.size = .id)
param.coef.summ <- param.coef.summ %>%
  separate(parameter, c("parameter", "pos"), sep = 2)

# plot
(grav.params.4.fish.plot <- ggplot(param.coef.summ,
              aes(x = pos, y = mean, color = fish.size)) +
  geom_jitter(width = 0.1) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean+sd),
                width = 0.1) +
  facet_grid(parameter~., scales = "free") +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black",
                                    fill=NA, size=1)))

saveRDS(grav.params.4.fish.plot, paste(getwd(), "/figs for MS/gravel params 4 fish sizes.rds", sep = ""))





# plot all fish values for a site
ggplot(training.list.list[[3]][[4]], 
       aes(x = logcon, y = logres)) +
  geom_point() +
  stat_quantile(quantiles = c(0.03, .97)) +
  ggtitle(names(training.list.list)[3])

ggplot(training.list.list[[1]][[4]], 
       aes(x = logcon, y = logres)) +
  geom_point() +
  stat_quantile(quantiles = c(0.03, .97)) +
  ggtitle(names(training.list.list)[1])

ggplot(training.list.list[[2]][[4]], 
       aes(x = logcon, y = logres)) +
  geom_point() +
  stat_quantile(quantiles = c(0.03, .97)) +
  ggtitle(names(training.list.list)[2])

ggplot(training.list.list[[4]][[4]], 
       aes(x = logcon, y = logres)) +
  geom_point() +
  stat_quantile(quantiles = c(0.03, .97)) +
  ggtitle(names(training.list.list)[4])



# scrap ####
# code to work out fish lenght -> weight formulas
# salmo trutta ####
# a = 2.294e-6
# ln(a) = -12.9852
#mean max saltru = 171.6028
#mean min saltru = 83.526
-12.985 + (3.05 * log(83.52, base = exp(1)))
#log(dw) = 0.5115127
exp(1)^0.5115127
# dw = 1.667443
log10(1.667812)
# 0.222

-12.985 + 3.05 * log(130)
#log(dw) = 1.86098
exp(1)^1.86098
# dw = 6.430035
log10(6.430035)
#0.8082133

# galax
# mean length ~ 70
-14.77496979 + (3.417 * log(70,	base = exp(1)))
# log dw = -0.2578615
exp(1)^-0.2578615
# dw = 0.7727
log10(0.7727)
# -0.1119891

# gobio 
# mean length ~60
# G. brevicepps, jellyman 2013
# a = 2.155e-07
# b = 3.616
log(2.155e-07)
# ln(a) = -15.3503
-15.3503 + 3.616 * log(60, base = exp(1))
# log(dw) = -0.5451501
exp(1)^-0.5451501
# dw = 0.5797548
log10(0.5797548)
# log10dw = -0.2367556
# Anguilla
# mean length ~ 500
# A dieffenbachii, jellyman 2013
# a = 1.390e-8
# b = 3.641
log(1.390e-8)
# ln(a) = -18.09138
-18.09138 + 3.641 * log(500)
# log(dw) = 4.536008
exp(1)^4.536008
# dw = 93.31753
log10(93.31753)
# log10dw = 1.969963
