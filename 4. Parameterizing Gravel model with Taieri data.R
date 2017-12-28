# parameterize gravel model with Taieri biomass data
# JPomz
# june 2017

# start with adj list corrected in "fixing taxa...R" script


# useful food web functions from Petchey
source("C:\\Users\\Justin\\Documents\\Data\\FW modelling Petchey Github\\ttl-resources-master\\food_web\\FoodWebFunctions.r")
# gravel functions
source("gravel_functions.R")

# read in estimated Taieri dry weights 
taieri <- readRDS("estimated dw taieri webs.rds")

# get all body sizes present at sites
# remove NAs
taieri <- taieri %>% filter(!is.na(dw))
# sort
taieri <- taieri %>% arrange(site, logdw)
# make into a list
taieri.list <- split(taieri, list(taieri$site))

# read in corrected adjaceny matrices
web.list <- readRDS("animal only adj list.rds")

# subset web.list to match only names in taieri.list
web.list <- web.list[names(web.list) %in% names(taieri.list)]
#saveRDS(web.list, "17 observed taieri food webs.rds")

# convert predation matrix to matrix of pairs
pairs.list <- web.list %>%
  llply(function (x){
    as.data.frame(Matrix.to.list(x))})

# name columns resource and consumer
pairs.list <- llply(pairs.list,
                    function (x){
                      colnames(x) <- c("resource", "consumer");x})

# make list of resource consumer pairs with dw
pairs.dw.list <- NULL
for (i in 1:length(pairs.list)){
  resource <- left_join( # add dw for resource column
    pairs.list[[i]], taieri.list[[i]][,c(2,6)],
    by = c("resource" = "taxa"))
  consumer <- left_join( # just consumer dw
    pairs.list[[i]], taieri.list[[i]][,c(2,6)],
    by = c("consumer" = "taxa"))[3]
  combined <- bind_cols(resource, consumer) # combine
  colnames(combined) <- c("resource", "consumer", "res.dw", "con.dw") # rename columns
  pairs.dw.list[[i]] <- combined # put in list
  rm("resource", "consumer", "combined") # remove objects
}
# should throw ~36 warnings
# this is bc output of Matrix.to.list() [line 24] is class == factor
# left_join() commands above [lines 34, 37] convert factor to character class and combine accordingly
# output is correct, can ignore warnings

# names elements in list to match site names
names(pairs.dw.list) <- names(pairs.list)

# all pairs ####
# #this is for chapter 3
# all.pairs <- ldply(pairs.dw.list)
# all.pairs <- all.pairs %>%
#   select(res.dw, con.dw) %>%
#   filter(!is.na(res.dw), !is.na(con.dw)) %>%
#   mutate(logres = log10(res.dw),
#          logcon = log10(con.dw)) %>%
#   select(logres, logcon)
# # # save this for use in inferring fw structure in mine streams
# # saveRDS(all.pairs, "all predation pairs.rds")
# 

# # differences by site?
# site.pairs <- ldply(pairs.dw.list)
# site.pairs <- site.pairs %>%
#   select(.id, res.dw, con.dw) %>%
#   filter(!is.na(res.dw), !is.na(con.dw)) %>%
#   mutate(logres = log10(res.dw),
#          logcon = log10(con.dw)) %>%
#   select(.id, logres, logcon)
# # plot of all pairs color = .id
# ggplot(site.pairs, aes(x = logcon,
#                       y = logres,
#                       color = .id)) +
#   geom_point() +
#   stat_quantile(quantiles = c(0.01, .97))
# 
# # sidebar looking at when resource > consumer
# all.pairs <- all.pairs %>% 
#   mutate(res.greater= logres > logcon)
# 
# all.pairs %>%
#   count(res.greater)
# # FALSE = 423
# # TRUE = 60
# # per cent res.greater = 0.141844
# 
# all.pairs %>%
#   filter(res.greater == TRUE) %>%
#   mutate(ratio = logres / logcon) %>%
#   summarize(mean = mean(ratio),
#             min = min(ratio), 
#             max = max(ratio))
# # smallest pred = 70% of the body size of prey
# # biggest = 100%
# # mean = 90%


# make list of res-con biomass pairs for model parameterizations
# each element in list contains all pairs EXCEPT those of
# the web they will be used to parameterize. 
# e.g. element named "blackrock" will contain all other
# pairs data [-"blackrock"], and be used to paramterize 
# the model in order to predict the blackrock web 
training.list <- NULL
for (i in 1:length(pairs.dw.list)){
  a <- ldply(pairs.dw.list[-i])
  training.list[[i]] <- a
  rm(a)
}

names(training.list) <- names(pairs.dw.list)

# select estimated dry weights
# filter out dryweight NA's
# log10 transform dry weights, test if res > con
# filter out res > con == TRUE
# commented out filter step to see if it changes results...
training.list <- training.list %>%
  llply(function (x){
    x %>%
      select(res.dw, con.dw) %>%
      filter(!is.na(res.dw), !is.na(con.dw)) %>%
      mutate(logres = log10(res.dw),
             logcon = log10(con.dw),
             res.greater = res.dw > con.dw) #%>%
#      filter(res.greater == FALSE)
  })
# saveRDS(training.list, "complete gravel training list.rds")

# # get info on training list for results ####
# n.vector <- numeric(length = length(training.list))
# for (i in 1:length(training.list)){
#   n.vector[i] <- nrow(training.list[[i]])
# }
# max(n.vector)
# # 471
# min(n.vector)
# #379

# # plot all training data sets ####
# llply(training.list, function (x){
#   ggplot(x, aes(x = logcon, y = logres)) +
#   geom_point() +
#   stat_quantile(quantiles = c(0.01, .97))}
# )
# # plot just one for power point outline
# ggplot(training.list[[1]], aes(x = logcon, y = logres)) +
#   geom_point(shape = 1, size = 3) +
#   stat_quantile(quantiles = c(0.01, .97))+
#   labs(y = "log prey body size", x = "log predator body size", size = 4) +
#   theme_classic() +
#   theme(text=element_text(size=20))

# parameterize gravel model ####
# parametrize model for each food web
# make a list of regression coefficients
pars.list <- NULL
for (i in 1:length(training.list)){
  Bprey <- training.list[[i]]$logres
  Bpred <- training.list[[i]]$logcon
  pars.list[[i]] <- reg_fn(
    Bprey,Bpred,quartil = c(0.01,0.97))
}
names(pars.list) <- names(training.list)

# make object for food web niche parameters
# 
web.pars <- NULL
for (i in 1:length(pars.list)){
  web.pars[[i]] <- get_pars_Niche(pars.list[[i]], Ball = taieri.list[[i]]$logdw)
}

# add taxa names to web.pars
for (i in 1:length(web.pars)){
  web.pars[[i]] <- as.data.frame(web.pars[[i]])
  web.pars[[i]]$taxa <- taieri.list[[i]]$taxa
}

# save taxa names as .Rds object for "truncating...R"
taxa.inf.links <- taieri.list %>%
  llply(function (x) {
    x$taxa
  })

#saveRDS(taxa.inf.links, "taxa names used to infer links.rds")

# # plot niche parameters for all training data sets
# llply(web.pars, function (x){
#   ggplot(as.data.frame(x), aes(x = n, y = n)) +
#     geom_point() +
#     geom_line(aes(x = n, y = low), lty = 4) +
#     geom_line(aes(x = n, y = high), lty = 3) +
#     geom_line(aes(x = n, y = c), lty = 1)}
# )

# # make a threshold for small bugs
# # threshold roughly determined in text file
# treshold.list <- llply(web.pars, function (x){
#   which(x[,"n"]<=-3.7)
# })
# 
# # set small bugs hi/low niche to 0
# for (i in 1:length(web.pars)){
#   web.pars[[i]][,"low"][treshold.list[[i]]] <- 0
#   web.pars[[i]][,"high"][treshold.list[[i]]] <- 0
# }


# calculate food web links for taieri
# predation matrix
web.links.inf <- llply(web.pars, function (x){
  L_fn(x[,"n"],x[,"c"],x[,"low"],x[,"high"])
})

# add taxa names to predation matrices
for (i in 1:length(web.links.inf)){
  dimnames(web.links.inf[[i]]) <- list(web.pars[[i]]$taxa,
                                   web.pars[[i]]$taxa)
}

# names elements in list to match site names
names(web.links.inf) <- names(pairs.list)

# forbidden links ####



# save inferred matrices as rds object
saveRDS(web.links.inf, file = "inferred adj matr gravel method.rds")


# make table of model information for manuscript / SI
# want site, n pairwise intxn, parameter coefficients
# na pairwise intxns
n.pairs <- training.list %>%
  ldply(function (x){
    data.frame(n = nrow(x))
  })
param.coef <- pars.list %>% 
  ldply(function (x){
    data.frame(B0hi = x[[1]][1],
               B0center = x[[2]][1],
               B0lo = x[[3]][1],
               B1hi = x[[1]][2],
               B1center = x[[2]][2],
               B1lo = x[[3]][2])
  })
# join two data frames
model.table <- left_join(n.pairs, param.coef, by = ".id")

# write csv of paramters
# modify in excel to make table for publication
write_csv(model.table, "C:/Users/Justin/Google Drive/Data/Predicting NZ Food Webs/figs for MS/gravel model parameters.csv")







