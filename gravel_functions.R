##############################################
##############################################
#
# R Code supplementing the paper 
# Inferring food web structure form predator-prey body size relationship
# by Gravel, Poisot, Albouy, Velez, Mouillot
# Methods in Ecology and Evolution
# PUTS THE VOLUME/ISSUE/PAGES HERE
# February 2013
# 
##############################################
##############################################

# 2. Useful functions
# from gravel et al. 2013


##############################################
# Get regression parameters
# Input arguments:
# Bprey = log10 biomass of the prey
# Bpred = log10 biomass of the predator
# Quartil = a vector of the inferior and the superio quartile c(0.03,0.97)
# Returns a list of regression objectis
# Requires the quantreg package

reg_fn = function(Bprey,Bpred,quartil) {
  
  library(quantreg)
  mean_reg = lm(Bprey~Bpred)			# For the n parameter
  qrsup = rq(Bprey~Bpred,tau = quartil[2])	# For the higher limit of the range
  qrinf = rq(Bprey~Bpred,tau = quartil[1])	# For the lower limit of the range
  
  return(list(mean_reg$coef,qrsup$coef,qrinf$coef))
  
}

##############################################
# Estimate the niche parameters for all species of a list 
# Input arguments: 
# pars = resulting parameters of the function reg_Niche
# Ball = list of body size
# Returns a matrix with four parameters for each species

get_pars_Niche = function(pars,Ball) {
  
  mean_reg = pars[[1]]
  qrsup = pars[[2]]
  qrinf = pars[[3]]
  
  # Estimate parameters for the allometric relationships
  delta = mean_reg[2]
  b1 = mean_reg[1]
  b2 = delta	
  
  # Estimate the parameters for the niche model
  n = Ball						# The niche n
  c = b1 + b2*Ball				# The centroid c
  low = qrinf[1] + qrinf[2]*Ball	# The lower limit of the range
  high = qrsup[1] + qrsup[2]*Ball	# The higher limit of the range
  
  return(cbind(n,c,low,high))	
}

##############################################
# Transform the parameters into an interaction matrix (the metaweb)
# Input:
# n = vector of size S with the parameter n for each of the S species
# c = vector of size S with the parameter c for each of the S species
# low = vector of size S with the parameter low for each of the S species
# high = vector of size S with the parameter high for each of the S species
# Returns a SxS matrix with 0 indicating absence of a link and 1 indicating the presence of a link
# Predators on columns, preys on rows

L_fn = function(n,c,low,high) {
  
  S = length(n)   	
  L = matrix(0,nr=S,nc=S)
  
  for(i in 1:S)
    for(j in 1:S)
      if(n[j]>low[i] && n[j]<high[i]) L[j,i] = 1
  
  return(L)	
  
}

##############################################
