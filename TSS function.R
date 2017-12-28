# True Skill Statistic function
# JPomz 
# 19 June 2017

# make a function which calculates the TSS for observed versus inferred adjacency matrices

tss <- function (observed, inferred){
  # make sure adjacency matrices are same dimensions and have same colnames
  stopifnot(dim(observed) == dim(inferred), 
            identical(colnames(observed), colnames(inferred)))
  
  # subtract inferred from observed
  minus <- observed - inferred
  # multiply observed and inferred
  multiply <- observed * inferred
  # change values of 1 in multiplied matrix to 2
  multiply[multiply == 1] <- 2
  # add minus and multiply matrices
  prediction <- minus + multiply
  # prediction outcome matrix now has all 4 possibilities repreented as different integer values
  # 2 = true positive (a); links both obserevd & predicted
  # -1 = false positive (b); predicted but not observed
  # 1 = false negative (c); observed but not predicted
  # 0 = true negative (d); not predicted, not observed
  tss.vars <- data.frame(
    a = length(prediction[prediction==2]),
    b = length(prediction[prediction==-1]), 
    c = length(prediction[prediction==1]), 
    d = length(prediction[prediction==0]),
    S = ncol(prediction)
  )
  # calculate TSS
  # TSS = (a*d - b*c)/((a+c)*(b+d))
  # also decompose TSS to see proportion of positive and negative predictions
  # abar = proportion of links correctly predicted
  # dbar = proportion of links correctly not predicted
  # bc = proportion of links incorrectly predicted
  tss <- tss.vars %>%
    mutate(tss =(a*d - b*c)/((a+c)*(b+d)),
           abar = a/S**2,
           dbar = d / S**2,
           bbar = b / S**2,
           cbar = c / S**2,
           bc = (b + c) / S**2)
  tss
}






