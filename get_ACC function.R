# get_ACC function
# J Pomz 2 dec 2017

# function to get accuracy (ACC) of observed and inferred adjacency matrices from confusion matrix derivations

get_ACC <- function(observed, inferred){
  # make sure adjacency matrices are same dimensions and have same colnames
  # identical() also checks that names are in same order
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
  # 2 = true positive (tp); links both obserevd & predicted
  # -1 = false positive (fp); predicted but not observed
  # 1 = false negative (fn); observed but not predicted
  # 0 = true negative (tn); not predicted, not observed
  
  tp = length(prediction[prediction==2])
  fp = length(prediction[prediction==-1]) 
  fn = length(prediction[prediction==1]) 
  tn = length(prediction[prediction==0])
  
  # ACC = accuracy = (tp + tn) / (tp + tn + fp +fn)
  ACC = (tp + tn) / (tp + tn + fp + fn)
  ACC
}