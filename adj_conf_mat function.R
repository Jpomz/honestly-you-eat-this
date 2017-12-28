# Confusion matrix function
# JPomz 
# 2 December 2017

# make a function which calculates the confusion matrix for observed versus inferred adjacency matrices
# for more information on confusion matrices, see:
# https://en.wikipedia.org/wiki/Confusion_matrix
# accessed 2 December 2017  

adj_conf_matrix <- function (observed, inferred){
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
  vars <- data.frame(
    tp = length(prediction[prediction==2]),
    fp = length(prediction[prediction==-1]), 
    fn = length(prediction[prediction==1]), 
    tn = length(prediction[prediction==0]),
    S = ncol(prediction)
  )
  # derivations
  # TPR = true positive rate = tp / (tp + fn)
  # TNR = true negative rate = tn / (tn + fp)
  # PPV = positive predictive value = tp / (tp + fp)
  # NPV = negative predictive value = tn / (tn + fn)
  # FNR = false negative rate = fn / (fn + tp)
  # FPR = false positive rate = fp /(fp + tn)
  # FDR = false discovery rate = fp / (fp + tp)
  # FOR = false omission rate = fn / (fn +tn)
  # ACC = accuracy = (tp + tn) / (tp + tn + fp +fn)
  
  tss <- data.frame(tp = vars$tp,
                    fp = vars$fp,
                    fn = vars$fn,
                    tn = vars$tn,
                    TPR = vars$tp / 
                      (vars$tp + vars$fn),
                    TNR = vars$tn / 
                      (vars$tn + vars$fp),
                    PPV = vars$tp /
                      (vars$tp + vars$fp),
                    NPV = vars$tn /
                      (vars$tn + vars$fn),
                    FNR = vars$fn / 
                      (vars$fn + vars$tp),
                    FPR = vars$fp /
                      (vars$fp + vars$tn),
                    FDR = vars$fp /
                      (vars$fp + vars$tp),
                    FOR = vars$fn /
                      (vars$fn +vars$tn),
                    ACC = (vars$tp + vars$tn) /
            (vars$tp + vars$tn + vars$fp +vars$fn))
  tss
}






