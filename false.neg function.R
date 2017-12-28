# false negative function
# JPOMZ October 2017

# function that extracts false negative (observed but not predicted) from inferred food web networks

false.neg <- function (observed, inferred){
  # make sure that observed and inferred matrices are same dimensions and have identical column names, and colnames in identical order
  stopifnot(dim(observed) == dim(inferred), 
            identical(colnames(observed),
                      colnames(inferred)))
  minus <- observed - inferred
  falseneg <- which(minus ==1, arr.ind = T)
  if(length(falseneg) == 0)
    return(matrix(0))
  else
    result <- NULL
    for (i in 1:length(falseneg[,1])){
      result <- rbind(result, names(minus[falseneg[i,]]))
    }
  colnames(result) <- c("resource", "consumer")
  result
}

