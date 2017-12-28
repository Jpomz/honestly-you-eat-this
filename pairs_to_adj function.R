# pairs_to_adj()
# J Pomz Aug 2017

# function to make adjacency matrix from dataframe of pairwise species interactions.  

pairs_to_adj <- function(taxa, pairs){
  
  # taxa is a data.frame which must have a "node" column and contains a single row for each taxa present
  # if duplicated nodes, could use unique(taxa$node)
  # pairs is a data frame that has all biotic interactions in paired "resource" and "consumer" columns
  
  S = nrow(taxa)
  A = matrix(0, nrow = S, ncol = S)
  dimnames(A) = list(taxa$node, taxa$node)
  for (i in 1:nrow(pairs))
      A[pairs$resource[i], pairs$consumer[i]] = 1
  A[match(taxa, rownames(A)),
    match(taxa, colnames(A))]
  return(A)
}