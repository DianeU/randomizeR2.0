condRefSet <- function(obs,refSet){
  #Check how many NA's in y and on which positions
  y <- obs@resp
  posNA <- which(is.na(y))
  #Which groups are on this position
  referenz<- (obs@rs@M)[posNA]
  #Sum of all B's on the NA positions in the responce vector
  B_missing <- sum(referenz)
  
  # NASeq <- t(as.matrix(obs@rs@M[!(is.na(y))]))
  
  hasSameBs <- apply(refSet@M, 1, function(x) sum(x[posNA]) == B_missing)
  newRefSet <- refSet@M[hasSameBs,]
  if (is.vector(newRefSet)) newRefSet <- t(as.matrix(newRefSet, ncol = length(newRefSet)))
  return(newRefSet)
}