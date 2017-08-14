library(data.table)
library(MASS)

################# BIAS CORRECTION ################# 
### Misclassification Method
misclass.rate <- function(confusion_matrix){
  return( t( t(confusion_matrix)/colSums(confusion_matrix) ) )
}
misclass.nx. <- function(confusion_matrix, target_n.y){
  return(as.vector( solve(misclass.rate(confusion_matrix)) %*% as.matrix(target_n.y) ))
}
misclass.nxy <- function(confusion_matrix, target_n.y){
  return(round( t( t(misclass.rate(confusion_matrix))*misclass.nx.(confusion_matrix, target_n.y)) ))
}

### Reclassification Method
reclass.rate <- function(confusion_matrix){
  return( confusion_matrix/rowSums(confusion_matrix) )
}
reclass.nx. <- function(confusion_matrix, target_n.y){
  return(as.vector( round(t(reclass.rate(confusion_matrix)) %*% as.matrix(target_n.y)) ))
}
reclass.nxy <- function(confusion_matrix, target_n.y){
  return(round( reclass.rate(confusion_matrix)*target_n.y ))
}

### Ratio-to-TP Method
rtp.rate <- function(confusion_matrix){
  return( t( t(confusion_matrix)/diag(confusion_matrix) ) )
}
rtp.nxx <- function(confusion_matrix, target_n.y){
  return(as.vector( round(solve(rtp.rate(confusion_matrix)) %*% as.matrix(target_n.y))  ))
}
rtp.nxy <- function(confusion_matrix, target_n.y){
  rtp_rate <- rtp.rate(confusion_matrix)
  return( rtp(rtp_rate, rtp.nxx(rtp_rate, target_n.y)) )
}
rtp.nx. <- function(confusion_matrix, target_n.y){
  return( colSums(rtp.nxy(confusion_matrix, target_n.y)) )
}

### Combinaison of Misclassification and Reclassification Methoreclass
combi.matrix <- function(confusion_matrix, target_n.y){
  test_nx. <- colSums(confusion_matrix)
  misclass_nx. <- misclass.nx.(confusion_matrix, target_n.y)
  scale <- misclass_nx./test_nx.
  reclass_matrix <- t(t(confusion_matrix)*scale)
  return( reclass_matrix)
}
combi.nx. <- function(confusion_matrix, target_n.y){
  return( reclass.nx.(combi.matrix(confusion_matrix, target_n.y), target_n.y) )
}
combi.nxy <- function(confusion_matrix, target_n.y){
  return( reclass.nxy(combi.matrix(confusion_matrix, target_n.y), target_n.y) )
}

