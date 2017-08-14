source("code_biasCorrection.R")

################# GENERIC FIELLER INTERVAL #################
### Confidence Interval for Ratio A/B
######### DEFINITION OF VARIABLES
# mean_A, mean_B : Expected value of numerator and denominator.
# var_A, var_B : Variance of numerator and denominator.
# cov_AB : Covariance of numerator and denominator.
# z : Intended coverage of confidence interval, i.e., 68% CI for z=1, 95% CI for z=1.96

fieller.interval <- function(z, mean_A, var_A, mean_B, var_B, cov_AB){
  # Interval's Denonimator
  den <- mean_B^2 - z^2*var_B 
  # Interval's Numerator
  num_comp1 <- mean_A*mean_B - z^2*cov_AB
  num_comp2 <- num_comp1^2 - (mean_A^2-z^2*var_A)*den 
  if(num_comp2<0){
    print("No Interval")
    return(rep(mean_A/mean_B, 2))
  }
  # Intervals' bounds
  low <- (num_comp1 - sqrt(num_comp2))/den
  high <- (num_comp1 + sqrt(num_comp2))/den
  return(c(low, high))
}

################# APPLIED FIELLER INTERVAL ################# 

######### CLASS SIZE nx.
fieller.interval.n0. <- function(z=1, t01, t10, test_n0., test_n1., target_n.0, target_n.1, misclass_n0., misclass_n1.){
  # Sample-to-Sample Variance
  # t01
  if(misclass_n0.==0) var_t01 <- t01*(1-t01)/test_n0. 
  else if(test_n0.==0) var_t01 <- t01*(1-t01)/misclass_n0.
  else var_t01 <- t01*(1-t01)/test_n0. + t01*(1-t01)/misclass_n0.
  # t10
  if(misclass_n1.==0) var_t10 <- t10*(1-t10)/test_n1. 
  else if(test_n1.==0) var_t10 <- t10*(1-t10)/misclass_n1. 
  else var_t10 <- t10*(1-t10)/test_n1. + t10*(1-t10)/misclass_n1.
  
  # Ratio's Numerator
  mean_A <- target_n.0 - t10*(target_n.0+target_n.1)
  var_A <- var_t10*(target_n.0+target_n.1)^2
  # Ratio's Denominator
  mean_B <- 1 - t01 - t10
  var_B <- var_t01 + var_t10
  # Covariance of Numerator & Denominator
  cov_AB <- (target_n.0+target_n.1)*var_t10
  
  ### Interval
  return(fieller.interval(z, mean_A, var_A, mean_B, var_B, cov_AB))
}
fieller.interval.n1. <- function(z=1, t01, t10, test_n0., test_n1., target_n.0, target_n.1, misclass_n0., misclass_n1.){
  ### Parameters
  # Sample-to-Sample Variance
  var_t01 <- t01*(1-t01)/test_n0. + t01*(1-t01)/misclass_n0.
  var_t10 <- t10*(1-t10)/test_n1. + t10*(1-t10)/misclass_n1.
  # Ratio's Numerator
  mean_A <- target_n.1 - t01*(target_n.0+target_n.1)
  var_A <- var_t01*(target_n.0+target_n.1)^2
  # Ratio's Denominator
  mean_B <- 1 - t01 - t10
  var_B <- var_t01 + var_t10
  # Covariance of Numerator & Denominator
  cov_AB <- (target_n.0+target_n.1)*var_t01
  
  ### Interval
  return(fieller.interval(z, mean_A, var_A, mean_B, var_B, cov_AB))
}

######### ERROR DECOMPOSITION nxy
fieller.interval.n01 <- function(z=1, t01, t10, test_n0., test_n1., target_n.0, target_n.1){
  ### Sample-to-Sample Variance
  # t01
  if(misclass_n0.==0) var_t01 <- t01*(1-t01)/test_n0. 
  else if(test_n0.==0) var_t01 <- t01*(1-t01)/misclass_n0.
  else var_t01 <- t01*(1-t01)/test_n0. + t01*(1-t01)/misclass_n0.
  # t10
  if(misclass_n1.==0) var_t10 <- t10*(1-t10)/test_n1. 
  else if(test_n1.==0) var_t10 <- t10*(1-t10)/misclass_n1. 
  else var_t10 <- t10*(1-t10)/test_n1. + t10*(1-t10)/misclass_n1.

  # Ratio's Numerator
  n_x._A <- target_n.0 - (target_n.0+target_n.1)*t10
  mean_A <- t01*(target_n.0 - t10*(target_n.0+target_n.1))
  var_A <- t01^2 * (target_n.0+target_n.1)^2*var_t10 + n_x._A^2*var_t01 + (target_n.0+target_n.1)^2*var_t01*var_t10

  # Ratio's Denominator
  mean_B <- 1 - t01 - t10
  var_B <- var_t01 + var_t10
  
  # Covariance of Numerator & Denominator
  cov_AB <- (target_n.0+target_n.1)* (t10*var_t01  + t01*var_t10) - target_n.0*var_t01
  
  ### Interval
  return(fieller.interval(z, mean_A, var_A, mean_B, var_B, cov_AB))
}

######### REPRODUCTION OF PRIOR WORK from Shieh 2009, Buonaccorsi 2010
### Prior work concerns intervals describing a general population from which target set is sample. 
### Prior work concerns estimates of class proportions, e.g., pi0 = n0./n..  
### It can easiliy be modified to target true class estimates, i.e., n0. as function fieller.interval.n1.()
fieller.interval.shieh <- function(z=1, t01, t10, test_n0., test_n1., target_n.0, target_n.1){
  ### Basic Components
  pi_naive <- target_n.0 / (target_n.0+target_n.1)
  mean_A <- pi_naive - t10 
  mean_B <- 1-t10-t01 
  
  ### Variance Components
  var_A <- pi_naive*(1-pi_naive)/(target_n.0+target_n.1) + t10*(1-t10)/test_n1.
  var_B <- t01*(1-t01)/test_n0. + t10*(1-t10)/test_n1.
  cov_AB <- t10*(1-t10)/test_n1.
  
  ### Interval
  return(fieller.interval(z, mean_A, var_A, mean_B, var_B, cov_AB))
}