source("code_biasCorrection.R")

################# FIELLER INTERVAL ################# 
fieller.interval.n0. <- function(z=1, t01, t10, test_n0., test_n1., target_n.0, target_n.1, misclass_n0., misclass_n1.){
  ### Parameters
  # Sample-to-Sample Variance
  var_t01 <- t01*(1-t01)/test_n0. + t01*(1-t01)/misclass_n0.
  var_t10 <- t10*(1-t10)/test_n1. + t10*(1-t10)/misclass_n1.
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
fieller.interval <- function(z, mean_A, var_A, mean_B, var_B, cov_AB){
  # Interval's Denonimator
  den <- mean_B^2 - z^2*var_B 
  # Interval's Numerator
  num_comp1 <- mean_A*mean_B - z^2*cov_AB
  num_comp2 <- sqrt( num_comp1^2 - (mean_A^2-z^2*var_A)*den )
  # Intervals' bounds
  low <- (num_comp1 - num_comp2)/den
  high <- (num_comp1 + num_comp2)/den
  return(c(low, high))
}

################# TEST INTERVAL ################# 
multi.test.fieller <- function(z=1, nTest=1000, nTarget=100, population_t01_t10=c(0.1,0.2), base_nx.=c(300,500,1000,2000)){
	# Init 
  res <- data.frame()
	
  # Evaluate Fieller's Interval for All Combination of test/target nx.
	for(test_n0. in base_nx.){
	  for(test_n1. in base_nx.){
		  r <- c()
		  for(target_n0. in base_nx.){
		    for(target_n1. in base_nx.){
  			  r <- c(r, test.fieller(z, nTest, nTarget, population_t01_t10, test_nx.=c(test_n0.,test_n1.), target_nx.=c(target_n0.,target_n1.)) )
		    }
  		}
		  res <- rbind(res, r)
	  }
	}
	
	# Format Results
  names <- c()
  for(i in 1:length(base_nx.)) {
    for(k in 1:length(base_nx.)) names <- c(names, paste(paste0("n0.", base_nx.[i]),paste0("n0.", base_nx.[k]), sep="_"))
  }
  colnames(res) <- rownames(res) <- names

  # Write Results
	spec <- paste0("t01t10=",paste(population_t01_t10, collapse="-"), "_nx.=",paste(base_nx., collapse="-")  )
	write.csv(res, paste0("Result_Fieller/",spec,".csv"))

	return(res)
}
test.fieller <- function(z=1, nTest=1000, nTarget=100, population_t01_t10=c(0.1,0.2), test_nx.=c(300,300), target_nx.=c(300,500)){
  # Init
  success <- 0
  
  ### Make Vector of Test Sets
  # Class 0
  test_n0. <- test_nx.[1]
  test_n01 <- rbinom(nTest, test_n0., population_t01_t10[1])
  test_n00 <- test_n0. - test_n01
  # CLass 1
  test_n1. <- test_nx.[2]
  test_n10 <- rbinom(nTest, test_n1., population_t01_t10[2])
  test_n11 <- test_n1. - test_n10
  
  ### Test Intervals
  for(i in 1:nTest){
    # Make Test Set's confusion matrix
    test_matrix <- matrix(c(test_n00[i], test_n01[i], test_n10[i], test_n11[i]), ncol=2, nrow=2)
    
    ### Make Vector of Target Sets
    # Class 0
    target_n0. <- target_nx.[1]
    target_n01 <- rbinom(nTarget, target_n0., population_t01_t10[1])
    target_n00 <- target_n0. - target_n01
    # CLass 1
    target_n1. <- target_nx.[2]
    target_n10 <- rbinom(nTarget, target_n1., population_t01_t10[2])
    target_n11 <- target_n1. - target_n10
    # target_n.y
    target_n.0 <- target_n00 + target_n10
    target_n.1 <- target_n01 + target_n11
    
	  ### Test Fieller
    for(k in 1:nTarget){
      # Apply Misclassification Method
      misclass_nx. <- misclass.nx.(test_matrix, c(target_n.0[k], target_n.1[k]))
      # Get Fieller's confidence interval (for Class 0)
      interval <- fieller.interval.n0.(z, t01=(test_n01[i]/test_n0.), t10=(test_n10[i]/test_n1.), test_n0., test_n1., 
                                      target_n.0[k], target_n.1[k], misclass_n0.=misclass_nx.[1], misclass_n1.=misclass_nx.[2])
      # Test interval
      if(interval[1]<target_n0. & interval[2]>target_n0.) success <- success + 1
    }
  }
  ### Return % of intervals containing target_n0.
  return(100*success/(nTest*nTarget))
}