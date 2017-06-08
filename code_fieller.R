source("code_biasCorrection.R")

################# FIELLER INTERVAL ################# 
fieller.interval <- function(z, mean_A, var_A, mean_B, var_B, cov_AB){
  #  print(paste("Interval Base", round(mean_A/mean_B) ))
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

fieller.interval.shieh2.corrected <- function(z=1, t01, t10, test_n0., test_n1., target_n.0, target_n.1, misclass_n0., misclass_n1.){
  # Sample-to-Sample Variance
  # t01
  if(misclass_n0.==0) var_t01 <- t01*(1-t01)/test_n0. 
  else if(test_n0.==0) var_t01 <- t01*(1-t01)/misclass_n0.
  else var_t01 <- t01*(1-t01)/test_n0. + t01*(1-t01)/misclass_n0.
  # t10
  if(misclass_n1.==0) var_t10 <- t10*(1-t10)/test_n1. 
  else if(test_n1.==0) var_t10 <- t10*(1-t10)/misclass_n1. 
  else var_t10 <- t10*(1-t10)/test_n1. + t10*(1-t10)/misclass_n1.
  
  ### Basic Components
  pi_naive <- target_n.0 / (target_n.0+target_n.1)
  mean_A <- pi_naive - t10 
  mean_B <- 1-t10-t01 
  
  ### Variance Components
  var_A <- var_t10
  #var_A <- pi_naive*(1-pi_naive)/(target_n.0+target_n.1) + var_t10
  var_B <- var_t01 + var_t10
  cov_AB <- var_t10
  
  ### Interval
  return(fieller.interval(z, mean_A, var_A, mean_B, var_B, cov_AB))
}
fieller.interval.shieh2 <- function(z=1, t01, t10, test_n0., test_n1., target_n.0, target_n.1){
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
fieller.interval.shieh <- function(z=1, t01, t10, test_n0., test_n1., target_n.0, target_n.1){
  ### Basic Components
  pi_naive <- target_n.0 / (target_n.0+target_n.1)
  N <- pi_naive - t10 
  D <- 1-t10-t01 
  #  N00 <- test_n0.*(1-t01)
  #  N11 <- test_n1.*(1-t10)
  
  ### Variance Components
  sigma11 <- pi_naive*(1-pi_naive)/(target_n.0+target_n.1) + t10*(1-t10)/test_n1.
  sigma22 <- t01*(1-t01)/test_n0. + t10*(1-t10)/test_n1.
  sigma12 <- t10*(1-t10)/test_n1.
  
  ### Advanced Components
  f0 <- N^2 - z^2*sigma11
  f1 <- D*N - z^2*sigma12
  f2 <- D^2 - z^2*sigma22
  C <- f1^2 - f2*f0
  
  ### Check type of interval & return results
  if((!C >= 0) | (!f2 >= 0)) print("Atypical")
  if(C<0) {
    print("C<0")
    return(rep(N/D,2))
  }
  
  ### Interval
  r1 <- (f1 + sqrt(C))/f2
  r2 <- (f1 - sqrt(C))/f2
  
  return(c(r2,r1))
}

fieller.interval.shieh.corrected <- function(z=1, t01, t10, test_n0., test_n1., target_n.0, target_n.1, misclass_n0., misclass_n1.){
  ### Sample-to-Sample Variance
  # t01
  if(misclass_n0.==0) var_t01 <- t01*(1-t01)/test_n0. 
  else if(test_n0.==0) var_t01 <- t01*(1-t01)/misclass_n0.
  else var_t01 <- t01*(1-t01)/test_n0. + t01*(1-t01)/misclass_n0.
  # t10
  if(misclass_n1.==0) var_t10 <- t10*(1-t10)/test_n1. 
  else if(test_n1.==0) var_t10 <- t10*(1-t10)/misclass_n1. 
  else var_t10 <- t10*(1-t10)/test_n1. + t10*(1-t10)/misclass_n1.
  
  ### Basic Components
  pi_naive <- target_n.0 / (target_n.0+target_n.1)
  N <- pi_naive - t10 
  D <- 1-t10-t01 
  
  ### Variance Components
  sigma11 <- pi_naive*(1-pi_naive)/(target_n.0+target_n.1) + var_t10
  sigma22 <- var_t01 + var_t10
  sigma12 <- var_t10
  
  ### Advanced Components
  f0 <- N^2 - z^2*sigma11
  f1 <- D*N - z^2*sigma12
  f2 <- D^2 - z^2*sigma22
  C <- f1^2 - f2*f0
  
  ### Check type of interval & return results
  if((!C >= 0) | (!f2 >= 0)) print("Atypical")
  if(C<0) {
    print("C<0")
    return(rep(N/D,2))
  }
  
  ### Interval
  r1 <- (f1 + sqrt(C))/f2
  r2 <- (f1 - sqrt(C))/f2
  
  return(c(r2,r1))
}

################# TEST INTERVAL ################# 
multi.test.fieller <- function(z=1, nTest=1000, nTarget=100, population_t01_t10=c(0.1,0.2), base_nx.=c(300,500,1000,2000), output="n01"){
	# Init 
  res <- data.frame()
	
  # Evaluate Fieller's Interval for All Combination of test/target nx.
	for(test_n0. in base_nx.){
	  for(test_n1. in base_nx.){
		  r <- c()
		  for(target_n0. in base_nx.){
		    for(target_n1. in base_nx.){
		     # r <- c(r, test.fieller(z, nTest, nTarget, population_t01_t10, test_nx.=c(test_n0.,test_n1.), target_nx.=c(target_n0.,target_n1.), output) )
		      r <- c(r, test.fieller.pop(z, nTest, nTarget, population_t01_t10, test_nx.=c(test_n0.,test_n1.), pop_nx.=c(target_n0.,target_n1.), output) )
		    }
  		}
		  res <- rbind(res, r)
	  }
	}
	
	# Format Results
  colnames <- c()
  for(i in 1:length(base_nx.)) {
    for(k in 1:length(base_nx.)) colnames <- c(colnames, paste(paste0("n0.", base_nx.[i]),paste0("n'0.", base_nx.[k]), sep="_"))
  }
  rownames <- gsub("n0", "n1", colnames)
  colnames(res) <- colnames
  rownames(res) <- rownames

  # Write Results
	spec <- paste0(output, "_z", z, "_",paste(population_t01_t10, collapse="-")  ) # , "_nx.=",paste(base_nx., collapse="-")
	try( write.csv(res, paste0("Result_Fieller/",spec,".csv")) )

	return(res)
}

test.fieller <- function(z=1, nTest=100, nTarget=100, population_t01_t10=c(0.1,0.2), test_nx.=c(300,300), target_nx.=c(300,500), output="n01"){
  # Init
  success <- count_small_misclass <- count_neg <- error <- 0
  if(output=="shieh") target_nx.<-2*target_nx.
  print(paste("n0.:", test_nx.[1], " n1.:", test_nx.[2], " -  n'0.:", target_nx.[1], " n'1.:", target_nx.[2]))
  
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
    count_neg <- 0
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
      ### Apply Misclassification Method
      misclass_nx. <- misclass.nx.(test_matrix, c(target_n.0[k], target_n.1[k]))
      # if(misclass_nx.[1]<=5) count_small_misclass <- count_small_misclass+1 
      # print(paste("Misclass.<=5:", misclass_nx.[1]))
      if(output=="n01") {
        misclass_n01 <- misclass_nx.[1]*(test_n01[i]/test_n0.) 
        # if(misclass_n01<=5) count_small_misclass <- count_small_misclass+1 
        # print(paste("Misclass.<=5:", misclass_n01))
      }
      
      ### Handle Negative Values
      if(misclass_nx.[1]<=0){
        #        print(paste("n0. <= 0", misclass_nx.[1], misclass_nx.[2]))
        misclass_nx.[2] <- misclass_nx.[2] + misclass_nx.[1]
        misclass_nx.[1] <- 0
        count_neg <- count_neg+1
      }
      if(misclass_nx.[2]<=0){
        #       print(paste("n1. <= 0", misclass_nx.[2], misclass_nx.[1]))
        misclass_nx.[1] <- misclass_nx.[2] + misclass_nx.[1]
        misclass_nx.[2] <- 0
        count_neg <- count_neg+1
      }
      
      ### Get Fieller's confidence interval
      # Rounding
      #      misclass_nx. <- round(misclass_nx.)
      # Make interval
      if(output=="n0.") interval <- fieller.interval.n0.(z, t01=(test_n01[i]/test_n0.), t10=(test_n10[i]/test_n1.), test_n0., test_n1., target_n.0[k], target_n.1[k], misclass_n0.=misclass_nx.[1], misclass_n1.=misclass_nx.[2]) 
      else if(output=="n1.") interval <- fieller.interval.n1.(z, t01=(test_n01[i]/test_n0.), t10=(test_n10[i]/test_n1.), test_n0., test_n1., target_n.0[k], target_n.1[k], misclass_n0.=misclass_nx.[1], misclass_n1.=misclass_nx.[2]) 
      else if(output=="n01") interval <- fieller.interval.n01(z, t01=(test_n01[i]/test_n0.), t10=(test_n10[i]/test_n1.), test_n0., test_n1., target_n.0[k], target_n.1[k], misclass_n0.=misclass_nx.[1], misclass_n1.=misclass_nx.[2]) 
      #else if(output=="shieh") interval <- (target_n.0[k]+target_n.1[k])*fieller.interval.shieh2(z, t01=(test_n01[i]/test_n0.), t10=(test_n10[i]/test_n1.), test_n0., test_n1., target_n.0[k], target_n.1[k])
      else if(output=="shieh") interval <- fieller.interval.shieh2(z, t01=(test_n01[i]/test_n0.), t10=(test_n10[i]/test_n1.), test_n0., test_n1., target_n.0[k], target_n.1[k])
     # else if(output=="shieh") interval <- fieller.interval.shieh2.corrected(z, t01=(test_n01[i]/test_n0.), t10=(test_n10[i]/test_n1.), test_n0., test_n1., target_n.0[k], target_n.1[k], misclass_n0.=misclass_nx.[1], misclass_n1.=misclass_nx.[2])
      # Rounding
      #      interval <- round(interval)
      #      if(output=="n01") misclass_n01 <- round(misclass_n01)
      
      ### Test interval
      if(output=="n01"){ 
        if(min(interval)<=target_n01[k] & max(interval)>=target_n01[k]) success <- success + 1
      } else if(output=="n0.") {
        if(min(interval)<=target_n0. & max(interval)>=target_n0.) success <- success + 1
      } else if(output=="shieh") {
        #if(interval[1]<=target_n0. & interval[2]>=target_n0.) success <- success + 1
        #if(interval[1]<=target_n0./(target_n.0[k]+target_n.1[k]) & interval[2]>=target_n0./(target_n.0[k]+target_n.1[k])) success <- success + 1
        if(min(interval)<=target_n0./(target_n.0[k]+target_n.1[k]) & max(interval)>=target_n0./(target_n.0[k]+target_n.1[k])) success <- success + 1
      } else {
        if(min(interval)<=target_n1. & max(interval)>=target_n1.) success <- success + 1
      }
    }
  }
  
  ### Print
  #      if(output=="n01") print(paste("Target:", target_n01[k], "Misclass:", misclass_n01, "Interval:", paste(interval, collapse=" - ")))
  #      print(paste("Target:", target_n0., "Output:", target_n.0[k], "Misclass:", misclass_nx.[1], "Interval:", paste(interval, collapse=" - ")))
  #      if(output=="n01") error <- c(error, misclass_n01-target_n01[k])
  #      print(paste("Reproduce Shieh (incorrect):", (target_n.0[k]/(target_n.0[k]+target_n.1[k])-target_n01/target_n0.)/(1-target_n01/target_n0.-target_n10/target_n1.) ))     
  #      print(paste("Corrected Shieh:", (target_n.0[k]/(target_n.0[k]+target_n.1[k])-target_n10/target_n1.)/(1-target_n01/target_n0.-target_n10/target_n1.) ))     
  #Verif Shieh:  print(paste("Target:", target_n0., "Output:", target_n.0[k], "Misclass:", misclass_nx.[1])) #, "Interval:", paste(interval, collapse=" - ")))
  #Verif Shieh:  print(paste("Corrected Shieh:", (target_n.0[k]+target_n.1[k])* (target_n.0[k]/(target_n.0[k]+target_n.1[k])-target_n10[k]/target_n1.)/(1-target_n01[k]/target_n0.-target_n10[k]/target_n1.) ))     
  #      print(paste("TRUE RATE:", target_n0./(target_n0.+target_n1.) ))
  #      print(paste("Misclass.<=5:", count_small_misclass))
  #      print(paste0("Neg.Values: ", count_neg, " /", 2*nTest*nTarget))
  #      print(mean(abs(error))) 
  #      print(paste("Interval:", interval))
  
  ### Return % of intervals containing target_n0. or target_n01
  return(100*success/(nTest*nTarget))
}

test.fieller.pop <- function(z=1, nTest=100, nTarget=100, population_t01_t10=c(0.1,0.2), test_nx.=c(300,300), pop_nx.=c(300,500), output="shieh"){
  # Init
  success <- count_small_misclass <- count_neg <- error <- 0
  pop_nx.<-2*pop_nx.
  print(paste("n0.:", test_nx.[1], " n1.:", test_nx.[2], " -  n*0.:", pop_nx.[1], " n*1.:", pop_nx.[2]))
  
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
    
    ### Test Fieller
    for(k in 1:nTarget){
      ### Make Target Sets
      # Class 0
      #    target_n0. <- target_nx.[1]
      target_n0. <- rbinom(1, pop_nx.[1]+pop_nx.[2], pop_nx.[1]/(pop_nx.[1]+pop_nx.[2]) )
      target_n01 <- rbinom(1, target_n0., population_t01_t10[1])
      target_n00 <- target_n0. - target_n01
      # CLass 1
      #   target_n1. <- target_nx.[2]
      target_n1. <- pop_nx.[1]+pop_nx.[2] - target_n0.
      target_n10 <- rbinom(1, target_n1., population_t01_t10[2])
      target_n11 <- target_n1. - target_n10
      # target_n.y
      target_n.0 <- target_n00 + target_n10
      target_n.1 <- target_n01 + target_n11
      
      ### Apply Misclassification Method
      misclass_nx. <- misclass.nx.(test_matrix, c(target_n.0, target_n.1))
      
      ### Handle Negative Values
      if(misclass_nx.[1]<=0){
        misclass_nx.[2] <- misclass_nx.[2] + misclass_nx.[1]
        misclass_nx.[1] <- 0
      }
      if(misclass_nx.[2]<=0){
        misclass_nx.[1] <- misclass_nx.[2] + misclass_nx.[1]
        misclass_nx.[2] <- 0
      }
      
      ### Get Fieller's confidence interval
      # Make interval
      interval <- fieller.interval.shieh2(z, t01=(test_n01[i]/test_n0.), t10=(test_n10[i]/test_n1.), test_n0., test_n1., target_n.0, target_n.1)
     # interval <- fieller.interval.shieh2.corrected(z, t01=(test_n01[i]/test_n0.), t10=(test_n10[i]/test_n1.), test_n0., test_n1., target_n.0, target_n.1, misclass_n0.=misclass_nx.[1], misclass_n1.=misclass_nx.[2])
      
      ### Test interval
      if(min(interval)<=pop_nx.[1]/(pop_nx.[1]+pop_nx.[2]) & max(interval)>=pop_nx.[1]/(pop_nx.[1]+pop_nx.[2]) ) success <- success + 1

    }
  }
  
  ### Return % of intervals containing target_n0. or target_n01
  return(100*success/(nTest*nTarget))
}