source("code_biasCorrection.R")
source("code_fieller.R")

################# TEST INTERVAL ################# 
### SPECIFICATION OF VARIABLE
# output="n0." to test intervals for class size estimates of class C0
# output="n1." to test intervals for class size estimates of class C1
# output="n01" to test intervals for detailed error estimates (error composition), i.e., for number of items belonging to class C0 and misclassified in C1.  
# output="shieh" to test intervals from prior work in Shieh 2009, describing a general population from which target set is sample, and describing estimates of class proportions, e.g., pi0 = n0./n.. 

######### EXECUTE MULTIPLE TESTS
multi.test.fieller <- function(z=1, nTest=1000, nTarget=100, population_t01_t10=c(0.1,0.2), base_nx.=c(300,500,1000,2000), output="n01", test_pop=FALSE){
	# Init 
  res <- data.frame()
	
  # Evaluate Fieller's Interval for All Combination of test/target nx.
	for(test_n0. in base_nx.){
	  for(test_n1. in base_nx.){
		  r <- c()
		  for(target_n0. in base_nx.){
		    for(target_n1. in base_nx.){
		      if(test_pop) r <- c(r, test.fieller.pop(z, nTest, nTarget, population_t01_t10, test_nx.=c(test_n0.,test_n1.), pop_nx.=c(target_n0.,target_n1.)) )
          else r <- c(r, test.fieller(z, nTest, nTarget, population_t01_t10, test_nx.=c(test_n0.,test_n1.), target_nx.=c(target_n0.,target_n1.), output) )
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
  dir.create("Result_Fieller")
	spec <- paste0(output, "_z", z, "_t=",paste(population_t01_t10, collapse="-"), "_nx.=",paste(base_nx., collapse="-"), "_pop=", test_pop) 
	try( write.csv(res, paste0("Result_Fieller/",spec,".csv")) )

	return(res)
}

######### TEST OF SAMPLE-TO-SAMPLE INTERVALS
test.fieller <- function(z=1, nTest=100, nTarget=100, population_t01_t10=c(0.1,0.2), test_nx.=c(300,300), target_nx.=c(300,500), output="n01"){
  # Init
  success <- count_small_misclass <- count_neg <- error <- 0
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
      ### Rounding (not recommended, it yields slight biases)
      # misclass_nx. <- round(misclass_nx.)
      ### Make interval
      if(output=="n0.") interval <- fieller.interval.n0.(z, t01=(test_n01[i]/test_n0.), t10=(test_n10[i]/test_n1.), test_n0., test_n1., target_n.0[k], target_n.1[k], misclass_n0.=misclass_nx.[1], misclass_n1.=misclass_nx.[2]) 
      else if(output=="n1.") interval <- fieller.interval.n1.(z, t01=(test_n01[i]/test_n0.), t10=(test_n10[i]/test_n1.), test_n0., test_n1., target_n.0[k], target_n.1[k], misclass_n0.=misclass_nx.[1], misclass_n1.=misclass_nx.[2]) 
      else if(output=="n01") interval <- fieller.interval.n01(z, t01=(test_n01[i]/test_n0.), t10=(test_n10[i]/test_n1.), test_n0., test_n1., target_n.0[k], target_n.1[k], misclass_n0.=misclass_nx.[1], misclass_n1.=misclass_nx.[2]) 
      else if(output=="shieh") interval <- fieller.interval.shieh(z, t01=(test_n01[i]/test_n0.), t10=(test_n10[i]/test_n1.), test_n0., test_n1., target_n.0[k], target_n.1[k])
      ### Rounding (not recommended, it yields slight biases)
      # interval <- round(interval)
      # if(output=="n01") misclass_n01 <- round(misclass_n01)
      
      ### Test interval
      if(output=="n01"){ 
        if(min(interval)<=target_n01[k] & max(interval)>=target_n01[k]) success <- success + 1
      } else if(output=="n0.") {
        if(min(interval)<=target_n0. & max(interval)>=target_n0.) success <- success + 1
      } else if(output=="shieh") {
        if(min(interval)<=target_n0./(target_n.0[k]+target_n.1[k]) & max(interval)>=target_n0./(target_n.0[k]+target_n.1[k])) success <- success + 1
      } else {
        if(min(interval)<=target_n1. & max(interval)>=target_n1.) success <- success + 1
      }
    }
  }
  
  ### Return % of intervals containing target_n0. or target_n01
  return(100*success/(nTest*nTarget))
}

######### REPRODUCTION OF PRIOR WORK, from Shieh 2009, Buonaccorsi 2010
### This function concerns intervals describing a general population from which target set is sample. 
### This function concerns estimates of class proportions, e.g., pi0 = n0./n..  
test.fieller.pop <- function(z=1, nTest=100, nTarget=100, population_t01_t10=c(0.1,0.2), test_nx.=c(300,300), pop_nx.=c(300,500)){
  # Init
  success <- count_small_misclass <- count_neg <- error <- 0
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
      interval <- fieller.interval.shieh(z, t01=(test_n01[i]/test_n0.), t10=(test_n10[i]/test_n1.), test_n0., test_n1., target_n.0, target_n.1)
      
      ### Test interval
      if(min(interval)<=pop_nx.[1]/(pop_nx.[1]+pop_nx.[2]) & max(interval)>=pop_nx.[1]/(pop_nx.[1]+pop_nx.[2]) ) success <- success + 1

    }
  }
  
  ### Return % of intervals containing target_n0. or target_n01
  return(100*success/(nTest*nTarget))
}