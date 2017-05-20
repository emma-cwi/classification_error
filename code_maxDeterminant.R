source("code_biasCorrection.R")

test.det.var <- function(dataset="Data/5_wave.csv", target_nx.=c(1092,753,455), test_nx.=c(50,50,50), nTest=1000, nTarget=100){
  data <- read.csv(dataset)
  data$true <- factor(data$true)
  classes <- levels(data$true)
  data <- data.table(data)
  nClass <- length(classes)
  
  ### Split Data
  split <- list()
  for(c in 1:nClass) split[[c]] <- data[ true==classes[c] ]
  
  # Make variable names
  colnames_nx. <- colnames_nxy <- c()
  for(x in 1:nClass){
    colnames_nx. <- c(colnames_nx., paste0("n", x, "."))
    for(y in 1:nClass) colnames_nxy <- c(colnames_nxy, paste0("n", x, y))
  }
  
  # Init Result
  test_matrix <- target_matrix <- var_nxy <- var_nx. <- det_misclass <- det_rtp <- c()
  
  ### Iterate Test Set
  for(i in 1:nTest){
    
    # Split Test and Target Sets
    test_set <- target_set <- list()
    for(c in 1:nClass){
      s <- sample(1:nrow( split[[c]] ), test_nx.[c])
      test_set[[c]] <- split[[c]][s]
      target_set[[c]] <- split[[c]][!s]
    }
    
    # Make Test Set Matrix
    test_matrix_temp <- c()
    for(trueClass in 1:nClass){
      for(outputClass in 1:nClass){
        nxy <- nrow(test_set[[trueClass]][ output==classes[outputClass] ])
        test_matrix_temp <- c(test_matrix_temp, nxy)
      }
    }
    test_matrix <- rbind(test_matrix, test_matrix_temp)
    test_matrix_temp <- matrix(test_matrix_temp, ncol=nClass, nrow=nClass)
    
    # Make Matrix Determinant
    det_misclass_temp <- det(misclass.rate(test_matrix_temp))
    det_rtp_temp <- det(rtp.rate(test_matrix_temp))
    det_misclass <- c(det_misclass, det_misclass_temp)
    det_rtp <- c(det_rtp, det_rtp_temp)
    
    if(det_misclass_temp == 0 | det_rtp_temp == 0){
      print("Determinant == 0")
      next
    }
    
    ### Iterate Target Set
    misclass_nxy <- misclass_nx. <- c()
    for(k in 1:nTarget){
      
      ### Sample Target Data
      target_set_temp <- target_set
      for(c in 1:nClass){
        if(target_nx.[c] > nrow(target_set[[c]])) {
          print(paste("Target set too large for class", c) )
          return()
        }
        s <- sample(1:nrow( target_set[[c]] ), target_nx.[c])
        target_set_temp[[c]] <- target_set_temp[[c]][s]
      }
      
      ### Make Target Set Matrix
      target_matrix_temp <- c()
      for(trueClass in 1:nClass){
        for(outputClass in 1:nClass){
          nxy <- nrow(target_set_temp[[trueClass]][ output==classes[outputClass] ])
          target_matrix_temp <- c(target_matrix_temp, nxy)
        }
      }
      target_matrix <- rbind(target_matrix, c(i, target_matrix_temp))
      target_matrix_temp <- matrix( target_matrix_temp, ncol=nClass, nrow=nClass)
      target_n.y <- rowSums(target_matrix_temp) 
      
      ### Get Results
      misclass_nxy <- rbind(misclass_nxy, as.vector( misclass.nxy(test_matrix_temp, target_n.y)) ) 
      misclass_nx. <- rbind(misclass_nx., misclass.nx.(test_matrix_temp, target_n.y))
    }
    
    ### Get Variance
    var_nxy <- rbind(var_nxy, colVar(misclass_nxy))
    var_nx. <- rbind(var_nx., colVar(misclass_nx.))
  }
  
  res_nxy <- cbind(det_misclass, det_rtp, var_nxy)
  colnames(res_nxy) <- c("Det.Misclass", "Det.RTP", colnames_nxy)
  
  res_nx. <- cbind(det_misclass, det_rtp, var_nx.)
  colnames(res_nx.) <- c("Det.Misclass", "Det.RTP", colnames_nx.)
  
  res_sum <- cbind(det_misclass, det_rtp, rowSums(var_nx.), rowSums(var_nxy)) 
  colnames(res_sum) <- c("Det.Misclass", "Det.RTP", "Sum.Var.nx.", "Sum.Var.nxy")
  
  colnames(test_matrix) <- colnames_nxy
  colnames(target_matrix) <- c("Iteration", colnames_nxy)
  res <- list( target_matrix=target_matrix, test_matrix=test_matrix, 
               res_nx.=res_nx., res_nxy=res_nxy, res_sum=res_sum )
  
  return(res)
  
}

colVar <- function(data){
  res <- c()
  for(i in 1:ncol(data)){
    res <- c(res, var(data[,i]))
  }
  return(res)
}

res_varplot <- function(res, name=deparse(substitute(res))){
  plot(res$res_sum[,c(1,3)], main=name )
  plot(res$res_sum[,c(1,4)], main=name )
  plot(res$res_sum[,c(2,3)], main=name )
  plot(res$res_sum[,c(2,4)], main=name )
}

res_cordet <- function(res){
  print("Correlation Misclass Determinant and Var(nx.), Var(nxy)")
  print( c(cor(res$res_sum[,c(1,3)])[2,1] , cor(res$res_sum[,c(1,4)])[2,1] ) )
  print("Correlation RTP Determinant and Var(nx.), Var(nxy)")
  print( c( cor(res$res_sum[,c(2,3)])[2,1] , cor(res$res_sum[,c(2,4)])[2,1] ) )
}


