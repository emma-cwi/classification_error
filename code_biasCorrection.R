library(data.table)
library(MASS)

################# BIAS CORRECTION ################# 
### Perle
perle.rate <- function(confusion_matrix){
  return( t( t(confusion_matrix)/colSums(confusion_matrix) ) )
}
perle.count <- function(perle_rate, output_count){
  return(as.vector( round(solve(perle_rate) %*% as.matrix(output_count)) ))
}
perle <- function(confusion_matrix, output_count){
  return( perle.count(perle.rate(confusion_matrix), output_count ) )
}
perle.all <- function(confusion_matrix, output_count){
  return(round( t( t(perle.rate(confusion_matrix))*perle(confusion_matrix, output_count)) ))
}

### Double Sampling
ds.rate <- function(confusion_matrix){
  return( confusion_matrix/rowSums(confusion_matrix) )
}
ds.count <- function(ds_rate, output_count){
  return(as.vector( round(t(ds_rate) %*% as.matrix(output_count)) ))
}
ds <- function(confusion_matrix, output_count){
  return( ds.count(ds.rate(confusion_matrix), output_count ) )
}
ds.all <- function(confusion_matrix, output_count){
  return(round( ds.rate(confusion_matrix)*output_count ))
}

### RTP
rtp.rate <- function(confusion_matrix){
  return( t( t(confusion_matrix)/diag(confusion_matrix) ) )
}
rtp.tp <- function(rtp_rate, output_count){
  return(as.vector( round(solve(rtp_rate) %*% as.matrix(output_count))  ))
}
rtp <- function(rtp_rate, rtp_tp){
  return( round( t( t(rtp_rate)*rtp_tp) ) )
}
rtp.all <- function(confusion_matrix, output_count){
  rtp_rate <- rtp.rate(confusion_matrix)
  return( rtp(rtp_rate, rtp.tp(rtp_rate, output_count)) )
}

### Combi
combi.all <- function(confusion_matrix, output_count){
  test_nx. <- colSums(confusion_matrix)
  perle_nx. <- perle(confusion_matrix, output_count)
  ratio <- perle_nx./test_nx.
  ds_matrix <- t(t(confusion_matrix)*ratio)
  return( ds.all(ds_matrix, output_count) )
}

################# REAL DATA TEST ################# 
test.real.data <- function(dataset="Data/5_wave.csv", target_nx.=c(1092,753,455), test_nx.=c(), nTest=100, 
                           do_perle=TRUE, do_ds=TRUE, do_combi=TRUE, do_rtp=TRUE, do_nx.=TRUE, do_nxy=TRUE, do_check=FALSE){
  data <- read.csv(dataset)
  data$true <- factor(data$true)
  classes <- levels(data$true)
  data <- data.table(data)
  nClass <- length(classes)
  
  ### Split Data
  split <- list()
  for(c in 1:nClass) split[[c]] <- data[ true==classes[c] ]
  
  # Make variable names
  colnames_nx. <- colnames_n.y <- colnames_nxy <- c()
  for(x in 1:nClass){
    colnames_nx. <- c(colnames_nx., paste0("n", x, "."))
    colnames_n.y <- c(colnames_n.y, paste0("n.", x))
    for(y in 1:nClass) colnames_nxy <- c(colnames_nxy, paste0("n", x, y))
  }
  
  ### Init Iterations
  target_matrix <- test_matrix <- target_n.y <- target_n.y_error <- c()
  if(do_perle) perle_nxy <- perle_nxy_error <- perle_nx. <- perle_nx._error <- perle_check <- c()
  if(do_ds) ds_nxy <- ds_nxy_error <- ds_nx. <- ds_nx._error <- ds_check <- c()
  if(do_combi) combi_nxy <- combi_nxy_error <- combi_nx. <- combi_nx._error <- combi_check <- c()
  if(do_rtp) rtp_nxy <- rtp_nxy_error <- rtp_nx. <- rtp_nx._error <- rtp_check <- c()
  
  ### Iterate Test
  for(i in 1:nTest){
    
    ### Sample Data
    test_set <- target_set <- list()
    for(c in 1:nClass){
      s <- sample(1:nrow( split[[c]] ), target_nx.[c])
      target_set[[c]] <- split[[c]][s]
      test_set[[c]] <- split[[c]][!s]
    }
    if(!is.null(test_nx.[c])){
      for(c in 1:nClass){
        if(test_nx.[c] > nrow(test_set[[c]])) {
          print(paste("Test set too large for class", c) )
          return()
        }
        s <- sample(1:nrow( test_set[[c]] ), test_nx.[c])
        test_set[[c]] <- test_set[[c]][s]
      }
    }
    
    ### Count Errors
    test_matrix_temp <- target_matrix_temp <- c()
    # Test set
    for(trueClass in 1:nClass){
      for(outputClass in 1:nClass){
        nxy <- nrow(test_set[[trueClass]][ output==classes[outputClass] ])
        test_matrix_temp <- c(test_matrix_temp, nxy)
      }
    }
    test_matrix <- rbind(test_matrix, test_matrix_temp)
    test_matrix_temp <- matrix(test_matrix_temp, ncol=nClass, nrow=nClass)
    # Target Set
    for(trueClass in 1:nClass){
      for(outputClass in 1:nClass){
        nxy <- nrow(target_set[[trueClass]][ output==classes[outputClass] ])
        target_matrix_temp <- c(target_matrix_temp, nxy)
      }
    }
    target_matrix <- rbind(target_matrix, target_matrix_temp)
    target_matrix_temp <- matrix( target_matrix_temp, ncol=nClass, nrow=nClass)
    target_n.y <- rbind(target_n.y, rowSums(target_matrix_temp) )
    
    # Get Results - nxy
    if(do_nxy){
      if(do_perle) {
        try( perle_nxy <- rbind(perle_nxy, as.vector( perle.all(test_matrix_temp, target_n.y[i,]) ) ) )
        try( perle_nxy_error <- rbind(perle_nxy_error, perle_nxy[i,] - target_matrix[i,] ) )
      }
      if(do_ds) {
        ds_nxy <- rbind(ds_nxy, as.vector( ds.all(test_matrix_temp, target_n.y[i,]) ) )
        ds_nxy_error <- rbind(ds_nxy_error, ds_nxy[i,] - target_matrix[i,] ) 
      }
      if(do_combi) {
        try( combi_nxy <- rbind(combi_nxy, as.vector( combi.all(test_matrix_temp, target_n.y[i,]) ) ) )
        try( combi_nxy_error <- rbind(combi_nxy_error, combi_nxy[i,] - target_matrix[i,]) )
      }
      if(do_rtp) {
        try( rtp_nxy <- rbind(rtp_nxy, as.vector( rtp.all(test_matrix_temp, target_n.y[i,]) ) ) )
        try( rtp_nxy_error <- rbind(rtp_nxy_error, rtp_nxy[i,] - target_matrix[i,] ) )
      }
    }
    
    # Get Results - nx.
    if(do_nx.){
      if(do_perle) {
        try( perle_nx. <- rbind(perle_nx., perle(test_matrix_temp, target_n.y[i,]) ) )
        try( perle_nx._error <- rbind(perle_nx._error, perle_nx.[i,] - target_nx. ) )
      }
      if(do_ds) {
        ds_nx. <- rbind(ds_nx., ds(test_matrix_temp, target_n.y[i,]) )
        ds_nx._error <- rbind(ds_nx._error, ds_nx.[i,] - target_nx. )
      }
      if(do_combi) {
        try( combi_nx. <- rbind(combi_nx., colSums(matrix( combi_nxy[i,], ncol=nClass, nrow=nClass) ) ) )
        try( combi_nx._error <- rbind(combi_nx._error, combi_nx.[i,] - target_nx. ) )
      }
      if(do_rtp) {
        try( rtp_nx. <- rbind(rtp_nx., colSums(matrix( rtp_nxy[i,], ncol=nClass, nrow=nClass) ) ) )
        try( rtp_nx._error <- rbind(rtp_nx._error, rtp_nx.[i,] - target_nx. ) )
      }
      target_n.y_error <- rbind(target_n.y_error, target_n.y[i,] - target_nx. ) 
    }
    
    # Check Code
    if(do_check){
      if(do_perle) try( perle_check <- rbind(perle_check, perle(target_matrix_temp, target_n.y[i,]) ) )
      if(do_ds) ds_check <- rbind(ds_check,  ds(target_matrix_temp, target_n.y[i,]) )
      if(do_combi) try( combi_check <- rbind(combi_check, colSums(combi.all(target_matrix_temp, target_n.y[i,])) ) )
      if(do_rtp) try( rtp_check <- rbind(rtp_check, colSums(rtp.all(target_matrix_temp, target_n.y[i,])) ) )
    }
  }
  
  # Colnames
  target_nx. <- rbind(target_nx.)
  colnames(target_nx.) <- colnames(target_n.y_error) <- colnames_nx.
  colnames(target_matrix) <- colnames(test_matrix) <- colnames_nxy
  colnames(target_n.y) <- colnames_n.y
  res <- list( target_matrix=target_matrix, test_matrix=test_matrix, 
               target_n.y=target_n.y, target_nx.=target_nx., target_n.y_error=target_n.y_error )
  
  if(do_perle) {
    if(do_nxy) {
      colnames(perle_nxy) <- colnames(perle_nxy_error) <- colnames_nxy
      res <- append(res, list( perle_nxy=perle_nxy, perle_nxy_error=perle_nxy_error))
    }
    if(do_nx.) {
      colnames(perle_nx.) <- colnames(perle_nx._error) <- colnames_nx.
      res <- append(res, list( perle_nx.=perle_nx., perle_nx._error=perle_nx._error))
    }
    if(do_check) {
      colnames(perle_check) <- colnames_nx.
      res <- append(res, list( perle_check=perle_check ))
    }
  }
  if(do_ds) {
    if(do_nxy) {
      colnames(ds_nxy) <- colnames(ds_nxy_error) <- colnames_nxy
      res <- append(res, list( ds_nxy=ds_nxy, ds_nxy_error=ds_nxy_error))
    }
    if(do_nx.) {
      colnames(ds_nx.) <- colnames(ds_nx._error) <- colnames_nx.
      res <- append(res, list( ds_nx.=ds_nx., ds_nx._error=ds_nx._error))
    }
    if(do_check) {
      colnames(ds_check) <- colnames_nx.
      res <- append(res, list( ds_check=ds_check ))
    }
  }
  if(do_combi) {
    if(do_nxy) {
      colnames(combi_nxy) <- colnames(combi_nxy_error) <- colnames_nxy
      res <- append(res, list( combi_nxy=combi_nxy, combi_nxy_error=combi_nxy_error))
    }
    if(do_nx.) {
      colnames(combi_nx.) <- colnames(combi_nx._error) <- colnames_nx.
      res <- append(res, list( combi_nx.=combi_nx., combi_nx._error=combi_nx._error))
    }
    if(do_check) {
      colnames(combi_check) <- colnames_nx.
      res <- append(res, list( combi_check=combi_check ))
    }
  }
  if(do_rtp) {
    if(do_nxy) {
      colnames(rtp_nxy) <- colnames(rtp_nxy_error) <- colnames_nxy
      res <- append(res, list( rtp_nxy=rtp_nxy, rtp_nxy_error=rtp_nxy_error))
    }
    if(do_nx.) {
      colnames(rtp_nx.) <- colnames(rtp_nx._error) <- colnames_nx.
      res <- append(res, list( rtp_nx.=rtp_nx., rtp_nx._error=rtp_nx._error))
    }
    if(do_check) {
      colnames(rtp_check) <- colnames_nx.
      res <- append(res, list( rtp_check=rtp_check ))
    }
  }
  return(res)
}



################ PLOT ##################
get_lim_nxy_error <- function(res){
  # Init
  max <- max(res$target_n.y_error)
  min <- min(res$target_n.y_error)
  # Browse
  if(!is.null(res$perle_nxy_error)) {
    max <- max(max, res$perle_nxy_error)
    min <- min(min, res$perle_nxy_error)
  }
  if(!is.null(res$ds_nxy_error)) {
    max <- max(max, res$ds_nxy_error)
    min <- min(min, res$ds_nxy_error)
  }
  if(!is.null(res$combi_nxy_error)) {
    max <- max(max, res$combi_nxy_error)
    min <- min(min, res$combi_nxy_error)
  }
  if(!is.null(res$rtp_nxy_error)) {
    max <- max(max, res$rtp_nxy_error)
    min <- min(min, res$rtp_nxy_error)
  }
  return(c(min, max))
}
get_lim_nxy <- function(res){
  # Init
  max <- max(res$target_n.y)
  min <- min(res$target_n.y)
  # Browse
  if(!is.null(res$perle_nxy)) {
    max <- max(max, res$perle_nxy)
    min <- min(min, res$perle_nxy)
  }
  if(!is.null(res$ds_nxy)) {
    max <- max(max, res$ds_nxy)
    min <- min(min, res$ds_nxy)
  }
  if(!is.null(res$combi_nxy)) {
    max <- max(max, res$combi_nxy)
    min <- min(min, res$combi_nxy)
  }
  if(!is.null(res$rtp_nxy)) {
    max <- max(max, res$rtp_nxy)
    min <- min(min, res$rtp_nxy)
  }
  return(c(min, max))
}
get_lim_nx._error <- function(res){
  # Init
  max <- max(res$target_n.y_error)
  min <- min(res$target_n.y_error)
  # Browse
  if(!is.null(res$perle_nx._error)) {
    max <- max(max, res$perle_nx._error)
    min <- min(min, res$perle_nx._error)
  }
  if(!is.null(res$ds_nx._error)) {
    max <- max(max, res$ds_nx._error)
    min <- min(min, res$ds_nx._error)
  }
  if(!is.null(res$combi_nx._error)) {
    max <- max(max, res$combi_nx._error)
    min <- min(min, res$combi_nx._error)
  }
  if(!is.null(res$rtp_nx._error)) {
    max <- max(max, res$rtp_nx._error)
    min <- min(min, res$rtp_nx._error)
  }
  return(c(min, max))
}
get_lim_nx. <- function(res){
  # Init
  max <- max(res$target_n.y)
  min <- min(res$target_n.y)
  # Browse
  if(!is.null(res$perle_nx.)) {
    max <- max(max, res$perle_nx.)
    min <- min(min, res$perle_nx.)
  }
  if(!is.null(res$ds_nx.)) {
    max <- max(max, res$ds_nx.)
    min <- min(min, res$ds_nx.)
  }
  if(!is.null(res$combi_nx.)) {
    max <- max(max, res$combi_nx.)
    min <- min(min, res$combi_nx.)
  }
  if(!is.null(res$rtp_nx.)) {
    max <- max(max, res$rtp_nx.)
    min <- min(min, res$rtp_nx.)
  }
  return(c(min, max))
}

res_boxplot <- function(res, name=deparse(substitute(res)), do_perle=TRUE, do_ds=TRUE, do_combi=FALSE, do_rtp=FALSE, do_nx.=TRUE, do_nxy=TRUE, do_nxx=FALSE, do_check=FALSE){
  if(do_nxy){
    lim <- get_lim_nxy(res)
    if(do_perle) boxplot(res$perle_nxy, ylim=lim, main=paste(name, '- PERLE nxy'))
    if(do_ds) boxplot(res$ds_nxy, ylim=lim, main=paste(name, '- DS nxy'))
    if(do_combi) boxplot(res$combi_nxy, ylim=lim, main=paste(name, '- COMBI nxy'))
    if(do_rtp) boxplot(res$rtp_nxy, ylim=lim, main=paste(name, '- RTP nxy'))
    boxplot(res$target_matrix, ylim=lim, main=paste(name, '- RAW nxy'))
  } 
  if(do_nx.){
    lim <- get_lim_nx.(res)
    if(do_perle) boxplot(res$perle_nx., ylim=lim, main=paste(name, '- PERLE nx.'))
    if(do_ds) boxplot(res$ds_nx., ylim=lim, main=paste(name, '- DS nx.'))
    if(do_combi) boxplot(res$combi_nx., ylim=lim, main=paste(name, '- COMBI nx.'))
    if(do_rtp) boxplot(res$rtp_nx., ylim=lim, main=paste(name, '- RTP nx.'))
    boxplot(res$target_n.y, ylim=lim, main=paste(name, '- RAW n.y'))
  } 
  if(do_nxx){
    lim <- get_lim_nxy(res)
    n_class <- ncol(res$target_nx.)
    index <- c(1:n_class)+n_class*(c(1:n_class)-1)
    if(do_perle) boxplot(res$perle_nxy[,index], ylim=lim, main=paste(name, '- PERLE nxx'))
    if(do_ds) boxplot(res$ds_nxy[,index], ylim=lim, main=paste(name, '- DS nxx'))
    if(do_combi) boxplot(res$combi_nxy[,index], ylim=lim, main=paste(name, '- COMBI nxx'))
    if(do_rtp) boxplot(res$rtp_nxy[,index], ylim=lim, main=paste(name, '- RTP nxx'))
    boxplot(res$target_matrix[,index], ylim=lim, main=paste(name, '- RAW nxx'))
  }  
  if(do_check){
    lim <- c(min(res$target_nx.), max(res$target_nx.))
    if(do_perle) boxplot(res$perle_check, ylim=lim, main=paste(name, '- PERLE CHECK nx.'))
    if(do_ds) boxplot(res$ds_check, ylim=lim, main=paste(name, '- DS CHECK nx.'))
    if(do_combi) boxplot(res$combi_check, ylim=lim, main=paste(name, '- COMBI CHECK nx.'))
    if(do_rtp) boxplot(res$rtp_check, ylim=lim, main=paste(name, '- RTP CHECK nx.'))
  } 
}
res_boxplot_error <- function(res, name=deparse(substitute(res)), do_perle=TRUE, do_ds=TRUE, do_combi=FALSE, do_rtp=FALSE, do_nx.=TRUE, do_nxy=TRUE, do_nxx=FALSE, do_check=FALSE){
  if(do_nxy){
    lim <- get_lim_nxy_error(res)
    if(do_perle) boxplot(res$perle_nxy_error, ylim=lim, main=paste(name, '- PERLE nxy Error'))
    if(do_ds) boxplot(res$ds_nxy_error, ylim=lim, main=paste(name, '- DS nxy Error'))
    if(do_combi) boxplot(res$combi_nxy_error, ylim=lim, main=paste(name, '- COMBI nxy Error'))
    if(do_rtp) boxplot(res$rtp_nxy_error, ylim=lim, main=paste(name, '- RTP nxy Error'))
  } 
  if(do_nx.){
    lim <- get_lim_nx._error(res)
    if(do_perle) boxplot(res$perle_nx._error, ylim=lim, main=paste(name, '- PERLE nx. Error'))
    if(do_ds) boxplot(res$ds_nx._error, ylim=lim, main=paste(name, '- DS nx. Error'))
    if(do_combi) boxplot(res$combi_nx._error, ylim=lim, main=paste(name, '- COMBI nx. Error'))
    if(do_rtp) boxplot(res$rtp_nx._error, ylim=lim, main=paste(name, '- RTP nx. Error'))
  } 
  if(do_nxx){
    lim <- get_lim_nxy_error(res)
    n_class <- ncol(res$target_nx.)
    index <- c(1:n_class)+n_class*(c(1:n_class)-1)
    if(do_perle) boxplot(res$perle_nxy_error[,index], ylim=lim, main=paste(name, '- PERLE nxx Error'))
    if(do_ds) boxplot(res$ds_nxy_error[,index], ylim=lim, main=paste(name, '- DS nxx Error'))
    if(do_combi) boxplot(res$combi_nxy_error[,index], ylim=lim, main=paste(name, '- COMBI nxx Error'))
    if(do_rtp) boxplot(res$rtp_nxy_error[,index], ylim=lim, main=paste(name, '- RTP nxx Error'))
  }  
  if(do_check){
    lim <- c(min(res$target_nx.), max(res$target_nx.))
    if(do_perle) boxplot(res$perle_check, ylim=lim, main=paste(name, '- PERLE CHECK nx.'))
    if(do_ds) boxplot(res$ds_check, ylim=lim, main=paste(name, '- DS CHECK nx.'))
    if(do_combi) boxplot(res$combi_check, ylim=lim, main=paste(name, '- COMBI CHECK nx.'))
    if(do_rtp) boxplot(res$rtp_check, ylim=lim, main=paste(name, '- RTP CHECK nx.'))
  } 
}

res_parplot_error <- function(res, name=deparse(substitute(res)), do_perle=TRUE, do_ds=TRUE, do_combi=TRUE, do_rtp=FALSE, do_nx.=TRUE, do_nxy=TRUE, do_nxx=FALSE, do_check=FALSE){
  if(do_nxy){
    if(do_perle) parcoord(res$perle_nxy_error, var.label=TRUE, main=paste(name, '- PERLE nxy Error'))
    if(do_ds) parcoord(res$ds_nxy_error, var.label=TRUE, main=paste(name, '- DS nxy Error'))
    if(do_combi) parcoord(res$combi_nxy_error, var.label=TRUE, main=paste(name, '- COMBI nxy Error'))
    if(do_rtp) parcoord(res$rtp_nxy_error, var.label=TRUE, main=paste(name, '- RTP nxy Error'))
  } 
  if(do_nx.){
    if(do_perle) parcoord(res$perle_nx._error, var.label=TRUE, main=paste(name, '- PERLE nx. Error'))
    if(do_ds) parcoord(res$ds_nx._error, var.label=TRUE, main=paste(name, '- DS nx. Error'))
    if(do_combi) parcoord(res$combi_nx._error, var.label=TRUE, main=paste(name, '- COMBI nx. Error'))
    if(do_rtp) parcoord(res$rtp_nx._error, var.label=TRUE, main=paste(name, '- RTP nx. Error'))
  } 
  if(do_nxx){
    n_class <- ncol(res$target_nx.)
    index <- c(1:n_class)+n_class*(c(1:n_class)-1)
    if(do_perle) parcoord(res$perle_nxy_error[,index], var.label=TRUE, main=paste(name, '- PERLE nxx Error'))
    if(do_ds) parcoord(res$ds_nxy_error[,index], var.label=TRUE, main=paste(name, '- DS nxx Error'))
    if(do_combi) parcoord(res$combi_nxy_error[,index], var.label=TRUE, main=paste(name, '- COMBI nxx Error'))
    if(dp_rtp) parcoord(res$rtp_nxy_error[,index], var.label=TRUE, main=paste(name, '- RTP nxx Error'))
  }  
  if(do_check){
    if(do_perle) parcoord(res$perle_check, var.label=TRUE, main=paste(name, '- PERLE CHECK nx.'))
    if(do_ds) parcoord(res$ds_check, var.label=TRUE, main=paste(name, '- DS CHECK nx.'))
    if(do_combi) parcoord(res$combi_check, var.label=TRUE, main=paste(name, '- COMBI CHECK nx.'))
    if(do_rtp) parcoord(res$rtp_check, var.label=TRUE, main=paste(name, '- RTP CHECK nx.'))
  } 
  parcoord(res$target_n.y_error, var.label=TRUE, main=paste(name, '- RAW n.y Error'))
}
res_parplot <- function(res, name=deparse(substitute(res)), do_perle=TRUE, do_ds=TRUE, do_combi=TRUE, do_rtp=FALSE, do_nx.=TRUE, do_nxy=TRUE, do_nxx=FALSE, do_check=FALSE){
  if(do_nxy){
    if(do_perle) parcoord(res$perle_nxy, col=rainbow(nrow(res$perle_nxy)), var.label=TRUE, main=paste(name, '- PERLE nxy'))
    if(do_ds) parcoord(res$ds_nxy, var.label=TRUE, main=paste(name, '- DS nxy'))
    if(do_combi) parcoord(res$combi_nxy, var.label=TRUE, main=paste(name, '- COMBI nxy'))
    if(do_rtp) parcoord(res$rtp_nxy, var.label=TRUE, main=paste(name, '- RTP nxy'))
  } 
  if(do_nx.){
    if(do_perle) parcoord(res$perle_nx., var.label=TRUE, main=paste(name, '- PERLE nx.'))
    if(do_ds) parcoord(res$ds_nx., var.label=TRUE, main=paste(name, '- DS nx.'))
    if(do_combi) parcoord(res$combi_nx., var.label=TRUE, main=paste(name, '- COMBI nx.'))
    if(do_rtp) parcoord(res$rtp_nx., var.label=TRUE, main=paste(name, '- RTP nx.'))
  } 
  if(do_nxx){
    n_class <- ncol(res$target_nx.)
    index <- c(1:n_class)+n_class*(c(1:n_class)-1)
    if(do_perle) parcoord(res$perle_nxy[,index], var.label=TRUE, main=paste(name, '- PERLE nxx'))
    if(do_ds) parcoord(res$ds_nxy[,index], var.label=TRUE, main=paste(name, '- DS nxx'))
    if(do_combi) parcoord(res$combi_nxy[,index], var.label=TRUE, main=paste(name, '- COMBI nxx'))
    if(dp_rtp) parcoord(res$rtp_nxy[,index], var.label=TRUE, main=paste(name, '- RTP nxx'))
  }  
  if(do_check){
    if(do_perle) parcoord(res$perle_check, var.label=TRUE, main=paste(name, '- PERLE CHECK nx.'))
    if(do_ds) parcoord(res$ds_check, var.label=TRUE, main=paste(name, '- DS CHECK nx.'))
    if(do_combi) parcoord(res$combi_check, var.label=TRUE, main=paste(name, '- COMBI CHECK nx.'))
    if(do_rtp) parcoord(res$rtp_check, var.label=TRUE, main=paste(name, '- RTP CHECK nx.'))
  } 
  parcoord(res$target_n.y, var.label=TRUE, main=paste(name, '- RAW n.y'))
}

