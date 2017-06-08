# p <- 0.5
# n <- 1000
# m <- 1000
# rbinom(10,n,p)/n

testDiff_cauchy <- function(p=0.5,
						#nmList=c(c(1:4)*20, c(1:10)*100, c(1:5)*10000), 
						nmList = c(1:5)*100,
						nRun=100){
	# Instantiate Population
	N <- 1000000
	n1 <- p*N
	n0 <- (1-p)*N
	data <- cbind(c(1:N),c(rep(0,n0),rep(1,n1)))

	# Fill in cells of results matrices
	for(i in 1:length(nmList)){
		m <- nmList[i]
		for(j in 1:length(nmList)){
			n <- nmList[j]

			# List of difference r1-r2 in confidence interval
			diff_list <- c()
			for(k1 in 1:nRun){
				# Sample r1 of n items
				sn <- sample(1:N, n, replace=TRUE)
				data_n <- data[sn,]
				data_n_1 <- sum(data_n[,2])
				data_n_0 <- n - sum(data_n[,2])
				# Draw rate
				r1 <- data_n_1/data_n_0

				for(k2 in 1:nRun){
					# Sample r2 of m items
					sm <- sample(1:N, m, replace=TRUE)
					data_m <- data[sm,]
					data_m_1 <- sum(data_m[,2])
					data_m_0 <- m - sum(data_m[,2])
					# Draw rate
					r2 <- data_m_1/data_m_0

					# Test Diff
					diff_list <- c(diff_list, r1-r2) 
				}			
			}

			pdf(paste("diff_cauchy","_p",p*100,"_n",nmList[j],"_m",nmList[i],".pdf", sep=""))
			hist(diff_list, ylab="Frequency", xlab="r1 - r2", breaks=20) # ,xlim=c(-0.005,0.005)) 
			dev.off()
		}
	}
}

nmMatrix <- function(   p=0.5,
						nmList=c(c(1:4)*20, c(1:10)*100, c(1:5)*10000), 
						nRun=100){
	# Instantiate Population
	n1 <- p*N
	n0 <- (1-p)*N
	data <- cbind(c(1:N),c(rep(0,n0),rep(1,n1)))

	# Instantiate matrices of results (filled with NA)
	# Row : Sample r2 of m items
	# Col : Sample r1 of n items
	matrix <- matrix(,ncol=19, nrow=19)
	rownames(m) <- colnames(m) <- nmList
	m_perle <- m_cauchy <- matrix

	# Fill in cells of results matrices
	for(i in 1:length(nmList)){
		m <- nmList[i]
		for(j in 1:length(nmList)){
			n <- nmList[j]

			# List of actual % of r2 in confidence interval
			eval_perle <- eval_cauchy <- eval_ps <- eval_sp <- c()
			for(k1 in 1:nRun){
				# Sample r1 of n items
				sn <- sample(1:N, n, replace=TRUE)
				data_n <- data[sn,]
				data_n_1 <- sum(data_n[,2])
				data_n_0 <- 1 - sum(data_n[,2])
				# Draw rates
				r1_perle <- data_n_1/n
				r1_cauchy <- data_n_1/data_n_0
				# Draw variance
				var1_perle <- r1_perle*(1-r1_perle)/n

				for(k2 in 1:nRun){
					# Sample r2 of m items
					sn <- sample(1:N, m, replace=TRUE)
					data_m <- data[sm,]
					data_m_1 <- sum(data_m[,2])
					data_m_0 <- 1 - sum(data_m[,2])
					# Draw rates
					r2_perle <- data_m_1/m
					r2_cauchy <- data_m_1/data_m_0
					# Draw variance
					var2_perle <- r1_perle*(1-r1_perle)/n

					# Test PERLE

				}
			}
		}
	}

	# Write matrices of results
	write.csv(m_perle, "CI_perle_simple.csv")
	write.csv(m_cauchy, "CI_perle_all.csv")
	write.csv(m_ds, "CI_double_sampling.csv")
	write.csv(m_ps, "CI_pop_to_sample.csv")
	write.csv(m_sp, "CI_sample_to_pop.csv")
}