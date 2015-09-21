##
## function to check if the csmf chains have converged
##
## csmf     : the csmf to be checked, either vector or list (with sub-population)
## conv.csmf: the minimum mean csmf to check
## test     : type of test
## verbose  : whether to return the test statistics, or simply outcome

csmf.diag <- function(csmf, conv.csmf = 0.02, test= c("gelman", "heidel")[2], verbose = TRUE, autoburnin = FALSE, external.sep = TRUE, ...){
	
	external.causes = seq(41, 51)
	
	# check if the input is insilico object or only csmf
	if(class(csmf) == "insilico"){
		if(external.sep){
			csmf <- csmf[, -external.causes]
		}else{
			csmf <- csmf$csmf
		}
	}  
	if(class(csmf) == "list" && class(csmf[[1]]) == "insilico"){
		csmf_only <- vector("list", length(csmf))
		for(i in 1:length(csmf)){
			if(external.sep){
				csmf_only[[i]] <- csmf[[i]]$csmf[, -external.causes]
			}else{
				csmf_only[[i]] <- csmf[[i]]$csmf
			}
		}
		csmf <- csmf_only
	}
	
	# check conv.csmf value
	if(conv.csmf >= 1 || conv.csmf <= 0) {
		stop("Invalid threshold.")
	}
	if(test == "gelman" && class(csmf) != "list"){
		stop("Need more than one chain to perform Gelman and Rubin's test.")
	}
	gelman.single <- function(list, conv.csmf){
		# initialize the mean csmf for all chains
		mean <- rep(0, dim(list[[1]])[2])
		# in case chains have different lengths, count the normalizing constant
		Z <- 0
		for(i in 1:length(list)){
			# add sums
			mean <- mean + apply(list[[i]], 2, sum)
			# add counts
			Z <- Z + dim(list[[i]])[1]
			# change to mcmc object
			list[[i]] <- mcmc(list[[i]])
		}
		# sort and find which causes to report
		mean <- mean / Z
		which <- which(mean > conv.csmf)
		order <- order(mean[which], decreasing = TRUE)
		for(i in 1:length(list)){
			list[[i]] <- list[[i]][, which[order]]
		}
		# get the mcmc list
		multi <- mcmc.list(list)
		test <- gelman.diag(multi, autoburnin=FALSE)
		return(test)
	}
	# function to perform heidel test on one chain
	heidel.single <- function(one, conv.csmf){
		one <- as.matrix(one)
		# sort chain by CSMF mean values
		mean <- apply(one, 2, mean)
		one <- one[, order(mean, decreasing = TRUE)]
		# remove CSMFs below conv.csmf
		mean <- apply(one, 2, mean)
		# use drop = FALSE to filter without making one an vector
		one <- one[, which(mean > conv.csmf), drop = FALSE]
		one <- mcmc(one)
		test <- heidel.diag(one)
		return(test)
	}
	# Heidel test
	if(test == "heidel"){
		testout <- NULL
		conv <- 1
		if(class(csmf) == "list"){
			testout <- vector("list", length(csmf))
			names(testout) <- names(csmf)
			for(i in 1:length(testout)){
				testout[[i]] <- heidel.single(csmf[[i]], conv.csmf)
				conv <- conv * prod(testout[[i]][, 1]) * prod(testout[[i]][, 4])
			}
		}else{
			testout <- heidel.single(csmf, conv.csmf)
			conv <- prod(testout[, 1]) * prod(testout[, 4])
		}
	}
	# Gelman test
	if(test == "gelman"){
		testout <- NULL
		conv <- 1
		if(class(csmf[[1]]) == "list"){
			testout <- vector("list", length(csmf[[1]]))
			names(testout) <- names(csmf[[1]])
			for(i in 1:length(testout)){
				# create new list for each sub population
				csmf_sub <- vector("list", length(csmf))
				for(j in 1:length(csmf)){
					csmf_sub[[j]] <- csmf[[j]][[i]]
				}
				testout[[i]] <- gelman.single(csmf_sub, conv.csmf)
			}
		}else{
			testout <- gelman.single(csmf, conv.csmf)
		}
	}
	
	if(is.na(conv)) conv <- 0
	if(!verbose && test == "heidel"){
		return(conv)
	}else{
		return(testout)
	}
}