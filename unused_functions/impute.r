
impute <- function(indic.w.missing, y.new, cond.prob, S, subbelong, missing.imp){
###########################################################
# function to impute missing data (NOT USED)
# @param:
#		indic.w.missing : matrix of indicators
#		y.new           : current cause vector
#		cond.prob       : cond prob matrix
#		S 				: number of symptoms
#		subbelong       : vector of sub-population membership
#		missing.imp     : as in the main function	
# @values:
# 		imputed indic matrix		
	## if no imputing
	if(missing.imp == "none"){
		return(indic.w.missing)

	## if impute for all missings
	}else if(missing.imp == "all"){
		## arrange indic into vectors, one death after another
		indic.w.missing <- t(indic.w.missing)
		## find where are the missing values
		where.miss <- which(indic.w.missing < 0)
	
	# if impute for all missings except sub-population specific total missing
	}else if(missing.imp == "sub"){
		## arrange indic into vectors, one death after another
		indic.w.missing <- t(indic.w.missing)
		## find where are the missing values
		where.miss <- which(indic.w.missing == -1)
	}

	if(length(where.miss)  == 0) return(indic.w.missing)
	## find the causes for each of the missing values
	what.cause <- y.new[trunc((where.miss-0.1)/S) + 1]
	## find the corresponding cond.prob
	getprob <- function(i, where.miss, what.cause, cond.prob, S){
		which.symp <- where.miss[i] %% S
		if(which.symp == 0) which.symp <- S
		cond.prob[which.symp, what.cause[i]]
	}
	prob <- sapply(seq(1:length(where.miss)), getprob, 
					where.miss, what.cause, cond.prob, S)
	## sample new indic.vec
	samp <- rbinom(n = length(where.miss), size = 1, 
					prob = prob)
	indic.w.missing[where.miss] <- samp
	return(t(indic.w.missing))
}