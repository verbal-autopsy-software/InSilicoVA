#' Get individual COD probabilities from InSilicoVA Model Fits
#' 
#' This function calculates individual probabilities for each death and provide posterior credible intervals for each estimates.
#' 
#' 
#' @param data Original data.
#' @param object Fitted \code{"insilico"} object.
#' @param CI Confidence interval for posterior estimates.
#' @param \dots Not used.
#' \item{mean}{ individual mean COD distribution matrix.} 
#' \item{median}{ individual median COD distribution matrix.} 
#' \item{lower}{ individual lower bound for each COD probability.} 
#' \item{upper}{ individual upper bound for each COD probability.} 

#' @author Zehang Li, Tyler McCormick, Sam Clark
#' 
#' Maintainer: Zehang Li <lizehang@@uw.edu>
#' @seealso \code{\link{insilico}}, \code{\link{plot.insilico}}
#' @references Tyler H. McCormick, Zehang R. Li, Clara Calvert, Amelia C.
#' Crampin, Kathleen Kahn and Samuel J. Clark Probabilistic cause-of-death
#' assignment using verbal autopsies, \emph{arXiv preprint arXiv:1411.3042}
#' \url{http://arxiv.org/abs/1411.3042} (2014)
#' @examples
#' \dontrun{
#' }
#' @export get.indiv
get.indiv <- function(object, CI = 0.95, java_option = "-Xmx1g", ...){
	if(is.null(java_option)) java_option = "-Xmx1g"
	options( java.parameters = java_option )
	id <- object$id
	obj <- .jnew("sampler/InsilicoSampler2")

	# 3d array of csmf to be Nsub * Nitr * C
	if(is.null(object$subpop)){
		csmf <- array(0, dim = c(1, dim(object$csmf)[1], dim(object$csmf)[2]))
		csmf[1, , ] <- object$csmf
		subpop <- rep(as.integer(1), length(object$id))
	}else{
		csmf <- array(0, dim = c(length(object$csmf), dim(object$csmf[[1]])[1], dim(object$csmf[[1]])[2]))
		for(i in 1:length(object$csmf)){
			csmf[i, , ] <- object$csmf[[i]]
		}
		subpop <- as.integer(match(object$subpop, names(object$csmf)))
	}	

	csmf.j <- .jarray(csmf, dispatch = TRUE)
	data.j <- .jarray(object$data, dispatch = TRUE)
	subpop.j <- .jarray(subpop, dispatch = TRUE)
	zero_mat.j <- .jarray(object$zero_mat, dispatch = TRUE)
	# get condprob to be Nitr * S * C array
	if(object$updateCondProb == FALSE){
		condprob <- array(0, dim = c(1, dim(object$probbase)[1], dim(object$probbase)[2]))
		condprob[1, , ] <- object$probbase
	}else{
		if(object$keepProbbase.level){
			condprob <- array(0, dim = c(dim(object$conditional.probs)[1], dim(object$probbase)[1], dim(object$probbase)[2]))
			for(i in 1:dim(condprob)[1]){
				#fix for interVA probbase
				object$conditional.probs[object$conditional.probs == "B -"] <- "B-"
				object$conditional.probs[object$conditional.probs == ""] <- "N"
				
				temp <- object$conditional.probs[i, match(object$probbase, colnames(object$conditional.probs))]
				condprob[i, , ] <- matrix(as.numeric(temp), 
											dim(object$probbase)[1], 
										    dim(object$probbase)[2])
			}
		}else{
			condprob <- object$conditional.probs
		}
	}
	condprob.j <- .jarray(condprob, dispatch = TRUE)

	indiv  <- .jcall(obj, "[[D", "IndivProb", 
					 data.j, csmf.j, subpop.j, condprob.j)

	indiv <- do.call(rbind, lapply(indiv, .jevalArray))
	K <- dim(indiv)[1] / 4
	mean <- indiv[1:K, ]
	median <- indiv[(K+1):(2*K), ]
	lower <- indiv[(2*K+1):(3*K), ]
	upper <- indiv[(3*K+1):(4*K), ]

	## add back all external cause death 41:51 in standard VA
	if(object$external){
		ext.flag <- apply(object$indiv.prob[, 41:51], 1, sum)
		ext.probs <- object$indiv.prob[which(ext.flag == 1), ]
		mean <- rbind(mean, ext.probs)
		median <- rbind(median, ext.probs)
		lower <- rbind(lower, ext.probs)
		upper <- rbind(upper, ext.probs)
	}
	return(list(mean = mean, median = median, 
				lower = lower, upper = upper))
}
