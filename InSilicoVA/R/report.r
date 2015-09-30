#' Summarizing InSilicoVA Model Fits
#' 
#' This function is the summary method for class \code{insilico}.
#' 
#' \code{summary.insilico} formats some basic information about the InSilicoVA
#' fitted object on screen and show the several top CSMFs of user's choice. See
#' below for more detail.
#' 
#' @aliases summary.insilico print.insilico_summary
#' @param object Fitted \code{"insilico"} object.
#' @param CI.csmf Confidence interval for CSMF estimates.
#' @param CI.cond Confidence interval for conditional probability estimates
#' @param file Optional .csv file to write to. If it is specified, individual
#' cause of death distribution will be saved to the file.
#' @param top Number of top causes to display on screen.
#' @param \dots Not used.
#' @return \item{id}{ the ID of the deaths.} \item{indiv}{ individual Cause of
#' Death distribution matrix.} \item{csmf}{ CSMF distribution and confidence
#' interval for each cause.} \item{csmf.ordered}{ CSMF distribution and
#' confidence interval for each cause, ordered by mean.} \item{condprob}{
#' Conditional probability matrix and confidence intervals.}
#' 
#' \item{useProbbase}{Component of \code{"insilico"} object.}
#' \item{keepProbbase.level}{Component of \code{"insilico"} object.}
#' \item{datacheck}{Component of \code{"insilico"} object.}
#' \item{length.sim}{Component of \code{"insilico"} object.}
#' \item{thin}{Component of \code{"insilico"} object.} \item{burnin}{Component
#' of \code{"insilico"} object.} % \item{HIV}{Component of \code{"insilico"}
#' object.} % \item{Malaria}{Component of \code{"insilico"} object.}
#' \item{jump.scale}{Component of \code{"insilico"} object.}
#' \item{levels.prior}{Component of \code{"insilico"} object.}
#' \item{levels.strength}{Component of \code{"insilico"} object.}
#' \item{trunc.min}{Component of \code{"insilico"} object.}
#' \item{trunc.max}{Component of \code{"insilico"} object.}
#' \item{subpop_counts}{Component of \code{"insilico"} object.}
#' \item{showTop}{Component of \code{"insilico"} object.}
#' @author Zehang Li, Tyler McCormick, Sam Clark
#' 
#' Maintainer: Zehang Li <lizehang@@uw.edu>
#' @seealso \code{\link{insilico}}, \code{\link{plot.insilico}}
#' @references Tyler H. McCormick, Zehang R. Li, Clara Calvert, Amelia C.
#' Crampin, Kathleen Kahn and Samuel J. Clark Probabilistic cause-of-death
#' assignment using verbal autopsies, \emph{arXiv preprint arXiv:1411.3042}
#' \url{http://arxiv.org/abs/1411.3042} (2014)
#' @examples
#' 
#' # load sample data together with sub-population list
#' data(SampleInput_insilico)
#' # extract INterVA style input data
#' data <- SampleInput_insilico$data
#' # extract sub-population information. 
#' # The groups are "HIV Positive", "HIV Negative" and "HIV status unknown".
#' subpop <- SampleInput_insilico$subpop
#' 
#' # run without subpopulation
#' fit1<- insilico( data, subpop = NULL, 
#'               length.sim = 400, burnin = 200, thin = 10 , seed = 1,
#'               external.sep = TRUE, keepProbbase.level = TRUE)
#' summary(fit1)
#' summary(fit1, top = 10)
#' 
#' @export summary.insilico
summary.insilico <- function(object, CI.csmf = 0.95, CI.cond = 0.95, 
				  file = NULL, top = 10, ...){
	id <- object$id
	prob <- object$indiv.prob
	csmf <- object$csmf
	if(is.null(object$conditional.probs)){
		cond.prob = FALSE
	}else{
		cond.prob = TRUE
	}
	if(cond.prob){
		cond <- object$conditional.probs
		if(length(dim(cond)) == 2){
			isLevel = TRUE
		}else{
			isLevel = FALSE
		}
	}

	indiv <- cbind(id, prob)
	## write individual COD distribution to file
	if(!is.null(file)){
		write.csv(indiv, file = file, row.names = FALSE)
	}

	# organize CSMF
	ci.low <- (1 - CI.csmf) / 2
	ci.up <- 1 - ci.low
	if(class(object$csmf) == "list"){
		csmf.out <- vector("list", length(csmf))
		csmf.out.ordered <- vector("list", length(csmf))
		for(i in 1:length(csmf)){
			mean <- apply(csmf[[i]], 2, mean)
			median <- apply(csmf[[i]], 2, median)
			sd <- apply(csmf[[i]], 2, sd)
			low <- apply(csmf[[i]], 2, function(object){quantile(object, prob = ci.low)})
			up <- apply(csmf[[i]], 2, function(object){quantile(object, prob = ci.up)})
			csmf.out[[i]]  <- cbind(mean, sd,  low, median, up)
			colnames(csmf.out[[i]]) <- cbind("Mean","Std.Error", "Lower", "Median", "Upper")
			csmf.out.ordered[[i]] <- csmf.out[[i]][order(csmf.out[[i]][,1],
												decreasing = TRUE), ]
			names(csmf.out)[i] <- names(csmf)[i]
			names(csmf.out.ordered)[i] <- names(csmf)[i]
		}
	}else{
		mean <- apply(csmf, 2, mean)
		median <- apply(csmf, 2, median)
		sd <- apply(csmf, 2, sd)
		low <- apply(csmf, 2, function(object){quantile(object, prob = ci.low)})
		up <- apply(csmf, 2, function(object){quantile(object, prob = ci.up)})
		csmf.out <- cbind(mean, sd, low, median, up)
		colnames(csmf.out) <- cbind("Mean","Std.Error","Lower",  "Median", "Upper")
		csmf.out.ordered <- csmf.out[order(csmf.out[,1], decreasing = TRUE), ]
	}
	
	

	# organize conditional probability matriobject
	ci.low <- (1 - CI.cond) / 2
	ci.up <- 1 - ci.low
	
	if(cond.prob){
		if(isLevel){
			mean <- apply(cond, 2, mean)
			median <- apply(cond, 2, median)
			sd <- apply(cond, 2, sd)
			low <- apply(cond, 2, function(object){quantile(object, prob = ci.low)})
			up <- apply(cond, 2, function(object){quantile(object, prob = ci.up)})
			cond.out <- cbind(mean, sd, low, median, up)
			colnames(cond.out) <- cbind("Mean","Std.Error",  "Lower", "Median","Upper")
			rownames(cond.out) <- colnames(cond)
		}else{
			mean <- t(apply(cond, c(2,3), mean))
			median <- t(apply(cond, c(2,3), median))
			sd <- t(apply(cond, c(2,3), sd))
			low <- t(apply(cond, c(2,3), function(object){quantile(object, probs =ci.low)}))
			up <- t(apply(cond, c(2,3), function(object){quantile(object, probs =ci.up)}))
			cond.out <- list(mean, sd, low,median, up)
			names(cond.out) <- c("Mean","Std.Error", "Lower", "Median", "Upper")
		}
	}else{
		cond.out <- NULL
	}

	# get subpopulation counts table
	if(!is.null(object$subpop)){
		subpop_counts <- table(object$subpop)
		subpop_counts <- subpop_counts[match(names(object$csmf), 
											 names(subpop_counts))]
	}else{
		subpop_counts <- NULL
	}

	out <- list( id = id, 
				 indiv = indiv, 
				 csmf = csmf.out, 
				 csmf.ordered = csmf.out.ordered,
				 condprob = cond.out, 				
				 useProbbase = object$useProbbase, 
				 keepProbbase.level = object$keepProbbase.level, 
				 datacheck = object$datacheck,
				 length.sim = object$length.sim, 
				 thin = object$thin, 
				 burnin = object$burnin, 
				 HIV = object$HIV, 
				 Malaria = object$Malaria, 
				 jump.scale = object$jump.scale, 
				 levels.prior = object$levels.prior, 
				 levels.strength = object$levels.strength, 
				 trunc.min = object$trunc.min, 
				 trunc.maobject = object$trunc.max, 
			     subpop_counts = subpop_counts,
				 showTop = top)
	class(out) <- "insilico_summary"
	return(out)
}

print.insilico <- function(x,...){
	cat("InSilicoVA fitted object:\n")
	cat(paste(length(x$id), "death processed\n"))
	cat(paste(x$length.sim, "iterations performed, with first", 
			  x$burnin, "iterations discarded\n",
			  trunc((x$length.sim - x$burnin)/x$thin), "iterations saved after thinning\n"))
		if(x$useProbbase){
			cat("Fitted with InterVA4 conditional probability matrix\n")
		}else if(x$keepProbbase.level){
			cat("Fitted with re-estimated InterVA4 conditional probability level table\n")	
		}else{
			cat("Fitted with re-estimating InterVA4 conditional probability matrix\n")	
		}    
}

print.insilico_summary <- function(x, ...) {
	cat("InSilicoVA Call: \n")
	cat(paste(length(x$id), "death processed\n"))
	cat(paste(x$length.sim, "iterations performed, with first", 
			  x$burnin, "iterations discarded\n",
			  trunc((x$length.sim - x$burnin)/x$thin), "iterations saved after thinning\n"))
	if(x$useProbbase){
			cat("Fitted with InterVA4 conditional probability matrix\n")
		}else if(x$keepProbbase.level){
			cat("Fitted with re-estimated InterVA4 conditional probability level table\n")	
		}else{
			cat("Fitted with re-estimating InterVA4 conditional probability matrix\n")	
		}   
	if(x$datacheck){
		cat("Data consistency check performed as in InterVA4 \n")
	}
	if(!is.null(x$HIV)){
		cat(paste("HIV level:     ", x$HIV, "\n"))
		cat(paste("Malaria level: ", x$Malaria, "\n"))
	}

	if(!is.null(x$subpop_counts)){
		cat("Sub population frequencies:")
		print(x$subpop_counts)
	}	
	
	if(class(x$csmf.ordered) != "list"){
		top <- min(x$showTop, dim(x$csmf.ordered)[1])
		cat("\n")
		cat(paste("Top", top,  "CSMFs:\n"))
		csmf.out.ordered <- x$csmf.ordered
	    out <- as.matrix(csmf.out.ordered[1:top, ])
			out[, 1] <- round(out[, 1], 4)
			out[, 2] <- round(out[, 2], 4)
			out[, 3] <- round(out[, 3], 4)
			out[, 4] <- round(out[, 4], 4)
			out[, 5] <- round(out[, 5], 4)
		print(out)		
	}else{
		for(i in 1:length(x$csmf.ordered)){
			top <- min(x$showTop, dim(x$csmf.ordered[[i]])[1])
			cat("\n")
			cat(paste(names(x$csmf.ordered)[i], "- Top", top,  "CSMFs:\n"))
			csmf.out.ordered <- x$csmf.ordered[[i]]
		    out <- as.matrix(csmf.out.ordered[1:top, ])
				out[, 1] <- round(out[, 1], 4)
				out[, 2] <- round(out[, 2], 4)
				out[, 3] <- round(out[, 3], 4)
				out[, 4] <- round(out[, 4], 4)
				out[, 5] <- round(out[, 5], 4)
			print(out)		
		}
	}
}
	

