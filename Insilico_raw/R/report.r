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
		csmf.out.ordered <- vector("list", length(csmf))
		for(i in 1:length(csmf)){
			mean <- apply(csmf[[i]], 2, mean)
			median <- apply(csmf[[i]], 2, median)
			sd <- apply(csmf[[i]], 2, sd)
			low <- apply(csmf[[i]], 2, function(object){quantile(object, prob = ci.low)})
			up <- apply(csmf[[i]], 2, function(object){quantile(object, prob = ci.up)})
			csmf.out <- cbind(mean, sd,  low, median, up)
			colnames(csmf.out) <- cbind("Mean","Std.Error", "Lower", "Median", "Upper")
			csmf.out.ordered[[i]] <- csmf.out[order(csmf.out[,1],
												decreasing = TRUE), ]
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
	

