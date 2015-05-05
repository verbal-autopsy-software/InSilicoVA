plot.insilico <- function(x, type = c("bar", "compare")[1], 
	top = 10, causelist = NULL, which.sub = NULL, 
	xlab = "Causes", ylab = "CSMF", title = "Top CSMF Distribution", 
	horiz = TRUE, angle = 60, fill = "lightblue", border = "black", bw = FALSE,
	...){
	
	sx <- summary(x)
	
	# to please R CMD Check with ggplot
	Group <- Mean <- Lower <- Upper <- Causes <- NULL

	if(type == "compare" && is.null(sx$subpop_counts)){
		type <- "bar"
		warning("No group detected, switch to simple bar plot")
	}
	##
	## for multiple populations plot 
	##
	if(type == "compare"){
		# get which causes to plot
		# if causelist is provided then fine, 
		# otherwise construct from top argument
		if(top > 0 && is.null(causelist)){
			for(i in 1:length(sx$csmf.ordered)){
				causelist <- c(causelist, 
						setdiff(rownames(sx$csmf.ordered[[i]])[1:round(top)], 
								causelist))
			}
		}
		
		csmf.toplot <- NULL
		csmf.sub <- NULL
		for(i in 1:length(sx$csmf.ordered)){
			vcauses <- rownames(sx$csmf.ordered[[i]])
			if(class(causelist) == "character"){
				toplot <- rep(1, length(causelist))
				for(j in 1:length(causelist)){
					full <- match.arg(tolower(causelist[j]), tolower(vcauses))
					toplot[j] <- match(full, tolower(vcauses))
				}
			}else{
				# if numeric cause list is given
				toplot <- match(rownames(sx$csmf[[1]])[causelist], 
								vcauses)
			}
			csmf.toplot.tmp <- sx$csmf.ordered[[i]][toplot, c(1, 3, 5)]
			csmf.toplot <- rbind(csmf.toplot, csmf.toplot.tmp)
			csmf.sub <- c(csmf.sub,  rep(names(sx$csmf.ordered)[i], 
									     dim(csmf.toplot.tmp)[1]))
		}
		rownames(csmf.toplot) <- NULL
		csmf.toplot <- data.frame(Causes=causelist, csmf.toplot, Group =csmf.sub)
	}
	##
	## for single population plot 
	##
	if(type == "bar"){
		# check which subpopulation to plot first.
		# if no subpopulation
		if(is.null(sx$subpop_counts)){
			csmf.ordered <- sx$csmf.ordered[, c(1, 3, 5)]
		}else{
			# if there are multiple subpopulation
			if(is.null(which.sub)) stop("More than one groups detected. Please specify which to plot")
			if(class(which.sub) == "character"){
				which.sub <- match(which.sub, names(sx$csmf.ordered))
				if(is.na(which.sub)) stop("Invalid 'which.sub', no match found.")
			}
			csmf.ordered <- sx$csmf.ordered[[which.sub]][, c(1, 3, 5)]	
		}
		# get cause names
		vcauses <- rownames(csmf.ordered)
		# get which causes to plot
		if(is.null(causelist)){
			if(top > 0){
				toplot <- 1:round(top)
			}else{
				stop("Invalid causes to plot. Please specify causelist or a positive top")
			}
		}else{
			if(class(causelist) == "character"){
				toplot <- rep(1, length(causelist))
				for(i in 1:length(causelist)){
					full <- match.arg(tolower(causelist[i]), tolower(vcauses))
					toplot[i] <- match(full, tolower(vcauses))
				}
			}else{
				# if numeric cause list is given
				toplot <- match(rownames(sx$csmf)[causelist], 
								vcauses)
			}
		}
		csmf.toplot <- data.frame(Causes = rownames(csmf.ordered[toplot, ]), 
								  csmf.ordered[toplot, ])
	}

	
	# making bar plot
	# require: number of top causes or which causes to plot;
	#		   color list
	if(type == "bar"){
		# initialize ggplot, force order of bars
		if(horiz){
			g <- ggplot(csmf.toplot, aes(x=reorder(Causes, seq(length(Causes), 1)),
									 y=Mean))	
		}else{
			g <- ggplot(csmf.toplot, aes(x=reorder(Causes, seq(1:length(Causes))),
									 y=Mean))
		}
		g <- g + geom_bar( stat="identity", 
				 colour=border, fill=fill,  size =.3)
		g <- g + geom_errorbar(aes(ymin = Lower, ymax = Upper), size = .3, width = .2, position = position_dodge(.9))
		g <- g + xlab(xlab) + ylab(ylab) 
		g <- g + ggtitle(title)
		if(horiz) g <- g + coord_flip()
		if(bw) g <- g + theme_bw()
		if(!horiz) g <- g + theme(axis.text.x = element_text(angle = angle, hjust = 1))	
		return(g)
	}

	# making bar plot for comparison
	if(type == "compare"){
		# initialize ggplot, force order of bars
		if(horiz){
			g <- ggplot(csmf.toplot, aes(x=reorder(Causes, seq(length(Causes), 1)),
									 y=Mean, 
									 fill = Group))
		}else{
			g <- ggplot(csmf.toplot, aes(x=reorder(Causes, seq(1:length(Causes))),
									 y=Mean, 
									 fill = Group))
		}
		g <- g + geom_bar( stat="identity", position = "dodge", 
				 colour=border,  size =.3)
		g <- g + geom_errorbar(aes(ymin = Lower, ymax = Upper), size = .3, width = .2, position = position_dodge(.9))
		g <- g + xlab(xlab) + ylab(ylab) 
		g <- g + ggtitle(title)
		g <- g + scale_y_continuous() 
		if(horiz) g <- g + coord_flip()
		if(bw) g <- g + theme_bw()
		if(!horiz) g <- g + theme(axis.text.x = element_text(angle = angle, hjust = 1))
		return(g)
	}   

}