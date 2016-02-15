
#' Implement InSilicoVA methods
#' 
#' This function implements InSilicoVA model with non-InterVA4 input data. 
#' 
#' Please see \code{insilico} for more details about choosing chain length and 
#' OS system differences. This function implements InSilico with customized
#' input format and training data. 
#' 
#' For more detail of model specification, see the paper on
#' \url{http://arxiv.org/abs/1411.3042}.
#' 
#' @param data The original data to be used. It is suggested to use similar
#' input as InterVA4, with the first column being death IDs and 245 symptoms. 
#' The only difference in input is InsilicoVA takes three levels: ``present'',
#' ``absent'', and ``missing (no data)''. Similar to InterVA software,
#' ``present'' symptoms takes value ``Y''; ``absent'' symptoms take take value
#' ``NA'' or ``''. For missing symptoms, e.g., questions not asked or answered
#' in the original interview, corrupted data, etc., the input should be coded
#' by ``.'' to distinguish from ``absent'' category. The order of the columns does
#' not matter as long as the column names are correct. It can also include more 
#' unused columns than the standard InterVA4 input. But the first column should be 
#' the death ID. For example input data format, see \code{RandomVA1} and 
#' \code{RandomVA2}.
#' @param train Training data, it should be in the same format as the testing data
#' and contains one additional column (see \code{cause} below) specifying known
#' cause of death. The first column is also assumed to be death ID. 
#' @param cause the name of the column in \code{train} that contains cause of death.
#' @param cause.table The list of causes of death used in training data.
#' @param thre a numerical value between 0 to 1. It specifies the maximum rate of
#' missing for any symptoms to be considered in the model. Default value is set to
#' 0.95, meaning if a symptom has more than 95\% missing in the training data, it
#' will be removed.
#' @param type Two types of learning conditional probabilities are provided: ``quantile''
#' or ``fixed''. Since InSilicoVA works with ranked conditional probabilities P(S|C), ``quantile''
#' means the rankings of the P(S|C) are obtained by matching the same quantile distributions
#' in the default InterVA P(S|C), and ``fixed'' means P(S|C) are matched to the closest values
#' in the default InterVA P(S|C) table. Empirically both types of rankings produce similar results.
#' @param isNumeric Indicator if the input is already in numeric form. If the
#' input is coded numerically such that 1 for ``present'', 0 for ``absent'',
#' and -1 for ``missing'', this indicator could be set to True to avoid
#' conversion to standard InterVA format.
#' @param updateCondProb Logical indicator. If FALSE, then fit InSilicoVA model without 
#' re-estimating conditional probabilities.
#' @param keepProbbase.level see \code{\link{insilico}} for more detail.
#' @param CondProb see \code{\link{insilico}} for more detail. 
#' @param CondProbNum see \code{\link{insilico}} for more detail. 
#' @param datacheck Not Implemented.
#' @param datacheck.missing Not Implemented.
#' @param warning.write Not Implemented.
#' @param external.sep Not Implemented.
#' @param length.sim see \code{\link{insilico}} for more detail.
#' @param thin see \code{\link{insilico}} for more detail.
#' @param burnin see \code{\link{insilico}} for more detail.
#' @param auto.length see \code{\link{insilico}} for more detail.
#' @param conv.csmf see \code{\link{insilico}} for more detail.
#' @param jump.scale see \code{\link{insilico}} for more detail.
#' @param levels.prior see \code{\link{insilico}} for more detail.
#' @param levels.strength see \code{\link{insilico}} for more detail.
#' @param trunc.min see \code{\link{insilico}} for more detail.
#' @param trunc.max see \code{\link{insilico}} for more detail.
#' @param subpop see \code{\link{insilico}} for more detail.
#' @param java_option see \code{\link{insilico}} for more detail.
#' @param seed see \code{\link{insilico}} for more detail.
#' @param phy.code see \code{\link{insilico}} for more detail.
#' @param phy.cat see \code{\link{insilico}} for more detail.
#' @param phy.unknown see \code{\link{insilico}} for more detail.
#' @param phy.external see \code{\link{insilico}} for more detail.
#' @param phy.debias see \code{\link{insilico}} for more detail.
#' @param exclude.impossible.cause Not Implemented
#' @param indiv.CI see \code{\link{insilico}} for more detail.
#'
#' @return \code{insilico} object
#' @export insilico.train
#'
#' @examples
#' \dontrun{
#' x <- 1
#' }
insilico.train <- function(data, train, cause, cause.table, thre = 0.95, type = c("quantile", "fixed")[1], isNumeric = FALSE, updateCondProb = TRUE, keepProbbase.level = TRUE,  CondProb = NULL, CondProbNum = NULL, datacheck = TRUE, datacheck.missing = TRUE, warning.write = FALSE, external.sep = TRUE, length.sim = 4000, thin = 10, burnin = 2000, auto.length = TRUE, conv.csmf = 0.02, jump.scale = 0.1, levels.prior = NULL, levels.strength = 1, trunc.min = 0.0001, trunc.max = 0.9999, subpop = NULL, java_option = "-Xmx1g", seed = 1, phy.code = NULL, phy.cat = NULL, phy.unknown = NULL, phy.external = NULL, phy.debias = NULL, exclude.impossible.cause = TRUE, indiv.CI = 0.95){ 
	
	prob.learn <- extract.prob(train = train, 
							  gs = cause, 
							  gstable = cause.table, 
							  thre = thre, 
							  type = type, 
							  isNumeric = isNumeric)
	# remove unused symptoms
	col.exist <- c(colnames(data)[1], cause, colnames(prob.learn$symps.train))
	remove <- which(colnames(data) %in% col.exist == FALSE)
	if(length(remove) > 0){
		warning(paste(length(remove), "symptoms deleted from testing data to match training data:", 
			paste(colnames(data)[remove], collapse = ", ")),
			immediate. = TRUE)	
		data <- data[, -remove]
	}

	if(updateCondProb){
		probbase.touse <- prob.learn$cond.prob.alpha
		CondProbNum <- NULL
	}else{
		probbase.touse <- prob.learn$cond.prob
		CondProbNum <- prob.learn$cond.prob
	}

	fit <- insilico.fit(data = data, 
						isNumeric = isNumeric, 
						updateCondProb = updateCondProb, 
						keepProbbase.level = keepProbbase.level, 
						CondProb = CondProb, 
						CondProbNum = CondProbNum, 
						datacheck = FALSE, 
						datacheck.missing = FALSE, 
						warning.write = FALSE, 
						external.sep = FALSE, 
						length.sim = length.sim, 
						thin = thin, 
						burnin = burnin, 
						auto.length = auto.length, 
						conv.csmf = conv.csmf, 
						jump.scale = jump.scale, 
						levels.prior = levels.prior, 
						levels.strength = levels.strength, 
						trunc.min = trunc.min, 
						trunc.max = trunc.max, 
						subpop = subpop, 
						java_option = java_option, 
						seed = seed, 
						phy.code = phy.code, 
						phy.cat = phy.cat, 
						phy.unknown = phy.unknown, 
						phy.external = phy.external, 
						phy.debias = phy.debias, 
						exclude.impossible.cause = FALSE, 
						indiv.CI = indiv.CI, 
						
						customization.dev = TRUE, 
						Probbase_by_symp.dev = FALSE, 
						probbase.dev = probbase.touse, 
						table.dev = prob.learn$table.alpha, 
						table.num.dev = prob.learn$table.num, 
						gstable.dev = colnames(prob.learn$cond.prob), 
						nlevel.dev = 15, 
						)
	return(fit)  	
} 
