\name{insilico}
\alias{insilico}
\alias{print.insilico}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Implement InSilicoVA methods
}
\description{ 
This function implements InSilicoVA model. The InSilicoVA model is fitted with MCMC implemented in Java. For more detail, see the paper on \url{http://arxiv.org/abs/1411.3042}.

For Windows user, this function will produce a popup window showing the progress. For Mac and Unix user, this function will print progress messages on the console. Special notice for users using default R GUI for mac, the output will not be printed on console while the function is running, and will only be printed out after it is completed. Thus if you use a Mac, we suggest using either RStudio for mac, or running R from terminal.
}

\usage{
  insilico(data, isNumeric = FALSE,useProbbase = FALSE, keepProbbase.level = TRUE,  
           cond.prob.touse = NULL, datacheck = TRUE, warning.write = FALSE, 
           external.sep = TRUE, length.sim = 4000, thin = 10, burnin = 2000, 
           auto.length = TRUE, conv.csmf = 0.02, jump.scale = 0.1, 
           levels.prior = NULL, levels.strength = 1, trunc.min = 0.0001, 
           trunc.max = 0.9999, subpop = NULL, java_option = "-Xmx1g", seed = 1)
}
\arguments{
  \item{data}{The original data to be used. It is suggested to use similar input as InterVA4, with the first column being death IDs. The only difference in input is InsilicoVA takes three levels: ``present'', ``absent'', and ``missing (no data)''. Similar to InterVA software, ``present'' symptoms takes value ``Y''; ``absent'' symptoms take take value ``NA'' or ``''. For missing symptoms, e.g., questions not asked or answered in the original interview, corrupted data, etc., the input should be coded by ``.'' to distinguish from ``absent'' category.}

  \item{isNumeric}{Indicator if the input is already in numeric form. If the input is coded numerically such that 1 for ``present'', 0 for ``absent'', and -1 for ``missing'', this indicator could be set to True to avoid conversion to standard InterVA format.}

  \item{useProbbase}{Logical indicator. If TRUE, then use InterVA conditional probability table without re-estimating.}

  \item{keepProbbase.level}{Logical indicator when \code{useProbbase} is FALSE. If TRUE, then only estimate the InterVA's conditional probability table interpretation table; if FALSE, estimate the whole conditional probability matrix. Default to TRUE.}

  \item{cond.prob.touse}{Customized conditional probability table to use. Currently only accepting the same configuration as InterVA-4 software. It should be a matrix of 245 rows of symptoms and 60 columns of causes, arranged in the same order as in InterVA-4 specification. The elements in the matrix should be the conditional probability of corresponding symptom given the corresponding cause, represented in alphabetic form indicating levels. For example input, see \code{\link{condprob}}}

  \item{datacheck}{Logical indicator for whether to check the data satisfying InterVA rules. Default set to be TRUE. If \code{warning.write} is set to true, the inconsistent input will be logged in file warnings.txt. It's strongly suggested to be set to TRUE.}

  \item{warning.write}{Logical indicator for whether to save the changes made to data input by \code{datacheck}. If set to TRUE, the changes will be logged in file warnings.txt in current working directory.}
%% 
%%
  \item{external.sep}{Logical indicator for whether to separate out external causes first. Default set to be TRUE. If set to TRUE, the algorithm will estimate external causes, e.g., traffic accident, accidental fall, suicide, etc., by checking the corresponding indicator only without considering other medical symptoms. It is strongly suggested to set to be TRUE.}

  
  \item{length.sim}{Number of iterations to run. Default to be 4000.}

  \item{thin}{Proportion of thinning for storing parameters. For example, if thin = k, the output parameters will only be saved every k iterations. Default to be 10}

  \item{burnin}{Number of iterations as burn-in period. Parameters sampled in burn-in period will not be saved.}

  \item{auto.length}{Logical indicator of whether to automatically increase chain length if convergence not reached.}

  \item{conv.csmf}{Minimum CSMF value to check for convergence if auto.length is set to TRUE. For example, under the default value 0.02, all causes with mean CSMF at least 0.02 will be checked for convergence.}
%   \item{alpha.scale}{
% Parameter for Dirichlet model. \eqn{\alpha} value used in Dirichlet prior. Larger \eqn{\alpha} gives less volatile CSMF.
% }
%   \item{InterVA.prior}{ Logical indicator of whether to use InterVA prior as prior on CSMF. Default to be TRUE.
%   }
%    \item{csmf.prior}{ If \code{UseInterVA} is set to 0, this allows a user defined distribution of CSMF to be used. If it is not provided, the algorithm will start from uniform distribution across all causes.
% }
 
%   \item{HIV}{ The HIV level used when \code{UseHIV} is set to 1 or 2.
% }
%   \item{Malaria}{ The Malaria level used when \code{UseHIV} is set to 1 or 2.
% }
  \item{jump.scale}{ The scale of Metropolis proposal in the Normal model. Default to be 0.1.
}

  \item{levels.prior}{ Vector of prior expectation of conditional probability levels. They do not have to be scaled. The algorithm internally calibrate the scale to the working scale through \code{levels.strength}. If NULL the algorithm will use InterVA table as prior.
}
  \item{levels.strength}{ Scaling factor for the strength of prior beliefs in the conditional probability levels. Larger value constrain the posterior estimates to be closer to prior expectation. Defult value 1 scales \code{levels.prior} to a suggested scale that works empirically.  
}
  \item{trunc.min}{Minimum possible value for estimated conditional probability table. Default to be 0.0001}

  \item{trunc.max}{Maximum possible value for estimated conditional probability table. Default to be 0.9999}
  
  \item{subpop}{A vector of sub-population assignments of the same length of death records. It could be numerical indicators or character vectors of names.}


  \item{java_option}{Option to initialize java JVM. Default to ``-Xmx1g'', which sets the maximum heap size to be 1GB. If R produces ``java.lang.OutOfMemoryError: Java heap space'' error message, consider increasing heap size using this option, or one of the following: (1) decreasing \code{length.sim}, (2) increasing \code{thin}, or (3) disabling \code{auto.length}.}

  \item{seed}{ Seed used for initializing sampler. The algorithm will produce the same outcome with the same seed in each machine.}

}
\details{
The chains could be set to run automatically longer. If set \code{auto.length} to be TRUE, the chain will assess convergence after finishing the length K chain input by user using Heidelberger and Welch's convergence diagnostic. If convergence is not reached, the chain will run another K iterations and use the first K iterations as burn-in. If the chain is still not converged after 2K iterations, it will proceed to another 2K iterations and again use the first 2K iterations as burn-in. If convergence is still not reached by the end, it will not double the length again to avoid heavy memory use. A warning will be given in that case. The extended chains will be thinned in the same way.

For more detail of model specification, see the paper on \url{http://arxiv.org/abs/1411.3042}.
}

\value{
\item{id}{A vector of death ID. Note the order of the ID is in general different from the input file. See \code{report} for organizing the report.}

\item{indiv.prob}{Matrix of individual mean cause of death distribution. Each row corresponds to one death with the corresponding ID.}

\item{csmf}{Matrix of CSMF vector at each iterations after burn-in and thinning. Each column corresponds to one cause.}

\item{conditional.probs}{If the model is estimated with \code{keepProbbase.level} = TRUE, this value gives a matrix of each conditional probability at each level at each iterations. Each column corresponds to one level of probability. If \code{keepProbbase.level} = FALSE, this value gives a three-dimensional array. If \code{UseProbbase} = TRUE, the value will be set to NULL. See \code{report} for more analysis.}

\item{missing.symptoms}{Vector of symptoms missing from all input data.}

\item{external}{Logical indicator of whether the model is fitted with external causes separated calculated.}

}
\author{
Zehang Li, Tyler McCormick, Sam Clark

Maintainer: Zehang Li <lizehang@uw.edu>

}
\references{
Tyler H. McCormick, Zehang R. Li, Clara Calvert, Amelia C. Crampin, Kathleen Kahn and Samuel J. Clark(2014) 
\emph{Probabilistic cause-of-death assignment using verbal autopsies}, \url{http://arxiv.org/abs/1411.3042}
\cr \emph{Working paper no. 147, Center for Statistics and the Social Sciences, University of Washington}
}
\keyword{ InSilicoVA }
% \note{
% }
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{plot.insilico}}, \code{\link{summary.insilico}} 
}
\examples{
# load sample data together with sub-population list
data(SampleInput_insilico)
# extract InterVA style input data
data <- SampleInput_insilico$data
# extract sub-population information. 
# The groups are "HIV Positive", "HIV Negative" and "HIV status unknown".
subpop <- SampleInput_insilico$subpop

# run without subpopulation
fit1<- insilico( data, subpop = NULL,  
              length.sim = 400, burnin = 200, thin = 10 , seed = 1,
              external.sep = TRUE, keepProbbase.level = TRUE)

# re-run with subpopulation
fit2<- insilico( data, subpop = subpop, 
              length.sim = 400, burnin = 200, thin = 10 , seed = 1,
              external.sep = TRUE, keepProbbase.level = TRUE)

# run without re-sampling conditional probabilities
fit3<- insilico( data, subpop = subpop, 
              length.sim = 400, burnin = 200, thin = 10 , seed = 1,
              external.sep = TRUE, useProbbase = TRUE)

}
