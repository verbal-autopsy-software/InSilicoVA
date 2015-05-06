\name{summary.insilico}
\alias{summary.insilico}
\alias{print.insilico_summary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summarizing InSilicoVA Model Fits
}
\description{
This function is the summary method for class \code{insilico}.}
\usage{
\method{summary}{insilico}(object, CI.csmf = 0.95, CI.cond = 0.95, 
    file = NULL, top = 10, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{Fitted \code{"insilico"} object.}
  \item{CI.csmf}{Confidence interval for CSMF estimates.}
  \item{CI.cond}{Confidence interval for conditional probability estimates}
  \item{file}{Optional .csv file to write to. If it is specified, individual cause of death distribution will be saved to the file.}
  \item{top}{Number of top causes to display on screen.}
  \item{\dots}{Not used.}
}
\details{
\code{summary.insilico} formats some basic information about the InSilicoVA fitted object on screen and show the several top CSMFs of user's choice. See below for more detail. 
}
\value{
\item{id}{ the ID of the deaths.} 
\item{indiv}{ individual Cause of Death distribution matrix.}
\item{csmf}{ CSMF distribution and confidence interval for each cause.}
\item{csmf.ordered}{ CSMF distribution and confidence interval for each cause, ordered by mean.}
\item{condprob}{ Conditional probability matrix and confidence intervals.}

\item{useProbbase}{Component of \code{"insilico"} object.}
\item{keepProbbase.level}{Component of \code{"insilico"} object.}
\item{datacheck}{Component of \code{"insilico"} object.}
\item{length.sim}{Component of \code{"insilico"} object.}
\item{thin}{Component of \code{"insilico"} object.}
\item{burnin}{Component of \code{"insilico"} object.}
% \item{HIV}{Component of \code{"insilico"} object.}
% \item{Malaria}{Component of \code{"insilico"} object.}
\item{jump.scale}{Component of \code{"insilico"} object.}
\item{levels.prior}{Component of \code{"insilico"} object.}
\item{levels.strength}{Component of \code{"insilico"} object.}
\item{trunc.min}{Component of \code{"insilico"} object.}
\item{trunc.max}{Component of \code{"insilico"} object.}
\item{subpop_counts}{Component of \code{"insilico"} object.}
\item{showTop}{Component of \code{"insilico"} object.}
}

\author{
Zehang Li, Tyler McCormick, Sam Clark

Maintainer: Zehang Li <lizehang@uw.edu>

}
\references{
Tyler H. McCormick, Zehang R. Li, Clara Calvert, Amelia C. Crampin, Kathleen Kahn and Samuel J. Clark Probabilistic cause-of-death assignment using verbal autopsies, \emph{arXiv preprint arXiv:1411.3042} \url{http://arxiv.org/abs/1411.3042} (2014)
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{insilico}}, \code{\link{plot.insilico}} 
}
\examples{
# load sample data together with sub-population list
data(SampleInput_insilico)
# extract INterVA style input data
data <- SampleInput_insilico$data
# extract sub-population information. 
# The groups are "HIV Positive", "HIV Negative" and "HIV status unknown".
subpop <- SampleInput_insilico$subpop

# run without subpopulation
fit1<- insilico( data, subpop = NULL, 
              length.sim = 400, burnin = 200, thin = 10 , seed = 1,
              external.sep = TRUE, keepProbbase.level = TRUE)
summary(fit1)
summary(fit1, top = 10)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

