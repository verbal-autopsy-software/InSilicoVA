pkgname <- "InSilicoVA"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('InSilicoVA')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("SampleInput_insilico")
### * SampleInput_insilico

flush(stderr()); flush(stdout())

### Name: SampleInput_insilico
### Title: 100 records of Sample Input
### Aliases: SampleInput_insilico
### Keywords: datasets

### ** Examples

data(SampleInput_insilico)



cleanEx()
nameEx("causetext")
### * causetext

flush(stderr()); flush(stdout())

### Name: causetext
### Title: Translation list of COD codes
### Aliases: causetext
### Keywords: datasets

### ** Examples

data(causetext)



cleanEx()
nameEx("condprob")
### * condprob

flush(stderr()); flush(stdout())

### Name: condprob
### Title: Conditional probability table used by InterVA-4
### Aliases: condprob
### Keywords: datasets

### ** Examples

data(condprob)



cleanEx()
nameEx("csmf.diag")
### * csmf.diag

flush(stderr()); flush(stdout())

### Name: csmf.diag
### Title: Convergence test for fitted InSilico model
### Aliases: csmf.diag
### Keywords: InSilicoVA

### ** Examples

# load sample data together with sub-population list
data(SampleInput_insilico)
## Not run: 
##D # extract INterVA style input data
##D data <- SampleInput_insilico$data
##D # extract sub-population information. 
##D # The groups are "HIV Positive", "HIV Negative" and "HIV status unknown".
##D subpop <- SampleInput_insilico$subpop
##D 
##D # run without sub-population
##D fit1a<- insilico( data, subpop = NULL, HIV = "h", Malaria = "h", 
##D               length.sim = 400, burnin = 200, thin = 10 , seed = 1, 
##D               auto.length = FALSE)
##D fit1b<- insilico( data, subpop = NULL, HIV = "h", Malaria = "h", 
##D               length.sim = 400, burnin = 200, thin = 10 , seed = 2, 
##D               auto.length = FALSE)
##D fit1c<- insilico( data, subpop = NULL, HIV = "h", Malaria = "h", 
##D               length.sim = 400, burnin = 200, thin = 10 , seed = 3, 
##D               auto.length = FALSE)
##D # single chain check
##D csmf.diag(fit1a)
##D # equivalent way of check one chain
##D csmf.diag(fit1a$csmf)
##D 
##D # multiple chains check
##D csmf.diag(list(fit1a, fit1b, fit1c), test = "gelman")
##D # equivalent way of check one chain
##D csmf.diag(list(fit1a$csmf, fit1b$csmf, fit1c$csmf), test = "gelman")
##D 
##D 
##D # with sub-populations
##D fit2a<- insilico( data, subpop = subpop, HIV = "h", Malaria = "h", 
##D               length.sim = 400, burnin = 200, thin = 10 , seed = 1, 
##D               auto.length = FALSE)
##D fit2b<- insilico( data, subpop = subpop, HIV = "h", Malaria = "h", 
##D               length.sim = 400, burnin = 200, thin = 10 , seed = 2, 
##D               auto.length = FALSE)
##D fit2c<- insilico( data, subpop = subpop, HIV = "h", Malaria = "h", 
##D               length.sim = 400, burnin = 200, thin = 10 , seed = 3, 
##D               auto.length = FALSE)
##D 
##D # single chain check
##D csmf.diag(fit2a)
##D # equivalent way of check one chain
##D csmf.diag(fit2a$csmf)
##D 
##D # multiple chains check
##D csmf.diag(list(fit2a, fit2b, fit2c), test = "gelman")
##D # equivalent way of check one chain
##D csmf.diag(list(fit2a$csmf, fit2b$csmf, fit2c$csmf), test = "gelman")
## End(Not run)




cleanEx()
nameEx("insilico")
### * insilico

flush(stderr()); flush(stdout())

### Name: insilico
### Title: Implement InSilicoVA methods
### Aliases: insilico print.insilico
### Keywords: InSilicoVA

### ** Examples

# load sample data together with sub-population list
data(SampleInput_insilico)
# extract INterVA style input data
data <- SampleInput_insilico$data
# extract sub-population information. 
# The groups are "HIV Positive", "HIV Negative" and "HIV status unknown".
subpop <- SampleInput_insilico$subpop

# run without subpopulation
fit1<- insilico( data, subpop = NULL, HIV = "h", Malaria = "h", 
              length.sim = 400, burnin = 200, thin = 10 , seed = 1,
              external.sep = TRUE, keepProbbase.level = TRUE)

# re-run with subpopulation
fit2<- insilico( data, subpop = subpop, HIV = "h", Malaria = "h", 
              length.sim = 400, burnin = 200, thin = 10 , seed = 1,
              external.sep = TRUE, keepProbbase.level = TRUE)




cleanEx()
nameEx("plot.insilico")
### * plot.insilico

flush(stderr()); flush(stdout())

### Name: plot.insilico
### Title: plot CSMF from a "insilico" object
### Aliases: plot.insilico
### Keywords: InSilicoVA

### ** Examples

# load sample data together with sub-population list
data(SampleInput_insilico)
# extract INterVA style input data
data <- SampleInput_insilico$data
# extract sub-population information. 
# The groups are "HIV Positive", "HIV Negative" and "HIV status unknown".
subpop <- SampleInput_insilico$subpop

# run without sub-population
fit1<- insilico( data, subpop = NULL, HIV = "h", Malaria = "h", 
              length.sim = 400, burnin = 200, thin = 10 , seed = 1,
              external.sep = TRUE, keepProbbase.level = TRUE)
# default plot
plot(fit1)
# customized plot
plot(fit1, top = 15, horiz = FALSE, fill = "gold", 
		   bw = TRUE, title = "Top 15 CSMFs", angle = 70)

# run with sub-populations
fit2<- insilico( data, subpop = subpop, HIV = "h", Malaria = "h", 
              length.sim = 400, burnin = 200, thin = 10 , seed = 1,
              external.sep = TRUE, keepProbbase.level = TRUE)
# default plot
plot(fit2, type = "compare", top = 5, title = "Top 5 causes comparison")
# customized single sub-population plot
plot(fit2, which.sub = "Unknown", 
	 title = "Top 10 causes in HIV status unknown population")
# customized plot with specified causes, with abbreviation here.
some_causes <- c("HIV", "Pulmonary", 
				"Other and unspecified infect dis")
plot(fit2, type = "compare", horiz = FALSE,  causelist = some_causes,
		   title = "Infectious diseases in three sub-populations", 
		   angle = 20)




cleanEx()
nameEx("probbase")
### * probbase

flush(stderr()); flush(stdout())

### Name: probbase
### Title: Conditional probability of InterVA4
### Aliases: probbase
### Keywords: datasets

### ** Examples

data(probbase)



cleanEx()
nameEx("summary.insilico")
### * summary.insilico

flush(stderr()); flush(stdout())

### Name: summary.insilico
### Title: Summarizing InSilicoVA Model Fits
### Aliases: summary.insilico print.insilico_summary

### ** Examples

# load sample data together with sub-population list
data(SampleInput_insilico)
# extract INterVA style input data
data <- SampleInput_insilico$data
# extract sub-population information. 
# The groups are "HIV Positive", "HIV Negative" and "HIV status unknown".
subpop <- SampleInput_insilico$subpop

# run without subpopulation
fit1<- insilico( data, subpop = NULL, HIV = "h", Malaria = "h", 
              length.sim = 400, burnin = 200, thin = 10 , seed = 1,
              external.sep = TRUE, keepProbbase.level = TRUE)
summary(fit1)
summary(fit1, top = 10)



### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
