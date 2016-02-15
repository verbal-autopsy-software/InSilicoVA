
# test script and all example functions in InSilicoVA package
install.packages("~/Bitbucket-repos/InSilicoVA-beta/InSilicoVA_1.1.tar.gz", repos = NULL, type = "source")
library(InSilicoVA)


#---------------------------------------------------------------------#
#-------------------------  insilico    ------------------------------#
#---------------------------------------------------------------------#
# Toy example with 1000 VA deaths
data(RandomVA1) 
fit0<- insilico(RandomVA1, subpop = NULL,  
                length.sim = 20, burnin = 10, thin = 1 , seed = 1,
			           auto.length = FALSE)
summary(fit0)
summary(fit0, id = "d199")

##
## Scenario 1: standard input without sub-population specification
##
fit1<- insilico(RandomVA1, subpop = NULL,  
              length.sim = 1000, burnin = 500, thin = 10 , seed = 1,
		          auto.length = FALSE)
summary(fit1)
plot(fit1)

##
## Scenario 2: standard input with sub-population specification
##
data(RandomVA2)
fit2<- insilico(RandomVA2, subpop = list("sex"),  
              length.sim = 1000, burnin = 500, thin = 10 , seed = 1,
		          auto.length = FALSE)
summary(fit2)
plot(fit2, type = "compare")
plot(fit2, which.sub = "Men")

##
## Scenario 3: standard input with multiple sub-population specification
##
fit3<- insilico(RandomVA2, subpop = list("sex", "age"),  
              length.sim = 1000, burnin = 500, thin = 10 , seed = 1,
		          auto.length = FALSE)
summary(fit3)

##
## Scenario 3: standard input with multiple sub-population specification
##
fit3<- insilico(RandomVA2, subpop = list("sex", "age"),  
              length.sim = 1000, burnin = 500, thin = 10 , seed = 1,
		          auto.length = FALSE)
summary(fit3)

##
## Scenario 5 - 7 are special situations rarely needed in practice,
##   but included here for completeness. 
##   The below examples use no sub-population or physician codes, 
##   but specifying sub-population is still possible as in Scenario 2 - 4.
## 

##
## Scenario 5: skipping re-estimation of conditional probabilities
##
# Though in practice the need for this situation is very unlikely, 
# use only the default conditional probabilities without re-estimation
fit5<- insilico(RandomVA1, subpop = NULL,  
              length.sim = 1000, burnin = 500, thin = 10 , seed = 1,
              updateCondProb = FALSE, 
		          auto.length = FALSE) 
summary(fit5)

##
## Scenario 6: modify default conditional probability matrix
##
# Load the default conditional probability matrix 
data(condprob)
# The conditional probabilities are given in levels such as I, A+, A, A-, etc.
condprob[1:5, 1:5]
# To modify certain cells 
new_cond_prob <- condprob
new_cond_prob["elder", "HIV/AIDS related death"] <- "C"
# or equivalently
new_cond_prob[1, 3] <- "C"

fit6<- insilico(RandomVA1, subpop = NULL,  
              length.sim = 1000, burnin = 500, thin = 10 , seed = 1,
              CondProb = new_cond_prob, 
		          auto.length = FALSE) 
# note: compare this with fit1 above to see the change 
# induced by changing Pr(elder | HIV) from "C+" to "C".
summary(fit6)

##
## Scenario 7: modify default numerical values in conditional probabilities directly
##
# Load the default conditional probability matrix 
data(condprobnum)
# The conditional probabilities are given in numerical values in this dataset
condprobnum[1:5, 1:5]
# To modify certain cells, into any numerical values you want 
new_cond_prob_num <- condprobnum
new_cond_prob_num["elder", "HIV/AIDS related death"] <- 0.004
# or equivalently
new_cond_prob_num[1, 3] <- 0.005

fit7<- insilico(RandomVA1, subpop = NULL,  
              length.sim = 1000, burnin = 500, thin = 10 , seed = 1,
              CondProbNum = new_cond_prob_num, 
		          auto.length = FALSE) 
# note: compare this with fit1, fit5, and fit6
summary(fit7)

##
## Scenario 8: physician coding
## see also the examples in physician_debias() function section
##
# Load sample input for physicians
data(RandomPhysician)
# The symptom section looks the same as standard input
head(RandomPhysician[, 1:5])
# At the end of file, including a few more columns of physician id and coded cause
head(RandomPhysician[, 245:250])

# load Cause Grouping (if physician-coded causes are in larger categories)
data(SampleCategory)
head(SampleCategory)

# existing doctor codes in the sample dataset
doctors <- paste0("doc", c(1:15))
causelist <- c("Communicable", "TB/AIDS", "Maternal",
               "NCD", "External", "Unknown")
phydebias <- physician_debias(RandomPhysician, 
  phy.id = c("rev1", "rev2"), phy.code = c("code1", "code2"), 
  phylist = doctors, causelist = causelist, 
  tol = 0.0001, max.itr = 100)

fit8 <- insilico(RandomVA1, subpop = NULL,  
              length.sim = 1000, burnin = 500, thin = 10 , seed = 1,
              phy.debias = phydebias,
              phy.cat = SampleCategory, 
              phy.external = "External", phy.unknown = "Unknown",
		          auto.length = FALSE) 
summary(fit8)

#########################################################################
###      More test                                                   ####
#########################################################################
fit2.p <- insilico(RandomVA2, subpop = list("sex"),  
              length.sim = 1000, burnin = 500, thin = 10 , seed = 1,
              phy.debias = phydebias,
              phy.cat = SampleCategory, 
              phy.external = "External", phy.unknown = "Unknown",
		   auto.length = FALSE) 
summary(fit2.p)
fit3.p <- insilico(RandomVA2, subpop = list("sex", "age"),  
              length.sim = 1000, burnin = 500, thin = 10 , seed = 1,
              phy.debias = phydebias,
              phy.cat = SampleCategory, 
              phy.external = "External", phy.unknown = "Unknown",
		   auto.length = FALSE) 
summary(fit3.p)
fit4.p <- insilico(RandomVA2, subpop = list("sex"),  
              length.sim = 1000, burnin = 500, thin = 10 , seed = 1,
              updateCondProb = FALSE, 
		   phy.debias = phydebias,
              phy.cat = SampleCategory, 
              phy.external = "External", phy.unknown = "Unknown",
		   auto.length = FALSE) 
summary(fit4.p)
fit5.p <- insilico(RandomVA2, subpop = list("sex"),  
              length.sim = 1000, burnin = 500, thin = 10 , seed = 1,
              CondProb = new_cond_prob, 
		   phy.debias = phydebias,
              phy.cat = SampleCategory, 
              phy.external = "External", phy.unknown = "Unknown",
		   auto.length = FALSE) 
summary(fit5.p)
fit6.p <- insilico(RandomVA2, subpop = list("sex"),  
              length.sim = 1000, burnin = 500, thin = 10 , seed = 1,
              CondProb = new_cond_prob, updateCondProb = FALSE,
		   phy.debias = phydebias,
              phy.cat = SampleCategory, 
              phy.external = "External", phy.unknown = "Unknown",
		   auto.length = FALSE) 
summary(fit6.p)
fit7.p <- insilico(RandomVA2, subpop = list("sex"),  
              length.sim = 1000, burnin = 500, thin = 10 , seed = 1,
              CondProbNum = new_cond_prob_num, 
		   phy.debias = phydebias,
              phy.cat = SampleCategory, 
              phy.external = "External", phy.unknown = "Unknown",
		   auto.length = FALSE) 
summary(fit7.p)

#---------------------------------------------------------------------#
#-----------------------  insilico.plot ------------------------------#
#---------------------------------------------------------------------#

data(RandomVA1) 
##
## Scenario 1: without sub-population specification
##
fit1<- insilico(RandomVA1, subpop = NULL,  
              length.sim = 1000, burnin = 500, thin = 10 , seed = 1,
              auto.length = FALSE)
# basic line plot
plot(fit1)
# basic bar plot
plot(fit1, type = "bar")
# line plot with customized look
plot(fit1, top = 15, horiz = FALSE, fill = "gold", 
           bw = TRUE, title = "Top 15 CSMFs", angle = 70, 
           err_width = .2, err_size = .6, point_size = 2)

##
## Scenario 2: with sub-population specification
##
data(RandomVA2)
fit2<- insilico(RandomVA2, subpop = list("sex"),  
              length.sim = 1000, burnin = 500, thin = 10 , seed = 1,
              auto.length = FALSE)
summary(fit2)
# basic side-by-side line plot for all sub-populations
plot(fit2, type = "compare", main = "Top 5 causes comparison")
# basic line plot for specific sub-population
plot(fit2, which.sub = "Women", main = "Top 5 causes for women")
# customized plot with only specified causes
# the cause names need not be exact as InterVA cause list
# substrings in InterVA cause list is enough for specification
# e.g. the following two specifications are the same
some_causes_1 <- c("HIV/AIDS related death", "Pulmonary tuberculosis")
some_causes_2 <- c("HIV", "Pulmonary")
plot(fit2, type = "compare", horiz = FALSE,  causelist = some_causes_1,
              title = "HIV and TB fractions in two sub-populations", 
              angle = 20)


#---------------------------------------------------------------------#
#-----------------------  groupplot     ------------------------------#
#---------------------------------------------------------------------#

data(RandomVA1) 
##
## Scenario 1: without sub-population specification
##
fit1<- insilico(RandomVA1, subpop = NULL,  
              length.sim = 1000, burnin = 500, thin = 10 , seed = 1,
              auto.length = FALSE)
# stack bar plot for grouped causes
# the default grouping could be seen from
data(SampleCategory)
stackplot(fit1, type = "dodge", xlab = "")

##
## Scenario 2: with sub-population specification
##
data(RandomVA2)
fit2<- insilico(RandomVA2, subpop = list("sex"),  
              length.sim = 1000, burnin = 500, thin = 10 , seed = 1,
              auto.length = FALSE)
stackplot(fit2, type = "stack", angle = 0)
stackplot(fit2, type = "dodge", angle = 0)
# Change the default grouping by separating TB from HIV
data(SampleCategory)
SampleCategory[c(3, 9), ]
SampleCategory[3, 2] <- "HIV/AIDS"
SampleCategory[9, 2] <- "TB"
stackplot(fit2, type = "stack", grouping = SampleCategory, 
          sample.size.print = TRUE, angle = 0)
stackplot(fit2, type = "dodge", grouping = SampleCategory,
          sample.size.print = TRUE, angle = 0)

# then change the order of display for sub-population and cause groups
groups <- c("HIV/AIDS", "TB", "Communicable", "NCD", "External",
            "Maternal", "causes specific to infancy") 
subpops <- c("Women", "Men")
stackplot(fit2, type = "stack", grouping = SampleCategory, 
          order.group = groups, order.sub = subpops, 
          sample.size.print = TRUE, angle = 0)

#---------------------------------------------------------------------#
#-------------------------  get.indiv    -----------------------------#
#-------------------------  indiv.plot   -----------------------------#
#---------------------------------------------------------------------#
# Toy example with 1000 VA deaths
data(RandomVA1)
fit1<- insilico(RandomVA1, subpop = NULL,  
              length.sim = 1000, burnin = 500, thin = 10 , seed = 1,
              auto.length = FALSE)
summary(fit1, id = "d199")

# update credible interval for individual probabilities to 90%
indiv.new <- get.indiv(fit1, CI = 0.9)
fit1$indiv.prob.lower <- indiv.new$lower
fit1$indiv.prob.upper <- indiv.new$upper
fit1$indiv.CI <- 0.9
summary(fit1, id = "d199")


# get empirical aggregated COD distribution 
agg.csmf <- get.indiv(data = RandomVA2, fit1, CI = 0.95, 
                      is.aggregate = TRUE, by = NULL)
head(agg.csmf)

# aggregate individual COD distribution by sex and age
# note the model was fitted assuming the same CSMF for all deaths
# this aggregation provides an approximate CSMF for each sub-groups
agg.by.sex.age <- get.indiv(data = RandomVA2, fit1, CI = 0.95, 
                            is.aggregate = TRUE, by = list("sex", "age"))
head(agg.by.sex.age$mean)

# plot of aggregated individual COD distribution
# 0. plot for all data
indivplot(agg.csmf, top = 10)
# 1. plot for specific one group
indivplot(agg.by.sex.age, which.plot = "Men 60-", top = 10)
# 2. comparing multiple groups
indivplot(agg.by.sex.age, which.plot = list("Men 60+", "Men 60-"), 
                          top = 5)
# 3. comparing multiple groups on selected causes
indivplot(agg.by.sex.age, which.plot = list("Men 60-", "Women 60-"), 
                          top = 0, causelist = c(
                            "HIV/AIDS related death", 
                            "Pulmonary tuberculosis", 
                            "Other and unspecified infect dis", 
                            "Other and unspecified NCD"))

 
# head(agg.by.sex.age$"Men 60+")
# head(agg.by.sex.age$"Men 60-")

# round(agg.by.sex.age$"Men 60-"[c(3, 9,20), ], 4)
# round(agg.by.sex.age$"Men 60+"[c(3, 9,20), ], 4)
# round(agg.by.sex.age$"Women 60-"[c(3, 9,20), ], 4)
# round(agg.by.sex.age$"Women 60+"[c(3, 9,20), ], 4)

# round(summary(fit3)$csmf[["Men 60-"]][c(3, 9,20), ], 4)
# round(summary(fit3)$csmf[["Men 60+"]][c(3, 9,20), ], 4)
# round(summary(fit3)$csmf[["Women 60-"]][c(3, 9,20), ], 4)
# round(summary(fit3)$csmf[["Women 60+"]][c(3, 9,20), ], 4)

# par(mfrow = c(2, 3))
# plot(apply(fit1$csmf, 2, mean), apply(fit1$indiv.prob, 2, mean))
# abline(a=0,b=1,col="red")
# plot(apply(fit1$csmf, 2, quantile, 0.025), apply(fit1$indiv.prob.low, 2, mean))
# abline(a=0,b=1,col="red")
# plot(apply(fit1$csmf, 2, quantile, 0.975), apply(fit1$indiv.prob.up, 2, mean))
# abline(a=0,b=1,col="red")

# hist(apply(fit1$csmf, 2, mean)- apply(fit1$indiv.prob, 2, mean))
# hist(apply(fit1$csmf, 2, quantile, 0.025)-apply(fit1$indiv.prob.low, 2, mean))
# hist(apply(fit1$csmf, 2, quantile, 0.975)-apply(fit1$indiv.prob.up, 2, mean))


