#mashr tutorial

library(ashr)
library(mashr)
set.seed(1)
simdata = simple_sims(500,5,1)


###############################################################################
#### Canonical covariances                                                 ####
###############################################################################

# Inputs are effect matrix (Bhat) and standard error matrix (Shat)
data = mash_set_data(simdata$Bhat, simdata$Shat)

# Make covariance matrix 

U.c = cov_canonical(data)  
print(names(U.c))
# .c stands for canonical. i.e it postulates a certain covariance structure
# and fits the model accordingly to a mixture over these matrices giving a weight wl
# to each one
# The alternative is to estimate the variance covariance structure from the data
# i.e data driven aproach
# identity -> single value for variance no covariance ammong conditions
# equal effects -> variance = covariance between all conditions
# het_1...het_3 single variance convarance multipled by different factors (1/4, 1/2. 3/4) 


# Fit the model
m.c = mash(data, U.c)

# Extract Posterior Summaries (get statistics)

# local false sign rate
head(get_lfsr(m.c))

# posterior mean
head(get_pm(m.c))

# (posteriore standard deviation)
head(get_psd(m.c))


# Use get_significant_results to find the indices of effects that are “significant”

head(get_significant_results(m.c))

print(length(get_significant_results(m.c)))

print(head(get_significant_results(m.c, conditions=1)))


#sharing of effects ammong conditions (same sign, factor up to .5)

print(get_pairwise_sharing(m.c)) 


# just same sign

print(get_pairwise_sharing(m.c, factor=0))

# log likelihood of teh fitted model
print(get_loglik(m.c))


# Estimated mixture proportions
print(get_estimated_pi(m.c))
barplot(get_estimated_pi(m.c),las = 2)


# metaplot for the most significant effect
par(bg = "white")
plot.new() 
mash_plot_meta(m.c,get_significant_results(m.c)[1])

###############################################################################
#### Data-driven covariances                                               ####                 
###############################################################################


# Step 1: select strong signals

# here we select the strong signals as those with lfsr<0.05 in any condition in the 1by1 analysis.

data   = mash_set_data(simdata$Bhat, simdata$Shat)
m.1by1 = mash_1by1(data)
strong = get_significant_results(m.1by1,0.05)

# Obtain initial data-driven covariance matrices
# this will serve as inicialization to esmatte B instead Bhat (?)
U.pca = cov_pca(data,5,subset=strong)
print(names(U.pca))


# Apply Extreme Deconvolution
# estimation of Betas and their covariance structure

U.ed = cov_ed(data, U.pca, subset=strong)

# Run mash
# Now that we have the covariance structure
m.ed = mash(data, U.ed)
print(get_loglik(m.ed),digits = 10)

# compare the likelihood with tha canonical covariance model
# using canonical and data driven covariance structure


U.c = cov_canonical(data)  
m   = mash(data, c(U.c,U.ed))
print(get_loglik(m),digits = 10)


# pick the best? (overfitting? AIC?)

###############################################################################
#### Accounting for correlations among measurements                        ####                 
###############################################################################

### Method 1: estimate_null_correlation_simple

# simulate data with correlations.
V = matrix(0.5,5,5)
diag(V) = 1
simdata$Bhat = simdata$B + mvtnorm::rmvnorm(2000, sigma = V)

# Read in the data, and estimate correlations:
data   = mash_set_data(simdata$Bhat, simdata$Shat)
V.simple = estimate_null_correlation_simple(data)
data.Vsimple = mash_update_data(data, V=V.simple) 


# mashr with simple canonical covariances as in the
U.c = cov_canonical(data.Vsimple) 
# fits with correlations because data.V includes correlation information 
m.Vsimple = mash(data.Vsimple, U.c) 
# log-likelihood of the fit with correlations set to V
print(get_loglik(m.Vsimple),digits=10) 

# fit without correlations because data object was set up without correlations
m.orig = mash(data, U.c) 

print(get_loglik(m.orig),digits=10)


## Compare model fit

loglik = c(get_loglik(m.orig), get_loglik(m.Vsimple))
significant = c(length(get_significant_results(m.orig)), length(get_significant_results(m.Vsimple)))
false_positive = c(sum(get_significant_results(m.orig) < 501), 
                   sum(get_significant_results(m.Vsimple) < 501))
tb = rbind(loglik, significant, false_positive)
colnames(tb) = c('without cor', 'V simple')
row.names(tb) = c('log likelihood', '# significance', '# False positive')
tb

### Method 2: mash_estimate_corr_em
# Initialize with simple canonical covariances (?)
#
# With details = TRUE in mash_estimate_corr_em, 
# it returns the estimates residual correlation matrix with the mash fit.
V.em = mash_estimate_corr_em(data, U.c, details = TRUE)

m.Vem = V.em$mash.model
print(get_loglik(m.Vem),digits=10) # log-likelihood of the fit


## Compare model fit

loglik = c(get_loglik(m.orig), get_loglik(m.Vsimple), get_loglik(m.Vem))
significant = c(length(get_significant_results(m.orig)), length(get_significant_results(m.Vsimple)),
                length(get_significant_results(m.Vem)))
false_positive = c(sum(get_significant_results(m.orig) < 501), 
                   sum(get_significant_results(m.Vsimple) < 501),
                   sum(get_significant_results(m.Vem) < 501))
tb = rbind(loglik, significant, false_positive)
colnames(tb) = c('without cor', 'V simple', 'V EM')
row.names(tb) = c('log likelihood', '# significance', '# False positive')
tb









