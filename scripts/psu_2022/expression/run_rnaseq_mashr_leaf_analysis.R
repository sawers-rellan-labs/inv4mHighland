
## ------------------------------------------------------------------
# Step 4: Run mash for different models
# using just the simple canonical covariances as in the initial introductory vignette.
# correcting for correlations between samples inside a "condition"
# 



tb <- data.frame()
results <- list()

for(coef in 2:4){
  
  effects = SEs = matrix(NA,nrow = nrow(y_filtered_bySample),ncol = 4)
  rownames(effects) = rownames(SEs) = rownames(y_filtered_bySample)
  
  for(x in 1:4) {
    y_filtered_by_leaf = y_filtered_bySample[,y_filtered_bySample$samples$leaf_tissue==x]
    y_filtered_by_leaf$group = interaction(y_filtered_by_leaf$samples$Genotype,y_filtered_by_leaf$samples$Treatment)
    # filter_genes = filterByExpr(y_filtered_by_leaf,group = y_filtered_by_leaf$group)
    # y_filtered_by_leaf = y_filtered_by_leaf[filter_genes,]
    d = y_filtered_by_leaf$samples
    design_interaction = model.matrix(~Treatment*Genotype,d)
    y_filtered_by_leaf = calcNormFactors(y_filtered_by_leaf)
    voomR = voom(y_filtered_by_leaf ,design=design_interaction,plot=T)
    fit = lmFit(voomR)
    ebfit = eBayes(fit)
    tt <- topTable(ebfit,coef = coef,sort.by = 'none',n=Inf)
    effects[,x] = tt$logFC
    SEs[,x] = tt$logFC/tt$t
  }
  
  
  
  # Step 2: Obtain initial data-driven covariance matrices
  
  data = mash_set_data(as.matrix(effects), as.matrix(SEs))
  U.c = cov_canonical(data)
  
  # Step 2: Obtain initial data-driven covariance matrices
  
  # select strong signals
  m.1by1 = mash_1by1(data)
  strong = get_significant_results(m.1by1,0.05)
  if(length(strong)<20) {
    strong = order(apply(m.1by1$result$lfsr,1,min))[1:min(20,nrow(effects))]
  }
  
  # Perform PCA on data and return list of candidate covariance matrices
  U.pca = cov_pca(data,npc=ncol(effects),subset=strong)
  # npc:	the number of PCs to use, should be less than or equal to n_conditions(data)
  # subset: indices of the subset of data to use (set to NULL for all data)
  # print(names(U.pca))
  
  ## ----------------------------------------------------------------
  # Step 3: prepare canonical/data-driven covariance matrices
  # Perform "extreme deconvolution" (Bovy et al) on a subset of the data
  U.ed = cov_ed(data, U.pca, subset=strong)
  # subset: a subset of data to be used when ED is run (set to NULL for all the data)
  # The function cov_ed is used to apply the ED algorithm from a specified initialization
  # (here U.pca) and to a specified subset of signals.
  
  
  
  # Why U.c and U.ed?  They are supposed to be different ways of estimating the conditions covariance matrix
  V.em_c_ed = mash_estimate_corr_em(data, c(U.c,U.ed), details = TRUE)
  m.Vem_c_ed = V.em_c_ed$mash.model
  m.Vem_c_ed$result$NAs = is.na(effects)
  m.Vem_c_ed$V = V.em_c_ed$V
  
  print(get_loglik(m.Vem_c_ed),digits=10) 
  length(m.Vem_c_ed$fitted_g$pi)
  aic(m.Vem_c_ed)
  
  results[[as.character(coef)]][["m.Vem_c_ed"]] <- m.Vem_c_ed
  
  
  ###
  V.em_ed = mash_estimate_corr_em(data, U.ed, details = TRUE)
  m.Vem_ed = V.em_ed$mash.model
  m.Vem_ed$result$NAs = is.na(effects)
  
  quartz()
  mash_plot_meta(m.Vem_ed,2)
  
  m.Vem_ed$V = V.em_ed$V
  
  print(get_loglik(m.Vem_ed),digits=10) 
  length(m.Vem_ed$fitted_g$pi)
  aic(m.Vem_ed)
  
  results[[as.character(coef)]][["m.Vem_ed"]] <- m.Vem_ed
  
  ###
  V.em_c = mash_estimate_corr_em(data, U.c, details = TRUE)
  m.Vem_c = V.em_c$mash.model
  m.Vem_c$result$NAs = is.na(effects)
  m.Vem_c$V = V.em_c$V
  
  print(get_loglik(m.Vem_c),digits=10) 
  length(m.Vem_c$fitted_g$pi)
  aic(m.Vem_c)
  
  results[[as.character(coef)]][["m.Vem_c"]] <- m.Vem_c
  ###
  m.c = mash(data,U.c)
  # Fitting model with 217 mixture components
  # -57582
  
  print(get_loglik(m.c),digits=10)
  length(m.c$fitted_g$pi)
  aic(m.c)
  
  
  ###
  m.c_ed  = mash(data, c(U.c,U.ed))
  
  # Fitting model with 337 mixture components.
  
  print(get_loglik(m.c_ed),digits=10)
  length(m.c_ed$fitted_g$pi)
  aic(m.c_ed)
  
  results[[as.character(coef)]][["m.c_ed"]] <- m.c_ed
  ####
  # Compare likelihood of the  models
  
  ## Compare model fit
  
  
  model =  c('canonical',
             'canonical + ed',
             'canonical + Vem',
             'ed + Vem', 
             'canonical + ed + Vem')
  
  significant = c(length(get_significant_results(m.c)), 
                  length(get_significant_results(m.c_ed)),
                  length(get_significant_results(m.Vem_c)),
                  length(get_significant_results(m.Vem_ed)),
                  length(get_significant_results(m.Vem_c_ed)))
  
  loglike = c(
    max(get_loglik(m.c)),
    max(get_loglik(m.c_ed)),
    max(get_loglik(m.Vem_c)),
    max(get_loglik(m.Vem_ed)),
    max(get_loglik(m.Vem_c_ed))
  )
  
  k = c(length(m.c$fitted_g$pi), 
        length(m.c_ed$fitted_g$pi), 
        length(m.Vem_c$fitted_g$pi),
        length(m.Vem_ed$fitted_g$pi), 
        length(m.Vem_c_ed$fitted_g$pi))
  
  AIC = c(aic(m.c), aic(m.c_ed), aic(m.Vem_c), aic(m.Vem_ed), aic(m.Vem_c_ed))
  
  tb <- rbind ( tb,
                data.frame(
                  coef = coef,
                  model = model,
                  significant = significant,
                  loglike = loglike,
                  k = k,
                  AIC = AIC)
  )
}

# write.csv(tb,"by_leaf_mashr_tb.csv")
# saveRDS(results, "by_leaf_results.RDS")

coef_names <- colnames(head(ebfit$coefficients))
tb$coef_name <- coef_names[tb$coef]

library(dplyr)
tb %>%
  dplyr::select(coef, coef_name, everything())

tb %>%
  dplyr::select(coef, coef_name, everything()) %>%
  dplyr::group_by(coef) %>%
  dplyr::arrange(AIC) %>%
  dplyr::slice(1)

hist(get_log10bf(results$`4`$m.Vem_ed))

length(get_significant_results(results$`4`$m.Vem_ed))


get_pairwise_sharing(results$`4`$m.Vem_ed)
# Get effects significant on all conditions (leafs)
# Signs and magnitude might be different
n_leaf <- get_n_significant_conditions(results$`4`$m.Vem_ed,)

bf[names(n_leaf[n_leaf >2])] %>% sort() %>% tail()

quartz()

bf <- get_log10bf(results$`4`$m.Vem_ed) %>% as.numeric()
names(bf) <-row.names(results$`4`$m.Vem_ed$result$lfsr)
sort(bf) %>% tail()
str(results$`4`$m.Vem_ed)
get_lfsr(results$`4`$m.Vem_ed)


deg_idx <- get_significant_results(results$`4`$m.Vem_ed)
deg <- names(deg_idx)
deg_1 <- deg
length(deg)
paste(deg, collapse =  ",")
deg_lfsr <- get_lfsr(results$`4`$m.Vem_e)[deg_idx,]
deg_lfsr[order(deg_lfsr[,2]),]

head(get_pm(results$`4`$m.Vem_ed)[deg_idx,])

quartz()
barplot(get_estimated_pi(results$`4`$m.Vem_ed),las = 2)

print(get_pairwise_sharing(results$`4`$m.Vem_ed)) 


deg_idx <- get_significant_results(results$`2`$m.Vem_ed)
deg <- names(deg_idx)


paste(dd, collapse =  ",")
deg_lfsr <- get_lfsr(results$`2`$m.Vem_e)[deg_idx,]
deg_lfsr[order(deg_lfsr[,2]),]

deg_idx <- get_significant_results(results$`2`$m.Vem_ed)
deg <- names(deg_idx)
deg_lfsr <- get_lfsr(results$`2`$m.Vem_ed)[deg_idx,]
deg_lfsr[order(deg_lfsr[,2]),]

