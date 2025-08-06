library(edgeR)
library(dplyr)
library(ggplot2)
library(limma)
library(ashr)
library(mashr)

counts <- read.csv('data/inv4mRNAseq_gene_sample_exp.csv')
sampleInfo <- read.csv('data/PSU-PHO22_Metadata.csv') 

tag <- sampleInfo$side_tag
names(tag) <-sampleInfo$library

genes <- data.frame(gene = counts[,2])
counts <- as.matrix(counts[,-c(1:2)])

rownames(counts) <- genes$gene

sampleNames <- tag[colnames(counts)]
colnames(counts) <- sampleNames 

# sampleNames = sapply(sampleNames,function(x) strsplit(x,'_')[[1]][1])
sampleNames %in% sampleInfo$side_tag
sampleInfo <- sampleInfo[match(sampleNames,sampleInfo$side_tag),]

y <- DGEList(counts = counts,samples = sampleInfo)
y$group <- interaction(y$samples$Treatment,y$samples$Genotype)
keep <- filterByExpr(y,group = y$group)
y$samples$lowCount <- y$samples$lib.size < 2e7


y_filtered <- y[keep,]

quartz()
mds <- plotMDS(y_filtered,pch=21)

quartz()
plot(mds$x,mds$y,pch=21,bg = (as.factor(y_filtered$samples$lowCount)==T)+1, 
     main = "Expression MDS")
y_filtered_bySample = y_filtered[,!y_filtered$samples$lowCount]

y_filtered_bySample$samples
table(y_filtered_bySample$samples$Treatment,
      y_filtered_bySample$samples$leaf_tissue)
table(y_filtered_bySample$samples$Genotype,
      y_filtered_bySample$samples$leaf_tissue)

table(y_filtered_bySample$samples$Treatment,
      y_filtered_bySample$samples$Genotype,
      y_filtered_bySample$samples$leaf_tissue)

quartz()
mds2 <- plotMDS(y_filtered_bySample,pch=21,
              label = y_filtered_bySample$samples$side_tag)


plot(mds2$x,mds2$y,pch=21,bg = factor(y_filtered_bySample$samples$Treatment),col=0)
plot(mds2$x,mds2$y,col = factor(y_filtered_bySample$samples$Genotype))

d = y_filtered_bySample$samples
d$x = mds2$x
d$y = mds2$y
d$Treatment <- factor(d$Treatment)
levels(d$Treatment) <- c("+P","-P")

quartz()
ggplot(d,aes(x=x,y=y)) + geom_point(aes(color = Treatment))

quartz()
ggplot(d,aes(x=x,y=y)) + geom_point(aes(color = decimal_time ))



quartz()
ggplot(d,aes(x=x,y=y)) + geom_point(aes(color = COLLECTOR))

quartz()
ggplot(d,aes(x=x,y=y)) + geom_point(
  aes(color = row,
      shape = Treatment)) 

quartz()
ggplot(d,aes(x=x,y=y)) + geom_point(aes(color = Genotype))



quartz()
d %>%
  mutate(leaf =factor(leaf_tissue)) %>%
ggplot(aes(x=x,y=y)) + 
  xlab(paste0("dim1 ","(",round(100*mds2$var.explained[1],),"%)")) +
  ylab(paste0("dim1 ","(",round(100*mds2$var.explained[2],),"%)"))+
  geom_point(aes(fill = leaf,shape = Treatment),size=4) +
  scale_fill_viridis_d()+
  scale_shape_manual(values=c(24,21))+
  guides(   shape = guide_legend(title = "Treatment",order = 1,
                                 override.aes = list(size=7)),
    fill = guide_legend( title = "Leaf", order = 2, 
      override.aes = list(geom = "point", shape = 22,size=7))) +
  ggpubr::theme_classic2(base_size = 25) +
  theme(legend.box = "horizontal", 
        legend.spacing =  unit(0,"line"),
        legend.box.spacing =unit(0, "in"),
        legend.position = c(0.75,0.17))



design = model.matrix(~ Plot_Column + Plot_Row + leaf_tissue + Treatment*Genotype,d)
y_filtered_bySample = calcNormFactors(y_filtered_bySample)

voomR = voom(y_filtered_bySample,design=design,plot=T)
fit = lmFit(voomR)
ebfit = eBayes(fit)

ebfit$coefficient
coef_cols <- ebfit$coefficients %>% colnames()
# of_interest <-  c("leaf_tissue","Treatment-P","GenotypeINV4","leaf_tissue:GenotypeINV4","GenotypeINV4:Treatment-P")

of_interest <-  c("leaf_tissue","Treatment-P","GenotypeINV4","Treatment-P:GenotypeINV4")

is_interesting <- which(coef_cols %in% of_interest )

results <-list()

for(x in of_interest){
  r <- cbind(topTable(ebfit,coef = x,sort.by = 'none',n=Inf),
             data.frame(predictor=x)) %>%
    tibble::rownames_to_column("Response")
  cr <- qt(0.975, ebfit$df.residual + ebfit$df.prior) * ebfit$stdev.unscaled[,x] * sqrt(ebfit$s2.post)
  # Calculating  confidence interval
  r$upper  <- r$logFC +cr
  r$lower  <- r$logFC -cr
  results[[x]]<-r
}



gene_symbol <- read.table("data/gene_symbol.tab",quote="",header=TRUE, sep ="\t", na.strings = "")
gene_symbol$locus_name
nrow(gene_symbol)

effect_order <-  c("leaf_tissue","Treatment-P","GenotypeINV4","Treatment-P:GenotypeINV4")

effects <- results %>% dplyr::bind_rows() %>%
  mutate(predictor=factor(predictor, levels=effect_order)) %>%
  mutate(is_significant = adj.P.Val < 0.05) %>%
  mutate(upregulated = logFC > 0 ) %>%
  mutate(downregulated = logFC < 0 ) %>%
  mutate( 
    regulation = case_when( (is_significant & upregulated)~ "up",
                            (is_significant & downregulated) ~ "down",
                            .default = NA
    )) %>%
   left_join(gene_symbol, by =c(Response="gene_model")) %>%
   select(predictor, locus_symbol, Response,  everything())
effects$predictor %>% as.factor() %>% levels()


effects %>%
  filter(adj.P.Val< 0.05)   %>%
  mutate(sign=sign(logFC) %>% as.factor()) %>%
  group_by(predictor) %>%
  arrange(adj.P.Val) %>%
  dplyr::slice(1:2) %>% select(1:5,8)

# Effect size plot ----

quartz()
effects %>%
  rename(Response="gene") %>%
  #filter(predictor=="Treatment-P") %>%
  filter(adj.P.Val< 0.05)  %>%
  arrange(adj.P.Val) %>% 
  ggplot(aes(y=abs(logFC), x= predictor)) +
  ggbeeswarm::geom_quasirandom() +
  scale_y_log10()


quartz()
effects %>%
  rename(Response="gene") %>%
  # filter(predictor=="Treatment-P") %>%
  filter(adj.P.Val< 0.05)  %>%
  arrange(adj.P.Val) %>% 
  ggplot(aes(y= AveExpr, x= predictor)) +
  ggbeeswarm::geom_quasirandom() 


# mir3
# effects %>%
#   filter(predictor=="leaf_tissue")  %>%
#   filter(Response=="Zm00001eb068400")


effects %>%
  rename(gene="Response") %>%
  filter(predictor=="Treatment-P:GenotypeINV4") %>% 
  filter(adj.P.Val< 0.05)  %>%
  arrange(adj.P.Val) %>%
  dplyr::select(chr,gene,logFC,adj.P.Val,enzyme_descriptors) %>% tibble()


inv4m_ft <- effects %>%
  rename(gene="Response") %>%
  filter(predictor=="GenotypeINV4") %>%
  inner_join(ft_genes) %>%
  filter(adj.P.Val< 0.05)  %>%
  arrange(adj.P.Val) %>%
  dplyr::select(gene,logFC,adj.P.Val, annotation) %>% tibble() %>% print(n=50)



phospho_ft <- effects %>%
  rename(gene="Response") %>%
  filter(predictor=="Treatment-P") %>%
  inner_join(ft_genes) %>%
  filter(adj.P.Val< 0.05)  %>%
  arrange(adj.P.Val) %>%
  dplyr::select(gene,logFC,adj.P.Val, annotation) %>% tibble() %>% print(n=50)




found_genes <- effects %>%
  rename(gene="Response") %>%
  filter(predictor=="Treatment-P") %>%
  filter(adj.P.Val< 0.05) %>%
  filter(gene %in% phos_genes)  %>%
  arrange(adj.P.Val) %>% pull(gene)

effects %>%
  rename(gene="Response") %>%
  filter(predictor=="Treatment-P") %>%
  inner_join(lipid_genes %>% 
               dplyr::select(gene,v4,chr,pos,enzyme_descriptors,pathway_ids)) %>%
  filter(adj.P.Val< 0.05)  %>%
  arrange(adj.P.Val) %>%
  dplyr::select(chr,gene,logFC,adj.P.Val,enzyme_descriptors) 

 

effects %>%
  rename(gene="Response") %>%
  filter(predictor=="leaf_tissue") %>%
  inner_join(lipid_genes %>% 
               dplyr::select(gene,v4,chr,pos,enzyme_descriptors,pathway_ids)) %>%
  filter(adj.P.Val< 0.05)  %>%
  arrange(adj.P.Val) %>%
  dplyr::select(chr,gene,logFC,adj.P.Val,enzyme_descriptors) %>% tibble() %>%
  print(n=220)


# 
# design_interaction = model.matrix(~leaf_tissue*Treatment*Genotype,d)
# fit = lmFit(voomR,design_interaction)
# ebfit = eBayes(fit)
# 
# nrow(ebfit$coefficients)
# r = list()
# for(x in 2:4){
# r[[as.character(x)]] <- topTable(ebfit,coef = x,sort.by = 'none',n=Inf)
# }
# r$`4`[order(r$`4`$adj.P.Val),]
# ?topTable
# 
# ?topTable
# head(ebfit$coefficients)
# 
# topTable(ebfit,coef = 3)
# topTable(ebfit,coef = 4) 
# 
# cat(row.names(topTable(ebfit,coef = 2,number = 5000)))
# 
# gene='Zm00001eb053360'
# # gene='Zm00001eb003820' # PILNCR1-miR399
# # gene='Zm00001eb191650' # PHOS2
# 
# # P x Inv4m limma
 gene='Zm00001eb189920' # aldh2

# # P x Inv4m mashr shared 4
# # gene='Zm00001eb349450' # glk3
# gene='Zm00001eb194380' # Exocyst complex component SEC6

gene='Zm00001eb010130'
gene='Zm00001eb154820'
gene='Zm00001eb380260'
d$y = voomR$E[gene,]
d$counts = y_filtered_bySample$counts[gene,]

phos_genes <- rownames(voomR$E)[rownames(voomR$E) %in% phosphorus_genes$gene]

voomR$E[phos_genes,]
quartz()
heatmap(voomR$E[phos_genes,], col= col)

quartz()
heatmap(voomR$E[found_genes,], col= col)

lg_in_exp<- rownames(voomR$E)[rownames(voomR$E) %in% lipid_genes$gene] 

point_col <- RColorBrewer::brewer.pal(11, "RdBu")

col = colorRampPalette(colors = point_col , space="Lab")(25)

quartz()
heatmap(voomR$E[lg_in_exp,], col= col)

quartz()
heatmap(voomR$E[inv4m_genes,], col= col)

quartz()
heatmap(voomR$E[leaf_genes,], col= col)

quartz()
heatmap(voomR$E[inv4m_ft$gene,], col= col)

quartz()
heatmap(voomR$E[phospho_ft$gene,], col= col)

lipid_ft


voomR$E[,]
quartz()
ggplot(d,aes(x=factor(leaf_tissue),y=y, color=Treatment)) +
  ggtitle(gene) +
  geom_boxplot() + 
  facet_wrap(~Genotype) 

# ggplot(d,aes(x=Treatment,y=counts)) + geom_boxplot(aes(color = Genotype,group = interaction(Treatment,Genotype))) + facet_wrap(~leaf_tissue)
ggplot(d,aes(x=Genotype,y=y))  +
  ggtitle(gene) +
  geom_boxplot(aes(color = Treatment,group = interaction(Treatment,Genotype))) + 
  facet_wrap(~leaf_tissue)




########################################################
# by leaf analysis #####################################
########################################################

# make design in each = ~Treatment*Genotype ############
# check effects for interaction




# Akaike information content calculation
aic <- function(x){
  # x is a mashr model object
  # I am assuming x$fitted_g$pi, the number of mixture components,
  # is the number of parameters 
  2*length(x$fitted_g$pi) - 2*max(get_loglik(x))
}


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


barplot(get_estimated_pi(results$`3`$m.Vem_ed),las = 2)
print(get_pairwise_sharing(results$`4`$m.Vem_ed))



regulated <- rownames(deg_lfsr)


universe <- rownames(get_lfsr(results$`2`$m.Vem_ed))


list_ego_results <- ego_analysis(regulated,universe)

lapply(list_ego_results, function(x) sum(x@result$p.adjust < 0.05))

quartz()
dotplot(list_ego_results$ego_BP, x="FoldEnrich",showCategory=15, title="BP")

quartz()
dotplot(list_ego_results$ego_MF,x="FoldEnrich", showCategory=10, title="MF",)

quartz()
dotplot(list_ego_results$ego_CC,x="FoldEnrich",showCategory=10, title="CC")


