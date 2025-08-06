library(dplyr)
library(janitor)
library(ggplot2)

library(edgeR)
library(limma)
library(ashr)
library(mashr)


pal <- c("chartreuse2", "darkblue")

pal2 <- c("darkorchid","chartreuse2", "darkblue")



# filter internal standards
to_exclude<- c("DG_12_0", "LPC_17_0", "LPE_17_1",
               "PC_25_0", 
               "Sphingosine_17_1","TG_17_0"
)


lipid_csv <- "data/PSU_RawData_MSDial_NewStdInt_240422.csv"
raw <- read.csv(lipid_csv,na.strings = c("","#N/A","NA","Inf")) %>%
  dplyr::select(-all_of(to_exclude)) %>%
  filter(!grepl("Methanol",sampleID)) %>%
  filter(!grepl("Pool",sampleID))
nrow(raw)

lipid_meta <- read.csv("data/lipid_metadata.csv")

plant_csv <- "../../data/22_NCS_PSU_LANGEBIO_FIELDS_PSU_P_field.csv"


psu <- read.csv(plant_csv) %>%
  rename(Genotype = "Who.What", rowid = "P22." ) %>%
  # dplyr::filter(Genotype %in% c("B73", "CTRL","INV4M",  "NUE")) %>%
  dplyr::filter(Genotype %in% c("CTRL","INV4M"))  %>%
  droplevels()

sampleInfo <- read.csv('data/PSU-PHO22_Metadata.csv') %>%
  rename(Genotype=genotype) %>% 
  rename(rowid=row)

sampleInfo <- psu  %>% 
  dplyr::select(rowid,Rep,Plot_Row, Plot_Column,DTA,DTS) %>%
  dplyr::inner_join(sampleInfo) %>% filter(!is.na(DTA))

sampleInfo$leaf_group <-NA

sampleInfo$leaf_group[sampleInfo$leaf_tissue <3] <- "apical"
sampleInfo$leaf_group[sampleInfo$leaf_tissue >2] <- "basal"

sampleInfo$tube  <- gsub("L|R","S",sampleInfo$tube, perl =TRUE)

sampleInfo$rowid  <-  sampleInfo$row

raw$tube <- gsub(".*_|.raw","",raw$sampleID, perl =TRUE)
raw$tube <- gsub("L|R","S",raw$tube, perl =TRUE)
colnames(raw)


pheno <- sampleInfo %>%
  dplyr::select(rowid,tube,Rep,Plot:Plot_Row,Treatment:leaf_tissue,leaf_group) %>%
  inner_join(
    raw %>%
      dplyr::select(tube,DGDG_34_0:TG_58_5)) %>% 
  mutate(block=as.factor(Rep)) %>%
  droplevels() %>%
  tidyr::pivot_wider(names_from = Rep, values_from = Rep, names_prefix = "block_",
                     values_fn = function(x) 1, values_fill = 0)  %>%
  dplyr::select(rowid, tube:Plot_Row,block,block_6:block_12,everything()) %>%
  mutate(Treatment = factor(Treatment,levels=c("Low_P","High_P")))
nrow(pheno)
colnames(pheno)
pheno$rowid

ratio_csv <- "data/rx2.csv"
rx <- read.csv(ratio_csv,na.strings = c("","#N/A","NA","Inf"))
ratios <- paste(rx[,1],rx[,2], sep =".")
nrow(rx)

for(num in rx[,1]){
  for(den in rx[,1]){
    pheno[paste(num,den, sep=".")] <- pheno[num]/pheno[den]
  }
}


m <- pheno%>%
  dplyr::select(DGDG_34_0:TG_58_5) %>%
  as.matrix() 
max_count <- max(m)
# m <- log2(pheno[,ratios] %>% as.matrix())

counts <-  pheno[,ratios] %>% as.matrix()%>%t()*max_count


counts[is.infinite(counts)] <- max_count
counts[is.na(counts)] <- max_count
colnames(counts)

dimnames(counts)
colnames(counts) <- pheno$tube

class(counts)
genes = data.frame(gene= rownames(counts))

sampleNames = pheno$tube
sampleNames %in% sampleInfo$tube
sampleInfo = sampleInfo[match(sampleNames,sampleInfo$tube),]

y = DGEList(counts = counts,samples = sampleInfo)
y$group = interaction(y$samples$Treatment,y$samples$Genotype)
y_filtered <-y
y_filtered_bySample <-y

d$y = voomR$E[gene,]
d$counts = y_filtered_bySample$counts[gene,]

design = model.matrix(~ Plot_Row + Plot_Column + Rep + leaf_tissue+ Treatment*Genotype,d)
voomR = voom(y_filtered_bySample,design=design,plot=T)


fit = lmFit(voomR)
ebfit = eBayes(fit)


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

for(coef in 5:7){
  
  effects = SEs = matrix(NA,nrow = nrow(y_filtered_bySample),ncol = 4)
  rownames(effects) = rownames(SEs) = rownames(y_filtered_bySample)
  
  for(x in 1:4) {
    y_filtered_by_leaf = y_filtered_bySample[,y_filtered_bySample$samples$leaf_tissue==x]
    y_filtered_by_leaf$group = interaction(y_filtered_by_leaf$samples$Genotype,y_filtered_by_leaf$samples$Treatment)
    # filter_genes = filterByExpr(y_filtered_by_leaf,group = y_filtered_by_leaf$group)
    # y_filtered_by_leaf = y_filtered_by_leaf[filter_genes,]
    d = y_filtered_by_leaf$samples
    design_interaction = model.matrix(~ Plot_Row + Plot_Column + Rep + Treatment*Genotype,d)
    y_filtered_by_leaf = calcNormFactors(y_filtered_by_leaf)
    voomR = voom(y_filtered_by_leaf ,design=design_interaction,plot=T)
    fit = lmFit(voomR)
    ebfit = eBayes(fit)
    ebfit$coefficients
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


coef_names <- colnames(head(ebfit$coefficients))
tb$coef_name <- coef_names[tb$coef]

tb %>%
  dplyr::select(coef, coef_name, everything())

tb %>%
  dplyr::select(coef, coef_name, everything()) %>%
  dplyr::group_by(coef) %>%
  dplyr::arrange(AIC) %>%
  dplyr::slice(1)

hist(get_log10bf(results$`7`$m.Vem_ed))



#------ plot_effects

coef="7"
conditions=1:2
coef_names 

get_ci_to_plot <- function(coef="5", conditions=1:4, coef_names =as.character(5:7)){
  coef_string <- coef_names[as.numeric(coef)]
  
  n_leaf <- get_n_significant_conditions(
    results[[coef]][["m.Vem_ed"]],
    thresh = 0.05,
    conditions = conditions)
  
  n_share <- length(conditions)
  
  bf <- get_log10bf(results[[coef]][["m.Vem_ed"]]) %>% as.numeric()
  names(bf) <-row.names(results[[coef]][["m.Vem_ed"]][["result"]][["lfsr"]])
  
  
  shared <- names(n_leaf[n_leaf == n_share])
  
  lfsr_leaf <- get_lfsr(results[[coef]][["m.Vem_ed"]])
  colnames (lfsr_leaf) <- paste0("leaf",1:4)
  
  # shared <- (lfsr_leaf[,"leaf1"] <0.05) &  (lfsr_leaf[,"leaf2"] < 0.05)
  
  lfsr <- lfsr_leaf[shared,]
  
  
  mn <- get_pm(results[[coef]][["m.Vem_ed"]])[shared,]
  is_shared <- (apply(sign(mn)[,conditions],1,sum) != 0 ) 
  if(any(is_shared )){ 
    shared <- shared[is_shared]
  }else{
    return(data.frame())
  }
  
  # some responses don't have variation :/
  
  colnames (mn) <- paste0("leaf",1:4)
  # drop false for nrow=1
  se <- get_psd(results[[coef]][["m.Vem_ed"]])[shared, ,drop = FALSE]/2 
  colnames (se) <- paste0("leaf",1:4)

  
  shared <- rownames(se)
  
  lfsr <- lfsr_leaf[shared, ,drop = FALSE]
  colnames (lfsr) <- paste0("leaf",1:4)
  mn <- get_pm(results[[coef]][["m.Vem_ed"]])[shared, ,drop = FALSE]/2
  colnames (mn) <- paste0("leaf",1:4)
  se <- get_psd(results[[coef]][["m.Vem_ed"]])[shared, ,drop = FALSE]/2
  colnames (se) <- paste0("leaf",1:4)
  
  conf.level <- 0.95
  ci.value <- -qnorm( ( 1 - conf.level ) / 2 )
  
  lower <- mn - ci.value  * se
  upper <- mn + ci.value * se
  
  colnames (lower) <- paste0("leaf",1:4)
  colnames (upper) <- paste0("leaf",1:4)
  
  cbind(
    mn %>% as.data.frame() %>%
      tibble::rownames_to_column("Response") %>% 
      tidyr::pivot_longer(
        cols =  paste0("leaf",1:4), 
        values_to = "logFC", names_to = "leaf"),
    lower %>% as.data.frame()  %>%
      tibble::rownames_to_column("Response") %>%
      tidyr::pivot_longer(
        cols =  paste0("leaf",1:4), 
        values_to = "lower", names_to = "leaf") %>%
      dplyr::select(lower),
    upper %>% as.data.frame()  %>%
      tibble::rownames_to_column("Response") %>%
      tidyr::pivot_longer(
        cols =  paste0("leaf",1:4), 
        values_to = "upper", names_to = "leaf")%>%
      dplyr::select(upper),
    lfsr %>% as.data.frame()  %>%
      tibble::rownames_to_column("Response") %>%
      tidyr::pivot_longer(
        cols =  paste0("leaf",1:4), 
        values_to = "lfsr", names_to = "leaf")%>%
      dplyr::select(lfsr)
  ) %>% mutate(predictor=coef_string) %>%
    filter(leaf  %in% paste0("leaf",conditions)) %>%
    # inner_join(lipid_meta, by=c(Response="col_name"))  %>%
    arrange(predictor,Response,leaf) 
}


my_coef <- c("Intercept","Plot_Row","Plot_Column","Rep",
             "-P","Inv4m","-P:Inv4m")
to_plot <- rbind(
  get_ci_to_plot("5",coef_names=my_coef, conditions=3:4),
  get_ci_to_plot("6",coef_names=my_coef,conditions=3:4),
  get_ci_to_plot("7",coef_names=my_coef,conditions=3:4)
) %>% ungroup %>%
  # mutate(head_group = factor(head_group, levels=c("glycolipid","phospholipid")) )%>%
  mutate(predictor = factor(predictor, levels=my_coef[5:7]))%>%
  mutate(alpha = as.numeric(factor(leaf))*0.5)%>%
  group_by(predictor,Response) %>%  
  filter(n()>1) %>% as.data.frame() %>%
  arrange(predictor, -abs(logFC)) %>% 
  mutate(.r = row_number()) %>% # Add a row number variable %>%
  mutate(Response = forcats::fct_reorder(Response,.r))  

pd = position_dodge(0.4)

quartz()
to_plot  %>%
  ggplot( aes(x     = logFC,
              y     = Response,
              group = leaf,
             # color = head_group,
            #  fill = head_group,
              shape= leaf)) +
  ggtitle("Differential lipid analysis older leaves (3 & 4)\n multiple hypothesis corrected lfsr < 0.05 ") +
  xlab("Effect (log2 Fold Change)") +
  geom_vline(xintercept = 0,lty =2)+
  geom_point(position= pd,
             size   = 2)+
  geom_errorbar(aes(xmin  =  upper,
                    xmax  =  lower),
                position= pd,
                width =  0.2,
                size  =  0.7) +
  guides( color = guide_legend(reverse = TRUE, title =""),
          shape = guide_legend(title =""),
          fill="none") +
  facet_wrap(.~predictor, ncol =4,scales="free_x")  +
  # scale_y_discrete(limits=rev) +
  scale_color_manual(values = c("blue","red"))+
  scale_fill_manual(values = c("blue","red"))+
  scale_shape_manual(values=c(17,25),
                     labels= c( "leaf 3","leaf 4")
  ) +
  theme_light(base_size =15) +
  theme(legend.position ="top",
        strip.background =element_rect(fill="white"),
        strip.text =  element_text(color = "black",size=15),
        axis.title.y=element_blank(),
        axis.text.y = element_text(hjust = 0, face = "bold"),
        plot.caption = element_text(hjust = 0))


to_plot <- rbind(
  get_ci_to_plot("5",coef_names=my_coef, conditions=1:2),
  get_ci_to_plot("6",coef_names=my_coef,conditions=1:2),
  get_ci_to_plot("7",coef_names=my_coef,conditions=1:2)
) %>% ungroup %>%  
  # mutate(head_group = factor(head_group, levels=c("glycolipid","phospholipid")) )%>%
  mutate(predictor = factor(predictor, levels=my_coef[5:7]))%>%
  mutate(alpha = as.numeric(lfsr <0.05))%>%
  group_by(predictor,Response) %>%  
  filter(n()>1) %>% as.data.frame() %>%
  arrange(predictor, -abs(logFC)) %>% 
  mutate(.r = row_number()) %>% # Add a row number variable %>%
  mutate(Response = forcats::fct_reorder(Response,.r))


to_plot %>% as.data.frame()
to_plot  %>%
  ggplot( aes(x     = logFC,
              y     = Response,
              group = leaf,
              shape= leaf)) +
  ggtitle("Differential lipid analysis younger leaves (1 & 2)\n multiple hypothesis corrected lfsr < 0.05 ") +
  xlab("Effect (log2 Fold Change)") +
  geom_vline(xintercept = 0,lty =2)+
  geom_point(position= pd,
             size   = 2)+
  geom_errorbar(aes(xmin  =  upper,
                    xmax  =  lower),
                position= pd,
                width =  0.2,
                size  =  0.7) +
  guides( color = guide_legend(reverse = TRUE, title =""),
          shape = guide_legend(title =""),
          fill="none") +
  facet_wrap(.~predictor, ncol =4,scales="free_x")  +
  # scale_y_discrete(limits=rev) +
  scale_color_manual(values = c("blue","red"))+
  scale_fill_manual(values = c("blue","red"))+
  scale_shape_manual(values=c(17,25),
                     labels= c( "leaf 1","leaf 2")
  ) +
  theme_light(base_size =15) +
  theme(legend.position ="top",
        strip.background =element_rect(fill="white"),
        strip.text =  element_text(color = "black",size=15),
        axis.title.y=element_blank(),
        axis.text.y = element_text(hjust = 0, face = "bold"),
        plot.caption = element_text(hjust = 0))


