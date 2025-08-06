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

internal_standards<- c(
  "CUDA",
  "Cholesterol",
  "CHOLESTEROL_D7_H20",
  "C17_CERAMIDE_H2O",
  "DG_12_0",
  "DG_18_1",
  "DG_18_1_2_0_0_0",
  "LPC_17_0",
  "LPE_17_1",
  "PC_25_0",
  "PE_34_0",
  "PG_34_0",
  "PG_17_0_17_0",
  "PG_17_0",
  "SM_35_1",
  "Sphingosine_17_1",
  "Sphingosine_D_17_1",
  "SPHINGOSINE_D_17_1",
  "TG_17_0",
  "TG_17_0_17_1_17_0_D5"
)


lipid_csv <- "data/PSU_RawData_MSDial_NewStdInt_240422.csv"
raw <- read.csv(lipid_csv,na.strings = c("","#N/A","NA","Inf")) 

to_exclude <-internal_standards[internal_standards %in% colnames(raw)]

raw  <- raw %>%  dplyr::select(-all_of(to_exclude)) %>%
  filter(!grepl("Methanol",sampleID)) %>%
  filter(!grepl("Pool",sampleID))
nrow(raw)
ncol(raw)


lipid_meta <- read.csv("data/lipid_metadata.csv") 

ms_order<- read.csv("data/PSU-PHO22 _ms_order.csv")
ms_order$tube <- gsub("L|R","S",ms_order$top_tag, perl =TRUE)



psu <- read.csv( "data/22_NCS_PSU_LANGEBIO_FIELDS_PSU_P_field.csv") %>%
  rename(Genotype = "Who.What", row = "P22." ) %>%
  # dplyr::filter(Genotype %in% c("B73", "CTRL","INV4M",  "NUE")) %>%
  dplyr::filter(Genotype %in% c("CTRL","INV4M"))  %>%
  droplevels()

metadata <- read.csv('data/PSU-PHO22_Metadata.csv')
   
metadata$tube <- gsub("L|R","S",metadata$tube, perl =TRUE)
colnames(metadata)


sampleInfo <- psu  %>% 
  dplyr::select(row,Rep,Plot_Row, Plot_Column,DTA,DTS) %>%
  dplyr::inner_join(metadata) %>% filter(!is.na(DTA)) %>%
  dplyr::right_join(ms_order) %>% distinct()
nrow(sampleInfo)
colnames(sampleInfo)


sampleInfo$leaf_node <- 2*sampleInfo$leaf_tissue -1 
sampleInfo$leaf_group <-NA
sampleInfo$leaf_group[sampleInfo$leaf_tissue <3] <- "apical"
sampleInfo$leaf_group[sampleInfo$leaf_tissue >2] <- "basal"


sampleInfo$rowid  <-  sampleInfo$row

raw$tube <- gsub(".*_|.raw","",raw$sampleID, perl =TRUE)
raw$tube <- gsub("L|R","S",raw$tube, perl =TRUE)
colnames(raw)


pheno <- sampleInfo %>%
  dplyr::select(tube,row,Rep,Plot:Plot_Row,Treatment:leaf_tissue,leaf_node,leaf_group) %>%
  inner_join(
    raw %>%
      dplyr::select(tube,DGDG_34_0:TG_58_5)) %>% 
  mutate(block=as.factor(Rep)) %>%
  droplevels() %>%
  tidyr::pivot_wider(names_from = Rep, values_from = Rep, names_prefix = "block_",
                     values_fn = function(x) 1, values_fill = 0)  %>%
  dplyr::select(tube:Plot_Row,block,block_6:block_12,everything()) %>%
  mutate(Treatment = factor(Treatment,levels=c("Low_P","High_P"))) %>% unique()


nrow(pheno)
colnames(pheno)

pheno %>%
  arrange(tube) %>% as.data.frame()




# MDS with raw data -----

colnames(pheno)

m <- pheno %>%
  dplyr::select(DGDG_34_0:TG_58_5) %>%
  as.matrix() 

m[is.na(m)] <- 0

dim(m)
diag(m)
rownames(m) <- pheno$tube

# hist(m)

# n_factor <- mean(m)/rowSums(m)
# dim(m)
# length(n_factor)
# 
# norm_counts <- round(diag(n_factor) %*% m,0) %>% t()
# colnames(norm_counts) <- pheno$tube
# 
# norm <-cbind(
#   data.frame(tube=rownames(m)),
#   t(norm_counts)) # maximum_signal/total_lipid_sum persample
# 
# 
# counts <-  norm %>%
#   dplyr::select(DGDG_34_0:TG_58_5) %>%
#   as.matrix() %>% t()
# colnames(counts) <- pheno$tube
# counts <- m

counts <- t(m)
colnames(counts) <- pheno$tube



class(counts)
genes = data.frame(gene= rownames(counts))

sampleNames = pheno$tube
sampleNames %in% sampleInfo$tube
sampleInfo = sampleInfo[match(sampleNames,sampleInfo$tube),]
#sampleInfo$LPC17_0_group <- sampleInfo$LPC17_0_group %>% as.numeric() %>% as.factor() %>% droplevels()
sampleInfo$Rep <- sampleInfo$Rep %>% as.numeric() %>% as.factor() %>% droplevels()
y = DGEList(counts = counts,samples = sampleInfo)
y$group = interaction(y$samples$Treatment,y$samples$Genotype)


quartz()
mds = plotMDS(y,pch=21,label = y$samples$side_tag)

d = y$samples
d$x = mds$x
d$y = mds$y


quartz()
d%>%
  # filter(sample_label !="Pool_6") %>%
  ggplot(aes(x=x,y=y, label=tube,color = Treatment )) +
  ggtitle( ggtitle("Multidimensional scaling\nMS-DIAL raw, Internal Standard QC,\nPSU Inv4m Phosphorus experiment") )+
  geom_point()+
  ggrepel::geom_label_repel(max.overlaps = 25,force_pull = 5) +
  scale_color_viridis_d()


# Modelling ---


design = model.matrix(~  injection_order + Plot_Column + Plot_Row + leaf_node+Treatment*Genotype,d)
voomR = voom(y,design=design,plot=T)

fit = lmFit(voomR)
ebfit= eBayes(fit)
# resids <- limma::removeBatchEffect(d, batch = "LPC17_0_group")
head(ebfit$coefficients)


quartz()
hist(ebfit$p.value[,"TreatmentLow_P:GenotypeINV4"])



classifyTestsF(ebfit,p.value=0.05)


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
ebfit$coefficients
of_interest <- c("TreatmentLow_P","GenotypeINV4","TreatmentLow_P:GenotypeINV4")

tb <- data.frame()
results <- list()

for(coef in of_interest){
  effects = SEs = matrix(NA,nrow = nrow(y),ncol = 4)
  rownames(effects) = rownames(SEs) = rownames(y)
  
  for(x in 1:4) {
    y_filtered_by_leaf = y[,y$samples$leaf_tissue==x]
    y_filtered_by_leaf$group = interaction(y_filtered_by_leaf$samples$Genotype,y_filtered_by_leaf$samples$Treatment)
    d = y_filtered_by_leaf$samples
    design_interaction = model.matrix(~  injection_order + Plot_Column + Plot_Row + Treatment*Genotype,d)
    y_filtered_by_leaf = calcNormFactors(y_filtered_by_leaf)
    voomR = voom(y_filtered_by_leaf ,design=design_interaction,plot=T)
    fit = lmFit(voomR)
    ebfit = eBayes(fit)
    colnames(ebfit$coefficients)
    tt <- topTable(ebfit,coef = coef,sort.by = 'none',n=Inf)
    effects[,x] = tt$logFC
    SEs[,x] = tt$logFC/tt$t
  }
  effects
  SEs

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
tb$coef_name <- tb$coef

tb %>%
  dplyr::select(coef, coef_name, everything())

tb %>%
  dplyr::select(coef, coef_name, everything()) %>%
  dplyr::group_by(coef) %>%
  dplyr::arrange(AIC) %>%
  dplyr::slice(1)

quartz()
hist(get_log10bf(results$TreatmentLow_P$m.Vem_ed))

# ed + Vem it is

length(get_significant_results(results$TreatmentLow_P$m.Vem_ed))

# get confidence intervals to plot effects

# Remember limma model you could calculate
# effects
# SEs
effects
# Get effects significant on at least 3 leaves
# Signs and magnitude might be different


get_n_significant_conditions(
  results[["TreatmentLow_P"]][["m.Vem_ed"]],
  thresh = 0.05,
  conditions = c(1,2),
  sig_fn = get_lfsr
) 



get_ci_to_plot <- function(coef="TreatmentLow_P", conditions=1:4, coef_names =as.character(5:7)){

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

lfsr <- lfsr_leaf[shared,,drop = FALSE]


mn <- get_pm(results[[coef]][["m.Vem_ed"]])[shared,,drop = FALSE]
is_shared <- (apply(sign(mn)[,conditions,drop = FALSE],1,sum) != 0 ) 
if(any(is_shared )){ 
  shared <- shared[is_shared]
 }else{
    return(data.frame())
}

# some responses don't have variation :/

colnames (mn) <- paste0("leaf",1:4)
se <- get_psd(results[[coef]][["m.Vem_ed"]])[shared,,drop = FALSE]
colnames (se) <- paste0("leaf",1:4)

shared <- rownames(se)

lfsr <- lfsr_leaf[shared,,drop = FALSE]
mn <- get_pm(results[[coef]][["m.Vem_ed"]])[shared,,drop = FALSE]
colnames (mn) <- paste0("leaf",1:4)
se <- get_psd(results[[coef]][["m.Vem_ed"]])[shared,,drop = FALSE]
colnames (se) <- paste0("leaf",1:4)

conf.level <- 0.95
ci.value <- -qnorm( ( 1 - conf.level ) / 2 )

lower <- mn - ci.value  * se
upper <- mn + ci.value * se

colnames (lower) <- paste0("leaf",1:4)
colnames (upper) <- paste0("leaf",1:4)

leaf_nodes <- c("node1","node3","node5","node7")
names(leaf_nodes) <- c("leaf1","leaf2","leaf3","leaf4")

leaf_groups <- c("top","top","bottom","bottom")
names(leaf_groups) <- c("leaf1","leaf2","leaf3","leaf4")

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
) %>% mutate(predictor=coef) %>%
  filter(leaf  %in% paste0("leaf",conditions)) %>%
  mutate(leaf_node=leaf_nodes[leaf]) %>%
  mutate(leaf_group=leaf_groups[leaf]) %>%
  inner_join(lipid_meta, by=c(Response="col_name"))  %>%
  arrange(head_group,predictor,Response,leaf) 
}

lipid_meta


to_plot  <-NULL
to_plot1 <- rbind(
  get_ci_to_plot("TreatmentLow_P",coef_names=coef_names, conditions=1:2),
  get_ci_to_plot("GenotypeINV4",coef_names=coef_names,conditions=1:2),
  get_ci_to_plot("TreatmentLow_P:GenotypeINV4",coef_names=coef_names,conditions=1:2)
) %>% ungroup %>%  
  mutate(head_group = factor(head_group, levels=c("phospholipid","glycolipid","neutral")) )%>%
  mutate(predictor = factor(predictor, levels=of_interest))%>%
  group_by(head_group,predictor,Response) %>%  
  filter(n()>1) %>% as.data.frame() %>%
  mutate(class = forcats::fct_reorder(class,-abs(logFC)))  %>%
  arrange(head_group,class,-abs(logFC)) %>% 
  mutate(.r = row_number()) %>% # Add a row number variable %>%
  mutate(Response = forcats::fct_reorder(Response,.r))  %>%
  mutate(class= forcats::fct_reorder(class,.r)) 


quartz()
to_plot  %>%
  ggplot( aes(x     = logFC,
              y     = Response,
              group = leaf,
              color = head_group,
              fill  = head_group,
              shape = leaf)) +
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
  facet_grid(class~predictor, scales="free", space="free_y")  +
  scale_color_manual(values = c("red","blue","orange"))+
  scale_fill_manual(values = c("red","blue","orange"))+
  scale_shape_manual(values=c(17,25),
                    labels= c( "leaf 1","leaf 2") )  +
  guides(shape=guide_legend(title="",override.aes=list(fill="black",col="black"))) +
  scale_y_discrete(limits=rev) +
  theme_light(base_size =15) +
  theme(legend.position ="top",
        strip.background =element_rect(fill="white"),
        strip.text =  element_text(color = "black",size=15),
        strip.text.y= element_text(color = "white",size=15),
        axis.title.y=element_blank(),
        axis.text.y = element_text(hjust = 0, face = "bold"),
        plot.caption = element_text(hjust = 0))


to_plot  <-NULL
to_plot2 <- rbind(
  get_ci_to_plot("TreatmentLow_P",coef_names=coef_names, conditions=3:4),
  get_ci_to_plot("GenotypeINV4",coef_names=coef_names,conditions=3:4),
  get_ci_to_plot("TreatmentLow_P:GenotypeINV4",coef_names=coef_names,conditions=3:4)
) %>% ungroup %>%
  mutate(head_group = factor(head_group, levels=c("phospholipid","glycolipid","neutral")) )%>%
  mutate(predictor = factor(predictor, levels=of_interest))%>%
  group_by(head_group,predictor,Response) %>%  
  filter(n()>1) %>% as.data.frame() %>%
  mutate(class = forcats::fct_reorder(class,-abs(logFC)))  %>%
  arrange(head_group,class,-abs(logFC)) %>% 
  mutate(.r = row_number()) %>% # Add a row number variable %>%
  mutate(Response = forcats::fct_reorder(Response,.r))  %>%
  mutate(class= forcats::fct_reorder(class,.r)) 


to_plot <- rbind(to_plot1,to_plot2) %>% ungroup %>%
  mutate(head_group = factor(head_group, levels=c("phospholipid","glycolipid","neutral")) )%>%
  mutate(predictor = factor(predictor, levels=of_interest))%>%
  group_by(head_group,predictor,Response) %>%  
  filter(n()>1) %>% as.data.frame() %>%
  mutate(class = forcats::fct_reorder(class,-abs(logFC)))  %>%
  arrange(head_group,class,-abs(logFC)) %>% 
  mutate(.r = row_number()) %>% # Add a row number variable %>%
  mutate(Response = forcats::fct_reorder(Response,.r))  %>%
  mutate(class= forcats::fct_reorder(class,.r)) 

  
  
  
  
levels(to_plot$head_group)
pd = position_dodge(1)

quartz()
 to_plot %>%
  ggplot( aes(x     = logFC,
              y     = Response,
              group = leaf_node,
              color = head_group,
              fill  = head_group,
              shape = leaf_group)) +
  ggtitle("Differential lipid analysis per leaf group \n multiple hypothesis corrected lfsr < 0.05 ") +
  xlab("Effect (log2 Fold Change)") +
    geom_vline(xintercept = 0,lty =2)+
    geom_point(position= pd,
               size   = 2)+
    geom_errorbar(aes(xmin  =  upper,
                      xmax  =  lower),
                  position= pd,
                  width =  0.2,
                  linewidth  =  0.7) +
   guides( color = guide_legend(reverse = TRUE, title =""),
           shape=guide_legend(override.aes=list(fill="black",col="black")),
   #         shape = guide_legend(title ="leaf group"),
             fill="none") +
    facet_grid(class~predictor, scales="free", space="free_y")  +
    scale_color_manual(values = c( "red","blue","orange"))+
    scale_fill_manual(values = c("red","blue","orange")) +
    #guides(shape=guide_legend(override.aes=list(shape=21))) +
    scale_shape_manual(values=c(25,17)) +
    #                  labels= c( "leaf 3","leaf 4")) +
    scale_y_discrete(limits=rev) +
    theme_light(base_size =15) +
    theme(legend.position ="top",
          strip.background =element_rect(fill="white"),
          strip.text =  element_text(color = "black",size=15),
          strip.text.y= element_text(color = "white",size=15),
          axis.title.y=element_blank(),
          axis.text.y = element_text(hjust = 0, face = "bold"),
          plot.caption = element_text(hjust = 0))





 
# deg_idx <- get_significant_results(results$`5`$m.Vem_ed)
# deg <- names(deg_idx)
# length(deg)
# paste(deg, collapse =  ",")
# deg_lfsr <- get_lfsr(results$`6`$m.Vem_e)[deg_idx,]
# deg_lfsr[order(deg_lfsr[,1]),]
# 
# head(get_pm(results$`6`$m.Vem_ed)[deg_idx,])
# 
# quartz()
# barplot(get_estimated_pi(results$`5`$m.Vem_ed),las = 2)
# 
# print(get_pairwise_sharing(results$`5`$m.Vem_ed)) 
# 
# 
# deg_idx <- get_significant_results(results$`6`$m.Vem_ed)
# deg <- names(deg_idx)
# quartz()
# mash_plot_meta(results$`6`$m.Vem_ed,deg)
# 
# 
# dd <- deg[deg %in% deg_1]
# 
# paste(dd, collapse =  ",")
# deg_lfsr <- get_lfsr(results$`6`$m.Vem_e)[deg_idx,]
# deg_lfsr[order(deg_lfsr[,2]),]
# 
# results$`6`$m.Vem_ed$result$PosteriorSD
# barplot(get_estimated_pi(results$`6`$m.Vem_ed),las = 2)
# print(get_pairwise_sharing(results$`6`$m.Vem_ed))
sampleInfo$rowid %>% sort %>% unique()



