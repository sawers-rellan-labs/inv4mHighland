library(dplyr)
library(ggplot2)

pal <-  c("gold","#4a0f82")



tail(psu)

plant_csv <- "/Users/fvrodriguez/Desktop/Desktop/22_NCS_PSU_LANGEBIO_FIELDS_PSU_P_field.csv"
# csv <- "22_NCS_PSU_LANGEBIO_FIELDS - PSU_P_field.csv"
ear_csv <- "/Users/fvrodriguez/Desktop/Desktop/22_NCS_PSU_LANGEBIO_FIELDS_PSU_P_field_ear_pheno.csv"
psu_ear <- read.csv(ear_csv, na.strings = c("","n/a","NA"), skip = 1) %>%
  dplyr::select(-description,-RK, -CC, -NIR, ear_rep = rep) %>%
  dplyr::rename(rowid = "row") %>%
  dplyr::arrange(rowid) %>%
  dplyr::group_by(rowid) %>%
  dplyr::select(-ear_rep) %>%
  dplyr::summarise_all(mean, na.rm = TRUE)  %>%
  droplevels()

ear_cols <-colnames(psu_ear)[-1]


psu <- read.csv(plant_csv) %>%
  dplyr::rename(Genotype = "Who.What", rowid = "P22." ) %>%
  dplyr::filter (rowid >= 3004, rowid <= 4192) %>%
  dplyr::mutate(Treatment = factor(Treatment, levels = c("HighP","LowP")))  %>%
  dplyr::filter(Genotype %in% c("CTRL","INV4M"))  %>%
  droplevels()


psu$Genotype <- factor(psu$Genotype, levels=c("CTRL","INV4M"))


# psu %>%  filter(Genotype == "NUE")
levels(psu$Treatment) <- c("+P","-P")
# levels(psu$Genotype) <- c("CTRL","Inv4m")
psu$Rep <-as.factor(psu$Rep)

pheno <- psu %>%
  dplyr::select(rowid, Rep,Plot_Row, Plot_Column,Treatment, Genotype, Height_Anthesis:DTS, X40_DAP_dw:harvest_dw) %>%
  inner_join(psu_ear) %>%
  dplyr::mutate(HI = TKW/harvest_dw)
colnames(pheno)






m <- pheno %>% dplyr::select(Height_Anthesis:HI)


get_mm_effects <- function(x){
  per_treatmemnt <- split(x,x$Treatment)
  vars <- colnames(x %>% dplyr::select(Height_Anthesis:HI))
  trt_names <- names(per_treatmemnt)
  effects <- data.frame(predictor=1:32)
  SE <-  data.frame(predictor=1:32)
  for(trt in trt_names){
    y <- per_treatmemnt[[trt]]
    per_trt <- lapply(vars, FUN=function(v){
    form <- paste( v, "~ (1|Rep) + Genotype") 
    model <- lmer(as.formula(form), data = y)
    out <- coef(summary(model)) %>% as.data.frame() %>% 
      tibble::rownames_to_column("response")
    out$predictor <- v
    out %>% dplyr::select(predictor,everything())
    }) %>% dplyr::bind_rows()
    colnames(per_trt)[4] <- "Std.Error"
    colnames(per_trt)[7] <- "p.value"
    effects$predictor <- per_trt$predictor 
    effects$response <- per_trt$response
    SE$predictor <- per_trt$predictor 
    SE$response <- per_trt$response
    effects <- cbind(effects,data.frame(trt=per_trt$Estimate))
    SE <- cbind(SE,data.frame(trt=per_trt$Std.Error))
  }
  colnames(effects)[3:4] <- trt_names
  colnames(SE)[3:4] <- trt_names
  list(effects = effects, SE=SE)
}

result <- get_mm_effects(pheno)

result$effects[,3:4]
result$SE[,3:4]

library(mashr)

# Step 2: Obtain initial data-driven covariance matrices

data = mash_set_data(as.matrix(result$effects[,3:4]), as.matrix(result$SE[,3:4]))
U.c = cov_canonical(data)

# Step 2: Obtain initial data-driven covariance matrices

# select strong signals
m.1by1 = mash_1by1(data)
strong = get_significant_results(m.1by1,0.05)
if(length(strong)<20) {
  strong = order(apply(m.1by1$result$lfsr,1,min))[1:min(20,nrow(effects))]
}

result$effects[c(24,4,20,2),]

quartz()
mash_plot_meta(m.1by1,c(4,2))

# Perform PCA on data and return list of candidate covariance matrices
U.pca = cov_pca(data,npc=2,subset=strong)
# npc:	the number of PCs to use, should be less than or equal to n_conditions(data)
# subset: indices of the subset of data to use (set to NULL for all data)
# print(names(U.pca))

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


V.em_c_ed = mash_estimate_corr_em(data, c(U.c,U.ed), details = TRUE)
m.Vem_c_ed = V.em_c_ed$mash.model
m.Vem_c_ed$result$NAs = is.na(effects)
m.Vem_c_ed$V = V.em_c_ed$V

print(get_loglik(m.Vem_c_ed),digits=10) 
length(m.Vem_c_ed$fitted_g$pi)
aic(m.Vem_c_ed)


V.em_ed = mash_estimate_corr_em(data, U.ed, details = TRUE)
m.Vem_ed = V.em_ed$mash.model
m.Vem_ed$result$NAs = is.na(effects)

get_significant_results(m.Vem_ed)

quartz()
mash_plot_meta(m.Vem_ed,24)
