# Check Rack effect 

library(dplyr)
library(janitor)
library(ggplot2)
library(FactoMineR)
library(factoextra)

pal <-  c("gold","#4a0f82")


mineral_names <-  c("Al", "B",  "Ca", "Fe",
                    "K",  "Mg", "Mn", "Mo",
                    "Ni", "P",  "S",  "Zn",
                    "Na", "Cu")


reliable <- c( "Ca", "Fe",  "K",
               "Mg", "Mn", "P", 
               "S",  "Zn")

unrelilable <- c("Al", "B", "Mo","Ni","Na","Cu")



mineral_csv <- "data/PSU_inv4m_ionome_all.csv"
raw <- read.csv(mineral_csv,na.strings = c("","#N/A","NA")) 
raw$tissue[raw$tissue=="stalk"] <- "stover"



plant_csv <- "../../data/22_NCS_PSU_LANGEBIO_FIELDS_PSU_P_field.csv"
# csv <- "22_NCS_PSU_LANGEBIO_FIELDS - PSU_P_field.csv"
ear_csv <- "../../data/22_NCS_PSU_LANGEBIO_FIELDS_PSU_P_field_ear_pheno.csv"
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
  dplyr::filter(Genotype %in% c("CTRL","INV4M"))



psu$Genotype <- factor(psu$Genotype, levels=c("CTRL","INV4M"))
levels(psu$Treatment) <- c("CTRL","Inv4m")

# psu %>%  filter(Genotype == "NUE")
levels(psu$Treatment) <- c("+P","-P")
# levels(psu$Genotype) <- c("CTRL","Inv4m")
psu$Rep <-as.factor(psu$Rep)

pheno <- psu %>%
  dplyr::select(rowid,Plot, Rep,Plot_Row, Plot_Column,Treatment, Genotype, Height_Anthesis:DTS, X40_DAP_dw:harvest_dw) %>%
  inner_join(psu_ear) %>%
  dplyr::mutate(HI = TKW/harvest_dw)   %>%
  dplyr::mutate(Plot_Column =case_when(
    Plot=="PlotVIII" ~ Plot_Column+10,
    Plot=="PlotVI" ~ Plot_Column,
    .default = NA
  ) )




short_names<-c("PH","STW40", "STW50","STW60","STWHV")

long_names <-c("Height_Anthesis","X40_DAP_dw","X50_DAP_dw","X60_DAP_dw","harvest_dw")

colnames(pheno)[colnames(pheno) %in% long_names] <- short_names

# Test for batch effect of the rack

ionome <- pheno %>%
  dplyr::select(rowid, Rep,Plot_Row, Plot_Column,Treatment, Genotype) %>%
  inner_join(raw %>% select(rowid,tissue,rack, reliable)) %>%
  droplevels()

# rack effect is confounded with tissue.
# stover were measured in racks 1-4 seed in racks 5-7

# futehermore tratments migh be nested within genotypes
# i.e. some racks 1-3 might be +P and 4 -P

with(ionome, table(Genotype,rack))
with(ionome, table(tissue,rack))



no_outliers <- apply(ionome[,reliable],2,FUN=function(x){
  Q1 <- quantile(x, .25)
  Q3 <- quantile(x, .75)
  IQR <- IQR(x)
  x[x < Q1 -3*IQR | x>Q3+3*IQR] <- NA
  x
})

ionome <-cbind( ionome %>% dplyr::select(rowid,tissue, rack, Rep,Plot_Row, Plot_Column,Treatment, Genotype),
                no_outliers)

ionome$Rep <- factor(ionome$Rep) 
ionome$rack <- factor(ionome$rack) 
ionome$tissue <- factor(ionome$tissue)

tissue <-NULL
numeric_predictors <- ionome  %>%
  dplyr::select(Plot_Column, Plot_Row, Rep, rack, tissue, Genotype, Treatment) %>%
  model.matrix(~ Plot_Column+Plot_Row+Rep+rack + tissue + Genotype + Treatment, .)

cor(numeric_predictors)

with(ionome,
     table(tissue,rack)
)

library(lme4)
library(lmerTest)

# nlme does not converge in the interaction model because singularities (not enough data)
# so I'll first test the interaction model with lme4 (boundary singularities are warnings, not errors)
# I can't specify the spatial  covariance structure in lme4
# So I use Plot_Row + Plot_Column as fixed effects
# If I run it as random effects  the model doesn't converge (not enough data)
# see https://idahoagstats.github.io/guide-to-field-trial-spatial-analysis/rcbd-r.html

quartz()

effects <- lapply(vars, FUN=function(x){
  form <- paste( x, "~ Plot_Row + Plot_Column + (1|Rep)  + rack +  tissue + Genotype*Treatment") 
  model <- lmer(as.formula(form), data = ionome) 
  out <- coef(summary(model)) %>% as.data.frame() %>% 
    tibble::rownames_to_column("predictor")
  colnames(out)[3] <- "Std.Error"
  colnames(out)[6] <- "p.value"
  out$response <- x
  out %>% dplyr::select(response,everything())
}) %>% dplyr::bind_rows()  %>%
  mutate(p.adjust = p.adjust(p.value, method="fdr")) %>% 
  arrange(p.adjust)

effects %>%
  filter(p.adjust < 0.05, 
         !grepl("Intercept|Row|Col",predictor)) 

# so everythung is rack 5-7 effect except for phosphorus, Mg and GenotypeINV4M:Treatment-P

# i'll ignore the rack effect




# lm(data= ionome %>%
#   filter(Genotype=="NUE", rack %in% 1:4),
#   P ~ Treatment+as.factor(rack)) %>% summary()
# 
# with(ionome %>%
#   filter(Genotype=="B73", rack %in% 5:7),
#   table(Treatment,rack))
# 
# lm(data= ionome %>%
#      filter(Genotype=="B73", rack %in% 1:4),
#    P ~ Treatment+as.factor(rack)) %>% summary()
# 
# lm(data= ionome %>%
#      filter(Genotype=="B73", rack %in% 5:7),
#    P ~ Treatment+as.factor(rack)) %>% summary()



