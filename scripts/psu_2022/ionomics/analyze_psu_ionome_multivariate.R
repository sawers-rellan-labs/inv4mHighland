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

# Read data ----

mineral_csv <- "../../data/PSU_inv4m_ionome_all.csv"
raw <- read.csv(mineral_csv,na.strings = c("","#N/A","NA")) 
raw$tissue[raw$tissue=="stalk"] <- "stover"

plant_csv <- "../../data/22_NCS_PSU_LANGEBIO_FIELDS_PSU_P_field.csv"

ear_csv <- "../../data/22_NCS_PSU_LANGEBIO_FIELDS_PSU_P_field_ear_pheno.csv"
psu_ear <- read.csv(ear_csv, na.strings = c("","n/a","NA"), skip = 1) %>%
  dplyr::select(-description,-RK, -CC, -NIR, ear_rep = rep) %>%
  #  dplyr::filter(EL > 13) %>%
  #  dplyr::filter(CD > 21) %>%
  dplyr::rename(rowid = "row") %>%
  dplyr::arrange(rowid) %>%
  dplyr::group_by(rowid) %>%
  dplyr::select(-ear_rep) %>%
  dplyr::summarise_all(mean, na.rm = TRUE)  %>%
  droplevels()

ear_cols <-colnames(psu_ear)[-1]


psu <- read.csv(plant_csv) %>%
  dplyr::rename(Genotype = "Who.What", rowid = "P22." ) %>%
  mutate(Treatment = factor(Treatment, levels = c("HighP","LowP")))  %>%
  dplyr::filter(Genotype %in% c("CTRL","INV4M"))  %>%
  droplevels()

levels(psu$Treatment) <- c("+P","-P")
colnames(psu)

pheno <- psu %>%
  dplyr::select(rowid, Rep,Plot_Row, Plot_Column,Treatment, Genotype, Height_Anthesis:DTS, X40_DAP_dw:harvest_dw) %>%
  inner_join(psu_ear) %>%
  dplyr::mutate(HI = TKW/harvest_dw)
colnames(pheno)

# Check batch effect for rack ----

ionome <- psu %>%
  dplyr::select(rowid, Rep,Plot_Row, Plot_Column,Treatment, Genotype) %>%
  inner_join( raw %>% dplyr::select(rowid,rack,tissue,all_of(reliable)))

nrow(ionome)    
no_outliers <- apply(ionome[,reliable],2,FUN=function(x){
  Q1 <- quantile(x, .25)
  Q3 <- quantile(x, .75)
  IQR <- IQR(x)
  x[x < Q1 -3*IQR | x>Q3+3*IQR] <- NA
  x
})
nrow(no_outliers)

ionome$Rep<- as.factor(ionome$Rep)
ionome$rack<- as.factor(ionome$rack)

ionome <- cbind( ionome %>%
  dplyr::select(rowid, Rep,Plot_Row, Plot_Column,Treatment, Genotype, tissue,rack),
  no_outliers)

  



# Multivariate Multiple regression

Y<- ionome  %>% dplyr::select(all_of(reliable)) 
scaled_Y <- Y %>% as.matrix() %>% scale()

# checking over different models there is no rack effect and no Rep effect
# so i just used plot row and column as covariates

Covariates <- ionome %>%
  dplyr::select(Plot_Column, Plot_Row, Rep, rack) %>%
  model.matrix(~  Plot_Column + Plot_Row + Rep + rack, .)

Genotype  <- matrix(as.factor(ionome$Genotype) %>% as.numeric())

Treatment  <- matrix(ionome$Treatment %>% as.numeric())

tissue  <- matrix(factor(ionome$tissue) %>% as.numeric())

mlm1 <- lm(scaled_Y ~ Covariates + tissue+Treatment*Genotype )
mlm2 <- lm(scaled_Y ~ Covariates + tissue*Treatment*Genotype )
mlm3 <- lm(scaled_Y ~ tissue*Treatment*Genotype)
mlm4 <- lm(scaled_Y ~ Covariates + tissue*Treatment )
mlm5 <- lm(scaled_Y ~ tissue*Treatment)

anova(mlm1,mlm2,mlm3,mlm4,mlm5)


result <- mlm3 %>% summary()
mlm3 %>% summary()
mlm5 %>% summary()

# No rack effect!


# Pivot wider -----

by_tissue <- psu %>%
  dplyr::select(rowid, Rep,Plot_Row, Plot_Column,Treatment, Genotype) %>%
  inner_join(
    raw %>%
      select(rowid,tissue,reliable) %>%
      tidyr::pivot_wider(
        names_from = tissue,
        values_from =all_of(reliable)
      )
  )

by_tissue$P_ratio <- by_tissue$P_seed/by_tissue$P_stover

no_outliers <- apply(by_tissue[,-(1:6)],2,FUN=function(x){
  Q1 <- quantile(x, .25, na.rm =TRUE)
  Q3 <- quantile(x, .75, na.rm =TRUE)
  IQR <- IQR(x, na.rm = TRUE)
  x[x < Q1 -3*IQR | x>Q3+3*IQR] <- NA
  x
})

by_tissue <- cbind( by_tissue %>%
  dplyr::select(rowid, Rep,Plot_Row, Plot_Column,Treatment, Genotype),
  no_outliers)

vars <- colnames(no_outliers )


Y<- by_tissue  %>% dplyr::select(all_of(vars)) 
scaled_Y <- Y %>% as.matrix() %>% scale()


Covariates <- by_tissue  %>%
  dplyr::select(Plot_Column, Plot_Row, Rep) %>%
  model.matrix(~  Plot_Column + Plot_Row + Rep, .)

Genotype  <- matrix(as.factor(by_tissue$Genotype) %>% as.numeric())

Treatment  <- matrix(by_tissue$Treatment %>% as.numeric())

mlm1 <- lm(scaled_Y ~ Covariates + Treatment*Genotype)

result <- mlm1 %>% summary()



effects <- lapply(names(result), FUN=function(response){
  r <- result[[response]]
  trait <- gsub(".* ", "", response, perl =TRUE)
  coeff <- r$coefficients
  data.frame(Response = trait, as.data.frame(coeff)) %>%
    tibble::rownames_to_column("predictor")
}) %>% dplyr::bind_rows()

colnames(effects)[4] <- "Std.Error"
colnames(effects)[6] <- "p.value"



effects$p.adjust <- p.adjust(effects$p.value, method="fdr")


effects$CI.upper <- effects$Estimate  + qnorm(0.975)*effects$Std.Error
effects$CI.lower <- effects$Estimate  - qnorm(0.975)*effects$Std.Error
effects %>% arrange(p.adjust) %>% filter(!grepl("Covariates|Intercept",predictor))

effects %>%
  arrange(p.adjust) %>%
  filter(!grepl("Covariates|Intercept",predictor, perl =TRUE)) %>%
  filter(p.adjust < 0.05) %>% tidyr::tibble()

significant_ionome <- effects %>%
  filter(!grepl("Covariates|Intercept",predictor, perl =TRUE)) %>%
  filter(p.adjust < 0.05) %>% 
  mutate(group = "Ionome") 
significant_ionome$Response <- paste("Stover",significant_ionome$Response)
significant_ionome$predictor <- "-P"


# get data toplot
to_plot <- effects %>%
  filter(p.adjust < 0.05) %>% # Select significant effects
  filter(!grepl("Intercept", predictor))  %>% # Just for this 2 conditions no intercepts
  mutate(predictor = gsub("stalk","", predictor)) %>%
  mutate(predictor = gsub("Treatment","P trt", predictor)) %>%
  ungroup() %>%
  group_by(Response) %>%
  mutate(max_effect = max(abs(Estimate))) %>%
  ungroup() %>%
  group_by(predictor) %>%   # As a precaution / handle in a separate .grouped_df method
  arrange(predictor) %>%   # arrange by facet variables and continuous values
  mutate(.r = row_number()) %>% # Add a row number variable %>%
  ungroup() %>%
  mutate(Response = forcats::fct_reorder(Response, max_effect))


pd = position_dodge(0.4)

plot_fx <- to_plot %>%
  ggplot( aes(x     = Estimate,
              y     = Response,
              shape  = predictor)) +
  xlab("Standardized Effect") +
  ylab("Trait") +
  geom_vline(xintercept = 0,lty =2)+
  geom_point(position= pd,
             size   = 4) +
  geom_errorbar(aes(xmin  =  CI.upper,
                    xmax  =  CI.lower),
                position= pd,
                width =  0.2,
                size  =  0.7) +
  guides(shape = guide_legend(
    title = NULL)) +
  scale_shape_manual(values= c(19,18,13)) +
  ggpubr::theme_classic2(base_size = 20) +
  theme(legend.position ="top",
        axis.title.y=element_blank(),
        axis.text.y = element_text(hjust = 0, face = "bold"),
        plot.caption = element_text(hjust = 0))


to_plot <- by_tissue %>%
  dplyr::select(rowid, Genotype,Treatment, all_of(vars)) %>%
  tidyr::pivot_longer(cols = c(vars,"P_ratio"), names_to = "trait", values_to = "ppm")

treatment.t.test <-to_plot %>%
  group_by(Genotype,trait) %>%
  rstatix::t_test(ppm~Treatment) %>%
  rstatix::adjust_pvalue(method="fdr") %>%
  rstatix::add_significance() %>%
  rstatix::add_y_position(scales="free") %>%
  arrange(p.adj)


treatment.t.test %>% 
  filter(p.adj <0.05)

ion <-NULL

library(ggfx)

t_test <-NULL

plot_ion <- function(df,t_test, ion){
  df %>%
    ggplot(aes(x=Treatment, y=ppm, color = Treatment)) +
    ggtitle(ion) +
    ggplot2::geom_boxplot(width = 0.25, linewidth=1, alpha = 0) %>% with_shadow(
      colour = "black",
      x_offset = 0,
      y_offset = 0,
      sigma = 1) +
    ggbeeswarm::geom_quasirandom(size=2) %>% with_shadow(
      colour = "black",
      x_offset = 0,
      y_offset = 0,
      sigma = 1) +
    # ggpubr::stat_pvalue_manual(
    #   t_test, # label = "p.adj",
    #   size = 10,
    #   bracket.size = 0.8,
    #   hide.ns = TRUE) +
    scale_y_continuous(expand = expand_scale(mult = c(0.05, 0.05)))+
    facet_wrap(.~Genotype)+
    ylab("mg/kg dry weight") +
    scale_fill_manual(values=pal)+
    scale_color_manual(values=pal)+
    pheno_theme2 
}


plot_ion


p01 <- plot_ion(to_plot %>% dplyr::filter(trait=="P_stover"),
         treatment.t.test %>% dplyr::filter(trait=="P_stover"),
         "Stover P")
p02 <- plot_ion(to_plot %>% dplyr::filter(trait=="P_seed"),
               treatment.t.test %>% dplyr::filter(trait=="P_seed"),
               "Seed P")

p03 <- plot_ion(to_plot %>% dplyr::filter(trait=="P_ratio"),
                treatment.t.test %>% dplyr::filter(trait=="P_ratio"),
                "Seed/Stover P ratio")

quartz(height = 14,width=14)
ggpubr::ggarrange(p01,p02,p03,
  nrow=1,
  ncol=3)

