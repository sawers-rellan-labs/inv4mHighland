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
  dplyr::filter (rowid >= 3004, rowid <= 4192) %>% # Inv4m experiment rows
  dplyr::mutate(Treatment = factor(Treatment, levels = c("HighP","LowP")))  %>%
  dplyr::filter(Genotype %in% c("CTRL","INV4M")) %>%
  dplyr::mutate(Genotype = factor(Genotype, levels = c("CTRL","INV4M")))  %>%
  droplevels()
psu$Rep <-as.factor(psu$Rep)
levels(psu$Genotype) <- c("CTRL","Inv4m")
levels(psu$Treatment) <- c("+P","-P")

short_names<-c("PH","STW40", "STW50","STW60","STWHV")
long_names <-c("Height_Anthesis","X40_DAP_dw","X50_DAP_dw","X60_DAP_dw","harvest_dw")
colnames(psu)[colnames(psu) %in% long_names] <- short_names



pheno <- psu %>%
  dplyr::select(rowid,Plot, Rep,Plot_Row, Plot_Column,Treatment, Genotype, PH:DTS, STW40:STWHV) %>%
  left_join(psu_ear) %>%
  dplyr::mutate(HI = TKW/STWHV)   %>%
  dplyr::mutate(Plot_Column =case_when(
    Plot=="PlotVIII" ~ Plot_Column+10,
    Plot=="PlotVI" ~ Plot_Column,
    .default = NA
  ) )%>% filter(!is.na(PH)) %>%
  droplevels() %>%
  left_join(
    raw %>%
      dplyr::select(rowid,tissue,reliable) %>%
      tidyr::pivot_wider(
        names_from = tissue,
        values_from =reliable
      )
  )

short_names<-c("PH","STW40", "STW50","STW60","STWHV")

long_names <-c("Height_Anthesis","X40_DAP_dw","X50_DAP_dw","X60_DAP_dw","harvest_dw")

colnames(pheno)[colnames(pheno) %in% long_names] <- short_names


# Test for batch effect of the rack

ionome <- pheno %>%
  dplyr::select(rowid, Rep,Plot_Row, Plot_Column,Treatment, Genotype) %>%
  inner_join(raw %>% dplyr::select(rowid,tissue,rack, reliable)) %>%
  droplevels()


numeric_predictors <- ionome  %>%
  dplyr::select(Plot_Column, Plot_Row, Rep, rack, tissue, Genotype, Treatment) %>%
  model.matrix(~ Plot_Column+Plot_Row+as.numeric(Rep)+rack + tissue + Genotype + Treatment, .)
cor(numeric_predictors)
with(ionome,
     table(tissue,rack)
)


# Pivot wider, split columns by tissue    

by_tissue <- pheno %>%
  # dplyr::select(rowid, Rep,Plot_Row, Plot_Column,Treatment, Genotype) %>%
  left_join(
    raw %>%
    dplyr::select(rowid,tissue,reliable) %>%
      tidyr::pivot_wider(
        names_from = tissue,
        values_from =all_of(reliable)
        )
  )



reliable_by_tissue <- gsub(".","_",interaction(reliable,c("seed","stover")) %>% levels(), fixed="TRUE")

# no_outliers <- apply(by_tissue[,reliable_by_tissue],2,FUN=function(x){
#   Q1 <- quantile(x, .25, na.rm =TRUE)
#   Q3 <- quantile(x, .75, na.rm =TRUE)
#   IQR <- IQR(x, na.rm = TRUE)
#   x[x < Q1 -3*IQR | x>Q3+3*IQR] <- NA
#   x
# })
# 
# 
# 
# ionome<- cbind(
#   by_tissue %>%
#     dplyr::select(rowid, Rep,Plot_Row, Plot_Column,Treatment, Genotype),
#     no_outliers)

ionome <- by_tissue

ionome$Rep<- as.factor(ionome$Rep)
ionome$Treatment<- as.factor(ionome$Treatment)
ionome$Genotype<- as.factor(ionome$Genotype)


for(ion in reliable){
ionome[,paste0(ion,"_ratio")] <-ionome[,paste0(ion,"_seed")]/ionome[,paste0(ion,"_stover")]
}

with(ionome,
     table(Plot_Row,Plot_Column)
)

vars <- c(reliable_by_tissue,paste0(reliable,"_ratio"))


# use spatial correlation
# the model converges and is nicely specified (I have enough data for this :)

library(nlme)

effects <- lapply(vars, FUN=function(x){
  form <- paste( x, "~ Treatment*Genotype") 
  model <- lme( as.formula(form),
                random = ~ 1 | Rep,
                corr = corSpher(form = ~ Plot_Row + Plot_Column),
                na.action = na.exclude, 
                data = ionome) 
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
         !grepl("Intercept",predictor)) %>% tibble()


significant  <- effects %>%
  filter(p.adjust< 0.05, !grepl("Intercept",predictor))

unscaled_significant <- significant

significant$predictor <-factor(significant$predictor, levels=c("Treatment-P", "GenotypeINV4M","GenotypeInv4m:Treatment-P"))

levels(significant$predictor) <- c("-P","Inv4m","Inv4m:-P")


#--- Plot effects

# get data toplot
scaled <- cbind(ionome%>% dplyr::select(rowid:Genotype),
                ionome %>% dplyr::select(all_of(vars)) %>% scale())
effects <- NULL
effects <- lapply(vars, FUN=function(x){
  form <- paste( x, "~ Genotype*Treatment") # because centering changes significance of  interactions, I removed them
  model <- lme( as.formula(form),
                random = ~ 1 | Rep,
                corr = corSpher(form = ~ Plot_Row + Plot_Column),
                na.action = na.exclude, 
                data = scaled) 
  out <- coef(summary(model)) %>% as.data.frame() %>% 
    tibble::rownames_to_column("predictor")
  colnames(out)[3] <- "Std.Error"
  colnames(out)[6] <- "p.value"
  out$response <- x
  out %>% dplyr::select(response,everything())
}) %>% dplyr::bind_rows()  %>%
  mutate(p.adjust = p.adjust(p.value, method="fdr")) %>% 
  arrange(p.adjust) %>%
  inner_join(unscaled_significant %>% dplyr::select(response,predictor)) # getting significant terms in unscaled model


effects$CI.upper <- effects$Value  + qnorm(0.975)*effects$Std.Error
effects$CI.lower <- effects$Value  - qnorm(0.975)*effects$Std.Error


significant  <- effects 

significant
levels(factor(significant$predictor))
significant$predictor <-factor(significant$predictor, levels=c("Treatment-P", "GenotypeInv4m", "GenotypeInv4m:Treatment-P"))

levels(significant$predictor) <- c("-P","Inv4m", "Inv4m:-P")

levels(factor(significant$predictor))

significant$trait_group <- if_else(grepl("_seed|_stover",significant$response, perl =TRUE),"Ionome","Direct")

# get data toplot
to_plot <- significant %>%
  mutate(response= gsub("_"," ", response)) %>%
  ungroup() %>%
  group_by(response) %>%
  mutate(max_effect = max(abs(Value))) %>%
  ungroup() %>%
  group_by(predictor) %>%   # As a precaution / handle in a separate .grouped_df method
  arrange(predictor) %>%   # arrange by facet variables and continuous values
  mutate(.r = row_number()) %>% # Add a row number variable %>%
  ungroup() %>%
  mutate(response = forcats::fct_reorder(response, max_effect))


pd = position_dodge(0.4)


p1 <- to_plot %>%
# filter(grepl("ratio",response)) %>%
  ggplot( aes(x     = Value,
              y     = response,
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
  ggpubr::grids()+
  theme(legend.position ="top",
        axis.title.y=element_blank(),
        axis.text.y = element_text(size=25,hjust = 0, face = "bold"),
        plot.caption = element_text(size=25, hjust = 0))



#--- Box plots  

# vars <- c(colnames(no_outliers))
vars <- c(colnames(no_outliers),paste0(reliable,"_ratio"))

to_plot <- ionome %>%
  dplyr::select(rowid, Genotype,Treatment, all_of(vars)) %>%
  tidyr::pivot_longer(cols = vars, names_to = "mineral_tissue", values_to = "ppm") %>%
  separate_wider_delim(cols = "mineral_tissue", names=c("mineral", "tissue"),delim = "_") 

genotype.t.test <-to_plot %>%
  dplyr::mutate(Genotype_label = gsub("Inv4m","<i> Inv4m </i>", Genotype)) %>%
  dplyr::mutate(Genotype_label = factor(Genotype_label, levels=c("CTRL","<i> Inv4m </i>"))) %>%
  group_by(tissue, Treatment,mineral) %>%
  rstatix::t_test(ppm~Genotype) %>%
  rstatix::adjust_pvalue(method="fdr") %>%
  rstatix::add_significance() %>%
  rstatix::add_y_position(scales="free")

genotype.t.test %>% 
  filter(p <0.05)
effects

treatment.t.test <-to_plot %>%
  dplyr::mutate(Genotype_label = gsub("Inv4m","<i> Inv4m </i>", Genotype)) %>%
  dplyr::mutate(Genotype_label = factor(Genotype_label, levels=c("CTRL","<i> Inv4m </i>"))) %>%
  group_by(tissue, Genotype_label,mineral) %>%
  rstatix::t_test(ppm~Treatment) %>%
#  rstatix::adjust_pvalue(method="fdr") %>%
  rstatix::add_significance() %>%
  rstatix::add_y_position(scales="free")

too_low <- treatment.t.test$tissue=="seed" & treatment.t.test$mineral=="P"
treatment.t.test$y.position[too_low]  <- treatment.t.test$y.position[too_low] + 200

too_low <- treatment.t.test$tissue=="stover" & treatment.t.test$mineral=="Zn"
treatment.t.test$y.position[too_low]  <- treatment.t.test$y.position[too_low] + 3

too_low <- treatment.t.test$tissue=="stover" & treatment.t.test$mineral=="S"
treatment.t.test$y.position[too_low]  <- treatment.t.test$y.position[too_low] + 50

treatment.t.test %>% 
  filter(p <0.05)

to_plot  %>%
  filter(tissue=="stover",mineral=="P") %>% arrange(rowid) %>% print(n=64)

treatment.t.test%>%
  filter(tissue=="stover",mineral=="P") 

library(ggfx)

plot_ion_t <- function(part,ion){
to_plot %>%
    dplyr::filter(tissue==part, mineral ==ion ) %>%
    dplyr::mutate(Genotype_label = gsub("Inv4m","<i> Inv4m </i>", Genotype)) %>%
    dplyr::mutate(Genotype_label = factor(Genotype_label, levels=c("CTRL","<i> Inv4m </i>"))) %>%
    droplevels() %>%
  ggplot(aes(x=Treatment, y=ppm, color = Treatment)) +
  ggtitle(paste(ion,part)) +
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
  ggpubr::stat_pvalue_manual(
    treatment.t.test%>% filter(tissue==part, mineral==ion), # label = "p.adj",
    size = 10,
    bracket.size = 0.8,
    hide.ns = TRUE) +
  scale_y_continuous(expand = expand_scale(mult = c(0.05, 0.075)))+
  facet_wrap(.~Genotype_label)+
  scale_fill_manual(values=pal)+
  scale_color_manual(values=pal)+
  pheno_theme2 
}


plot_ion_g <- function(part,ion){
  to_plot %>%
    dplyr::filter(tissue==part, mineral ==ion ) %>%
    dplyr::mutate(Genotype_label = gsub("Inv4m","<i> Inv4m </i>", Genotype)) %>%
    dplyr::mutate(Genotype_label = factor(Genotype_label, levels=c("CTRL","<i> Inv4m </i>"))) %>%
    droplevels() %>%
    ggplot(aes(x=Genotype_label, y=ppm, color = Genotype)) +
    ggtitle(paste(ion,part)) +
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
    ggpubr::stat_pvalue_manual(
      genotype.t.test%>% filter(tissue==part, mineral==ion), # label = "p.adj",
      size = 10,
      bracket.size = 0.8,
      hide.ns = TRUE) +
    scale_y_continuous(expand = expand_scale(mult = c(0.05, 0.075)))+
    facet_wrap(.~Treatment)+
    scale_fill_manual(values=pal)+
    scale_color_manual(values=pal)+
    pheno_theme2 +
    theme(axis.text.x = element_markdown(face="bold", color= "black"))
}



stover_P <- plot_ion_t("stover","P") + ylab("mg/kg dry weight") +ylim(500,4000)
significant
quartz()
plot_ion_t("seed","Ca")

quartz()
plot_ion_t("seed","Mg")

significant

p_P <- ggpubr::ggarrange(
  plot_ion_t("seed","P")+ ylab("mg/kg dry weight") + ylim(500,4000),
  plot_ion_t("stover","P") + ylab("mg/kg dry weight") +ylim(500,4000),
  # plot_ion_t("ratio","P")+ ylab(""),
  ncol=3
)

plot_ion_t("ratio","P")+ ylab("")

p2 <- ggpubr::ggarrange(
plot_ion_t("ratio","P")+ ylab(""),
# plot_ion_t("stover","P") + ylab("mg/kg dry weight") +ylim(500,4000),
# plot_ion_t("seed","P")+ ylab("mg/kg dry weight") + ylim(500,4000),
plot_ion_t("stover","Zn")+ ylab("ppm [mg/kg]") ,
plot_ion_t("ratio","Zn")+ ylab("") ,
plot_ion_t("stover","S")+ ylab("ppm [mg/kg]"),
plot_ion_t("seed","Ca")+ ylab("ppm [mg/kg]"),
plot_ion_t("seed","Mg")+ ylab("ppm [mg/kg]"),
align="hv",
nrow=2,
ncol=3)

quartz(height = 12,width=12)
ggpubr::ggarrange(p1,p2,ncol=1, heights = c(0.5,1))



