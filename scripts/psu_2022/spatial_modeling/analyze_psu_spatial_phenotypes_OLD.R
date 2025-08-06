library(dplyr)
library(ggplot2)

pal <-  c("gold","#4a0f82")


mineral_names <-  c("Al", "B",  "Ca", "Fe",
                    "K",  "Mg", "Mn", "Mo",
                    "Ni", "P",  "S",  "Zn",
                    "Na", "Cu")


reliable <- c( "Ca", "Fe",  "K",
               "Mg", "Mn", "P", 
               "S",  "Zn")

unrelilable <- c("Al", "B", "Mo","Ni","Na","Cu")

# ionomics data
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


with(pheno,
table(Plot_Row,Plot_Column)
)

nrow(pheno)

sd(pheno$DTA[pheno$Genotype=="CTRL"])
sd(pheno$DTA[pheno$Genotype=="Inv4m"])

# vars <- colnames(pheno %>% dplyr::select(PH:HI))
vars <- colnames(pheno %>% dplyr::select(PH:Zn_seed))

# I have enough data for GxE analysis
# use spatial correlation
# the model converges and is nicely specified (I have enough data for this :)

library(nlme)

effects <- lapply(vars, FUN=function(x){
  form <- paste( x, "~ Genotype*Treatment") 
  model <- lme( as.formula(form),
                # Random effect for block
                random = ~ 1 | Rep, 
                # Residual spatial correlation
                corr = corSpher(form = ~ Plot_Row + Plot_Column), 
                na.action = na.exclude, 
                data = pheno) 
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
         !grepl("Intercept",predictor)) 


significant  <- effects %>%
  filter(p.adjust< 0.05, !grepl("Intercept",predictor)) %>% tibble()

unscaled_significant <- significant

significant$predictor <-factor(significant$predictor, levels=c("Treatment-P", "GenotypeINV4M","GenotypeInv4m:Treatment-P"))

levels(significant$predictor) <- c("-P","Inv4m","Inv4m:-P")



write.csv(effects,"~/Desktop/phenotypic_effects.csv'")

## plot effects ----
# Redo the addditive model with scaled data
scaled <- cbind(pheno %>% dplyr::select(rowid:Genotype),
      pheno %>% dplyr::select(all_of(vars)) %>% scale())
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
  mutate(trait_group= forcats::fct_reorder(factor(trait_group),p.adjust)) %>%
  dplyr::select(trait_group, everything()) %>%
  ungroup() %>%
  group_by(response) %>%
  mutate(max_effect = max(abs(Value))) %>%
  ungroup() %>%
  group_by(trait_group,predictor) %>%   # As a precaution / handle in a separate .grouped_df method
  arrange(trait_group,max_effect) %>%   # arrange by facet variables and continuous values
  ungroup() %>%
  mutate(.r = row_number()) %>% # Add a row number variable %>%
  ungroup() %>%
  mutate(response = forcats::fct_reorder(response, .r))

signif_traits <- c(
  Ca_seed = "Seed Calcium",
  Mg_seed = "Seed Magnesium",
  S_stover = "Stover Sulphur",
  Ca_stover = "Stover Calcium",
  Zn_stover = "Stover Zinc",
  P_seed = "Seed Phosphorus",
  P_stover = "Stover Phosphorus",
  KW50 = "50 Kernel Weight", 
  CD = "<b>Cob Diameter</b>",
  PH = "<b>Plant Height at Anthesis</b>",
  DTA = "<b>Days to Anthesis</b>",
  STWHV = "Stover Weight at Harvest",
  DTS = "<b>Days to Silking</b>",
#  STW40 = "<b>Stover Weight at 40 DAP</b>",
  STW40 = "Stover Weight at 40 DAP",
  STW60 = "Stover Weight at 60 DAP",
  STW50 = "Stover Weight at 50 DAP"
)

#names(signif_traits) <- levels(to_plot$response)
to_plot$trait_label <- factor(signif_traits[as.character(to_plot$response)])
to_plot$trait_label <- forcats::fct_reorder(to_plot$trait_label, to_plot$.r)



plot_table <- data.frame(
  x =0,
  y = 1:length(signif_traits),
  Trait = factor(c(names(signif_traits)),
                 levels= names(signif_traits)),
  label = signif_traits
)

# left <- plot_table %>%
#   ggplot(aes(y=y, x=x)) +
#   geom_text(aes(label = label), hjust = 0, size= 5) +
#   coord_cartesian()+
#   theme_void()
# 
# left <- to_plot %>%
#   ggplot( aes(x     = Value,
#               y     = response,
#               label  = predictor)) +
#   geom_blank()+
#   ylab("Trait") +
#   facet_wrap(.~predictor)+
#   scale_y_discrete(labels= signif_traits) +
#   ggpubr::theme_classic2(base_size = 25) +
#   theme(axis.line.y=element_blank(),
#         axis.ticks.y=element_blank(),
#         axis.text.y=element_text(hjust=0),
#         axis.title.y=element_text(face="bold"),
#         axis.line.x=element_line(color="white"),
#         axis.text.x=element_text( color = "white"),
#         axis.ticks.x=element_blank(),
#         axis.title.x=element_text( color = "white"),
#         panel.grid.minor.x=element_blank(),
#         panel.grid.major.x=element_blank(),
#         strip.background = element_blank(),
#         strip.text = element_blank()
#   )

to_plot

quartz(height = 7 , width =12)

pd = position_dodge(0.4)

to_plot %>%
  dplyr::mutate(predictor =gsub("Inv4m", "<i> Inv4m </i>", predictor)) %>%
  dplyr::mutate(predictor =gsub("<i> Inv4m </i>:", "<i>Inv4m</i> &#215; ", predictor)) %>%
  ggplot( aes(x     = Value,
              y     = trait_label,
              group = trait_group)) +
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
  facet_wrap(. ~ predictor, ncol = 3)+
  guides(shape = guide_legend(
    title = NULL, reverse= TRUE)) +
  ggpubr::theme_classic2(base_size = 25) +
  theme(strip.background = element_blank(),
        strip.text   = ggtext::element_markdown( face = "bold", size=25),
        axis.title.y = element_blank(),
        axis.text.y  = ggtext::element_markdown(hjust = 1, color = "black"),
        panel.grid.major.y = element_line(color = "grey90"),
        plot.caption = element_text(hjust = 0))


# Phenotype PCA ----

library(FactoMineR)
library(factoextra)
library(MASS)
pca1 <- PCA(Y, ncp=10, graph = FALSE)
pca1$var$coord[,1]
X <- as.matrix(pca1$call$X)
psych::cortest.bartlett(X)
psych::cortest.bartlett(cor(X), n =120)


# pca2 <- psych::principal(X, nfactors=10, rotate="oblimin", scores=TRUE)
# 
# scores <- pca1$ind$coord
# colnames(scores)
# pheno2<- cbind(pheno,pca1$ind$coord)
# plots <- list()
# 
# plots$dim12 <- with(pheno,{
#   fviz_pca_biplot(pca1,
#                   geom ="blank",
#                   col.var = "grey75",
#                   title= "Trait PCA",
#                   select.var = list(contrib = 8),
#                   labelsize = 6,
#                   repel =TRUE) +
#     xlim(-5,5) +
#     ylim(-5,5)+
#     geom_point(aes(shape = Genotype, color = Treatment),
#                size= 4) +
#     scale_color_manual(values = rev(pal)) +
#     scale_shape_manual(values=c(21,17)) +
#     coord_equal() +
#     ggpubr::theme_classic2(base_size = 20)
#   
# })
# 
# quartz()
# plots$dim12
# 
# plots$dim34 <- with(pheno,{
#   fviz_pca_biplot(pca1,axes = c(3,4), geom ="blank",
#                   col.var = "grey75",
#                   title= "Trait PCA",
#                   select.var = list(contrib = 8),
#                   labelsize = 6,
#                   repel =TRUE) +
#     xlim(-3.5,2.5) +
#     geom_point(aes(shape = Genotype, color = Treatment),
#                size= 4) +
#     scale_color_manual(values = rev(pal)) +
#     scale_shape_manual(values=c(21,17)) +
#     coord_equal() +
#     ggpubr::theme_classic2(base_size = 20)
# }
# )
# 
# quartz(height = 7, width = 14)
# ggpubr::ggarrange(plots$dim12,plots$dim34,
#                   ncol = 2,
#                   align = "hv",
#                   common.legend = TRUE)
# 
# plots$dim56 <-  with(pheno,{
#   fviz_pca_ind(pca1,axes = c(5,6), geom ="point", alpha=0) +
#     ggtitle("Morphology PCA") +
#     geom_point(aes(shape = Genotype, color = Treatment)) +
#     scale_color_manual(values = rev(pal)) +
#     scale_shape_manual(values=c(21,19,17,24))
# }
# )
# 
# 
# 
# plot_con <- function(pca, axes = c(1,2)){
#   l <- sqrt(length(axes))
#   to_grid <- lapply( axes,
#                      function(pc){
#                        fviz_contrib(pca, choice="var", axes = pc, top = ncol(pca$var$coord)) +
#                          ggplot2::ggtitle(paste0("PC",pc))
#                      }
#   )
#   ggpubr::ggarrange(plotlist = to_grid,
#                     ncol = ceiling(l),
#                     nrow = floor(l))
# }

# 
# pca1$var$contrib
# plots$con <- plot_con(pca1,1:6)
# 
# plots$cor_pc12 <- fviz_pca_var(pca1,axes = c(1,2))
# 
# plots$cor_pc13 <- fviz_pca_var(pca1,axes = c(1,3))
# 
# plots$cor_pc34 <- fviz_pca_var(pca1,axes = c(3,4))
# 
# plots$cor_pc56 <- fviz_pca_var(pca1,axes = c(5,6))

# 
# pdf(file = "~/Desktop/PSU_2022_P_field_Morphology_PCA_inv4m.pdf")
# lapply(plots,print)
# dev.off()
# 
# pheno$Genotype <- factor(pheno$Genotype)
# pheno$label <- if_else(pheno$Genotype=="INV4M","INV4M","OTHER")

# 
# summary(factor(psu$Genotype))
to_t_test <- pheno %>%
  dplyr::select(rowid, Genotype,Treatment, all_of(vars)) %>%
  tidyr::pivot_longer(cols = vars, names_to = "trait", values_to = "value") %>% 
  dplyr::group_by(Treatment) %>%
  filter(trait %in% significant$response)


pheno
# genotype.t.test$y.position <- c(189,187,76,81,77,82,30,30)
genotype_significant <- significant$response[significant$predictor=="Inv4m"]
#genotype_significant <- c("P_stover","P_seed")
genotype.t.test <- to_t_test %>%
  filter(trait %in% genotype_significant) %>%
  dplyr::group_by(Treatment,trait) %>%
  rstatix::t_test( value ~ Genotype, ref.group = "CTRL") %>%
  rstatix::add_significance() %>% 
  rstatix::add_x_position(x = "Genotype") %>%
  rstatix::add_y_position(scales = "free") %>%
  mutate(group2 = gsub("Inv4m","<i>Inv4m</i>",group2))


genotype.t.test

treatment_significant <- significant$response[significant$predictor=="-P"]

treatment.t.test <- to_t_test %>%
  filter(trait %in% treatment_significant) %>%
  dplyr::group_by(Genotype,trait) %>%
  rstatix::t_test( value ~ Treatment, ref.group = "+P") %>%
  rstatix::add_significance() %>% 
  rstatix::add_x_position(x = "Treatment") %>%
  rstatix::add_y_position(scales = "free") %>%
  dplyr::mutate(y.position=1.05*y.position) %>%
  mutate(Genotype = gsub("Inv4m","<i>Inv4m</i>",Genotype))

treatment_significant <- significant$response[significant$predictor=="Inv4m:-P"]

interaction.t.test <- to_t_test %>%
  filter(trait %in% treatment_significant) %>%
  dplyr::group_by(Treatment,trait) %>%
  rstatix::t_test( value ~ Genotype, ref.group = "CTRL") %>%
  rstatix::add_significance() %>% 
  rstatix::add_x_position(x = "Genotype") %>%
  rstatix::add_y_position(scales = "free")

# treatment.t.test <- rbind(
#   pheno %>%
#   dplyr::group_by(Genotype) %>%
#   rstatix::t_test( KW50 ~ Treatment, ref.group = "+P"),
#   pheno %>%
#     dplyr::group_by(Genotype) %>%
#     rstatix::t_test( STWHV ~ Treatment, ref.group = "+P"),
#   pheno %>% 
#     dplyr::group_by(Genotype) %>%
#     rstatix::t_test( P_stover ~ Treatment, ref.group = "+P")
# ) %>%
#   rstatix::adjust_pvalue(method = "BH") %>%
#   rstatix::add_significance() %>% 
#   rstatix::add_x_position(x = "Treatment") 

# treatment.t.test$y.position <- c(16,14,155,125,3800,3700)

pheno_theme <- ggpubr::theme_classic2(base_size = 20) +
  ggplot2::theme(
    plot.title = ggtext::element_markdown(hjust = 0.5, face ="bold"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(face="bold", color= "black", size=20, angle =45, hjust = 1),
    strip.background = element_blank(),
    strip.text = element_text(size=25),
    legend.position = "none")

library(ggtext)
library(ggfx)
pheno_theme2 <- ggpubr::theme_classic2(base_size = 20) +
  ggplot2::theme(
    plot.title = ggtext::element_markdown(hjust = 0.5, face ="bold"),
    axis.title.y = ggtext::element_markdown(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size=25, face="bold", color= "black"),
    strip.background = element_blank(),
    strip.text = ggtext::element_markdown(size=20),
    legend.position = "none")

 pheno_theme3 <- ggpubr::theme_classic2(base_size = 20) +
  ggplot2::theme(
    plot.title = ggtext::element_markdown(hjust = 0.5, face ="bold"),
    axis.text.x = element_text(color= "black", size=20),
    strip.background = element_blank(),
    strip.text = element_text(size=20),
    legend.position = "none")

height <- pheno %>%
  mutate(Genotype = gsub("Inv4m","<i>Inv4m</i>",Genotype)) %>%
  mutate(Genotype = factor(Genotype, levels= c("CTRL","<i>Inv4m</i>"))) %>%
  #filter(Genotype!="B73") %>%
  droplevels() %>%
  ggplot2::ggplot(aes( x = Genotype, y = PH, col = Genotype )) +
  ggplot2::ggtitle("Plant Height") +
  ggplot2::ylab("cm") +
  ggplot2::geom_boxplot(width = 0.25, linewidth=1.5,alpha=0) %>% with_shadow(
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
    genotype.t.test %>% filter(trait=="PH") %>% mutate(y.position= y.position+2), 
    size = 10,
    bracket.size = 0.8,
    hide.ns = TRUE) +
  scale_color_manual(values = pal) +
  coord_cartesian(ylim = c(145,190)) +
  ggplot2::facet_wrap(~Treatment) +
  pheno_theme  +
  theme(axis.text.x =  element_markdown())

DTA <- pheno %>%
  mutate(Genotype = gsub("Inv4m","<i>Inv4m</i>",Genotype)) %>%
  mutate(Genotype = factor(Genotype, levels= c("CTRL","<i>Inv4m</i>"))) %>%
  filter(Genotype!="B73") %>%
  droplevels() %>%
  ggplot2::ggplot(aes( x = Genotype, y = DTA, col = Genotype )) +
  ggplot2::ggtitle("Anthesis")  +
  ggplot2::ylab("Days") +
  ylim(68,83)+
  ggplot2::geom_boxplot(width = 0.25, linewidth=1,alpha=0) %>% with_shadow(
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
    genotype.t.test %>% filter(trait=="DTA") %>% mutate(y.position= y.position+1), 
    size = 10,
    bracket.size = 0.8,
    hide.ns = TRUE) +
  scale_color_manual(values = pal) +
  scale_y_continuous(breaks =  68:83)+
  coord_cartesian(ylim = c(68,83)) +
  ggplot2::facet_wrap(~Treatment ) +
  pheno_theme +
  theme(axis.text.x =  element_markdown())



DTS <- pheno %>%
  filter(Genotype!="B73") %>%
  mutate(Genotype = gsub("Inv4m","<i>Inv4m</i>",Genotype)) %>%
  mutate(Genotype = factor(Genotype, levels= c("CTRL","<i>Inv4m</i>"))) %>%
  droplevels() %>%
  ggplot2::ggplot(aes( x = Genotype, y = DTS, col = Genotype )) +
  ggplot2::ggtitle("Silking")  +
  ggplot2::ylab("Days")  +
  ylim(68,83)+
  ggplot2::geom_boxplot(width = 0.25, linewidth=1,alpha=0) %>% with_shadow(
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
    genotype.t.test %>% filter(trait=="DTS") %>% mutate(y.position = y.position+1), 
    size = 10,
    bracket.size = 0.8,
    hide.ns = TRUE) +
  scale_color_manual(values = pal) +
  scale_y_continuous(breaks =  68:83)+
  coord_cartesian(ylim = c(68,83)) +
  ggplot2::facet_wrap(~Treatment ) +
  pheno_theme  +
  theme(axis.text.x =  element_markdown())


CD <- pheno %>%
  mutate(Genotype = gsub("Inv4m","<i>Inv4m</i>",Genotype)) %>%
  mutate(Genotype = factor(Genotype, levels= c("CTRL","<i>Inv4m</i>"))) %>%
  filter(Genotype!="B73") %>%
  droplevels() %>%
  ggplot2::ggplot(aes( x = Genotype, y = CD, col = Genotype )) +
  ggplot2::ggtitle("Cob Diameter")  +
  ggplot2::ylab("cm")  +
 # ylim(68,83)+
  ggplot2::geom_boxplot(width = 0.25, linewidth=1,alpha=0) %>% with_shadow(
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
    interaction.t.test %>% filter(trait=="CD") %>% mutate(y.position= y.position+0.5), 
    size = 10,
    bracket.size = 0.8,
    hide.ns = TRUE) +
  scale_color_manual(values = pal) +
   scale_y_continuous(breaks =  21:31)+
  coord_cartesian(ylim = c(21,31)) +
  ggplot2::facet_wrap(~Treatment ) +
  pheno_theme  +
  theme(axis.text.x =  element_markdown())

P_seed <- pheno %>%
  filter(Genotype!="B73") %>%
  droplevels() %>%
  ggplot2::ggplot(aes( x = Genotype, y = P_seed, col = Genotype )) +
  ggplot2::ggtitle("P seed")  +
  ggplot2::ylab("ppm [mg/kg]")  +
  # ylim(68,83)+
  ggplot2::geom_boxplot(width = 0.25, linewidth=1,alpha=0) %>% with_shadow(
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
    genotype.t.test %>% filter(trait=="P_seed") %>% mutate(y.position= y.position+100), 
    size = 5,
    bracket.size = 0.8,
    hide.ns = FALSE,label = "p") +
  scale_color_manual(values = pal) +
  #scale_y_continuous(breaks =  1500:4000)+
  coord_cartesian(ylim = c(500,4000)) +
  ggplot2::facet_wrap(~Treatment ) +
  pheno_theme 

# 
# P_stover <- pheno %>%
#   filter(Genotype!="B73") %>%
#   droplevels() %>%
#   ggplot2::ggplot(aes( x = Genotype, y = P_stover, col = Genotype )) +
#   ggplot2::ggtitle("P stover")  +
#   ggplot2::ylab("ppm [mg/kg]")  +
#   # ylim(68,83)+
#   ggplot2::geom_boxplot(width = 0.25, linewidth=1,alpha=0) %>% with_shadow(
#     colour = "black",
#     x_offset = 0,
#     y_offset = 0,
#     sigma = 1) +
#   ggbeeswarm::geom_quasirandom(size=2) %>% with_shadow(
#     colour = "black",
#     x_offset = 0,
#     y_offset = 0,
#     sigma = 1) +
#   ggpubr::stat_pvalue_manual(
#     genotype.t.test %>% filter(trait=="P_stover") %>% mutate(y.position= y.position+100), 
#     size = 5,
#     bracket.size = 0.8,
#     tip.length = 0.015, 
#     hide.ns = FALSE,label = "p") +
#   scale_color_manual(values = pal) +
#   #scale_y_continuous(breaks =  1500:4000)+
#   coord_cartesian(ylim = c(500,4000)) +
#   ggplot2::facet_wrap(~Treatment ) +
#   pheno_theme
# 
# 
# #quartz(height=12, width=15)
# quartz()
# ggpubr::ggarrange(P_seed, P_stover,
#                    ncol=2,nrow =1,
#                    font.label = list(size=30))


library(ggfx)
KW50 <- pheno %>%
  filter(Genotype!="B73") %>%
  mutate(Genotype = gsub("Inv4m","<i>Inv4m</i>",Genotype)) %>%
  mutate(Genotype = factor(Genotype, levels= c("CTRL","<i>Inv4m</i>"))) %>%
  droplevels() %>%
  ggplot2::ggplot(aes( x = Treatment, y = KW50, col = Treatment )) +
  ggplot2::ggtitle("50 Kernel <br> Weight")  +
  ggplot2::ylab("g")  +
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
  treatment.t.test %>% filter(trait=="KW50") , 
    size = 10,
    bracket.size = 0.8,
    hide.ns = TRUE) +
  scale_color_manual(values = pal) +
  scale_y_continuous(breaks =  5:17)+
  coord_cartesian(ylim = c(5,17)) +
  ggplot2::facet_wrap(. ~ factor(Genotype, levels= c("CTRL","<i>Inv4m</i>") )) +
  pheno_theme2 +
  theme(strip.text = element_markdown())

quartz()
KW50 


STWHV <- pheno %>%
  filter(Genotype!="B73") %>%
  mutate(Genotype = gsub("Inv4m","<i>Inv4m</i>",Genotype)) %>%
  mutate(Genotype = factor(Genotype, levels= c("CTRL","<i>Inv4m</i>"))) %>%
  droplevels() %>%
  ggplot2::ggplot(aes( x = Treatment, y = STWHV, col = Treatment )) +
  ggplot2::ggtitle("Stover Dry Weight <br> at Harvest")  +
  ggplot2::ylab("g")  +
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
    treatment.t.test %>% filter(trait=="STWHV"), 
    size = 10,
    bracket.size = 0.8,
    hide.ns = TRUE) +
  scale_color_manual(values = pal) +
  # scale_y_continuous(breaks =  5:17)+
  coord_cartesian(ylim = c(50,175)) +
  ggplot2::facet_wrap(. ~ factor(Genotype, levels= c("CTRL","<i>Inv4m</i>") )) +
  pheno_theme2


quartz()
STWHV 


 P_seed <- pheno %>%
   filter(Genotype!="B73") %>%
   mutate(Genotype = gsub("Inv4m","<i>Inv4m</i>",Genotype)) %>%
   mutate(Genotype = factor(Genotype, levels= c("CTRL","<i>Inv4m</i>"))) %>%
  droplevels() %>%
  ggplot2::ggplot(aes( x = Treatment, y = P_seed , col = Treatment )) +
  ggplot2::ggtitle("Seed <br> Phosphorus")  +
  ggplot2::ylab("ppm [mg/kg]")  +
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
    treatment.t.test %>% filter(trait=="P_seed"), 
    size = 10,
    bracket.size = 0.8,
    hide.ns = TRUE) +
  scale_color_manual(values = pal) +
  # scale_y_continuous(breaks =  5:17)+
  coord_cartesian(ylim = c(500,4200)) +
   ggplot2::facet_wrap(. ~ factor(Genotype, levels= c("CTRL","<i>Inv4m</i>") )) +
   pheno_theme2


P_stover <- pheno %>%
  filter(Genotype!="B73") %>%
  mutate(Genotype = gsub("Inv4m","<i>Inv4m</i>",Genotype)) %>%
  mutate(Genotype = factor(Genotype, levels= c("CTRL","<i>Inv4m</i>"))) %>%
  droplevels() %>%
  ggplot2::ggplot(aes( x = Treatment, y = P_stover , col = Treatment )) +
  ggplot2::ggtitle("Stover <br> Phosphorus")  +
  ggplot2::ylab("ppm [mg/kg]")  +
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
    treatment.t.test %>% filter(trait=="P_stover"), 
    size = 10,
    bracket.size = 0.8,
    hide.ns = TRUE) +
  scale_color_manual(values = pal) +
  # scale_y_continuous(breaks =  5:17)+
  coord_cartesian(ylim = c(500,4200)) +
  ggplot2::facet_wrap(. ~ factor(Genotype, levels= c("CTRL","<i>Inv4m</i>") )) +
  pheno_theme2


quartz(height=12, width=15)
ggpubr::ggarrange( DTA, DTS, height,CD, KW50, STWHV, P_seed, P_stover,
                   ncol=4,nrow =2,
                   labels = LETTERS[2:9],
                   font.label = list(size=30))


quartz(height=12, width=15)
ggpubr::ggarrange( DTA, DTS, height,CD, KW50, STWHV, P_seed, P_stover,
                   ncol=4,nrow =2,
                   labels = LETTERS[2:9],
                   font.label = list(size=30))


dw_raw <- pheno %>%
  group_by(rowid ) %>%
  mutate(STW121 = STWHV) %>%
  dplyr::select(rowid =  rowid, Rep, Plot_Row,Plot_Column,Rep,
                Treatment,Genotype, starts_with("STW")) %>%
  tidyr::pivot_longer(
    cols = c("STW40","STW50","STW60","STW121"),
    names_to = "DAP_string",
    values_to = "STW") %>%
  mutate(DAP = gsub("STW","",DAP_string, perl =TRUE) %>% as.integer())  %>%
  filter(!is.na(STW)) %>%
  droplevels()


dw_data <- pheno %>%
  group_by(rowid ) %>%
  mutate(STW121 = pmax(STW40,STW50,STW60,STWHV)) %>% # replace harvest weight with maximum
  mutate(STW01 = 0) %>%
  dplyr::select(rowid =  rowid, Rep, Plot_Row,Plot_Column,Rep,
                Treatment,Genotype, starts_with("STW")) %>%
  tidyr::pivot_longer(
    cols = c("STW01","STW40","STW50","STW60","STW121"),
    names_to = "DAP_string",
    values_to = "STW") %>%
  mutate(DAP = gsub("STW","",DAP_string, perl =TRUE) %>% as.integer())  %>%
  filter(!is.na(STW)) %>%
  droplevels()


# I don't have
dw_data$DW1 <- dw_data$STW+1

#  Drawing growth curves for each replicate----


# Using growthcurver to draw curves and estimate params
# for complete data  replicates(n=12)

library(growthcurver)

library(tidyr)

morethan4 <- dw_data %>%
  ungroup() %>%
  group_by(Genotype,Treatment,rowid) %>%
  summarise(n=sum(!is.na(DW1))) %>%
  filter(n >3) %>%
  arrange(rowid) %>% pull(rowid)

for_curves <- paste0("row_", morethan4)

dw_data$group <- factor(paste(dw_data$Rep,dw_data$Genotype))

without_na <- dw_data %>%
  ungroup() %>%
  dplyr::mutate(row_id = paste0("row_",rowid)) %>%
  filter(row_id%in% for_curves) %>%
  dplyr::select(row_id,time = DAP,DW1) %>%
  tidyr::pivot_wider(names_from ="row_id", values_from = DW1 )

treatment.t.test <- dw_raw %>%
  ungroup() %>%
  group_by(Genotype) %>%
  filter(Genotype!="B73") %>%
  mutate(DAP = as.factor(DAP)) %>%
  group_by(Genotype, DAP)  %>%
  droplevels() %>%
  rstatix::t_test( STW ~ Treatment, ref.group = "+P") %>%
  rstatix::adjust_pvalue(method="fdr") %>%
  rstatix::add_significance() %>%
  rstatix::add_y_position(scales = "free_y") %>%
  rstatix::add_x_position(x = "DAP", group = NULL, dodge = 1)

treatment.t.test$y.position <- treatment.t.test$y.position+5

library(ggfx)

pd <- position_dodge(1)

p_dw <- dw_raw %>%
  ungroup() %>%
  filter(Genotype!="B73" & DAP > 1) %>%
  droplevels() %>%
  ggplot2::ggplot(aes( x = as.factor(DAP), y = STW, col = Treatment )) +
  ggplot2::ggtitle("Stover Dry Weight")  +
  ggplot2::ylab("g")  +
  ggplot2::xlab("Days After Planting")+
  ggplot2::geom_boxplot(width = 0.25, linewidth=1, alpha = 0, position = pd) %>% with_shadow(
    colour = "black",
    x_offset = 0,
    y_offset = 0,
    sigma = 1) +
  ggbeeswarm::geom_quasirandom(size=2, width =0.25, dodge.width = 1) %>% with_shadow(
    colour = "black",
    x_offset = 0,
    y_offset = 0,
    sigma = 1) +
  facet_wrap(.~Genotype)+
  ggpubr::stat_pvalue_manual(
    treatment.t.test ,
    size = 10,
    bracket.size = 0.8,
    tip.length = 0.01,
    hide.ns = TRUE) +
  scale_color_manual(values = pal) +
  scale_y_continuous(expand = expansion(mult=c(0.05, 0.1))) +
  pheno_theme3 +
  theme(legend.position = c(0.1,0.9))


quartz()
p_dw 

# Imputing NAs with row mean
# m <- as.matrix(without_na)
#
# k <- which(is.na(m), arr.ind=TRUE)
# m[k] <- rowMeans(m, na.rm=TRUE)[k[,1]]

# dry_weight <- as.data.frame(m)

dry_weight <- without_na

growth_fit <- lapply(colnames(dry_weight)[-1],
                     FUN=function(x){
                       growth <- SummarizeGrowth(dry_weight$time, dry_weight[,x])
                       t <- 40:121
                       dw <- predict(growth$model,newdata = data.frame(t=t))
                       data.frame(sample=x, time =t, fittedDW=dw)
                     }) %>%
  dplyr::bind_rows()


growth_sum <- SummarizeGrowthByPlate(dry_weight, bg_correct = "none" )


# Plot the results



sample <- psu %>%
  dplyr::select(sample= "rowid", Rep, Plot_Row, Plot_Column, Treatment, Genotype) %>%
  #dplyr::mutate(Treatment= factor(Treatment, levels =c("LowP","HighP"))) %>%
  dplyr::mutate(sample = paste0("row_",sample)) %>%
  dplyr::inner_join(
    data.frame(sample = colnames(dry_weight) )
  ) %>%
  dplyr::inner_join(growth_sum) %>%
  dplyr::rename(AUCE = auc_e, AUCL= auc_l, STWMAX = k, TMID = t_mid)



selected_cols <- c("AUCE","AUCL","STWMAX","TMID")
# ylabels <-   c("AUC<sub>emp</sub>","AUC<sub>fit</sub>","STW<sub>max</sub>","T<sub>1/2</sub>")
ylabels <-   c("AUC empirical",
               "AUC logistic",
               "STW<sub>max</sub>",
               "T<sub>1/2</sub>")

names(ylabels) <- selected_cols
yunit_labels <-   c( AUCE ="kg &times; day",
                     AUCL ="kg &times; day",
                     STWMAX ="g",
                     TMID = "Days")

to_plot <- sample  %>%
  mutate(AUCE = AUCE/1000, AUCL = AUCL/1000) %>%
  ungroup() %>%
  dplyr::select(Genotype,Treatment, all_of(selected_cols) )%>%
  tidyr::pivot_longer(
    cols=all_of(selected_cols),
    names_to = "pheno",
    values_to = "value") %>%
  mutate(pheno = factor(pheno, levels= selected_cols)) %>%
  group_by(pheno) %>%
  mutate(z_score = scale(value))

t_test <-to_plot %>%
  group_by(Genotype,pheno) %>%
  rstatix::t_test( value~Treatment, ref.group="+P",
                   var.equal = FALSE,
                   detailed = FALSE,
  )  %>%
  rstatix::adjust_pvalue(method = "fdr") %>%
  rstatix::add_significance() %>%
  rstatix::add_y_position(scales = "free_y") %>%
  rstatix::add_x_position(dodge = 1,group = "Treatment")

t_test$y.position <- 1.05*t_test$y.position


pd <- position_dodge(1)

plot_pheno <- function(x) {
  to_plot %>%
    filter(pheno==x) %>%
    ggplot2::ggplot(aes( x = Treatment, y = value, col = Treatment )) +
    ggtitle(ylabels[x]) +
    ylab(yunit_labels[x]) +
    ggplot2::geom_boxplot(width = 0.25, linewidth=1, alpha=0) %>% with_shadow(
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
      t_test %>% filter(pheno==x), label = "p.adj.signif",
      size = 10,
      bracket.size = 0.8,
      hide.ns = TRUE) +
    facet_wrap(.~Genotype) +
    scale_color_manual(values = pal) +
    scale_y_continuous(expand = expansion(mult=c(0.05, 0.1))) +
    pheno_theme2 
    # theme( plot.margin = unit(c(1,-1,1,-1), 'lines'))
}

all_plots <- lapply(selected_cols, FUN = plot_pheno)
class(all_plots)


p_gc <- ggpubr::ggarrange(plotlist = all_plots[c(1,2,4,3)],
                  nrow = 1, ncol=4)

out_plot <- ggpubr::ggarrange(p_dw,p_gc, nrow=3,ncol=1)


quartz(width =12,height =15) 
out_plot

ggsave(out_plot,
       filename = "~/Desktop/growth_curves_2.svg",
       height=10, width = 14, units = "in")
# sample$PC1 <- pca.res$x[,1]
# sample$PC2 <- pca.res$x[,2]
# 
# quartz()
# sample %>%
#   ggplot2::ggplot(aes(x=PC1,y=PC2, color  = Treatment, group = Treatment, shape=Genotype)) +
#   ggplot2::geom_point()+
#   ggplot2::scale_color_manual(values = rev(pal))


table(growth_fit$time)

gc_fitted <- sample %>%
  dplyr::select(sample:Genotype) %>%
  dplyr::right_join(growth_fit) %>%
  arrange(Treatment, sample) %>%
  mutate(plot_order = 1:n()) %>% dplyr::select(plot_order, everything()) %>%
  mutate(sample = forcats::fct_reorder(as.factor(sample), plot_order, .fun=first))

gc <- sample %>%
  dplyr::select(sample:Genotype) %>%
  dplyr::right_join(
    dry_weight %>%
      pivot_longer(cols = c(starts_with('row_')),
                   names_to = "sample", values_to = "DW1")
  )


quartz()

gc %>%
  ggplot2::ggplot(aes( x = time, y = DW1, color=Treatment, group = sample)) +
  #ggplot2::ggtitle("Shoot Dry Mass")  +
  ggplot2::ylab("Stover Dry Weight [g]")  +
  ggplot2::xlab("Days After Planting")  +
  ggplot2::xlim(c(40,130))+
  ggplot2::ylim(c(0,150))+
  ggplot2::geom_line() +
  ggplot2::scale_color_manual(values = rev(pal)) +
  ggplot2::facet_wrap(~Genotype) +
  ggpubr::theme_classic2(base_size = 25) +
  ggplot2::theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5),
    strip.background = element_blank())

pd <- position_dodge(0.5)

p_00 <-  sample %>%
  ggplot2::ggplot(aes( x = Genotype,
                       y = STWMAX,
                       fill=Genotype,
                       # shape=Genotype,
                       color= Genotype)) +
  ylab("Maximum Stover Dry Weight [g]")+
  geom_boxplot(fill="white", width=0.25,outliers = FALSE, position = pd) +
  # stat_summary(fun.y = mean,
  #              geom = "point", color="white",size= 5,
  #              position= pd) +
  # stat_summary(fun.data = mean_cl_normal,
  #              key_glyph = "rect",
  #              geom = "errorbar",
  #              size=1,
  #              width=0.125,
  #              position= pd)+
  ggbeeswarm::geom_quasirandom(size=2,
                               shape = 21,
                               key_glyph = "rect",
                               color= "white",
                               dodge.width =0.5) +
  facet_wrap(~Treatment)+
  scale_shape_manual(values=c(21,24)) +
  scale_color_manual(values=rev(pal)) +
  scale_fill_manual(values=rev(pal)) +
  guides(shape = FALSE,fill=FALSE) +
  ggpubr::theme_classic2(base_size = 25) +
  ggplot2::theme(
    plot.title = element_text(hjust = 0.5),
    strip.text = element_text(face="bold"),
    strip.background = element_blank(),
    legend.position = "none")

quartz()
p_00

library(grid)
draw_path_with_shadow <- function (data, params, size) 
{
  if (is.null(data$linetype)) {
    data$linetype <- 0
  }
  else {
    data$linetype[is.na(data$linetype)] <- 0
  }
  grob <- segmentsGrob(0.1, 0.5, 0.9, 0.5, gp = gpar(col = alpha(data$colour %||% 
                                                                   data$fill %||% "black", data$alpha), fill = alpha(params$arrow.fill %||% 
                                                                                                                       data$colour %||% data$fill %||% "black", data$alpha), 
                                                     lwd = 1.5* .pt, lty = data$linetype %||% 
                                                       1, lineend = params$lineend %||% "butt"), arrow = params$arrow)
  if (!is.null(params$arrow)) {
    angle <- deg2rad(params$arrow$angle)
    length <- convertUnit(params$arrow$length, "cm", valueOnly = TRUE)
    attr(grob, "width") <- cos(angle) * length * 1.25
    attr(grob, "height") <- sin(angle) * length * 2
  }
  grob %>%   with_shadow(
    colour = "black",
    x_offset = 0,
    y_offset = 0,
    sigma = 1)
}


p_01 <- gc_fitted %>%
  ggplot2::ggplot(aes( x = time, y = fittedDW, color=Treatment, group = sample)) +
  ggtitle("Fitted Stover Weight") +
  ggplot2::ylab("g")  +
  ggplot2::xlab("Days After Planting")  +
  ggplot2::xlim(c(40,130))+
  ggplot2::ylim(c(0,155))+
  ggplot2::geom_line(key_glyph=draw_path_with_shadow) %>% with_shadow(
    colour = "black",
    x_offset = 0,
    y_offset = 0,
    sigma = 1) +
  ggplot2::scale_color_manual(values = pal) +
  ggplot2::guides(color=guide_legend(override.aes=list(fill=NA)))+
  ggplot2::labs(color="")+
  ggplot2::facet_wrap(~Genotype) +
  ggpubr::theme_classic2(base_size = 20) +
  ggplot2::theme(
    
    legend.position = c(0.8,0.9),
    plot.title = element_text(hjust = 0.5, face="bold",size=25),
    strip.text = element_text( size=20),
    legend.text = element_text(face="bold", size=20),
    plot.background = element_rect(fill = "transparent",
                                   colour = NA_character_), # necessary to avoid drawing plot outline
    legend.background = element_rect(fill = "transparent",color=NA),
    legend.box.background = element_rect(fill = "transparent", color=NA),
    legend.key = element_rect(fill = "transparent"),
    strip.background = element_blank())

out_plot <- ggpubr::ggarrange(p_dw,p_01, p_gc, nrow=3,ncol=1)


quartz(width =12,height =15) 
out_plot



quartz(height=10, width=15)
p_biomasss <- ggpubr::ggarrange(KW50,p_dw, widths = 1:2)

p_inv4m <- ggpubr::ggarrange(DTA, DTS, height, ncol=3)

quartz(height=14, width=14)
ggpubr::ggarrange( p_inv4m, p_biomasss, p_P, nrow = 3)
                     
effects_fig <- ggpubr::ggarrange( 
  DTA, DTS, height,
  KW50,p_01,stover_P,
  ncol=3,nrow =2,
  labels = LETTERS[2:7],
  font.label = list(size=30)
  )
quartz(height=10, width=15)
effects_fig 
ggsave(effects_fig,file="~/Desktop/effects_fig_1.svg",height = 10, width=15)
ggsave(effects_fig,file="~/Desktop/effects_fig_2.svg",height = 14, width=14)
# Only the harvest_dry weight
#
# mod_HDM1  <- lm(data = sample,
#    STWmax ~ Rep + Plot_Row + Plot_Column + Treatment * Genotype)
#
# mod_HDM2 <-lm(data = sample,
#              STWmax ~ Rep + Plot_Row + Plot_Column + Treatment +Genotype)
# anova(mod_HDM2,mod_HDM1)

# mod_HDM <- mod_HDM2
# anova(mod_HDM)

# comparing marginal means and confidence intervals
# library(multcomp)
# library(emmeans)
#
# em_all <- data.frame(
#   emmeans(mod_HDM, ~ Treatment) %>% as.data.frame(),
#   Genotype = "all") %>% dplyr::select( Genotype,Treatment, everything())
#
#
# marginal_HDM <- rbind(
#   em_all,
#   emmeans(mod_HDM, ~  Genotype + Treatment) %>% as.data.frame()
# )


#
#
# pd = position_dodge(0.4)    ### How much to jitter the points on the plot
#
#
# marginal_HDM
#
# p_01 <-marginal_HDM  %>%
#   filter(Genotype != "all") %>%
#   droplevels() %>%
#   ggplot( aes(y     = Genotype,
#               x     = emmean,
#               color = Treatment)) +
#   xlab("Maximum Stover Weight [g]") +
#   geom_point(shape  = 15,
#              size   = 4,
#              position = pd) +
#
#   geom_errorbar(aes(xmin  =  lower.CL,
#                     xmax  =  upper.CL),
#                 width =  0.2,
#                 size  =  0.7,
#                 position = pd) +
#   ggplot2::scale_color_manual(values = rev(pal)) +
#   ggpubr::theme_classic2(base_size = 25) +
#   theme(legend.position ="top",
#         axis.title   = element_text(face = "bold"),
#         axis.text    = element_text(face = "bold"),
#         plot.caption = element_text(hjust = 0)) +
#   coord_flip()

quartz(width=12)
ggpubr::ggarrange(p_01,p_00, ncol=2,
                  align="hv",
                  widths=c(0.7,0.3),
                  common.legend = TRUE
)

# Plot effects ----

vars <- colnames(sample %>% dplyr::select(STWMAX:AUCE))

Y <- pheno %>%  dplyr::select( all_of(vars))


scaled_Y <- Y %>% as.matrix() %>% scale()

Y <- pheno %>%  dplyr::select( all_of(vars))

Covariates <- pheno %>%
  dplyr::select(Plot_Column, Plot_Row, Rep) %>%
  model.matrix(~  Plot_Column + Plot_Row+Rep,.)

Genotype  <- matrix(pheno$Genotype %>% as.numeric())

model <- lm(Y ~ Genotype + Treatment + Covariates)


Treatment  <- matrix(sample$Treatment %>% as.numeric())

mlm1 <- lm(scaled_Y ~ Genotype + Treatment + Covariates)

mlm4 <- lm(scaled_Y ~ Genotype*Treatment + Covariates)

anova(mlm4,mlm3)

mlm3 <- lm(as.matrix(Y) ~ Genotype + Treatment + Covariates)

result <- mlm3 %>% summary()



effects <- lapply(names(result), FUN=function(response){
  r <- result[[response]]
  trait <- gsub(".* ", "", response, perl =TRUE)
  coeff <- r$coefficients
  data.frame(Response = trait, as.data.frame(coeff)) %>%
    tibble::rownames_to_column("predictor")
}) %>% dplyr::bind_rows()

colnames(effects)[4] <- "Std.Error"
colnames(effects)[6] <- "p.value"

effects %>%
  filter(p.value < 0.05)

effects$p.adjust <- p.adjust(effects$p.value, method="fdr")


effects$CI.upper <- effects$Estimate  + qnorm(0.975)*effects$Std.Error
effects$CI.lower <- effects$Estimate  - qnorm(0.975)*effects$Std.Error

effects %>%
  filter(p.adjust< 0.05)


significant_gc  <- effects %>%
  filter(p.adjust< 0.05, !grepl("Intercept",predictor)) %>%
  arrange(-abs(Estimate))

significant_gc$predictor <- "-P"
significant$group <- "Field"
significant_gc$group <- "Growth Curves"
out <- rbind(significant, significant_gc, significant_ionome) %>% dplyr::select(group, everything())


out$group<- factor(out$group,levels=c("Ionome","Growth Curves","Field"))
out$predictor <- factor(out$predictor, levels=c("-P","Inv4m"))
#write.csv(out,"~/Desktop/phenotypic_effects.csv",row.names = FALSE)

out <-read.csv("~/Desktop/phenotypic_effects.csv")
## plot effects ----


signif_traits <- c(
  STW50 = "Stover Weight at 50 DAP",
  DTA ="Days to Anthesis",
  DTS = "Days to Silking",
  PH ="Plant Height at Anthesis",
  STW60 = "Stover Weight at 60 DAP",
  STW40 = "Stover Weight at 40 DAP",
  KW50 = "Kernel Weight 50 Seeds",
  STWHV = "Stover Weight at Harvest",
  CD = "Cob Diameter",
  AUCL = "Logistic Growth AUC",
  AUCE = "Emprical Growth AUC",
  STWMAX = "Maximum Stover Weight",
  TMID = "Time to Half Maximum STW",
  P_stover = "Stover P content",
  P_seed = "Stover Zn content",
  Zn_stover = "Stover Zn content",
  S_stover = "Stover S content"
)

# get data toplot
to_plot <- out %>%
  #filter(p_adj < 0.05) %>% # Select significant effects
  #filter(predictor %in% c("Genotype","Treatment"))  %>% # Just for this 2 conditions no intercepts
  ungroup() %>%
  group_by(Response) %>%
  mutate(max_effect = max(abs(Estimate))) %>%
  ungroup() %>%
  group_by(group) %>%   
  mutate(source = forcats::fct_reorder(group,max_effect )) %>%
  ungroup() %>%
  group_by(group) %>%
  arrange(group,desc(max_effect)) %>%# arrange by facet variables and continuous values
  mutate(.r = row_number()) %>% # Add a row number variable %>%
  mutate(Response = forcats::fct_reorder(Response,.r)) 
  


trait_order <- levels(to_plot$Response)
trait_order 

signif_traits <- signif_traits[trait_order]

plot_table <- data.frame(
  x = 0,
  y = 1:length(signif_traits),
  Trait = factor(c(names(signif_traits)),
                 levels= names(signif_traits)),
  label = signif_traits
)

# 
# left <- to_plot %>%
#   ggplot( aes(x     = Estimate,
#               y     = Response,
#               label  = predictor)) +
#   geom_blank()+
#   ylab("Trait") +
#   scale_y_discrete(labels= signif_traits) +
#   ggpubr::theme_classic2(base_size = 25) +
#   theme(axis.line.y=element_blank(),
#         axis.ticks.y=element_blank(),
#         axis.text.y=element_text(hjust=0),
#         axis.title.y=element_text(face="bold"),
#         axis.line.x=element_line(color="white"),
#         axis.text.x=element_text( color = "white"),
#         axis.ticks.x=element_blank(),
#         axis.title.x=element_text( color = "white"),
#         panel.grid.minor.x=element_blank(),
#         panel.grid.major.x=element_blank()
#   )


left <- to_plot %>%
  ggplot( aes(x     = Estimate,
              y     = Response,
              label  = predictor)) +
  geom_blank()+
  ylab("Trait") +
  facet_wrap(.~predictor)+
  scale_y_discrete(limits = rev,labels= rev(signif_traits))+
  ggpubr::theme_classic2(base_size = 25) +
  theme(axis.line.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_text(hjust=0),
        axis.title.y=element_text(face="bold"),
        axis.line.x=element_line(color="white"),
        axis.text.x=element_text( color = "white"),
        axis.ticks.x=element_blank(),
        axis.title.x=element_text( color = "white"),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.x=element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank()
  )



pd = position_dodge(0.4)


right <- to_plot %>%
  ggplot( aes(x     = Estimate,
              y     = Response)) +
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
  scale_y_discrete(limits = rev,labels= rev(signif_traits))+
  facet_wrap(.~predictor)+
  guides(shape = guide_legend(
    title = NULL, reverse= TRUE)) +
  scale_shape_manual(values=c(17,21),
                     labels= c( expression(italic("inv4m")), "+P")) +
  ggpubr::theme_classic2(base_size = 25) +
  ggpubr::grids()+
  theme(strip.background = element_blank(),
        strip.text = element_text(face="bold", size=25),
        axis.text.x = element_text( size= 20),
        axis.title = element_text(size= 20),
        axis.title.y=element_blank(),
        axis.text.y = element_text(hjust = 0),
        plot.caption = element_text(hjust = 0))

quartz(width = 10, height=6)
right

ggsave(right, file="~/Desktop/inv4m_effect_ci.svg",width = 10, height=6)

quartz(width = 12, height=5.5)
ggpubr::ggarrange(left,right, nrow = 1,align = "h",
                  common.legend = TRUE,
                  widths = c(0.35,0.65)
)


quartz(width = 10, height=6)
right


quartz(width = 10, height=6)
ggpubr::ggarrange(left,right, nrow = 1,align = "h",
                  common.legend = TRUE,
                  legend = "top",
                  widths = c(0.35,0.65)
)



#####----



quartz()
curves %>%
  ggplot2::ggplot(aes(x=PC1,y=PC2, color  = Treatment, group = Treatment, shape=Genotype)) +
  ggplot2::geom_point()+
  ggplot2::scale_color_manual(values = pal)



to_plot <- dry_weight  %>%
  dplyr::select(sample:Genotype) %>%
  inner_join(
    dry_weight %>%
      pivot_longer(cols = c(starts_with('row_')),
                   names_to = "sample", values_to = "DW1")
    
  )



quartz()
hist(to_plot$DW1)

quartz()
to_plot  %>%
  ggplot2::ggplot(aes( x = time, y = DW1, color=Treatment, group = sample))+
  #ggplot2::ggtitle("Shoot Dry Mass")  +
  ggplot2::ylab("Shoot Dry Mass (g)")  +
  ggplot2::xlab("Days After Planting")  +
  ggplot2::xlim(c(40,130))+
  ggplot2::ylim(c(0,160))+
  ggplot2::geom_line() +
  ggplot2::scale_color_manual(values = pal) +
  ggplot2::facet_wrap(~Genotype) +
  ggpubr::theme_classic2(base_size = 12) +
  ggplot2::theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5),
    strip.background = element_blank())


table(to_plot$Rep)


#######

quartz()
dw_data  %>%
  dplyr::mutate(row_id = paste0("row_",rowid))%>%
  dplyr::filter(row_id %in% curves$sample) %>%
  ggplot2::ggplot(aes( x = DAP, y = DW, col = Treatment, group = rowid)) +
  #ggplot2::ggtitle("Shoot Dry Mass")  +
  ggplot2::ylab("Shoot Dry Mass (g)")  +
  ggplot2::xlab("Days After Planting")  +
  ggplot2::xlim(c(40,130))+
  ggplot2::ylim(c(0,155))+
  ggplot2::geom_line(stat="smooth", method = "nls",
                     formula = y ~ theta1/(1 + exp(-(theta2 + theta3*x))),
                     method.args = list(start=list(theta1 = maxw, theta2 = -5.4, theta3 = 0.0689575)),
                     alpha =0.8,
                     se = FALSE) +
  ggplot2::scale_color_manual(values = rev(pal)) +
  ggplot2::scale_x_continuous(breaks = 10*(4:13)) +
  ggbeeswarm::geom_beeswarm(
    aes(fill = Treatment), pch =21,
    col = "white",
    dodge.width = 3) +
  ggplot2::scale_fill_manual(values = rev(pal)) +
  ggplot2::facet_wrap(~Genotype) +
  ggpubr::theme_classic2(base_size = 12) +
  ggplot2::theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5),
    strip.background = element_blank())


DMH_t <- dw_data  %>%
  dplyr::mutate(row_id = paste0("row_",rowid))%>%
  dplyr::filter(row_id %in% curves$sample) %>%
  dplyr::filter(DAP ==121) %>%
  rstatix::t_test(formula =  DW ~ txg, detailed = TRUE)

DMH_t[c(2,5),]

dw_data  %>%
  dplyr::mutate(row_id = paste0("row_",rowid))%>%
  dplyr::filter(row_id %in% curves$sample) %>%
  dplyr::filter(DAP ==40) %>%
  rstatix::t_test(formula =  DW ~ txg)









#######

to_test <-  curves %>%
  dplyr::select(Rep,Treatment,Genotype,r) %>%
  pivot_wider(names_from = Treatment, values_from = r) %>%as.data.frame() %>%
  dplyr::filter(!is.na(LowP) & !is.na(HighP))


curves %>% rstatix::t_test( k ~ Treatment, paired = TRUE)

#filter for rowids in better_fit, replace
curves$group <- factor(paste(curves$Rep,curves$Genotype))
curves$rtg <- factor(paste(curves$Rep,curves$Treatment,curves$Genotype))
table(curves$rtg)
paired <- names(table(curves$group)[table(curves$group) ==2])
curves$txg <- factor(interaction(curves$Treatment,curves$Genotype))
comp <- list(c("LowP.INV4M","HighP.INV4M"),
             c("LowP.CTRL","HighP.CTRL"))

test1 <- curves %>%
  filter(group %in% paired) %>%
  dplyr::arrange(Rep,Genotype,Treatment) %>%
  dplyr::group_by(txg) %>%
  rstatix::t_test( r~1, mu=0,
                   # alternative = "greater" ,
                   #paired = TRUE,
                   # comparisons = comp,
                   detailed = TRUE,
                   p.adjust.method = "fdr")

test2 <-   curves %>%
  filter(group %in% paired) %>%
  dplyr::arrange(Rep,Genotype,Treatment) %>%
  dplyr::group_by(txg) %>%
  rstatix::t_test( k~1, mu=0,
                   # alternative = "greater" ,
                   #paired = TRUE,
                   # comparisons = comp,
                   detailed = TRUE,
                   p.adjust.method = "fdr")


curves %>%
  filter(group %in% paired) %>%
  dplyr::arrange(Rep,Genotype,Treatment) %>%
  dplyr::select(Rep,Genotype,Treatment,k) %>%
  pivot_wider(names_from = Treatment,values_from = k) %>%
  dplyr::rename(High_P = "HighP", Low_P = "LowP") %>%
  dplyr::mutate(P_response = High_P-Low_P) %>%
  rstatix::t_test( P_response ~ Genotype,
                   #alternative = "less",
                   paired = TRUE,
                   comparisons = comp,
                   detailed = TRUE,
                   p.adjust.method = "fdr")

sample %>%
  filter(group %in% paired) %>%
  dplyr::arrange(Rep,Genotype,Treatment) %>%
  dplyr::select(Rep,Genotype,Treatment,r) %>%
  pivot_wider(names_from = Treatment,values_from = r) %>%
  dplyr::rename(High_P = "HighP", Low_P = "LowP") %>%
  dplyr::mutate(P_response = High_P-Low_P) %>%
  rstatix::t_test( P_response ~ Genotype,
                   alternative = "less" ,
                   paired = TRUE,
                   comparisons = comp,
                   detailed = TRUE,
                   p.adjust.method = "fdr")

dry_weight %>%
  filter(group %in% paired) %>%
  dplyr::arrange(Rep,Genotype,Treatment) %>%
  #dplyr::select(Rep,Genotype,Treatment,t_mid) %>%
  lm(data=., r~Treatment*Genotype) %>%
  anova()

sample %>%
  filter(group %in% paired) %>%
  dplyr::arrange(Rep,Genotype,Treatment) %>%
  dplyr::select(Rep,Genotype,Treatment,t_mid) %>%
  pivot_wider(names_from = Treatment,values_from = t_mid) %>%
  dplyr::rename(High_P = "HighP", Low_P = "LowP") %>%
  dplyr::mutate(P_response = Low_P-High_P) %>%
  rstatix::t_test( P_response ~ Genotype,
                   alternative = "greater" ,
                   paired = TRUE,
                   comparisons = comp,
                   detailed = TRUE,
                   p.adjust.method = "fdr")


sample %>%
  filter(group %in% paired) %>%
  dplyr::arrange(Rep,Genotype,Treatment) %>%
  # dplyr::group_by(txg) %>%
  rstatix::t_test( k~txg,
                   # alternative = "greater" ,
                   #paired = TRUE,
                   comparisons = comp,
                   detailed = TRUE,
                   p.adjust.method = "fdr")

sample %>%
  #filter(group %in% paired) %>%
  dplyr::arrange(Rep,Genotype,Treatment) %>%
  # dplyr::group_by(txg) %>%
  rstatix::t_test( t_mid~txg,
                   # alternative = "greater" ,
                   #paired = TRUE,
                   comparisons = comp,
                   detailed = TRUE,
                   p.adjust.method = "fdr")


sample %>%
  filter(group %in% paired) %>%
  dplyr::arrange(Rep,Genotype,Treatment) %>%
  #dplyr::group_by(txg) %>%
  rstatix::t_test( t_mid~txg,
                   # alternative = "greater" ,
                   paired = TRUE,
                   comparisons = comp,
                   detailed = TRUE,
                   p.adjust.method = "fdr")






test3 <-  sample %>%
  # filter(group %in% paired) %>%
  arrange(Rep,Genotype,Treatment) %>%
  rstatix::t_test( t_mid ~ group1,
                   # alternative = "greater" ,
                   #paired = TRUE,
                   comparisons = comp,
                   detailed = TRUE,
                   p.adjust.method = "fdr")


?rstatix::t_test
?t.test
to_test <-  better_fit%>% dplyr::select(Rep,Treatment,Genotype,r) %>%
  dplyr::filter(Genotype == "INV4M") %>%
  pivot_wider(names_from = Treatment, values_from = r) %>%as.data.frame() %>%
  dplyr::rename(Low_P ="LowP",
                High_P ="HighP")
t.test( to_test$Low_P, to_test$High_P, paired = TRUE)

