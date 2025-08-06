library(dplyr)
library(janitor)
library(ggplot2)

pal <- c("chartreuse2", "darkblue")

pal2 <- c("darkorchid","chartreuse2", "darkblue")


lipid_csv <- "data/PSU_inv4m_lipids.csv"
raw <- read.csv(lipid_csv,na.strings = c("","#N/A","NA")) %>%
  janitor::clean_names()
colnames(raw) <- gsub("_sum_normalized_area","",colnames(raw))

plant_csv <- "../../data/22_NCS_PSU_LANGEBIO_FIELDS_PSU_P_field.csv"


# qfiltered <- gsub("L|R","S",y_filtered_bySample$samples$tube, perl =TRUE)


psu <- read.csv(plant_csv) %>%
  rename(Genotype = "Who.What", rowid = "P22." ) %>%
  mutate(Treatment = factor(Treatment, levels = c("LowP", "HighP")))  %>%
  # dplyr::filter(Genotype %in% c("B73", "CTRL","INV4M",  "NUE")) %>%
  dplyr::filter(Genotype %in% c("CTRL","INV4M"))  %>%
  droplevels()

sampleInfo <- read.csv('data/inv4mRNAseq_metadata.csv') %>%
  rename(Genotype=genotype)

sampleInfo$tube  <- gsub("L|R","S",sampleInfo$tube, perl =TRUE)

raw$tube <- gsub(".*_|.raw","",raw$file_name, perl =TRUE)
raw$tube <- gsub("L|R","S",raw$tube, perl =TRUE)
colnames(raw)

colnames(m)[grepl("^lpc",colnames(m), perl =TRUE)]

sampleInfo$tube
pheno <- sampleInfo %>%
  dplyr::select(tube,Plot:Plot_Column,Treatment:leaf_tissue) %>%
  inner_join(
    raw %>%
      dplyr::select(tube,sphingosine_d_17_1:lpc_18_3)) %>%
  mutate(block=as.factor(Rep)) %>%
  droplevels() %>%
  tidyr::pivot_wider(names_from = Rep, values_from = Rep, names_prefix = "block_",
                     values_fn = function(x) 1, values_fill = 0)  %>%
  dplyr::select(tube:Plot_Column,block,block_6:block_12,everything()) %>%
  mutate(Treatment = factor(Treatment,levels=c("Low_P","High_P")))
colnames(pheno)
#  filter(tube %in% qfiltered)


b <- pheno%>%
  dplyr::select(sphingosine_d_17_1:lpc_18_3) %>%
  dplyr::select(-cholesterol_d7_h20,-tg_17_0_17_1_17_0_d5, -pc_40_3) %>%
  mutate(y=NA) %>%
  as.matrix() +1

b[,"y"] <- b[,"pc_36_6"]/b[,"lpc_18_3"]

quartz()

m <- log10(b)


library(FactoMineR)
library(factoextra)
library(MASS)
pca1 <- PCA(m, ncp=10, graph = FALSE)

pca1$var$coord[,1]
X <- as.matrix(pca1$call$X)

psych::cortest.bartlett(X)
psych::cortest.bartlett(cor(X), n =120)


pca2 <- psych::principal(X, nfactors=10, rotate="oblimin", scores=TRUE)

scores <- pca1$ind$coord
colnames(scores)
pheno2<- cbind(pheno,pca1$ind$coord)


quartz()
pheno2%>%
  ggplot2::ggplot(aes(x=Dim.1,y=Dim.2, color  = Treatment, group = Treatment, shape=Genotype)) +
  ggplot2::geom_point()+
  ggplot2::scale_color_manual(values = pal)


table(pheno$leaf_tissue,pheno$Genotype)
table(pheno$leaf_tissue,pheno$Treatment)
which.min(X[,1])

plots<-list()

plots$dim12 <- with(pheno,{
  fviz_pca_ind(pca1,axes = c(1,2), geom ="point", alpha=0) +
    geom_point(aes(shape = Genotype, color = Treatment)) +
    ggtitle("Morphology PCA") +
    scale_color_manual(values = rev(pal)) +
    scale_shape_manual(values=c(21,19,17,24))
})




plots$dim34 <- with(pheno,{
  fviz_pca_ind(pca1,axes = c(3,4), geom ="point", alpha=0) +
    ggtitle("Morphology PCA") +
    geom_point(aes(shape = Genotype, color = Treatment)) +
    scale_color_manual(values = rev(pal)) +
    scale_shape_manual(values=c(21,19,17,24))
}
)


plots$dim56 <-  with(pheno,{
  fviz_pca_ind(pca1,axes = c(5,6), geom ="point", alpha=0) +
    ggtitle("Morphology PCA") +
    geom_point(aes(shape = Genotype, color = Treatment)) +
    scale_color_manual(values = rev(pal)) +
    scale_shape_manual(values=c(21,19,17,24))
}
)



plot_con <- function(pca, axes = c(1,2)){
  l <- sqrt(length(axes))
  to_grid <- lapply( axes,
                     function(pc){
                       fviz_contrib(pca, choice="var", axes = pc, top = ncol(pca$var$coord)) +
                         ggplot2::ggtitle(paste0("PC",pc))
                     }
  )
  ggpubr::ggarrange(plotlist = to_grid,
                    ncol = ceiling(l),
                    nrow = floor(l))
}

pca1$var$contrib
plots$con <- plot_con(pca1,1:6)

plots$cor_pc12 <- fviz_pca_var(pca1,axes = c(1,2))

plots$cor_pc13 <- fviz_pca_var(pca1,axes = c(1,3))

plots$cor_pc34 <- fviz_pca_var(pca1,axes = c(3,4))

plots$cor_pc56 <- fviz_pca_var(pca1,axes = c(5,6))

pdf(file = "results_psu_2022/PSU_2022_P_lipid_PCA.pdf")
lapply(plots,print)
dev.off()

pheno$Genotype <- factor(pheno$Genotype)
pheno$label <- if_else(pheno$Genotype=="INV4M","INV4M","OTHER")


#for random forest
X <- cbind(Genotype = factor(pheno$Genotype), pca1$call$X, Treatment =factor(pheno$Treatment))

X <- cbind(Genotype = factor(pheno$Genotype), pca1$call$X)

X <- cbind(Genotype = factor(pheno$Genotype), as.data.frame(pca1$ind$coord))



X <- cbind(Genotype = interaction(pheno$Genotype,pheno$Treatment), pca1$call$X)

colnames(X)
emod <-lda(Genotype ~ ., data = X )

quartz()
plot(emod)


emod.val <-NULL
emod.val <- predict(emod,X[,-1])

pheno_LD <- cbind(pheno,emod.val$x)

quartz()
barplot(sort(emod$scaling[,"LD1"]), las =2)

quartz()
barplot(sort(emod$scaling[,"LD2"]), las =2)

quartz()
barplot(sort(-abs(emod$scaling[,"LD1"])), las =2)


quartz()
barplot(sort(-abs(emod$scaling[,"LD2"])), las =2)


quartz()
pheno %>%
  ggplot(aes(x=log10(pc_36_6), y=log10(pc_36_2), label= leaf_tissue,
             color=Treatment)) +
  geom_text()








quartz()
pheno_LD %>%
  ggplot(aes(x = LD1, y=LD2)) +
  ggtitle("Skyline Normalized") +
  geom_point(aes(shape = Genotype, color = Treatment)) +
  scale_color_manual(values = rev(pal)) +
  scale_shape_manual(values=c(21,19,17,24)) +
  ggpubr::theme_classic2(base_size = 20)



# Multitrait multiple regression --------------------
# Probably I should use mashr with this too!
# Multivariable multiple regression
colnames(pheno)
vars <- colnames(m)

Y <- m

scaled_Y <- Y %>% as.matrix() %>% scale()

Covariates <- pheno %>%
  dplyr::select(Plot_Row,Plot_Column,block_6:block_12) %>%
  as.matrix()

Genotype  <- matrix(as.factor(pheno$Genotype) %>% as.numeric())

Treatment  <- matrix(pheno$Treatment %>% as.numeric())

leaf_tissue <- matrix(as.factor(pheno$leaf_tissue) %>% as.numeric())

mlm1 <- lm(scaled_Y ~ Covariates + leaf_tissue + Genotype + Treatment)
mlm2 <- lm(scaled_Y ~ Covariates  + leaf_tissue + Genotype*Treatment )
mlm1 <- lm(scaled_Y ~ Covariates + leaf_tissue + Genotype + Treatment)
mlm3 <- lm(scaled_Y ~ Covariates  + Genotype*Treatment )


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
  filter(p.value < 0.1) %>%
  filter(!grepl("Covariate",predictor)) %>%
  filter(!grepl("Intercept",predictor)) %>% arrange(p.value)


effects$p_adj <- p.adjust(effects$p.value, method="fdr")

quartz()
hist(effects$p.value)

effects$CI.upper <- effects$Estimate  + qnorm(0.975)*effects$Std.Error
effects$CI.lower <- effects$Estimate  - qnorm(0.975)*effects$Std.Error


effects %>%
  filter(p.value < 0.1) %>%
  filter(!grepl("Covariate",predictor)) %>%
  filter(!grepl("Intercept",predictor)) %>% arrange(p.value)

pheno$y <- log10(pheno$pc_36_6/pheno$lpc_18_3)
lm(data=pheno, y~  Plot_Row + Plot_Column + block + leaf_tissue +Genotype*Treatment) %>% summary()

effects %>%
  filter(p.value < 0.1) %>%
  filter(!grepl("Covariate",predictor)) %>%
  filter(!grepl("Intercept",predictor)) %>% arrange(p.value)

pheno$y <- log10(pheno$pc_36_6/pheno$lpc_18_3)
hist(pheno$y)
# pheno$y <-log(pheno$lpc_18_3)



mod00 <- lm(data=pheno, y ~  Plot_Row + Plot_Column + block + leaf_tissue + Genotype*Treatment) 

mod01 <- lm(data=pheno, y~ Plot_Row + Plot_Column + block + leaf_tissue + Genotype+Treatment)

mod02 <- lm(data=pheno, y ~  block + leaf_tissue + Genotype+ Treatment)

mod03 <- lm(data=pheno, y ~ leaf_tissue + Genotype*Treatment)

mod04 <- lm(data=pheno, y ~  Plot_Row + Plot_Column + block + leaf_tissue*Genotype*Treatment) 

anova(mod00)
anova(mod02)

anova(mod04,mod03, mod02,mod01,mod00)
AIC(mod03, mod02,mod01,mod00)

lm(data=pheno, log10(pheno$pc_36_6) ~  Plot_Row + Plot_Column + block + leaf_tissue + Genotype*Treatment) %>% summary()

lm(data=pheno, log10(pheno$pc_36_6) ~  Plot_Row + Plot_Column + block + leaf_tissue + Genotype*Treatment) %>% summary()

anova(mod00)
summary(mod00)


quartz()
pheno %>%
#  filter(Genotype== "INV4") %>%
  mutate(leaf_tissue=as.factor(leaf_tissue)) %>%
  ggplot(aes(x=leaf_tissue,
             y=log10(pc_36_6/lpc_18_3),
             col = Treatment,
             fill= Treatment)) +
  geom_boxplot(position=position_dodge(width=1), fill ="white",outlier.shape = NA )+
  ggbeeswarm::geom_quasirandom(dodge.width = 1, col="white",shape=21,width = 0.2) +
  scale_fill_manual(values=pal)+
  scale_color_manual(values=pal)+
  ggpubr::theme_classic2(base_size = 20)


quartz()
 pheno %>%
  mutate(leaf_tissue=as.factor(leaf_tissue)) %>%
  ggplot(aes(x=leaf_tissue,
             y=log10(pc_36_6/lpc_18_3),
             col = Treatment,
             fill= Treatment)) +
  geom_boxplot(position=position_dodge(width=1), fill ="white",outlier.shape = NA )+
  ggbeeswarm::geom_quasirandom(dodge.width = 1, col="white",shape=21) +
  scale_fill_manual(values=pal)+
  scale_color_manual(values=pal)+
  facet_wrap(.~Genotype, ncol=2) +
  ggpubr::theme_classic2(base_size = 20)+
  theme(strip.background = element_blank())





# quartz()
# pheno %>%
#   ggplot(aes(x=pc_36_6,
#              y=DTS,
#              group=leaf_tissue) +
#            facet_wrap(.~Treatment))


p_pc <- pheno %>%
  mutate(leaf_tissue=as.factor(leaf_tissue)) %>%
  ggplot(aes(x=leaf_tissue,
             y=log10(pc_36_6),
             col = Treatment,
             fill= Treatment)) +
  geom_boxplot(position=position_dodge(width=1), fill ="white",outlier.shape = NA )+
  ggbeeswarm::geom_quasirandom(dodge.width = 1, col="white",shape=21) +
  scale_fill_manual(values=pal)+
  scale_color_manual(values=pal)+
  ggpubr::theme_classic2(base_size = 20)

p_lpc <- pheno %>%
  mutate(leaf_tissue=as.factor(leaf_tissue)) %>%
  ggplot(aes(x=leaf_tissue,
             y=log10(lpc_18_3),
             col = Treatment,
             fill= Treatment)) +
  geom_boxplot(position=position_dodge(width=1), fill ="white",outlier.shape = NA )+
  ggbeeswarm::geom_quasirandom(dodge.width = 1, col="white",shape=21) +
  scale_fill_manual(values=pal)+
  scale_color_manual(values=pal)+
  ggpubr::theme_classic2(base_size = 20)




p_ratio <- pheno %>%
  mutate(leaf_tissue=as.factor(leaf_tissue)) %>%
  ggplot(aes(x=leaf_tissue,
             y=log10(pc_36_6/lpc_18_3),
             col = Treatment,
             fill= Treatment)) +
  geom_boxplot(position=position_dodge(width=1), fill ="white",outlier.shape = NA )+
  ggbeeswarm::geom_quasirandom(dodge.width = 1, col="white",shape=21) +
  scale_fill_manual(values=pal)+
  scale_color_manual(values=pal)+
  ggpubr::theme_classic2(base_size = 20)

quartz()
ggpubr::ggarrange( p_pc,p_lpc,p_ratio, 
                   align = "hv", 
                   common.legend = TRUE,
                   ncol=3)
## plot effects ----

# get data toplot
to_plot <- effects %>%
  filter(p.value < 0.05) %>%
  filter(!grepl("Covariate",predictor)) %>%
  filter(!grepl("Intercept",predictor)) %>% arrange(p.value) %>%
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

quartz()
 to_plot %>%
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
    title = NULL, reverse= TRUE)) +
  scale_shape_manual(values=c(13,5,21),
                     labels= c( expression(italic("Inv4m") %*% P), "leaf","-P")
                     ) +
  ggpubr::theme_classic2(base_size = 25) +
  theme(legend.position ="top",
        axis.title.y=element_blank(),
        axis.text.y = element_text(hjust = 0, face = "bold"),
        plot.caption = element_text(hjust = 0))

quartz()
 pheno %>%
   mutate(leaf_tissue=as.factor(leaf_tissue)) %>%
   ggplot(aes(x=leaf_tissue,
              y=log10(lpc_18_1),
              col = Treatment,
              fill= Treatment)) +
   geom_boxplot(position=position_dodge(width=1), fill ="white",outlier.shape = NA )+
   ggbeeswarm::geom_quasirandom(dodge.width = 1, col="white",shape=21) +
   scale_fill_manual(values=pal)+
   scale_color_manual(values=pal)+
   ggpubr::theme_classic2(base_size = 20)
 
 
quartz()
 pheno %>%
   mutate(leaf_tissue=as.factor(leaf_tissue)) %>%
   ggplot(aes(x=leaf_tissue,
              y=log10(lpc_18_3),
              col = Treatment,
              fill= Treatment)) +
   geom_boxplot(position=position_dodge(width=1), fill ="white",outlier.shape = NA )+
   ggbeeswarm::geom_quasirandom(dodge.width = 1, col="white",shape=21) +
   scale_fill_manual(values=pal)+
   scale_color_manual(values=pal)+
   ggpubr::theme_classic2(base_size = 20)
 
quartz()
 pheno %>%
   mutate(leaf_tissue=as.factor(leaf_tissue)) %>%
   ggplot(aes(x=leaf_tissue,
              y=log10(lpc_18_2_a),
              col = Treatment,
              fill= Treatment)) +
   geom_boxplot(position=position_dodge(width=1), fill ="white",outlier.shape = NA )+
   ggbeeswarm::geom_quasirandom(dodge.width = 1, col="white",shape=21) +
   scale_fill_manual(values=pal)+
   scale_color_manual(values=pal)+
   ggpubr::theme_classic2(base_size = 20)
 
 quartz()
 pheno %>%
   mutate(leaf_tissue=as.factor(leaf_tissue)) %>%
   ggplot(aes(x=leaf_tissue,
              y=log10(lpc_18_2_b),
              col = Treatment,
              fill= Treatment)) +
   geom_boxplot(position=position_dodge(width=1), fill ="white",outlier.shape = NA )+
   ggbeeswarm::geom_quasirandom(dodge.width = 1, col="white",shape=21) +
   scale_fill_manual(values=pal)+
   scale_color_manual(values=pal)+
   ggpubr::theme_classic2(base_size = 20)

 
quartz()
 pheno %>%
   mutate(leaf_tissue=as.factor(leaf_tissue)) %>%
   ggplot(aes(x=leaf_tissue,
              y=log10(pc_36_6/lpc_18_3),
              col = Treatment,
              fill= Treatment)) +
   geom_boxplot(position=position_dodge(width=1), fill ="white",outlier.shape = NA )+
   ggbeeswarm::geom_quasirandom(dodge.width = 1, col="white",shape=21) +
   scale_fill_manual(values=pal)+
   scale_color_manual(values=pal)+
   ggpubr::theme_classic2(base_size = 20)
 