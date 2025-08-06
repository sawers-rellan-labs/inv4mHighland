library(dplyr)
library(janitor)
library(ggplot2)

pal <- c("chartreuse2", "darkblue")

pal2 <- c("darkorchid","chartreuse2", "darkblue")


internal_standards<- c(
  "CUDA",
  "Cholesterol",
  "CHOLESTEROL_D7_H20",
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
  "TG_17_0",
  "TG_17_0_17_1_17_0_D5"
)


# filter internal standards
to_exclude<- c(internal_standards, "TG_57_6" #exclude the odd number lipid
               )

lipid_csv <- "data/PSU_Normalized_MSDial_240416.csv"
lipid_csv <- "data/PSU_RawData_MSDial_NewStdInt_240422.csv"
raw <- read.csv(lipid_csv,na.strings = c("","#N/A","NA","Inf"))
to_exclude <- to_exclude[to_exclude %in% colnames(raw)]
raw <- raw %>% dplyr::select(-all_of(to_exclude)) %>%
  filter(!grepl("Methanol",sampleID)) %>%
  filter(!grepl("Pool",sampleID))
nrow(raw)


plant_csv <- "../../data/22_NCS_PSU_LANGEBIO_FIELDS_PSU_P_field.csv"


# qfiltered <- gsub("L|R","S",y_filtered_bySample$samples$tube, perl =TRUE)


psu <- read.csv(plant_csv) %>%
  rename(Genotype = "Who.What", rowid = "P22." ) %>%
  # dplyr::filter(Genotype %in% c("B73", "CTRL","INV4M",  "NUE")) %>%
  dplyr::filter(Genotype %in% c("CTRL","INV4M"))  %>%
  droplevels()

sampleInfo <- read.csv('data/inv4mRNAseq_metadata.csv') %>%
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


pheno$rowid
pheno$x <- log10(pheno$PC_36_6)
pheno$x <- log10((pheno$PC_36_6)/(pheno$LPC_18_3))
pheno$x <- log10((pheno$PE_36_6)/(pheno$PC_36_6))
colnames(pheno)

lm( DTS ~ block + x ,data=pheno[pheno$leaf_tissue==1,]) %>% summary()
lm( DTS ~ block + x ,data=pheno[pheno$leaf_tissue==2,]) %>% summary()
lm( DTS ~ block + x ,data=pheno[pheno$leaf_tissue==3,]) %>% summary()
lm( DTS ~ block + x ,data=pheno[pheno$leaf_tissue==4,]) %>% summary()
lm( DTS ~ block + x ,data=pheno[pheno$leaf_tissue==2,]) %>% summary()
lm( DTA ~ block + x ,data=pheno[pheno$leaf_tissue==1,]) %>% summary()


pheno$x <- log10(pheno$PC_36_6)

quartz()
pheno %>%
  ggplot(aes(x= x, y= DTA, 
             group=factor(leaf_tissue), 
             color=leaf_tissue)) +
  ggtitle("log10(PC_36_6)") +
  geom_smooth(method = "lm") +
  geom_point()

pheno$x <- log10((pheno$PC_36_6)/(pheno$LPC_18_3))

quartz()
pheno %>%
  ggplot(aes(x= x, y= DTA, 
             group=factor(leaf_tissue), 
             color=leaf_tissue)) +
  
  ggtitle("log10([(PC_36_6)/(LPC_18_3)]") +
  geom_smooth(method = "lm") +
  geom_point()

pheno$x <- log10((pheno$PE_36_6)/(pheno$PC_36_6))


to_plot <- pheno %>% ungroup()%>%
  group_by(rowid,Plot_Row,Plot_Column,block,DTA,DTS,Genotype, Treatment, leaf_group) %>%
  filter(leaf_group =="basal") %>%
  summarise(PE = mean(PE_36_6),
            PC = mean(PC_36_6),
            PC_PE = mean((PC_36_6)/(PE_36_6))) %>%
  mutate(PC_PE_ratio = log10(PC_PE)) 

lm(DTA~log10(PC_PE)+ Genotype +Treatment, data = to_plot) %>% summary()

quartz()
to_plot %>%
  ggplot(aes(x= PC_PE_ratio, y= DTA, group=Genotype, color= Genotype)) +
  xlab("log10[(PC_18_3+1)/(PE_36_6+1)]") +
  geom_point() +
  geom_smooth(method = "lm") +
  ggpubr::theme_classic2(base_size = 20)+
  theme(legend.position = "top")


pheno %>%
  group_by(rowid,Plot_Row,Plot_Column,block,DTA,DTS,Treatment,leaf_group) %>%
  filter(leaf_group =="basal") %>%
  summarise(PC_36_6 = log10(mean(PC_36_6))) %>%
  ggplot(aes(x= PC_36_6 , y= DTA, group=factor(Treatment), color= Treatment)) +
  ggtitle("log10(PC_36_6)") +
  geom_smooth(method = "lm") +
  geom_point()

to_plot <- pheno %>%
  filter(leaf_group=="basal") %>%
  group_by(rowid,Plot_Row,Plot_Column,block,DTA,DTS,Genotype,Treatment) %>%
  summarise(PC_36_6 = log10(mean(PC_36_6)))

to_plot$rowid
to_plot %>% nrow()


lm(DTA~Treatment*Genotype+ PC_36_6 , data = to_plot) %>% summary()

to_plot %>%
  ggplot(aes(x= PC_36_6 , y= DTA, group=factor(Genotype), color= Genotype)) +
  ggtitle("log10(PC_36_6)") +
  geom_smooth(method = "lm") +
  geom_point()



pheno %>%
  group_by(rowid,Plot_Row,Plot_Column,block,DTA,DTS,Genotype) %>%
  filter(leaf_group=="basal") %>%
  summarise(PC_36_6 = log10(mean(PC_36_6+1))) %>%
  ggplot(aes(x= PC_36_6 , y= DTS, group=factor(Genotype), color= Genotype)) +
  ggtitle("log10(PC_36_6)") +
  geom_smooth(method = "lm") +
  geom_point()


pheno$x <- log10((pheno$PC_36_6)/(pheno$LPC_18_3))
  
lm( DTA ~  Genotype +x,data=pheno) %>% summary()


  
b <- pheno%>%
  dplyr::select(DGDG_34_0:TG_58_5) %>%
  as.matrix() 

#b[,"y"] <- b[,"PC_36_6"]/b[,"LPC_18_3"]


m <- log10(b)

m[m==-Inf] <-NA
m[m==Inf] <-NA
colnames(m)



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



X <- cbind(exp_group = interaction(pheno$Genotype,pheno$Treatment), pca1$call$X)

emod <-lda(exp_group ~ ., data = X )

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
pheno_LD %>%
  ggplot(aes(x = LD1, y=LD2)) +
  ggtitle("MS-DIAL normalized") +
  geom_point(aes(shape = Genotype, color = Treatment), size =5) +
  scale_color_manual(values = rev(pal)) +
  scale_shape_manual(values=c(21,19,17,24)) +
  ggpubr::theme_classic2(base_size = 20)



# Multitrait multiple regression --------------------
# Probably I should use mashr with this too!
# Multivariable multiple regression
colnames(pheno)
vars <- colnames(m)
summary(m)
dim(m)
Y <- m

scaled_Y <- Y %>% as.matrix() %>% scale()


Covariates <- pheno%>%
  dplyr::select(Plot_Row,Plot_Column,block_6:block_12) %>%
  as.matrix()

Genotype  <- matrix(as.factor(pheno$Genotype) %>% as.numeric())

Treatment  <- matrix(pheno$Treatment %>% as.numeric())

leaf_tissue <- matrix(as.factor(pheno$leaf_tissue) %>% as.numeric())
leaf_group <- matrix(as.factor(pheno$leaf_group) %>% as.numeric())

mlm1 <- lm(scaled_Y ~ Covariates + leaf_tissue + Treatment + Genotype )
mlm2 <- lm(scaled_Y ~ Covariates  + leaf_tissue + Treatment*Genotype)
mlm3 <- lm(scaled_Y ~ Covariates  + leaf_tissue*Treatment*Genotype )
mlm4 <- lm(scaled_Y ~ Covariates + Treatment*Genotype )
mlm5 <- lm(scaled_Y ~ Covariates + leaf_group + Treatment*Genotype )




result <- mlm5 %>% summary()


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
  filter(p.value < 0.05) %>%
  filter(!grepl("Covariate",predictor)) %>%
  filter(!grepl("Intercept",predictor)) %>% arrange(p.value)


effects$p_adj <- p.adjust(effects$p.value, method="fdr")

quartz()
hist(effects$p.value, breaks=40)

effects$CI.upper <- effects$Estimate  + qnorm(0.975)*effects$Std.Error
effects$CI.lower <- effects$Estimate  - qnorm(0.975)*effects$Std.Error


## plot effects ----

# get data toplot

to_plot <- effects %>%
  filter(p_adj< 0.05) %>%
  filter(!grepl("Covariate",predictor)) %>%
  filter(!grepl("Intercept",predictor)) %>% 
  mutate(sign= sign(Estimate)) %>%
  ungroup() %>%
  group_by(Response) %>%
  mutate(max_effect = max(abs(Estimate))) %>%
  ungroup() %>%
  group_by(predictor) %>%   # As a precaution / handle in a separate .grouped_df method
  arrange(sign, abs(Estimate)) %>%   # arrange by facet variables and continuous values
  mutate(.r = row_number()) %>% # Add a row number variable %>%
  ungroup()


pd = position_dodge(0.4)

colnames(to_plot)
quartz()
to_plot %>%
  mutate(predictor =gsub("_tissue","", predictor)) %>%
  mutate(predictor=factor(predictor, levels= c("leaf_group","Treatment","Genotype", "Genotype:Treatment"))) %>%
  droplevels() %>%
  ggplot( aes(x     = Estimate,
              y     = Response,
              shape  = predictor)) +
  ggtitle("FDR < 0.05") +
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
  facet_wrap(.~predictor, scales="free_x", ncol=4)  +
  # scale_shape_manual(values=c(13,5,21),
  #                   labels= c( expression(italic("Inv4m") %*% P), "leaf","P")
  #                   ) +
  #ggpubr::theme_classic2(base_size = 25) +
  theme(legend.position ="none",
        axis.title.y=element_blank(),
        axis.text.y = element_text(hjust = 0, face = "bold"),
        plot.caption = element_text(hjust = 0))



# get data toplot
to_plot <- effects %>%
  mutate(sign= sign(Estimate)) %>%
  filter(p.value< 0.05) %>%
  filter(predictor != "leaf_tissue") %>%
  filter(!grepl("Covariate",predictor)) %>%
  filter(!grepl("Intercept",predictor)) %>%
  droplevels() %>%
  ungroup() %>%
  group_by(Response) %>%
  mutate(max_effect = max(abs(Estimate))) %>%
  ungroup() %>%
  group_by(predictor) %>%   # As a precaution / handle in a separate .grouped_df method
  arrange(sign, abs(Estimate)) %>%   # arrange by facet variables and continuous values
  mutate(.r = row_number()) %>% # Add a row number variable %>%
  ungroup() 


quartz()
to_plot %>%
  mutate(predictor =gsub("_tissue","", predictor)) %>%
  droplevels() %>%
  mutate(predictor=factor(predictor, levels= c("Treatment","Genotype", "Genotype:Treatment"))) %>%
  ggplot( aes(x     = Estimate,
              y     = Response,
              shape  = predictor)) +
  ggtitle("Not leaf, p < 0.05") +
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
  facet_wrap(.~predictor)  +
  # scale_shape_manual(values=c(13,5,21),
  #                   labels= c( expression(italic("Inv4m") %*% P), "leaf","P")
  #                   ) +
  #ggpubr::theme_classic2(base_size = 25) +
  theme(legend.position ="none",
        axis.title.y=element_blank(),
        axis.text.y = element_text(hjust = 0, face = "bold"),
        plot.caption = element_text(hjust = 0))

write.csv(
  to_plot %>%
    mutate(predictor =gsub("_tissue","", predictor)) %>%
    droplevels() %>% arrange(predictor),
  file="results_psu_2022/inv4m_+P_effects_in_lipids.csv")

colnames(pheno)

# pheno$y <-log10(pheno$TG_52_6)
#pheno$y <-log10(pheno$LPE_18_3+1)
#pheno$y <-log10((pheno$LPE_18_3+1)/(pheno$LPC_18_3+1))
# hist(pheno$y)
# pheno$y <-log10(pheno$DGDG_40_2+1)
# pheno$y <-log10(pheno$SQDG_36_6+1)

pheno$y <- log10((pheno$PC_36_6)/(pheno$LPC_18_3))
# pheno$y <- log10((pheno$LPE_18_3+1)/(pheno$PE_34_3+1))
# pheno$y <- log10((pheno$PE_34_3+1)/(pheno$LPE_18_3+1))


lm(data=pheno, y ~ block + leaf_tissue + Genotype + Treatment)  %>% summary()
lm(data=pheno, y ~ leaf_tissue + Genotype*Treatment)  %>% summary()
lm(data=pheno, y ~ leaf_group+Genotype*Treatment) %>% summary()
lm(data=pheno, y ~  Plot_Row + Plot_Column + block + leaf_group*Genotype*Treatment) %>% summary()

effects %>%
  filter(p_adj < 0.05) %>%
  arrange(p.value)




pheno$y <-  log10(pheno$LPE)
mod00 <- lm(data=pheno, y ~ Plot_Row + Plot_Column + block + leaf_tissue*Treatment*Genotype)
mod01 <- lm(data=pheno, y ~ Plot_Row + Plot_Column + block + leaf_tissue + Treatment*Genotype) 
mod02 <- lm(data=pheno, y ~ Plot_Row + Plot_Column + block + leaf_tissue + Treatment*Genotype)
mod03 <- lm(data=pheno, y ~ Plot_Row + block + leaf_tissue + Genotype + Treatment)
mod04 <- lm(data=pheno, y ~ Plot_Row + leaf_tissue*Treatment*Genotype)
mod05 <- lm(data=pheno, y ~ block + leaf_tissue + Treatment + Genotype )
mod06 <- lm(data=pheno, y ~ leaf_tissue + Treatment*Genotype)
mod07 <- lm(data=pheno, y ~ leaf_tissue*Treatment*Genotype)
mod08 <- lm(data=pheno, y ~ leaf_tissue*Treatment+ Genotype)
mod09 <- lm(data=pheno, y ~ leaf_tissue*Genotype+Treatment)
AIC(mod09,mod08,mod07,mod06,mod05,mod04,mod03, mod02,mod01,mod00) %>% arrange(AIC)
mod02 %>% summary() 
mod00 %>% summary() 
#pheno$y = predict(mod04,pheno)

pheno$y <-  predict(mod00,pheno)
emmeans::emmeans(mod04, specs = "Treatment")

mod00 %>% summary() 
mod04 %>% summary() 


mod11 <- lm(data=pheno, y ~  Plot_Row + Plot_Column + block + leaf_group*Genotype*Treatment) 
mod12 <- lm(data=pheno, y ~  Plot_Row + Plot_Column + block + leaf_group + Genotype*Treatment)
mod13 <- lm(data=pheno, y ~  block + leaf_group*Genotype*Treatment)
mod14 <- lm(data=pheno, y ~  block + leaf_group + Genotype*Treatment)
mod15 <- lm(data=pheno, y ~ leaf_group*Genotype*Treatment)
mod16 <- lm(data=pheno, y ~ leaf_group + Genotype*Treatment)
AIC(mod11,mod12,mod13,mod14, mod15, mod16) %>% arrange(AIC)

mod12 %>% summary()
mod16 %>% summary()




quartz()
pheno %>%
  mutate(leaf_tissue=as.factor(leaf_tissue)) %>%
  ggplot(aes(x=leaf_tissue,
             y=log10(DGGA_36_5),
             col = Genotype,
             fill= Genotype)) +
  geom_boxplot(position=position_dodge(width=1), fill ="white",outlier.shape = NA , width = 0.25)+
  ggbeeswarm::geom_quasirandom(dodge.width = 1, col="white",shape=21, width = 0.25, size =2) +
  facet_wrap(.~Treatment, ncol=2) +
  scale_fill_manual(values=pal)+
  scale_color_manual(values=pal)+
  ggpubr::theme_classic2(base_size = 20) +
  theme(legend.position ="top",
        strip.background = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5))

quartz()

plot_a <- pheno %>%
  filter(leaf_group == "apical") %>%
  # mutate(leaf_group=as.factor(leaf_group)) %>%
  ggplot(aes(x=Genotype,
             y=log10(PC_36_6),
             col = Genotype,
             fill= Genotype)) +
  ggtitle("Apical leaves") +
  geom_boxplot(position=position_dodge(width=1), fill ="white",outlier.shape = NA , width = 0.25)+
  ggbeeswarm::geom_quasirandom(dodge.width = 1, col="white",shape=21, width = 0.25, size =2) +
  facet_wrap(.~Treatment, ncol=2) +
  scale_fill_manual(values=pal)+
  scale_color_manual(values=pal)+
#  ylim(5.5,8.5) +
  ggpubr::theme_classic2(base_size = 20) +
  theme(strip.background = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5))

plot_b <- pheno %>%
  filter(leaf_group == "basal") %>%
  # mutate(leaf_group=as.factor(leaf_group)) %>%
  ggplot(aes(x=Genotype,
             y=log10(LPC_18_3),
             col = Genotype,
             fill= Genotype)) +
  ggtitle("Basal leaves") +
  geom_boxplot(position=position_dodge(width=1), fill ="white",outlier.shape = NA , width = 0.25)+
  ggbeeswarm::geom_quasirandom(dodge.width = 1, col="white",shape=21, width = 0.25, size =2) +
  facet_wrap(.~Treatment, ncol=2) +
  scale_fill_manual(values=pal)+
  scale_color_manual(values=pal)+
#  ylim(5.5,8.5) +
  ggpubr::theme_classic2(base_size = 20) +
  theme(strip.background = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5))







quartz()
ggpubr::ggarrange(plot_a,plot_b, ncol=2, align = "hv",common.legend = TRUE)

quartz()
pheno %>%
  #  filter(Genotype== "INV4") %>%
  mutate(leaf_tissue=as.factor(leaf_tissue)) %>%
  ggplot(aes(x=Treatment,
             y=log10(SQDG_36_6),
             col = Genotype,
             fill= Genotype)) +
  geom_boxplot(position=position_dodge(width=1), fill ="white",outlier.shape = NA )+
  ggbeeswarm::geom_quasirandom(dodge.width = 1, col="white",shape=21,width = 0.2) +
  # facet_wrap(.~Genotype, ncol=2) +
  scale_fill_manual(values=pal)+
  scale_color_manual(values=pal)+
  ggpubr::theme_classic2(base_size = 20) +
  theme(strip.background = element_blank())

quartz()
 pheno %>%
  mutate(leaf_tissue=as.factor(leaf_tissue)) %>%
  ggplot(aes(x=leaf_tissue,
             y=log10(PC_36_6/LPC_18_3),
             col = Treatment,
             fill= Treatment)) +
  geom_boxplot(position=position_dodge(width=1), fill ="white",outlier.shape = NA, width =0.25 )+
  ggbeeswarm::geom_quasirandom(dodge.width = 1, col="white",shape=21, width = 0.25) +
  scale_fill_manual(values=pal)+
  scale_color_manual(values=pal)+
  facet_wrap(.~Genotype, ncol=2) +
  ggpubr::theme_classic2(base_size = 20)+
  theme(strip.background = element_blank())





# quartz()
# pheno %>%
#   ggplot(aes(x=PC_36_6,
#              y=DTS,
#              group=leaf_tissue) +
#            facet_wrap(.~Treatment))


p_pc <- pheno %>%
  mutate(leaf_tissue=as.factor(leaf_tissue)) %>%
  ggplot(aes(x=leaf_tissue,
             y=log10(PC_36_6),
             col = Treatment,
             fill= Treatment)) +
  geom_boxplot(position=position_dodge(width=1), fill ="white",outlier.shape = NA, width =0.25 )+
  ggbeeswarm::geom_quasirandom(dodge.width = 1, col="white",shape=21, width = 0.25) +
  scale_fill_manual(values=pal)+
  scale_color_manual(values=pal)+
  ggpubr::theme_classic2(base_size = 20)

p_lpc <- pheno %>%
  mutate(leaf_tissue=as.factor(leaf_tissue)) %>%
  ggplot(aes(x=leaf_tissue,
             y=log10(LPC_18_3),
             col = Treatment,
             fill= Treatment)) +
  geom_boxplot(position=position_dodge(width=1), fill ="white",outlier.shape = NA, width =0.25 )+
  ggbeeswarm::geom_quasirandom(dodge.width = 1, col="white",shape=21, width = 0.25) +
  scale_fill_manual(values=pal)+
  scale_color_manual(values=pal)+
  ggpubr::theme_classic2(base_size = 20)




p_ratio <- pheno %>%
  mutate(leaf_tissue=as.factor(leaf_tissue)) %>%
  ggplot(aes(x=leaf_tissue,
             y=log10(PC_36_6/LPC_18_3),
             col = Treatment,
             fill= Treatment)) +
  geom_boxplot(position=position_dodge(width=1), fill ="white",outlier.shape = NA, width =0.25 )+
  ggbeeswarm::geom_quasirandom(dodge.width = 1, col="white",shape=21, width = 0.25) +
  scale_fill_manual(values=pal)+
  scale_color_manual(values=pal)+
  ggpubr::theme_classic2(base_size = 20)

quartz()
ggpubr::ggarrange( p_pc,p_lpc,p_ratio, 
                   align = "hv", 
                   common.legend = TRUE,
                   ncol=3)



quartz()
 pheno %>%
   mutate(leaf_tissue=as.factor(leaf_tissue)) %>%
   ggplot(aes(x=leaf_tissue,
              y=log10(LPC_18_1),
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
              y=log10(LPC_18_3),
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
              y=log10(LPC_18_2_a),
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
              y=log10(LPC_18_2_b),
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
              y=log10(PC_36_6/LPC_18_3),
              col = Treatment,
              fill= Treatment)) +
   geom_boxplot(position=position_dodge(width=1), fill ="white",outlier.shape = NA )+
   ggbeeswarm::geom_quasirandom(dodge.width = 1, col="white",shape=21) +
   scale_fill_manual(values=pal)+
   scale_color_manual(values=pal)+
   ggpubr::theme_classic2(base_size = 20)
 