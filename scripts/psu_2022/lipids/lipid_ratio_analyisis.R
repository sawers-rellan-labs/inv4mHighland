
pal <- c("chartreuse2", "darkblue")

pal2 <- c("darkorchid","chartreuse2", "darkblue")



# filter internal standards
to_exclude<- c("DG_12_0", "LPC_17_0", "LPE_17_1",
               "PC_25_0", "PG_17_0","SM_35_1",
               "Sphingosine_17_1","TG_17_0",
               "TG_57_6" #exclude the odd number lipid
)



head_group <- c(
  DG ="glycolipid", 
  DGDG="glycolipid",
  DGGA="glycolipid",
  LPC="phospholipid",
  LPE="phospholipid",
  MGDG="glycolipid",
  PC="phospholipid",
  PE="phospholipid",
  PG="phospholipid",
  PI="phospholipid",
  SQDG="glycolipid",
  TG="glycolipid")

sp <- data.frame(
  colname = colnames(m),
  class = gsub("_.*","", colnames(m), perl =TRUE)
)

sp$head_group <-head_group[sp$class]

glyco <- sp$colname[sp$head_group == "glycolipid"]
phospho <- sp$colname[sp$head_group == "phospholipid"]

lipid_csv <- "data/PSU_Normalized_MSDial_240416.csv"
lipid_csv <- "data/PSU_RawData_MSDial_240416.csv"
raw <- read.csv(lipid_csv,na.strings = c("","#N/A","NA","Inf")) %>%
  dplyr::select(-all_of(to_exclude)) %>%
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



# with sums

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
  mutate(Treatment = factor(Treatment,levels=c("High_P","Low_P")))
m <- pheno %>%
  dplyr::select(DGDG_34_0:TG_58_5) %>% as.matrix()

# ananlysis by GROUP
class <- colnames(m) 
class <- gsub("_.*","", group, perl =TRUE) %>% sort() %>% unique()
n <- length(class)
pairs <- combn(n,2) %>% t
# ratios <- paste(group[pairs[,1]],group[pairs[,2]], sep = "_")

for (c in class){
  pheno[c] <- pheno  %>%
    dplyr::select(dplyr::starts_with(paste0(g,"_"))) %>%
    rowSums(na.rm = TRUE)
}

pheno$glyco <-   pheno  %>%
  dplyr::select(all_of(glyco)) %>%
  rowSums(na.rm = TRUE)

ratios <- paste(sp$colname,"glyco", sep="_")
for(num in sp$colname){
    pheno[paste(num,"glyco", sep="_")] <- (pheno[num]+0.5)/(pheno["glyco"]+0.5)
}



for(num in class[pairs[,1]]){
  for(den in class[pairs[,2]]){
  pheno[paste(num,den, sep="_")] <- pheno[num]/pheno[den]
}
}


m <- log2(pheno[,ratios] %>% as.matrix())

m[m==-Inf] <-NA
m[m==Inf] <-NA
colnames(m)







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

mlm0 <- lm(scaled_Y ~ Covariates + leaf_tissue*Treatment*Genotype )
mlm1 <- lm(scaled_Y ~ Covariates  + leaf_tissue + Treatment*Genotype)
mlm2 <- lm(scaled_Y ~ Covariates + leaf_tissue + Treatment + Genotype )
mlm3 <- lm(scaled_Y ~ Covariates  + leaf_tissue*Treatment*Genotype )
mlm4 <- lm(scaled_Y ~ Covariates + Treatment*Genotype )
mlm5 <- lm(scaled_Y ~ Covariates + leaf_group + Treatment*Genotype )
mlm6 <- lm(scaled_Y ~ Covariates + leaf_group*Treatment*Genotype )



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

effects$p_adj <- p.adjust(effects$p.value, method="fdr")

quartz()
hist(effects$p.value, breaks=40)

effects$CI.upper <- effects$Estimate  + qnorm(0.975)*effects$Std.Error
effects$CI.lower <- effects$Estimate  - qnorm(0.975)*effects$Std.Error

effects %>%
  filter(p.value< 0.05) %>%
  filter(!grepl("Covariate",predictor)) %>%
  filter(!grepl("Intercept",predictor)) %>%
  arrange(p.value)

## plot effects ----

# get data toplot

to_plot <- effects %>%
  mutate(Response = gsub("_glyco","",Response)) %>%
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

write.csv(
  to_plot %>%
    mutate(predictor =gsub("_tissue","", predictor)) %>%
    droplevels() %>% arrange(predictor),
  file="results_psu_2022/inv4m_+P_effects_in_lipid_sums_FDR.csv")


diff <- strsplit(to_plot$Response, split = "_") %>% unlist() %>% sort() %>% unique()
#diff <- to_plot$Response

quartz()
pheno %>%
  dplyr::select(leaf_tissue,Genotype,all_of(diff)) %>%
  tidyr::pivot_longer(cols = all_of(diff), values_to = "signal", names_to = "Response") %>%
  mutate(Response=forcats::fct_reorder(Response, signal)) %>%
  mutate(signal= signal) %>%
  ggplot(aes(x=Response,
             y=signal,
             col = Genotype,
             fill=  Genotype)) +
  geom_boxplot(position=position_dodge(width=1), fill ="white",outlier.shape = NA , width = 0.25) +
  ggbeeswarm::geom_quasirandom(dodge.width = 1, col="white",shape=21, width = 0.25, size =2) +
#  facet_wrap(.~Genotype, ncol=2) +
  scale_y_continuous(trans='log10') +
  scale_fill_manual(values=pal)+
  scale_color_manual(values=pal)+
  ggpubr::theme_classic2(base_size = 20) +
  theme(legend.position ="top",
        strip.background = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=90),
        plot.title = element_text(hjust = 0.5))



pd = position_dodge(0.4)

colnames(to_plot)
quartz()
to_plot %>%
  mutate(predictor =gsub("_group","", predictor)) %>%
  mutate(predictor =gsub("_tissue","", predictor)) %>%
  mutate(predictor =gsub("Treatment","-P", predictor)) %>%
  mutate(predictor =gsub("Genotype","Inv4m", predictor)) %>%
  mutate(predictor=factor(predictor, levels= c("leaf","-P","Inv4m", "-P:Inv4m"))) %>%
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




quartz()

pheno %>%
  mutate(leaf_tissue=as.factor(leaf_tissue)) %>%
  ggplot(aes(x=leaf_tissue,
             y=log10(DGDG_LPE),
             col = Treatment,
             fill=  Treatment)) +
  geom_boxplot(position=position_dodge(width=1), fill ="white",outlier.shape = NA , width = 0.25) +
  ggbeeswarm::geom_quasirandom(dodge.width = 1, col="white",shape=21, width = 0.25, size =2) +
  facet_wrap(.~Genotype, ncol=2) +
  scale_fill_manual(values=pal)+
  scale_color_manual(values=pal)+
  ggpubr::theme_classic2(base_size = 20) +
  theme(legend.position ="top",
        strip.background = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5))




quartz()

pheno %>%
  mutate(leaf_tissue=as.factor(leaf_tissue)) %>%
  ggplot(aes(x=leaf_tissue,
             y=log10(LPE_TG),
             col = Genotype,
             fill=  Genotype)) +
  geom_boxplot(position=position_dodge(width=1), fill ="white",outlier.shape = NA , width = 0.25) +
  ggbeeswarm::geom_quasirandom(dodge.width = 1, col="white",shape=21, width = 0.25, size =2) +
  facet_wrap(.~Treatment, ncol=2) +
  scale_fill_manual(values=pal)+
  scale_color_manual(values=pal)+
  ggpubr::theme_classic2(base_size = 20) +
  theme(legend.position ="top",
        strip.background = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5))




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

nrow(to_plot)

quartz()
to_plot %>%
  mutate(predictor =gsub("_tissue","", predictor)) %>%
  mutate(predictor =gsub("Treatment","-P", predictor)) %>%
  mutate(predictor =gsub("Genotype","Inv4m", predictor)) %>%
  mutate(predictor=factor(predictor, levels= c("-P","Inv4m", "-P:Inv4m"))) %>%
  droplevels() %>%
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
  file="results_psu_2022/inv4m_+P_effects_in_lipid_sums_raw_pvalue.csv")


y<- log10(pheno$LPE_SQDG)

mod00 <- lm(data=pheno, y ~ Plot_Row + Plot_Column + block + leaf_tissue + Genotype*Treatment) 
mod01 <- lm(data=pheno, y ~ Plot_Row + Plot_Column + block + leaf_tissue + Genotype+Treatment)
mod02 <- lm(data=pheno, y ~ Plot_Row + block + leaf_tissue + Genotype + Treatment)
mod03 <- lm(data=pheno, y ~ Plot_Row ++ leaf_tissue + Genotype + Treatment)
mod04 <- lm(data=pheno, y ~ Plot_Row + leaf_tissue*Genotype*Treatment)
mod05 <- lm(data=pheno, y ~ block + leaf_tissue + Genotype + Treatment)
mod06 <- lm(data=pheno, y ~ leaf_tissue + Genotype*Treatment)
mod07 <- lm(data=pheno, y ~ leaf_tissue*Genotype*Treatment)
mod08 <- lm(data=pheno, y ~ leaf_tissue*Genotype+Treatment)
AIC(mod08,mod07,mod06,mod05,mod04,mod03, mod02,mod01,mod00) %>% arrange(AIC)

mod00 %>% summary()


# individual species all_ratios ----

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

# ananlysis by GGROUP
spp <- colnames(pheno %>% dplyr::select(DGDG_34_0:TG_58_5))
n <- length(spp)
pairs <- combn(n,2) %>% t()

ratios <- paste(spp[pairs[,1]],spp[pairs[,2]], sep = "X")

for (r in ratios){
  pheno[r] <- pheno[gsub("X.*","",r,perl=TRUE)]/pheno[gsub(".*X","",r,perl=TRUE)]
}


m <- log10(pheno[,ratios] %>% as.matrix())

m[m==-Inf] <-NA
m[m==Inf] <-NA
colnames(m)


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

mlm0 <- lm(scaled_Y ~ Covariates + leaf_tissue*Treatment*Genotype )
mlm1 <- lm(scaled_Y ~ Covariates  + leaf_tissue + Treatment*Genotype)
mlm2 <- lm(scaled_Y ~ Covariates + leaf_tissue + Treatment + Genotype )
mlm3 <- lm(scaled_Y ~ Covariates + Treatment*Genotype )
mlm4 <- lm(scaled_Y ~ Covariates + leaf_group + Treatment*Genotype )




result <- mlm1 %>% summary()


effects <- lapply(names(result), FUN=function(response){
  r <- result[[response]]
  trait <- gsub(".* ", "", response, perl =TRUE)
  coeff <- r$coefficients
  data.frame(Response = trait, as.data.frame(coeff)) %>%
    tibble::rownames_to_column("predictor")
}) %>% dplyr::bind_rows()

summary(effects$p.value)
colnames(effects)[4] <- "Std.Error"
colnames(effects)[6] <- "p.value"

effects$p_adj <- p.adjust(effects$p.value, method="fdr")

quartz()
hist(effects$p.value, breaks=100)

quartz()
hist(effects$p_adj, breaks=100)

effects$CI.upper <- effects$Estimate  + qnorm(0.975)*effects$Std.Error
effects$CI.lower <- effects$Estimate  - qnorm(0.975)*effects$Std.Error


effects %>%
  filter(p.value < 0.05) %>%
  filter(!grepl("Covariate",predictor)) %>%
  filter(!grepl("Intercept",predictor)) %>% 
  filter(grepl("MGDG",Response, perl = TRUE)) %>% 
#  filter(!grepl("DGDG",Response)) %>% 
#  filter(grepl("SQDG_",Response)) %>%
  arrange(p.value)


quartz()
pheno %>%
  ggplot(aes(x=Treatment,
             y=log10(PC_36_1XTG_54_1),
             col = Genotype,
             fill=  Genotype)) +
  geom_boxplot(position=position_dodge(width=1), fill ="white",outlier.shape = NA , width = 0.25) 


quartz()
pheno %>%
  mutate(leaf_tissue=as.factor(leaf_tissue)) %>%
  ggplot(aes(x=leaf_tissue,
             y=log10(PC_TG),
             col = Treatment,
             fill=  Treatment)) +
  geom_boxplot(position=position_dodge(width=1), fill ="white",outlier.shape = NA , width = 0.25) +
  ggbeeswarm::geom_quasirandom(dodge.width = 1, col="white",shape=21, width = 0.25, size =2) +
 facet_wrap(.~Genotype, ncol=2) +
  scale_fill_manual(values=pal)+
  scale_color_manual(values=pal)+
  ggpubr::theme_classic2(base_size = 20) +
  theme(legend.position ="top",
        strip.background = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5))


