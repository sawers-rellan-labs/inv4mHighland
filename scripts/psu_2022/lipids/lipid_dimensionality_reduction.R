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


plant_csv <- "../../data/22_NCS_PSU_LANGEBIO_FIELDS_PSU_P_field.csv"


psu <- read.csv(plant_csv) %>%
  rename(Genotype = "Who.What", rowid = "P22." ) %>%
  # dplyr::filter(Genotype %in% c("B73", "CTRL","INV4M",  "NUE")) %>%
  dplyr::filter(Genotype %in% c("CTRL","INV4M"))  %>%
  droplevels()

sampleInfo <- read.csv('data/PSU-PHO22_Metadata.csv') %>%
  #rename(Genotype=genotype) %>% 
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



# MDS with raw data -----



m <- pheno%>%
  dplyr::select(DGDG_34_0:TG_58_5) %>%
  as.matrix() 

counts <-  pheno%>%
  dplyr::select(DGDG_34_0:TG_58_5) %>%
  as.matrix() %>% 
  t() 

hist(counts)

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

quartz()
mds = plotMDS(main ="Lipid MDS", y_filtered,pch=21,label = y$samples$side_tag)


d= y_filtered$samples
d$x = mds$x
d$y = mds$y

quartz()
ggplot(d,aes(x=x,y=y)) + geom_point(aes(color = Treatment))

quartz()
ggplot(d,aes(x=x,y=y)) + geom_point(aes(color = decimal_time ))



quartz()
ggplot(d,aes(x=x,y=y)) + geom_point(aes(color = COLLECTOR))

quartz()
ggplot(d,aes(x=x,y=y)) + geom_point(aes(color = Genotype))

d$Treatment <- factor(d$Treatment)
d$Rep <- factor(d$Rep)

quartz()
ggplot(d,aes(x=x,y=y)) + geom_point(aes(color = factor(leaf_tissue),shape = Treatment),size=3) 

# PCA
pheno <- with_dup %>% 
  filter(!grepl("Pool", sampleID)) %>%
  inner_join(sampleInfo)
  
m <- pheno%>%
  dplyr::select(DGDG_34_0:TG_58_5) %>%
  as.matrix() %>% log2()

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


pheno2<- cbind(pheno,pca1$ind$coord)


quartz()
pheno2%>%
  ggplot2::ggplot(aes(x=Dim.1,y=Dim.2, color  = Treatment, group = Treatment, shape=Genotype)) +
  ggplot2::geom_point()+
  ggplot2::scale_color_manual(values = pal)


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


# LDA -----
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
  ggtitle("MS-DIAL raw integration corrected") +
  geom_point(aes(shape = Genotype, color = Treatment), size =5) +
  scale_color_manual(values = rev(pal)) +
  scale_shape_manual(values=c(21,19,17,24)) +
  ggpubr::theme_classic2(base_size = 20)


