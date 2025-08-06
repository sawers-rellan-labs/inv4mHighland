library(dplyr)
library(janitor)
library(ggplot2)


# CHECK internal standards
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


#Skyline----

lipid_csv <- "data/PSU_inv4m_lipids.csv"
raw <- read.csv(lipid_csv,na.strings = c("","#N/A","NA")) %>%
janitor::clean_names(case="all_caps") %>%
dplyr::filter(!grepl("ZCN8",FILE_NAME)) %>%
  filter(!grepl("Control",FILE_NAME))   


raw$sampleID <- gsub(".*_|.raw","",raw$FILE_NAME) 
raw$sample_group <- gsub("_.*","",raw$FILE_NAME)
raw$tube <- gsub("L|R","S",raw$sampleID, perl =TRUE)
raw$tube
colnames(raw)
raw$sampleID


colnames(raw) <- gsub("_SUM_NORMALIZED_AREA","",colnames(raw))
INTERNAL_STANDARDS <- stringr::str_to_upper(internal_standards)
col_std <- sapply(INTERNAL_STANDARDS, function(x){
       grepl(x,colnames(raw))
  })

rownames(col_std ) <- colnames(raw)

is_std <- apply(col_std ,1, function(x) Reduce("|", x))
is_std[2] <-TRUE
in_raw <- colnames(raw)[is_std]


in_raw <- INTERNAL_STANDARDS[INTERNAL_STANDARDS %in% colnames(raw)] 

int_std <- raw %>%  dplyr::select(all_of(in_raw))
#dplyr::select(-all_of(to_exclude)) %>%
# filter(!grepl("Methanol",sampleID)) %>%




stat_summmary <- raw %>%
  dplyr::group_by(sample_group) %>%
  summarise(across(
    .cols = in_raw, 
    .fns = list(Mean = function(x){mean(x, na.rm = TRUE)},
                SD = function(x){sd(x, na.rm = TRUE)},
                CV_pct = function(x){100*sd(x, na.rm = TRUE)/mean(x, na.rm = TRUE)}), 
    .names = "{col}.{fn}"
  )) %>%
  tidyr::pivot_longer(cols=-"sample_group")%>%
  tidyr::pivot_wider(names_from="sample_group") %>%
  tidyr::separate(name,into = c("standard","metric"),sep = "\\.") %>%
  tidyr::pivot_longer(cols=c(3),names_to ="sample_group") %>% as.data.frame()
stat_summmary 


quartz()
stat_summmary %>%
  filter(sample_group=="230113") %>%
  group_by(metric) %>%
  mutate(standard=forcats::fct_reorder(standard, value,.na_rm = FALSE)) %>%
  droplevels() %>%
  ggplot(aes(y=value, x=standard,dodge=sample_group, fill=sample_group)) +
  ggtitle("Skyline Normalized, Internal Standard QC,\nPSU Inv4m Phosphorus experiment") +
  geom_col(position=position_dodge()) +
  facet_grid(metric~standard, scales= "free")



m <- raw %>%
  dplyr::select(-all_of(in_raw)) %>%
  dplyr::select(PC_12_0_13_0:LPC_18_3) %>%
  as.matrix() 

m[is.na(m)] <- 0

dim(m)
diag(m)
colnames(raw)

hist(m)

n_factor <- max(m)/rowSums(m)

norm <-cbind(
  data.frame(sampleID =raw$sampleID,
             sample_group=raw$sample_group,
             tube=raw$tube),
  diag(n_factor) %*%  m ) # maximum_signal/total_lipid_sum persample

colnames(norm)




# MSDIAL----

ms_order<- read.csv("data/PSU-PHO22 _ms_order.csv")
ms_order$tube <- gsub("L|R","S",ms_order$top_tag, perl =TRUE)

metadata <- read.csv('data/PSU-PHO22_Metadata.csv')
metadata$tube <- gsub("L|R","S",metadata$tube, perl =TRUE)



lipid_csv <- "data/PSU_RawData_MSDial_NewStdInt_240422.csv"


raw <- read.csv(lipid_csv,na.strings = c("","#N/A","NA","Inf")) %>%
  filter(!grepl("Pool_7",sampleID)) %>%
  filter(!grepl("Pool_6",sampleID)) %>%
  filter(sampleID !="221018_Methanol_1") %>%
  filter(sampleID !="230113_Methanol_6")

cat(raw$sampleID)
raw$sample_group <- gsub("_.*","",raw$sampleID) 

weird <- c(
"230113_1_R04",
"230113_14_R53",
"230113_15_R57",
"230113_2_R06",
"230113_5_R20")


pools <- grepl("Pool",raw$sampleID)
pools 
raw$sampleID[pools]
raw$sample_group[pools] <- "Pool1_5"

raw$sample_group[raw$sampleID %in% weird ] <- "Weird"

raw$sampleID <- gsub("Pool_","Pool",raw_dup$sampleID, perl =TRUE)
raw$sampleID <- gsub(".*_|.raw","",raw_dup$sampleID, perl =TRUE)
raw$tube <- gsub("L|R","S",raw$sampleID, perl =TRUE)
raw$tube
colnames(raw)
raw$sampleID

# normalize to lipid sum

m <- raw %>%
  dplyr::select(CUDA:TG_58_5) %>%
  # dplyr::select(-LPC_17_0) %>%
  as.matrix() 
dim(m)
diag(m)
colnames(raw)

hist(m)
 
n_factor <- max(m)/rowSums(m)

norm <-cbind(
  data.frame(sampleID =raw$sampleID,
             sample_group=raw$sample_group),
  diag(n_factor) %*%  m ) # maximum_signal/total_lipid_sum persample

quartz()
hist(norm[,-(1:2)] %>% as.matrix())
tail(norm)
norm$sampleID
# pool1_5 <-norm[norm$sample_group=="Pool1_6",]
# nrow(pool1_5)
# pool1_5 <- pool1_5[pool1_5$sampleID != "230113_Pool_6",]
# nrow(pool1_5)
# pool1_5$sample_group <- "Pool1_5"
# with_dup <- rbind(norm, pool1_5)
with_dup <-norm

with_dup$sampleID <- gsub("Pool_","Pool",with_dup$sampleID, perl =TRUE)
with_dup$sampleID <- gsub(".*_|.raw","",with_dup$sampleID, perl =TRUE)
with_dup$tube <- gsub("L|R","S",with_dup$sampleID, perl =TRUE)
with_dup$tube
colnames(with_dup)
with_dup$sampleID


in_norm <- internal_standards[internal_standards %in% colnames(norm)] 

internal_std <- norm %>%  dplyr::select(all_of(in_norm))
  #dplyr::select(-all_of(to_exclude)) %>%
  # filter(!grepl("Methanol",sampleID)) %>%


with_dup$sample_group[with_dup$sample_group =="230113"] <- "230113_other"
with_dup$sampleID

stat_summmary <- with_dup %>%
  dplyr::group_by(sample_group) %>%
  summarise(across(
    .cols = in_norm , 
    .fns = list(Mean = function(x){mean(x, na.rm = TRUE)},
                SD = function(x){sd(x, na.rm = TRUE)},
                CV_pct = function(x){100*sd(x, na.rm = TRUE)/mean(x, na.rm = TRUE)}), 
    .names = "{col}.{fn}"
  )) %>%
  tidyr::pivot_longer(cols=-"sample_group")%>%
  tidyr::pivot_wider(names_from="sample_group") %>%
  tidyr::separate(name,into = c("standard","metric"),sep = "\\.") %>%
  tidyr::pivot_longer(cols=-c(1,2),names_to ="sample_group") %>%
  group_by(sample_group,metric) %>%
  mutate(standard=forcats::fct_reorder(standard, value,.na_rm = FALSE))%>%
  droplevels()

LPCs <- colnames(with_dup)[grepl("LPC_", colnames(with_dup))]
TGs <- colnames(with_dup)[grepl("TG_", colnames(with_dup))]

stat_summmary <- with_dup %>%
  dplyr::group_by(sample_group) %>%
  summarise(across(
    .cols =   in_norm, 
    .fns = list(Mean = function(x){mean(x, na.rm = TRUE)},
                SD = function(x){sd(x, na.rm = TRUE)},
                CV_pct = function(x){100*sd(x, na.rm = TRUE)/mean(x, na.rm = TRUE)}), 
    .names = "{col}.{fn}"
  )) %>%
  tidyr::pivot_longer(cols=-"sample_group")%>%
  tidyr::pivot_wider(names_from="sample_group") %>%
  tidyr::separate(name,into = c("standard","metric"),sep = "\\.") %>%
  tidyr::pivot_longer(cols=-c(1,2),names_to ="sample_group") %>%
  group_by(sample_group,metric) %>%
  mutate(standard=forcats::fct_reorder(standard, value,.na_rm = FALSE))%>%
  droplevels()

quartz()
stat_summmary  %>%
  ggplot(aes(y=value, x=standard,dodge=sample_group, fill=sample_group)) +
  ggtitle("MS-DIAL raw, Internal Standard QC,\nPSU Inv4m Phosphorus experiment") +
  geom_col(position=position_dodge()) +
  facet_grid(metric~standard, scales= "free")

quartz()
with_dup %>%
  dplyr::select(sample_group,all_of(TGs)) %>%
  tidyr::pivot_longer(cols=-"sample_group",names_to = "standard") %>%
  dplyr::mutate(standard=factor(standard,levels = levels(stat_summmary$standard))) %>%
  ggplot(aes(y=value+1, x=standard,dodge=sample_group, color=sample_group)) +
  ggtitle("MS-DIAL raw, Internal Standard QC,\nPSU Inv4m Phosphorus experiment") +
  ggbeeswarm::geom_quasirandom(dodge.width=4, size=1) +
  facet_grid(.~standard, scales= "free") +
  scale_y_log10()

quartz()
with_dup %>%
  dplyr::select(sample_group,all_of(in_norm)) %>%
  tidyr::pivot_longer(cols=-"sample_group",names_to = "standard") %>%
  dplyr::mutate(standard=factor(standard,levels = levels(stat_summmary$standard))) %>%
  ggplot(aes(y=value+1, x=standard,dodge=sample_group, color=sample_group)) +
  ggtitle("MS-DIAL raw, Internal Standard QC,\nPSU Inv4m Phosphorus experiment") +
  ggbeeswarm::geom_quasirandom(dodge.width=4, size=1) +
  facet_grid(.~standard, scales= "free") +
  scale_y_log10()

quartz()
hist(with_dup$LPC_17_0, breaks=100) 
LPC_17_0_thresh <- 3.25E8
abline(v= LPC_17_0_thresh, col ="red")
abline(v= 2.0E8, col ="red")
abline(v= 4.5E8, col ="red")


raw$LPC17_0_group <- cut(
  raw$LPC_17_0,
  breaks= c(0,4.5E8,1.0e10),
  labels = c("<4.5e8",">4.5e8"))
levels(raw$LPC17_0_group )
raw$LPC17_0_group
colnames(raw)    

to_check <- raw %>%
  dplyr::select(sampleID, tube, LPC17_0_group) %>%
  filter(LPC17_0_group == ">4.5e8")  %>% pull(sampleID)
cat(to_check)
# MDS with raw data -----

# add sample_info

plant_csv <- "../../data/22_NCS_PSU_LANGEBIO_FIELDS_PSU_P_field.csv"


psu <- read.csv(plant_csv) %>%
  rename(Genotype = "Who.What", row = "P22." ) %>%
  # dplyr::filter(Genotype %in% c("B73", "CTRL","INV4M",  "NUE")) %>%
  dplyr::filter(Genotype %in% c("CTRL","INV4M"))  %>%
  droplevels()
metadata <- read.csv('data/PSU-PHO22_Metadata.csv')

metadata$tube <- gsub("L|R","S",metadata$tube, perl =TRUE)
colnames(metadata)


sampleInfo <- psu  %>% 
  dplyr::select(row,Rep,Plot_Row, Plot_Column,DTA,DTS) %>%
  dplyr::inner_join(sampleInfo) %>% filter(!is.na(DTA)) %>%
  right_join(
    with_dup %>% 
      dplyr::select(tube,sampleID,sample_group,
                    starts_with("LPC_"), 
                    starts_with("TG_"))) %>%
  mutate(sample_label = gsub("230113_","",sampleID)) 

nrow(sampleInfo)
colnames(sampleInfo)

sampleInfo$leaf_group <-NA

sampleInfo$leaf_group[sampleInfo$leaf_tissue <3] <- "apical"
sampleInfo$leaf_group[sampleInfo$leaf_tissue >2] <- "basal"




library(edgeR)
library(limma)

cols <- with_dup%>%
  dplyr::select(CUDA:TG_58_5) %>% colnames()
length(cols)

counts <-  with_dup%>%
  dplyr::select(all_of(cols)) %>%
  as.matrix() %>% 
  t() 
nrow(counts )
ncol(counts)
nrow(with_dup)
quartz()
hist(counts)

dimnames(counts)
colnames(counts) <-with_dup$sampleID

colnames(counts)
genes = data.frame(gene= cols)

sampleNames = paste(with_dup$sampleID,with_dup$sample_group)

nrow(sampleInfo)

y = DGEList(counts = counts,samples = sampleInfo)
y$group =y$samples$sample_group

quartz()
mds = plotMDS(main ="Lipid MDS", y)


d = y$samples
d$x = mds$x
d$y = mds$y

quartz()
ggplot(d,aes(x=x,y=y)) + geom_point(aes(color = sample_group))

quartz()
ggplot(d,aes(x=x,y=y)) + geom_point(aes(color = factor(leaf_tissue),shape = Treatment),size=3) 


colnames(d)

d %>% dplyr::select(sampleID,sample_group,lib.size)

#distance <- sqrt((d$x+4)^2+(d$y+1)^2)
d$Treatment <- as.factor(d$Treatment) %>% as.numeric
d$Genotype <- as.factor(d$Genotype)   %>% as.numeric


correlates <-cor(d$x ,d %>% dplyr::select_if(is.numeric),
    use="pairwise.complete.obs")  
nc <- colnames(correlates)

correlates <- data.frame(id=nc, cor=correlates %>% as.numeric())

correlates  %>%
  mutate_if(is.numeric, round, 2) %>%
  filter(id !="x") %>% arrange(-cor) %>% tibble()  %>% print(n=50)

correlates <-cor(d$y ,d %>% dplyr::select_if(is.numeric),
                 use="pairwise.complete.obs")  
nc <- colnames(correlates)

correlates <- data.frame(id=nc, cor=correlates %>% as.numeric())

correlates  %>%
  mutate_if(is.numeric, round, 2) %>%
  filter(id !="y") %>% 
  arrange(-cor) %>% tibble() %>%
  filter(abs(cor) >0.1) %>% print(n=50)

quartz(height=7,width =7)
d %>%
  filter(x>-5 & y> -5) %>%
ggplot(aes(x=x,y=y, label=sample_label,color = lib.size)) +
  ggtitle( ggtitle("Multidimensional scaling\nMS-DIAL raw, Internal Standard QC,\nPSU Inv4m Phosphorus experiment") )+
  geom_point()+
  # xlim(-6.5,6.5)+
  # ylim(-1.25,1.25)+
  ggrepel::geom_label_repel(max.overlaps = 25,force_pull = 5) +
  labs(color="Total\nLipid\nSignal")+
 scale_color_viridis_c()
d$row

quartz(height=7,width =7)
d %>%
  filter(x>-5 & y> -5) %>%
  ggplot(aes(x=x,y=y, label=sample_label,color = row)) +
  ggtitle( ggtitle("Multidimensional scaling\nMS-DIAL raw, Internal Standard QC,\nPSU Inv4m Phosphorus experiment") )+
  geom_point()+
  ggrepel::geom_label_repel(max.overlaps = 25,force_pull = 5) +
  scale_color_viridis_c()


quartz(height=7,width =7)
d %>%
  filter(x>-5 & y> -5) %>%
  ggplot(aes(x=x,y=y, label=sample_label,color = injection_order)) +
  ggtitle( ggtitle("Multidimensional scaling\nMS-DIAL raw, Internal Standard QC,\nPSU Inv4m Phosphorus experiment") )+
  geom_point()+
  ggrepel::geom_label_repel(max.overlaps = 25,force_pull = 5) 


quartz()
d %>%
  filter(x>-5 & y> -5) %>%
  filter(!is.na(row)) %>%
  droplevels() %>%
  ggplot(aes(y=x,x=row, label=sample_label,color = row)) +
  ggtitle( ggtitle("Multidimensional scaling\nMS-DIAL raw, Internal Standard QC,\nPSU Inv4m Phosphorus experiment") )+
  facet_grid(Treatment~COLLECTOR,scales ="free") +
  ylab("Dim1")  +
  geom_point()


quartz()
d %>%
  filter(x>-5 & y> -5) %>%
  filter(!is.na(row)) %>%
  droplevels() %>%
  ggplot(aes(y=y,x=row, label=sample_label,color = row)) +
  ggtitle( ggtitle("Multidimensional scaling\nMS-DIAL raw, Internal Standard QC,\nPSU Inv4m Phosphorus experiment") )+
  #facet_grid(Treatment~COLLECTOR,scales ="free") +
  ylab("Dim2")  +
  geom_point()


d %>%
  filter(row==3095)

d %>% dplyr::select_if(is.numeric)

quartz()
hist(d$lib.size, breaks=40)


quartz(height=7,width =7)
d%>%
  # filter(sample_label !="Pool_6") %>%
  ggplot(aes(x=x,y=y, label=sample_label,color = sample_group)) +
  ggtitle( ggtitle("Multidimensional scaling\nMS-DIAL raw, Internal Standard QC,\nPSU Inv4m Phosphorus experiment") )+
  geom_point()+
  ggrepel::geom_label_repel(max.overlaps = 25,force_pull = 5)



d$LPC_17_0

d$LPC17_0_group
quartz(height=7,width =7)
d%>%
 # filter(sample_label !="Pool_6") %>%
  ggplot(aes(x=x,y=y, label=sample_label,color = LPC17_0_group)) +
  ggtitle( ggtitle("Multidimensional scaling\nMS-DIAL raw, Internal Standard QC,\nPSU Inv4m Phosphorus experiment") )+
  geom_point()+
  ggrepel::geom_label_repel(max.overlaps = 25,force_pull = 5) +
  scale_color_viridis_d()

quartz(height=7,width =7)
d%>%
  # filter(sample_label !="Pool_6") %>%
  ggplot(aes(x=x,y=y, label=sample_label,color = Treatment )) +
  ggtitle( ggtitle("Multidimensional scaling\nMS-DIAL raw, Internal Standard QC,\nPSU Inv4m Phosphorus experiment") )+
  geom_point()+
  ggrepel::geom_label_repel(max.overlaps = 25,force_pull = 5) +
  scale_color_viridis_d()

correlates %>% arrange(-abs(cor))
quartz(height=7,width =7)
d%>%
  # filter(sample_label !="Pool_6") %>%
  ggplot(aes(x=x,y=y, label=sample_label,color = rna_batch1)) +
  ggtitle( ggtitle("Multidimensional scaling\nMS-DIAL raw, Internal Standard QC,\nPSU Inv4m Phosphorus experiment") )+
  geom_point()+
  ggrepel::geom_label_repel(max.overlaps = 25,force_pull = 5) +
  scale_color_viridis_c()

quartz(height=7,width =7)
d%>%
  # filter(sample_label !="Pool_6") %>%
  ggplot(aes(x=x,y=y, label=sample_label,color = Plot_Column )) +
  ggtitle( ggtitle("Multidimensional scaling\nMS-DIAL raw, Internal Standard QC,\nPSU Inv4m Phosphorus experiment") )+
  geom_point()+
  ggrepel::geom_label_repel(max.overlaps = 25,force_pull = 5) +
  scale_color_viridis_c()

quartz(height=7,width =7)
d%>%
  # filter(sample_label !="Pool_6") %>%
  ggplot(aes(x=x,y=y, label=sample_label,color = Treatment )) +
  ggtitle( ggtitle("Multidimensional scaling\nMS-DIAL raw, Internal Standard QC,\nPSU Inv4m Phosphorus experiment") )+
  geom_point()+
  ggrepel::geom_label_repel(max.overlaps = 25,force_pull = 5) +
  scale_color_viridis_d()

quartz(height=7,width =7)
d%>%
  # filter(sample_label !="Pool_6") %>%
  ggplot(aes(x=x,y=y, label=sample_label,color = Genotype )) +
  ggtitle( ggtitle("Multidimensional scaling\nMS-DIAL raw, Internal Standard QC,\nPSU Inv4m Phosphorus experiment") )+
  geom_point()+
  ggrepel::geom_label_repel(max.overlaps = 25,force_pull = 5) +
  scale_color_viridis_d()

quartz(height=7,width =7)
d%>%
  # filter(sample_label !="Pool_6") %>%
  ggplot(aes(x=x,y=y, label=sample_label,color = decimal_time )) +
  ggtitle( ggtitle("Multidimensional scaling\nMS-DIAL raw, Internal Standard QC,\nPSU Inv4m Phosphorus experiment") )+
  geom_point()+
  ggrepel::geom_label_repel(max.overlaps = 25,force_pull = 5) +
  scale_color_viridis_c()

with_dup

# rescaled---- 

# factor_per_spp <- apply(m,2,mean)
# length(factor_per_spp)
# with_dup$LPC17_0_group <- d$LPC17_0_group
# g1<- m[as.numeric(d$LPC17_0_group)==1,]
# dim(g1)
# 
# sc1 <- apply(g1,2,scale)
# sc1[is.na(sc1)] <- 0
# resc1<- apply(sc1,2,function(x){x - min(x)}) %*% diag(factor_per_spp) +1
# colnames(resc1) <-colnames(m)
# rownames(resc1) <- d$sampleID[as.numeric(d$LPC17_0_group)==1]
# 
# quartz()
# hist(resc1)
# 
# g2<- m[as.numeric(d$LPC17_0_group)==3,]
# 
# sc2 <- apply(g2,2,scale)  
# 
# sc2[is.na(sc2)] <- 0
# resc2<- apply(sc2,2,function(x){x - min(x)})%*%diag(max_per_spp) +1
# 
# colnames(resc2) <-colnames(m)
# rownames(resc2) <- d$sampleID[as.numeric(d$LPC17_0_group)==3]
# 
# m2 <- rbind(resc1,resc2)[d$sampleID,]


# 
# counts <- m2 %>%
#   as.matrix() %>% 
#   t() 
# nrow(counts )
# ncol(counts)
# quartz()
# hist(counts)
# 
# 
# genes = data.frame(gene= rownames(counts))
# 
# sampleNames = paste(with_dup$sampleID,with_dup$sample_group)
# 
# 
# y = DGEList(counts = counts,samples = sampleInfo)
# y$group =y$samples$sample_group
# 
# quartz()
# mds = plotMDS(main ="Lipid MDS", y)
# 
# 
# d = y$samples
# d$x = mds$x
# d$y = mds$y
# 
# quartz()
# ggplot(d,aes(x=x,y=y)) + geom_point(aes(color = sample_group))
