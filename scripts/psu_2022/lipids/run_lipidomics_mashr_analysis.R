library(edgeR)
library(ggplot2)
library(limma)
library(ashr)
library(mashr)

# counts = fread('../data/inv4mRNAseq_gene_sample_exp.csv',data.table=F)
#sampleInfo = fread('../data/inv4mRNAseq_metadata.csv',data.table=F) 
# counts = fread('inv4mRNAseq_gene_sample_exp.csv',data.table=F)


pal <- c("chartreuse2", "darkblue")

pal2 <- c("darkorchid","chartreuse2", "darkblue")



# filter internal standards
to_exclude<- c("DG_12_0", "LPC_17_0", "LPE_17_1",
               "PC_25_0", "PG_17_0","SM_35_1",
               "Sphingosine_17_1","TG_17_0",
               "TG_57_6" #exclude the odd number lipid
)

# lipid_csv <- "data/PSU_Normalized_MSDial_240416.csv"
# lipid_csv <- "data/PSU_RawData_MSDial_240416.csv"

lipid_csv <- "data/PSU_RawData_MSDial_NewStdInt_240422.csv"

raw <- read.csv(lipid_csv,na.strings = c("","#N/A","NA","Inf"))
to_exclude  <- to_exclude[to_exclude %in% colnames(raw)]
                          
raw <- raw %>%  dplyr::select(-all_of(to_exclude )) %>%
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
levels(pheno$Treatment) <- c("+P","-P") 

m <- pheno %>%
  dplyr::select(DGDG_34_0:TG_58_5) %>% as.matrix()

# analysis by GROUP

head_group <- c(
  DG ="neutral", 
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
  TG="neutral")

sp <- data.frame(
  colname = colnames(m),
  class = gsub("_.*","", colnames(m), perl =TRUE)
)

sp$head_group <-head_group[sp$class]

with(sp,
     table(head_group,class))
#group <- colnames(m) 
#group <- gsub("_.*","", group, perl =TRUE) %>% sort() %>% unique()
#n <- length(group)
#pairs <- combn(n,2) %>% t
# ratios <- paste(group[pairs[,1]],group[pairs[,2]], sep = "_")
# 
# for (g in group){
#   pheno[[g]] <- pheno  %>%
#     dplyr::select(dplyr::starts_with(paste0(g,"_"))) %>%
#     rowSums(na.rm = TRUE)
# }
# 
# 
# for (r in ratios){
#   pheno[r] <- pheno[gsub("_.*","",r,perl=TRUE)]/pheno[gsub(".*_","",r,perl=TRUE)]
# }
# ncol(pheno)

pheno <- pheno[-29,] # what's wrong with this column?

counts <- pheno[,colnames(m)] %>% as.matrix() %>% t() 

# counts <- pheno[,colnames(m)] %>% as.matrix() %>% t() * 4e7
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


y$samples$lowCount = y$samples$lib.size < 1
summary(y$samples$lib.size)


mds = plotMDS(y,pch=21,label = y$samples$side_tag,)

quartz()
plot(mds$x,mds$y,pch=21, main = "Lipid MDS")
# keep = mds$y > -1


y_filtered = y

nrow(y_filtered)


y_filtered_bySample = y_filtered[,!y_filtered$samples$lowCount]

y_filtered_bySample$samples
table(y_filtered_bySample$samples$Treatment,
      y_filtered_bySample$samples$leaf_tissue)
table(y_filtered_bySample$samples$Genotype,
      y_filtered_bySample$samples$leaf_tissue)

table(y_filtered_bySample$samples$Treatment,
      y_filtered_bySample$samples$Genotype,
      y_filtered_bySample$samples$leaf_tissue)

quartz()
mds2 = plotMDS(y_filtered_bySample,pch=21,
               label = y_filtered_bySample$samples$side_tag,bg = (as.factor(y_filtered_bySample$samples$lowCount)==T)+1)


d = y_filtered$samples
d$x = mds$x
d$y = mds$y

quartz()
ggplot(d,aes(x=x,y=y)) + geom_point(aes(color = Treatment))

quartz()
ggplot(d,aes(x=x,y=y)) + geom_point(aes(color = decimal_time ))



quartz()
ggplot(d,aes(x=x,y=y)) + geom_point(aes(color = COLLECTOR))

quartz()
ggplot(d,aes(x=x,y=y)) + geom_point(
  aes(color = as.factor(rowid),
      shape = Treatment)) 

quartz()
ggplot(d,aes(x=x,y=y)) + geom_point(aes(color = Genotype))

d$Treatment <- factor(d$Treatment)
d$Rep <- factor(d$Rep)

quartz()
ggplot(d,aes(x=x,y=y)) + geom_point(aes(color = factor(leaf_tissue),shape = Treatment),size=3) 


plot(mds2$x,mds2$y,pch=21,bg = factor(y_filtered_bySample$samples$Treatment),col=0)
plot(mds2$x,mds2$y,col = factor(y_filtered_bySample$samples$Genotype))

sampleInfo$Plot_Row
d$Treatment <- factor(d$Treatment,levels = c("High_P","Low_P"))
y_filtered_bySample <-y
design = model.matrix(~ Plot_Row + Plot_Column + Rep + leaf_tissue+ Treatment*Genotype,d)
y_filtered_bySample = calcNormFactors(y_filtered_bySample)
quartz()

voomR = voom(y_filtered_bySample,design=design,plot=T)

fit = lmFit(voomR)
ebfit = eBayes(fit)


head(ebfit$coefficients)
# 
# 
# design_interaction = model.matrix(~ Rep + leaf_tissue*Treatment*Genotype,d)
# fit = lmFit(voomR,design)
# ebfit = eBayes(fit)

nrow(ebfit$coefficients)
ebfit$coefficients
r = list()
for(x in 8:11){
r[[as.character(x)]] <- topTable(ebfit,coef = x,sort.by = 'none',n=Inf)
cr <- qt(0.975, ebfit$df.residual + ebfit$df.prior) * ebfit$stdev.unscaled[,x] * sqrt(ebfit$s2.post)
# Calculating  confidence interval
r[[as.character(x)]][["upper"]]  <- r[[as.character(x)]][["logFC"]] +cr
r[[as.character(x)]][["lower"]]  <- r[[as.character(x)]][["logFC"]] -cr
}
r$`8`[order(r$`8`$adj.P.Val),]
r$`9`[order(r$`9`$adj.P.Val),]
r$`10`[order(r$`10`$adj.P.Val),]
r$`11`[order(r$`11`$adj.P.Val),]


r$`8`["predictor"] <- "leaf"
r$`8`["Response"]   <- rownames(counts)
r$`9`["predictor"] <- "-P"
r$`9`["Response"]   <- rownames(counts)
r$`10`["predictor"] <- "Inv4m"
r$`10`["Response"]   <- rownames(counts)
r$`11`["predictor"] <- "-P:Inv4m"
r$`11`["Response"]   <- rownames(counts)

effects <- rbind(r$`8`,r$`9`,r$`10`,r$`11`) 

quartz()
effects %>%
  filter(predictor=="leaf") %>%
  ggplot(aes(y=-log10(adj.P.Val),x=logFC))+
  geom_point()

quartz()
effects %>%
  filter(predictor=="-P") %>%
  ggplot(aes(y=-log10(adj.P.Val),x=logFC))+
  geom_point()
  

quartz()
effects %>%
  filter(predictor=="Inv4m") %>%
  ggplot(aes(y=-log10(adj.P.Val),x=logFC))+
  geom_point()

quartz()
effects %>%
  filter(predictor=="-P:Inv4m") %>%
  ggplot(aes(y=-log10(adj.P.Val),x=logFC))+
  geom_point()

effect_order <- c("leaf","-P","Inv4m","-P:Inv4m")

effects <- effects %>%
  mutate(predictor = factor(predictor, levels=effect_order))%>%
  mutate(neglogP = -log10(adj.P.Val)) %>%
  mutate(is_significant = adj.P.Val<0.05) %>%
  mutate(upregulated = case_when(
    (predictor !="leaf" ) &   (logFC >= 2)    & is_significant ~ TRUE,
    (predictor =="leaf" ) &   (logFC >= 0.7)  & is_significant ~ TRUE,
    .default = FALSE)) %>%
  mutate(downregulated = case_when(
    (predictor !="leaf" ) & (logFC <= -2)   & is_significant ~ TRUE,
    (predictor =="leaf" ) & (logFC <= -0.7) & is_significant ~ TRUE,
    .default = FALSE))  %>%
  mutate(regulation = case_when(
    is_significant & upregulated ~ "Upregulated",
    is_significant & downregulated ~ "Downregulated",
    .default = "Unregulated"),
    is_DEG= is_significant & regulation != "Unregulated",
    label = ifelse(is_DEG, Response, ""))

effects$label <- sub("_","",effects$label, perl =TRUE)
effects$label <- sub("_",":",effects$label, perl =TRUE)

draw_key_custom <- function (data, params, size) {
  colour <- data$colour 
  data <- data[data$colour == colour,]
  grid::segmentsGrob(0, 1/2, 1, 1/2,
                     gp = grid::gpar(col = data$colour, 
                                     lwd = data$linewidth * 2))
}


quartz( height=6, width =12)
effects %>%
  filter(predictor !="-P:Inv4m") %>%
  droplevels() %>%
  inner_join(sp, by=c(Response="colname")) %>%
  ggplot(aes(x=logFC, y=-log10(adj.P.Val),color=head_group, fill= regulation, label= label))+
  ylab(expression(-log[10](italic(FDR)))) +
  xlab(expression(-log[2]("Fold Change"))) +
  ylim(0,10) +
  geom_point(alpha=0.7, shape=21, color="grey25") +
  scale_fill_manual( name="", values = c("#00AFBB", "grey","#bb0c00"), # to set the colours of our variable
                      labels = c("Downregulated","Not significant","Upregulated")) +
  ggrepel::geom_text_repel(force=3, ,
                           #fontface="italic",
                           min.segment.length = 0,
                           max.overlaps = 20, 
                           key_glyph = draw_key_custom) +
  scale_color_manual( name="", values = c("darkblue", "orange2","darkred"), # to set the colours of our variable
                      labels = c("Glycolipid","Neutral lipid","Phospholipid"))+
  # stat_ellipse(col= "grey75",level = 0.99)+
  facet_wrap(.~predictor,scales = "free_x", ncol=3) +
  guides( color = guide_legend(override.aes = list(linewidth=1.5),),
         fill="none")+
  ggpubr::theme_classic2(base_size = 25) +
  theme(plot.title = element_blank(),
        strip.text = ggtext::element_markdown(),
        legend.position = c(0.5,0.95),
        legend.spacing =  unit(0,"line"),
        legend.box.spacing =unit(0, "line"),
        legend.text = element_text(size=15),
        legend.direction = "horizontal",
        strip.background = element_blank())







  
  
# plot effects ----
to_plot <- effects %>%
  mutate(Response= gsub("_glyco","",Response)) %>%
  mutate(predictor=factor(predictor, levels= c("leaf","-P","Inv4m","-P:Inv4m"))) %>%
  filter(P.Value< 0.05) %>%  inner_join(sp, by=c(Response="colname"))   %>%
  mutate(sign= sign(logFC)) %>%
  ungroup() %>%
  group_by(Response) %>%
  mutate(max_effect = max(abs(logFC))) %>%
  ungroup() %>%
  group_by(predictor,head_group,class) %>%   
  mutate(head_group = forcats::fct_reorder(head_group,adj.P.Val)) %>%
  ungroup() %>%
  group_by(predictor) %>% # As a precaution / handle in a separate .grouped_df method
  arrange(predictor,head_group,class,adj.P.Val) %>%   # arrange by facet variables and continuous values
  mutate(.r = row_number()) %>% # Add a row number variable %>%
  ungroup() 

to_plot %>% as.data.frame() %>%
  filter(Response=="LPC_18_3")



pd = position_dodge(0.4)

colnames(to_plot)
quartz()
to_plot %>%
  droplevels() %>%
  mutate(Response = forcats::fct_reorder(Response,.r)) %>%
  ungroup() %>%
  ggplot( aes(x     = logFC,
              y     = Response,
              col = head_group)) +
  xlab("Effect (log2 Fold Change)") +
  geom_vline(xintercept = 0,lty =2)+
  geom_point(position= pd,
             size   = 3) +
  geom_errorbar(aes(xmin  =  upper,
                    xmax  =  lower),
                position= pd,
                width =  0.2,
                size  =  0.7) +
  guides(shape = guide_legend(
    title = NULL), color = guide_legend(
      reverse = TRUE)) +
  facet_wrap(.~predictor, scales="free", ncol =4)  +
  scale_y_discrete(limits=rev) +
  scale_color_manual(values = c("blue","red","chartreuse"))+
  # scale_shape_manual(values=c(13,5,21),
  #                   labels= c( expression(italic("Inv4m") %*% P), "leaf","P")
  #                   ) +
  theme_light(base_size =15) +
  theme(legend.position ="top",
        strip.background =element_rect(fill="white"),
        strip.text =  element_text(color = "black",size=15),
        axis.title.y=element_blank(),
        axis.text.y = element_text(hjust = 0, face = "bold"),
        plot.caption = element_text(hjust = 0))



pd = position_dodge(0.4)

colnames(to_plot)
quartz()
response_p <-  to_plot %>%
  dplyr::filter(predictor=="-P") %>%
  droplevels() %>%
  mutate(Response = forcats::fct_reorder(Response,.r)) %>%
  ungroup() %>%
  ggplot( aes(x     = logFC,
              y     = Response,
              col = head_group)) +
  xlab("Effect (log2 Fold Change)") +
  geom_vline(xintercept = 0,lty =2)+
  geom_point(position= pd,
             size   = 3) +
  geom_errorbar(aes(xmin  =  upper,
                    xmax  =  lower),
                position= pd,
                width =  0.2,
                size  =  0.7) +
  guides(shape = guide_legend(
    title = NULL), color = guide_legend(
      reverse = TRUE)) +
  facet_wrap(.~predictor, scales="free")  +
  scale_y_discrete(limits=rev) +
  scale_color_manual(values = c("blue","red","chartreuse"))+
  # scale_shape_manual(values=c(13,5,21),
  #                   labels= c( expression(italic("Inv4m") %*% P), "leaf","P")
  #                   ) +
  theme_light(base_size =15) +
  theme(legend.position ="top",
        strip.background =element_rect(fill="white"),
        strip.text =  element_text(color = "black",size=15),
        axis.title.y=element_blank(),
        axis.text.y = element_text(hjust = 0, face = "bold"),
        plot.caption = element_text(hjust = 0))


response_leaf <-  to_plot %>%
  dplyr::filter(predictor=="leaf") %>%
  droplevels() %>%
  mutate(Response = forcats::fct_reorder(Response,.r)) %>%
  ungroup() %>%
  ggplot( aes(x     = logFC,
              y     = Response,
              col = head_group)) +
  xlab("Effect (log2 Fold Change)") +
  geom_vline(xintercept = 0,lty =2)+
  geom_point(position= pd,
             size   = 3) +
  geom_errorbar(aes(xmin  =  upper,
                    xmax  =  lower),
                position= pd,
                width =  0.2,
                size  =  0.7) +
  guides(shape = guide_legend(
    title = NULL), color = guide_legend(
      reverse = TRUE)) +
  facet_wrap(.~predictor, scales="free")  +
  scale_y_discrete(limits=rev) +
  scale_color_manual(values = c("blue","red","chartreuse"))+
  # scale_shape_manual(values=c(13,5,21),
  #                   labels= c( expression(italic("Inv4m") %*% P), "leaf","P")
  #                   ) +
  theme_light(base_size =15) +
  theme(legend.position ="top",
        strip.background =element_rect(fill="white"),
        strip.text =  element_text(color = "black",size=15),
        axis.title.y=element_blank(),
        axis.text.y = element_text(hjust = 0, face = "bold"),
        plot.caption = element_text(hjust = 0))


quartz()
ggpubr::ggarrange(response_leaf, response_p, ncol=2,common.legend = TRUE)

# Just p reponsive genes ----

sp_up <- to_plot %>%
  filter( predictor =="-P") %>%
  filter( sign =="1") %>% pull(Response)

p_down <- to_plot %>%
  filter( predictor =="-P") %>%
  filter( sign =="-1") %>% pull(Response)

p_up <- sp$class[sp$colname %in% sp_up] %>% sort %>% unique()


ratio <- expand.grid(p_down,p_up)
ratio$Var1 <- as.character(ratio$Var1)
ratio$Var2 <- as.character(ratio$Var2)
ratios <- paste(ratio$Var1,ratio$Var2, sep="X")

for (i in 1:nrow(ratio)){
   r <- paste(ratio[i,],collapse = 'X')
   pheno[r] <- (pheno[,ratio[i,1]]+0.5)/(pheno[,ratio[i,2]]+0.5)
 }
colnames(pheno)

counts <- pheno[,ratios] %>% as.matrix() %>% t()
#counts[is.na(counts)] <-0
colnames(counts) <- pheno$tube



# run the limma model model from above 


to_plot <- effects %>%
  mutate(predictor=factor(predictor, levels= c("leaf","-P","Inv4m","-P:Inv4m"))) %>%
  filter(adj.P.Val< 0.05) %>%
  filter(predictor=="-P:Inv4m") %>%
  # inner_join(sp, by=c(Response="colname"))   %>%
  mutate(sign= sign(logFC)) %>%
  ungroup() %>%
  group_by(Response) %>%
  mutate(max_effect = max(abs(logFC))) %>%
  ungroup() %>%
  # group_by(predictor,head_group,class) %>%   
  #mutate(head_group = forcats::fct_reorder(head_group,adj.P.Val)) %>%
  # ungroup() %>%
  group_by(predictor) %>% # As a precaution / handle in a separate .grouped_df method
  arrange(predictor,adj.P.Val) %>%
  #arrange(predictor,head_group,class,adj.P.Val) %>%   # arrange by facet variables and continuous values
  mutate(.r = row_number()) %>% # Add a row number variable %>%
  ungroup() 

to_plot %>% as.data.frame() %>%
  filter(Response=="LPC_18_3")

pd = position_dodge(0.4)

colnames(to_plot)
quartz()
to_plot %>%
  dplyr::filter(predictor=="-P:Inv4m") %>%
  droplevels() %>%
  mutate(Response = forcats::fct_reorder(Response,.r)) %>%
  ungroup() %>%
  ggplot( aes(x     = logFC,
              y     = Response,
              #col = head_group
  )) +
  xlab("Effect (log2 Fold Change)") +
  geom_vline(xintercept = 0,lty =2)+
  geom_point(position= pd,
             size   = 3) +
  geom_errorbar(aes(xmin  =  upper,
                    xmax  =  lower),
                position= pd,
                width =  0.2,
                size  =  0.7) +
  guides(shape = guide_legend(
    title = NULL), color = guide_legend(
      reverse = TRUE)) +
  facet_wrap(.~predictor, scales="free")  +
  scale_y_discrete(limits=rev) +
  scale_color_manual(values = c("blue","red"))+
  # scale_shape_manual(values=c(13,5,21),
  #                   labels= c( expression(italic("Inv4m") %*% P), "leaf","P")
  #                   ) +
  theme_light(base_size =15) +
  theme(legend.position ="top",
        strip.background =element_rect(fill="white"),
        strip.text =  element_text(color = "black",size=15),
        axis.title.y=element_blank(),
        axis.text.y = element_text(hjust = 0, face = "bold"),
        plot.caption = element_text(hjust = 0))

# run the limma model model from above 


to_plot <- effects %>%
  mutate(predictor=factor(predictor, levels= c("leaf","-P","Inv4m","-P:Inv4m"))) %>%
  filter(adj.P.Val< 0.05) %>%
  filter(predictor=="Inv4m") %>%
  # inner_join(sp, by=c(Response="colname"))   %>%
  mutate(sign= sign(logFC)) %>%
  ungroup() %>%
  group_by(Response) %>%
  mutate(max_effect = max(abs(logFC))) %>%
  ungroup() %>%
  # group_by(predictor,head_group,class) %>%   
  #mutate(head_group = forcats::fct_reorder(head_group,adj.P.Val)) %>%
  # ungroup() %>%
  group_by(predictor) %>% # As a precaution / handle in a separate .grouped_df method
  arrange(predictor,adj.P.Val) %>%
  #arrange(predictor,head_group,class,adj.P.Val) %>%   # arrange by facet variables and continuous values
  mutate(.r = row_number()) %>% # Add a row number variable %>%
  ungroup() 


pd = position_dodge(0.4)

colnames(to_plot)
quartz()
to_plot %>%
  dplyr::filter(predictor=="Inv4m") %>%
  droplevels() %>%
  mutate(Response = forcats::fct_reorder(Response,.r)) %>%
  ungroup() %>%
  ggplot( aes(x     = logFC,
              y     = Response,
              #col = head_group
  )) +
  xlab("Effect (log2 Fold Change)") +
  geom_vline(xintercept = 0,lty =2)+
  geom_point(position= pd,
             size   = 3) +
  geom_errorbar(aes(xmin  =  upper,
                    xmax  =  lower),
                position= pd,
                width =  0.2,
                size  =  0.7) +
  guides(shape = guide_legend(
    title = NULL), color = guide_legend(
      reverse = TRUE)) +
  facet_wrap(.~predictor, scales="free")  +
  scale_y_discrete(limits=rev) +
  scale_color_manual(values = c("blue","red"))+
  # scale_shape_manual(values=c(13,5,21),
  #                   labels= c( expression(italic("Inv4m") %*% P), "leaf","P")
  #                   ) +
  theme_light(base_size =15) +
  theme(legend.position ="top",
        strip.background =element_rect(fill="white"),
        strip.text =  element_text(color = "black",size=15),
        axis.title.y=element_blank(),
        axis.text.y = element_text(hjust = 0, face = "bold"),
        plot.caption = element_text(hjust = 0))



#---
#gene='Zm00001eb053360'
gene='DGDG_LPE' # PILNCR1-miR399

d$y = voomR$E[gene,]
d$counts = y_filtered_bySample$counts[gene,]
# ggplot(d,aes(x=Treatment,y=y)) + geom_point(aes(color = Genotype,group = interaction(Treatment,Genotype)),position = position_jitterdodge()) + facet_wrap(~leaf_tissue)

quartz()
ggplot(d,aes(x=factor(Treatment, levels= c("Low_P","High_P")),y=y)) +
  ggtitle(gene) +
  geom_boxplot(aes(color = Genotype,group = interaction(Treatment,Genotype))) + 
  facet_wrap(~leaf_tissue) 

# ggplot(d,aes(x=Treatment,y=counts)) + geom_boxplot(aes(color = Genotype,group = interaction(Treatment,Genotype))) + facet_wrap(~leaf_tissue)
ggplot(d,aes(x=Genotype,y=y))  +
  ggtitle(gene) +
  geom_boxplot(aes(color = Treatment,group = interaction(Treatment,Genotype))) + 
  facet_wrap(~leaf_tissue)




########################################################
# by leaf analysis #####################################
########################################################

# make design in each = ~Treatment*Genotype ############
# check effects for interaction




# Akaike information content calculation
aic <- function(x){
  # x is a mashr model object
  # I am assuming x$fitted_g$pi, the number of mixture components,
  # is the number of parameters 
  2*length(x$fitted_g$pi) - 2*max(get_loglik(x))
}


## ------------------------------------------------------------------
# Step 4: Run mash for different models
# using just the simple canonical covariances as in the initial introductory vignette.
# correcting for correlations between samples inside a "condition"
# 



tb <- data.frame()
results <- list()

for(coef in 2:4){

  
  effects = SEs = matrix(NA,nrow = nrow(y_filtered_bySample),ncol = 4)
  rownames(effects) = rownames(SEs) = rownames(y_filtered_bySample)
  
  for(x in 1:4) {
    y_filtered_by_leaf = y_filtered_bySample[,y_filtered_bySample$samples$leaf_tissue==x]
    y_filtered_by_leaf$group = interaction(y_filtered_by_leaf$samples$Genotype,y_filtered_by_leaf$samples$Treatment)
    # filter_genes = filterByExpr(y_filtered_by_leaf,group = y_filtered_by_leaf$group)
    # y_filtered_by_leaf = y_filtered_by_leaf[filter_genes,]
    d = y_filtered_by_leaf$samples
    design_interaction = model.matrix(~Treatment*Genotype,d)
    y_filtered_by_leaf = calcNormFactors(y_filtered_by_leaf)
    voomR = voom(y_filtered_by_leaf ,design=design_interaction,plot=T)
    fit = lmFit(voomR)
    ebfit = eBayes(fit)
    tt <- topTable(ebfit,coef = coef,sort.by = 'none',n=Inf)
    effects[,x] = tt$logFC
    SEs[,x] = tt$logFC/tt$t
  }
  

  
  # Step 2: Obtain initial data-driven covariance matrices
  
  data = mash_set_data(as.matrix(effects), as.matrix(SEs))
  U.c = cov_canonical(data)
  
  # Step 2: Obtain initial data-driven covariance matrices
  
  # select strong signals
  m.1by1 = mash_1by1(data)
  strong = get_significant_results(m.1by1,0.05)
  if(length(strong)<20) {
    strong = order(apply(m.1by1$result$lfsr,1,min))[1:min(20,nrow(effects))]
  }
  
  # Perform PCA on data and return list of candidate covariance matrices
  U.pca = cov_pca(data,npc=ncol(effects),subset=strong)
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
  
  
  
  # Why U.c and U.ed?  They are supposed to be different ways of estimating the conditions covariance matrix
  V.em_c_ed = mash_estimate_corr_em(data, c(U.c,U.ed), details = TRUE)
  m.Vem_c_ed = V.em_c_ed$mash.model
  m.Vem_c_ed$result$NAs = is.na(effects)
  m.Vem_c_ed$V = V.em_c_ed$V
  
  print(get_loglik(m.Vem_c_ed),digits=10) 
  length(m.Vem_c_ed$fitted_g$pi)
  aic(m.Vem_c_ed)
  
  results[[as.character(coef)]][["m.Vem_c_ed"]] <- m.Vem_c_ed
  
  
  ###
  V.em_ed = mash_estimate_corr_em(data, U.ed, details = TRUE)
  m.Vem_ed = V.em_ed$mash.model
  m.Vem_ed$result$NAs = is.na(effects)
  
  quartz()
  mash_plot_meta(m.Vem_ed,2)
  
  m.Vem_ed$V = V.em_ed$V
  
  print(get_loglik(m.Vem_ed),digits=10) 
  length(m.Vem_ed$fitted_g$pi)
  aic(m.Vem_ed)
  
  results[[as.character(coef)]][["m.Vem_ed"]] <- m.Vem_ed
  
  ###
  V.em_c = mash_estimate_corr_em(data, U.c, details = TRUE)
  m.Vem_c = V.em_c$mash.model
  m.Vem_c$result$NAs = is.na(effects)
  m.Vem_c$V = V.em_c$V
  
  print(get_loglik(m.Vem_c),digits=10) 
  length(m.Vem_c$fitted_g$pi)
  aic(m.Vem_c)
  
  results[[as.character(coef)]][["m.Vem_c"]] <- m.Vem_c
  ###
  m.c = mash(data,U.c)
  # Fitting model with 217 mixture components
  # -57582
  
  print(get_loglik(m.c),digits=10)
  length(m.c$fitted_g$pi)
  aic(m.c)
  
  
  ###
  m.c_ed  = mash(data, c(U.c,U.ed))
  
  # Fitting model with 337 mixture components.
  
  print(get_loglik(m.c_ed),digits=10)
  length(m.c_ed$fitted_g$pi)
  aic(m.c_ed)
  
  results[[as.character(coef)]][["m.c_ed"]] <- m.c_ed
  ####
  # Compare likelihood of the  models
  
  ## Compare model fit
  
  
  model =  c('canonical',
             'canonical + ed',
             'canonical + Vem',
             'ed + Vem', 
             'canonical + ed + Vem')
  
  significant = c(length(get_significant_results(m.c)), 
                  length(get_significant_results(m.c_ed)),
                  length(get_significant_results(m.Vem_c)),
                  length(get_significant_results(m.Vem_ed)),
                  length(get_significant_results(m.Vem_c_ed)))
  
  loglike = c(
    max(get_loglik(m.c)),
    max(get_loglik(m.c_ed)),
    max(get_loglik(m.Vem_c)),
    max(get_loglik(m.Vem_ed)),
    max(get_loglik(m.Vem_c_ed))
  )
  
  k = c(length(m.c$fitted_g$pi), 
        length(m.c_ed$fitted_g$pi), 
        length(m.Vem_c$fitted_g$pi),
        length(m.Vem_ed$fitted_g$pi), 
        length(m.Vem_c_ed$fitted_g$pi))
  
  AIC = c(aic(m.c), aic(m.c_ed), aic(m.Vem_c), aic(m.Vem_ed), aic(m.Vem_c_ed))
  
  tb <- rbind ( tb,
    data.frame(
      coef = coef,
      model = model,
      significant = significant,
      loglike = loglike,
      k = k,
      AIC = AIC)
    )
}

# write.csv(tb,"by_leaf_mashr_tb.csv")
# saveRDS(results, "by_leaf_results.RDS")

coef_names <- colnames(head(ebfit$coefficients))
tb$coef_name <- coef_names[tb$coef]

library(dplyr)
tb %>%
  dplyr::select(coef, coef_name, everything())

tb %>%
  dplyr::select(coef, coef_name, everything()) %>%
  dplyr::group_by(coef) %>%
  dplyr::arrange(AIC) %>%
  dplyr::slice(1)

hist(get_log10bf(results$`4`$m.Vem_ed))

length(get_significant_results(results$`4`$m.Vem_ed))


get_pairwise_sharing(results$`4`$m.Vem_ed)
# Get effects significant on all conditions (leafs)
# Signs and magnitude might be different
n_leaf <- get_n_significant_conditions(results$`4`$m.Vem_ed,)

bf[names(n_leaf[n_leaf >2])] %>% sort() %>% tail()

quartz()

bf <- get_log10bf(results$`4`$m.Vem_ed) %>% as.numeric()
names(bf) <-row.names(results$`4`$m.Vem_ed$result$lfsr)
sort(bf) %>% tail()
str(results$`4`$m.Vem_ed)
get_lfsr(results$`4`$m.Vem_ed)


deg_idx <- get_significant_results(results$`4`$m.Vem_ed)
deg <- names(deg_idx)
deg_1 <- deg
length(deg)
paste(deg, collapse =  ",")
deg_lfsr <- get_lfsr(results$`4`$m.Vem_e)[deg_idx,]
deg_lfsr[order(deg_lfsr[,2]),]

head(get_pm(results$`4`$m.Vem_ed)[deg_idx,])

quartz()
barplot(get_estimated_pi(results$`4`$m.Vem_ed),las = 2)

print(get_pairwise_sharing(results$`4`$m.Vem_ed)) 



deg_idx <- get_significant_results(results$`2`$m.Vem_ed)
deg <- names(deg_idx)
dd <- deg[deg %in% deg_1]

paste(dd, collapse =  ",")
deg_lfsr <- get_lfsr(results$`2`$m.Vem_e)[deg_idx,]
deg_lfsr[order(deg_lfsr[,2]),]

deg_lfsr[order(deg_lfsr[,4]),]


barplot(get_estimated_pi(results$`3`$m.Vem_ed),las = 2)
print(get_pairwise_sharing(results$`2`$m.Vem_ed))



