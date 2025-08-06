# DEG characterization.
library(dplyr)
library(scales)
library(ggplot2)
library(rtracklayer)
library(GenomicRanges)

myGFF <- "~/ref/zea/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.59.chr.gff3"
B73<- rtracklayer::import(myGFF)  %>%
  subset(type=="gene"  & seqnames %in% 1:10)


genes <- as.data.frame(B73) 

genes$ID <- gsub("gene:","",genes$ID)

# Extract these coordinates from the MCScan data


inv4m_start <- genes[genes$ID=="Zm00001eb190470","start"]
inv4m_start

inv4m_end <- genes[genes$ID=="Zm00001eb194780","end"]
inv4m_end

inv4m <- GRanges(seqnames="4", 
                 ranges=IRanges(start=inv4m_start, end=inv4m_end),
                 strand = "+")

shared_introgression <- GRanges(seqnames="4", 
                                ranges=IRanges(start=157012149, end=195900523),
                                strand = "+")

introgression_limit <- round(c(157012149,195900523)/1000000,2) # checked manually


olap <-findOverlaps(inv4m,B73)
inv4m_gene_ids <- genes$ID[subjectHits(olap)]
length(inv4m_gene_ids)



olap2 <- findOverlaps(shared_introgression, B73)
shared_introgression_gene_ids <- genes$ID[subjectHits(olap2)]
length(shared_introgression_gene_ids)


flanking_introgression_gene_ids <- shared_introgression_gene_ids[!(shared_introgression_gene_ids %in% inv4m_gene_ids)]

length(flanking_introgression_gene_ids)

inv4m_start
# 172883675

inv4m_end
# 188082081

inversion_limit <-round(c(inv4m_start,inv4m_end)/1000000,2)
inversion_limit[2] - inversion_limit[1]
# 15.2

round(introgression_limit[2]-introgression_limit[1],1)

effects <- read.csv( file= "~/Desktop/predictor_effects.csv")
effects %>%
  filter(gene == "Zm00001eb191830") %>%
  select(predictor:logFC,adj.P.Val,locus_name, regulation)
effects$gene[grepl("Zm00001eb1918", effects$gene)] %>% sort() %>% unique()

#fisher.test(contingency)
#library(R.devices)

#quartz()
#fastqq(effect_inv4m$P[effect_inv4m$in_inv4m], speedup = FALSE)

# nulldev()
# lambda <- c(inv4m = fastqq(effect_inv4m$P[effect_inv4m$in_inv4m], speedup = FALSE),
#             flanking = fastqq(effect_inv4m$P[effect_inv4m$in_flanking],speedup = FALSE),
#             outside = fastqq(effect_inv4m$P[!effect_inv4m$in_shared], speedup = FALSE)
# ) %>% round( 2)
# dev.off()

# calculating Standard errors of lambda
# CIs overlap between inv4m and flanking

pvalue <- effect_inv4m$P[effect_inv4m$in_flanking]


chisq <- qchisq(1 - pvalue, 1)
Lambda <- median(chisq) / qchisq(0.5, 1)


sum(!is.infinite(chisq))

SE.median <- qchisq(0.975, 1) * (1.253 * ( sd(chisq[!is.infinite(chisq)]) / sqrt(sum(!is.infinite(chisq)) ) ) )

c(Lambda - (SE.median / qchisq(0.5, 1)),Lambda + (SE.median / qchisq(0.5, 1)))



# This is it!!!!!
# It is an statistical test of the difference between the two Pvalues distributions
# The alternative hypothesis must be related
#  to the number of outliers / in inv4m

ks.test(effect_inv4m$P.Value[effect_inv4m$in_inv4m], 
        effect_inv4m$P.Value[effect_inv4m$in_flanking])

ks.test(effect_inv4m$P.Value[effect_inv4m$in_inv4m], 
        effect_inv4m$P.Value[effect_inv4m$in_flanking],
        alternative = "greater")

ks.test(effect_inv4m$P.Value[effect_inv4m$in_flanking],
        effect_inv4m$P.Value[effect_inv4m$region=="outside"], 
        alternative = "greater")

ks.test(effect_inv4m$P.Value[effect_inv4m$in_shared], 
        effect_inv4m$P.Value[effect_inv4m$region=="outside"],
        alternative = "greater")



inv4m_dens <- qqplot(runif(n=sum(10000)),effect_inv4m$P[effect_inv4m$in_inv4m], plot.it = FALSE)
flanking_dens <- qqplot(runif(n=10000),effect_inv4m$P[effect_inv4m$in_flanking], plot.it = FALSE)
outside_dens <- qqplot(runif(n=10000),effect_inv4m$P[!effect_inv4m$in_shared], plot.it = FALSE)


qq_data <- rbind(
inv4m_dens %>% as.data.frame() %>%
  mutate(region = "inv4m") %>%
  mutate(expected = -log10(x), observed =-log10(y)),
flanking_dens %>% as.data.frame() %>%
  mutate(region = "flanking") %>%
  mutate(expected = -log10(x), observed =-log10(y)),
outside_dens %>% as.data.frame() %>%
  mutate(region = "outside") %>%
  mutate(expected = -log10(x), observed =-log10(y))
)

qq_data$region <- factor(qq_data$region,levels=c("inv4m","flanking","outside"))
levels(qq_data$region) <-c("*Inv4m*","flanking","outside")


key <- paste0(names(lambda), " $\\lambda = ", round(lambda,2),"$")
key_order <- c("*Inv4m*","flanking","outside")

pal0 <-c("purple4","darkgreen","gold")
pal1 <- c("black","grey50","white")
pal2 <- viridis_pal(option = "plasma")(3)
pal <- pal2

key_order <- c("*Inv4m*","flanking","outside")
names(pal) <- key_order


# QQplot ----
qqplot <- qq_data %>%
  ggplot(aes(y=observed, x=expected, color=region)) +
  ylab(expression("Observed "-log[10](italic(P))))+
  xlab(expression("Expected "-log[10](italic(P))))+
  geom_abline(slope=1, color="black") +
  geom_point(size=3) +
  # annotate('text', x = 2.5, y = c(31,22,13),
  #          hjust= 0,
  #          label = c(expression(lambda == "27.80"), 
  #                    expression(lambda == "22.30"),
  #                    expression(lambda == "1.81")),
  #         size = 5) +
  annotate('text', x = 2.3, y = c(25,15),
           hjust= 0,
           label = c(expression(italic(p) == "0.29"), 
                     expression(italic(p) == "2.2e-16")),
           size = 5) +
  scale_color_manual(values=pal) +
  ggpubr::theme_classic2(base_size = 25) +
  theme(legend.position = c(0.17,0.83),
        legend.margin=margin(0,0,0,0),
        legend.box.spacing = unit(0, "pt"),
        legend.text = ggtext::element_markdown())


quartz()
qqplot


key_order <- c("inv4m","flanking","outside")
names(pal) <- key_order


levels(factor(effects$predictor))
table(effects$predictor)
table(effect_inv4m$predictor)




# PLot manhattan ------


get_plot_data <- function(stats) {
  stats %>%
    dplyr::select(gene,SNP="locus_symbol", CHR, BP, P) %>%
      mutate(SNP=paste("  ",SNP)) %>%
    distinct() %>%
    arrange(CHR,BP)
}

plot_manhattan <-function(colticks=NA, highlight, annotateHighlight=FALSE, annotateTop=TRUE, annotateN=3) {
  fastman( m = to_plot, maxP = maxnegLogP,
           #ylab=expression(-log[10](italic(FDR))),
           ylab=expression(-log[10](italic(P))),
           cex.axis = 1.5, cex.lab=1.5,
           cex.text=0.8,
           highlight=highlight,
           #annotateHighlight = annotateHighlight,
           #annotationCol = "black",
           #annotationAngle=90,
           #annotateTop = annotateTop,
           #annotateN = annotateN,
           suggestiveline = BH_threshold,
           genomewideline = bonferroni_threshold,
           chrlabs=rep("",10),
           xlab ="", colticks=colticks) 
}



# 

# leaf_tissue #################

levels(factor(effects$predictor))

r <- effects %>% filter(predictor=="leaf_tissue")

p <- r$P.Value

names(p) <- r$gene


maxnegLogP <- -log10(min(p[p>0]))

names(p)[p<0.05]

sum(p<0.05)

to_plot <- get_plot_data(r)


BH_threshold <-RAINBOWR::CalcThreshold(
  to_plot %>% dplyr::select(gene,CHR,BP,P) %>%
    mutate(neglogP = -log10(P))  %>%
    dplyr::select(-P)
)

leaf.plot <- ggplotify::as.ggplot(
  ~{
    par(mgp=c(2.5,1,0))
    par(mar = c(5,4,4,2) + 0.1) ## default is c(5,4,4,2) + 0.1
    fastman( m = to_plot, maxP = maxnegLogP,
             #ylab=expression(-log[10](italic(P))),
             cex.axis = 1.5, cex.lab=1.5,
             cex.text=0.8,
             #annotationCol = "black",
             annotationAngle=90,
             annotateN = 3,
             annotateTop = TRUE,
             suggestiveline = BH_threshold,
             genomewideline = bonferroni_threshold,
             #chrlabs=rep("",10)
             ylab =""
    )
    title(ylab = expression(-log[10](italic(P))), cex.lab = 1.5)
  }
)


quartz()
leaf.plot


# lowP.plot ----

#  "Treatment-P" -----
factor(effects$predictor) %>% levels()

r <- effects %>% filter(predictor=="Treatment-P")

p <- r$P.Value
names(p) <- r$gene

maxnegLogP <- -log10(min(p[p>0]))

names(p)[p<0.05]

sum(r$neglogP > bonferroni_threshold)

to_plot <- get_plot_data(r)

# To  get other genes to highlight
# not useful for now
#
# highlight <- r %>%
#   filter(!is.na(locus_name)) %>% 
#   arrange(-neglogP) %>% 
#   pull(locus_symbol) %>%
#   paste("  ",.)

BH_threshold <-RAINBOWR::CalcThreshold(
  to_plot %>% dplyr::select(gene,CHR,BP,P) %>%
    mutate(neglogP = -log10(P))  %>%
    dplyr::select(-P)
)

lowP.plot <- ggplotify::as.ggplot(
  ~{par(mgp=c(2.5,1,0))
    par(mar = c(5,4,4,2) + 0.1) ## default is c(5,4,4,2) + 0.1
    plot_manhattan()}
)


quartz()
lowP.plot

###############################################################################
# inv4m.plot ----

levels(factor(effects$predictor))
r <-effect_inv4m%>% filter(predictor=="GenotypeINV4")
r$is_significant <- r$adj.P.Val < 0.05
r %>% arrange(P)
p <- r$adj.P.Val
names(p) <- r$gene
table(r$is_significant)

# top hit with no interaction model
# Zm00001eb194380 
# SEC6 exocyst complex protein
# secretory pathway of the cell plate
# growth in Arabidopsis and psychomitrella
# expression pattern in B73v4
# https://www.maizegdb.org/gene_center/gene/Zm00001d052245


maxnegLogP <- -log10(min(p[p>0]))

names(p)[p<0.05]

sum(p<0.05)

sum(r$neglogP<bonferroni_threshold)


to_plot <- get_plot_data(r) 


table(effect_inv4m$region)
table(effects$predictor)
highlight <- r %>%
  filter(!is.na(locus_name)) %>% 
  arrange(-neglogP) %>% 
  pull(locus_symbol) %>%
  paste("  ",.)


BH_threshold <-RAINBOWR::CalcThreshold(
  to_plot %>% dplyr::select(gene,CHR,BP,P) %>%
    mutate(neglogP = -log10(P))  %>%
    dplyr::select(-P)
)

# inv4m.plot <- ggplotify::as.ggplot(
#   ~{plot_manhattan()}
# )

inv4m.plot <- ggplotify::as.ggplot(
  ~{par(mgp=c(2.5,1,0))
    par(mar = c(5,4,4,2) + 0.1) ## default is c(5,4,4,2) + 0.1
    plot_manhattan()}
)

# Align Manhattan plots ------

png(file="~/Desktop/Manhattan_by_predictor.png", width =9, height=11, units = "in", res=300)

cowplot::plot_grid(inv4m.plot,lowP.plot,leaf.plot,
                   #inv4m.by.lowP.plot, 
                   ncol = 1, nrow = 3,
                   rel_heights = c(1,1,1),
                   align = 'v', axis = 'lr') 
dev.off()

effects$P.Value
factor(effects$predictor) %>% levels()
effects %>%
  ungroup() %>%
  filter(predictor=="GenotypeINV4",CHR==4) %>%
  filter(BP>150e6 & BP<200e6) %>%
  arrange(BP) %>% tibble() %>%
  select(-locus_symbol,-(logFC:t), -(adj.P.Val:lower),-B,  
         -(upregulated:is_mh_outlier) )%>%
  print(n=100)

pal <-c("#C84B79","#0F1D87","#F1FA41")
FDR_thresh = -log10(0.05)
p1 <- r  %>%
  filter(predictor=="GenotypeINV4",CHR==4) %>%
  mutate(x = BP/1e6) %>%
  ggplot(aes( x= x, y= -log10(adj.P.Val), group=region,color=region))+
  # ggtitle("inv4m") +
  ylab(expression(-log[10](italic(FDR)))) +
  # geom_vline(xintercept = inversion_limit, linetype="dashed", linewidth=1) +
  # geom_vline(xintercept = introgression_limit, linetype="dashed",linewidth=1) +
  annotate('text', x = c(inversion_limit,introgression_limit) , y = -0.5,
          label = rep("|",4), size = 10) +
  geom_point() + 
  scale_color_manual(values =pal)+
  scale_y_continuous(expand= c(0,0.5)) +
  geom_hline(yintercept = FDR_thresh, col= "red") +
  #geom_hline(yintercept = BH_threshold, col= "blue") +
  coord_cartesian(xlim=c(150,200)) +
  #xlab(NULL)+
  #ggtitle("Chromosome 4 Position [Mb]") +
  xlab("Chromosome 4 Position [Mb]")+
  ggpubr::theme_classic2(base_size = 25) +
  theme(legend.position="none", 
        # axis.title.x = element_blank(),
        # axis.line.x = element_blank(),
        # axis.text.x = element_blank(),
        # axis.ticks.x = element_blank(),
        plot.margin = unit(c(0.5, 0, 0.5, 0), "cm"))


contingency <- with(
  r %>% mutate( is_significant = P.Value <0.05/nrow(r)),
  table(region, is_significant))
contingency



contingency <- with(r,
     table(region, is_significant))

chisq.test(contingency)
contingency

contingency <- with(r %>% filter(region !="outside"),
     table(region, is_significant))

chisq.test(contingency)
contingency

# genotype x treatment effect

# trt_deg <- names(get_significant_results(results$`4`$m.Vem_ed))
# de_gr <-   genes %>% dplyr::filter(gene_id %in% trt_deg)

#:) this look good
# lfsr <- as.data.frame(get_lfsr(results$`4`$m.Vem_ed))
# p <- pmin( lfsr$V1,lfsr$V2,lfsr$V3,lfsr$V4)

# -P x inv4m effect

###############################################################################
# 
# levels(factor(effects$predictor))
# 
# r <- effects %>% filter(predictor=="Treatment-P:GenotypeINV4")
# 
# p <- r$P.Value
# names(p) <- r$gene
# 
# maxnegLogP <- -log10(min(p[p>0]))
# 
# 
# ?p.adjust()
# # inv4m.by.P.plot ----
# to_plot <- get_plot_data(r,genes)
# 
# 
# inv4m.by.lowP.plot <- ggplotify::as.ggplot(
#   ~{
#     fastman( m = to_plot, maxP = maxnegLogP,
#              ylab=expression(-log[10](italic(P))),
#              cex.axis = 1.5, cex.lab=1.5,
#              cex.text=0.8,
#              #annotationCol = "black",
#              annotationAngle=90,
#              annotateN = 3,
#              annotateTop = TRUE,
#              suggestiveline = BH_threshold,
#              genomewideline = bonferroni_threshold,
#              #chrlabs=rep("",10),
#              xlab ="",  yaxt="n") 
#     axis(side = 2, at=0:3,
#          cex.lab=1.5,
#          cex.axis=1.5,
#          las =2)
#   }
# )
# 


# Align zoom

# p2 <- to_plot %>%
#   filter(CHR==4) %>%
#   mutate(span = BP/1e6) %>%
#   ggplot(aes( x= span, y= -log10(P)))+
#   # ggtitle(label = "-P x inv4m")+
#   ylab(expression(-log[10](italic(FDR)))) +
#   geom_vline(xintercept = inversion_limit, col = "purple4", linetype="dashed", size=1) +
#   geom_vline(xintercept = introgression_limit, col= "purple4",linetype="dashed",size=1) +
#   geom_point(col= "#60a500") +
#   geom_hline(yintercept = bonferroni_threshold, col= "red") +
#   geom_hline(yintercept = BH_threshold, col= "blue") +
#   xlab("Chr4")+
#   coord_cartesian(xlim=c(150,200)) +
#   # scale_x_continuous(labels =  unit_format(unit = "",scale = 1e-6))+
#   ggpubr::theme_classic2(base_size = 25) +
#   theme( # plot.title = element_text(face = "italic", hjust = 1),
#         # axis.title.x = element.text(color=NA)
#         )

# png(file="~/Desktop/Chr04_expression.png", width =7, height=10, units = "in", res=300)

# Align left panel ------


load(file="~/Desktop/gt3_plot.rda")

detail <- cowplot::plot_grid(qqplot,gt3_plot,p1, 
                             ncol = 1, nrow = 3,
                             rel_heights = c(1.5,0.5,1),
                             align = 'v', axis = 'lr')


pal
png(file="~/Desktop/Chr04_expression.png", width =8, height=12, units = "in", res=300)
print(detail)
dev.off()

detail <- cowplot::plot_grid(gt3_plot,p1, qqplot,
                             ncol = 1, nrow = 3,
                             rel_heights = c(0.5,1,1.5),
                             align = 'v', axis = 'lr')

png(file="~/Desktop/Chr04_expression2.png", width =8, height=12, units = "in", res=300)
print(detail)
dev.off()
lambda


# annotate("path", 
#          x = bracket1$x, 
#          y = bracket1$y,
#          color = bracket1$color,
#          size=2) +
# annotate("path", 
#          x = bracket2$x, 
#          y = bracket2$y,
#          color = bracket2$color,
#          size=2) +
#   annotate("path", 
#            x = bracket3$x, 
#            y = bracket3$y,
#            color = bracket3$color,
#            size=2) +
#   annotate('text', x = mean(inversion_limit) ,
#            y = -0.6, fontface = 'italic',
#            label = "inv4m", size = 10, color=pal[1]) +
  
  
# png(file="~/Desktop/fig2.png", width =8, height=12, units = "in", res=300)
# 
# cowplot::plot_grid(byrow = FALSE, 
#                    leaf.plot,lowP.plot,inv4m.plot, 
#                    qqplot,gt3_plot,p1, 
#                    ncol = 2, nrow = 3,
#                    rel_heights = c(1,1,1),
#                    rel_widths = c(0.8,1),
#                    align = 'v', axis = 'lr') 
# dev.off()


