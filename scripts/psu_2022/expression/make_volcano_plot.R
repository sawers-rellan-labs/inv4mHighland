library(dplyr)
library(ggplot2)
library(kableExtra)
 
effects <- read.csv("~/Desktop/DEGs_FDR.csv")  %>%
  dplyr::filter(predictor %in% c("leaf_tissue","Treatment-P" ,"GenotypeINV4")) %>%
  droplevels() 
  
effects$predictor <-factor(effects$predictor)

levels(effects$predictor) 

levels(effects$predictor) <- c("inv4m","Leaf","-P")

predictor_order <- c("Leaf", "-P","inv4m")

effects$predictor <- factor(effects$predictor, levels=predictor_order )


effects$locus_name[effects$gene=="Zm00001eb034810"] <-"monogalactosyldiacylglycerol synthase 2"
effects[effects$gene=="Zm00001eb191790",]
n_top <- 10

top_labelled <- effects %>%
  filter(regulation!="Unregulated") %>%
  group_by(predictor,regulation) %>%
  filter(!is.na(locus_symbol) & !is.na(locus_name)) %>%
  arrange(-mahalanobis, .by_group = TRUE) %>%
  select(predictor,locus_symbol,gene) %>%
  dplyr::slice(1:n_top) %>%
  mutate(label = locus_symbol)  %>% print(n=60)

top_labelled$label[top_labelled$gene=="Zm00001eb034810"] <-"mgd2"

nrow(top_labelled)

effects <-effects %>% left_join(top_labelled)
levels(effects$predictor)


#ggpubr::ggarrange(v2,v1,v3,ncol=3,common.legend = TRUE,align = "hv")


scaleFUN <- function(x) sprintf("%.1f", x)

png(file="~/Desktop/volcanoplot_A.png",height=6, width=12, res=300, units = "in")
effects %>%
  ggplot(aes(x=logFC, y=-log10(adj.P.Val), col= regulation, label= label))+
  ylab(expression(-log[10](italic(FDR)))) +
  xlab(expression(-log[2]("Fold Change"))) +
  scale_y_continuous(labels=scaleFUN) +
  geom_point(alpha=0.3) +
  # stat_ellipse(col= "grey75",level = 0.99)+
  ggrepel::geom_text_repel(color ="black",
                           force=1,
                           fontface="italic",
                           min.segment.length = 0,
                           max.overlaps = 50) +
  scale_color_manual( name="", values = c("#00AFBB", "grey","#bb0c00"), # to set the colours of our variable  
                      labels = c("Downregulated","Not significant","Upregulated")) +
  facet_wrap(.~predictor,scales = "free_x") +
  guides(color = guide_legend(override.aes = list(size=7)))+
  ggpubr::theme_classic2(base_size = 25) +
  theme(plot.title = element_blank(),
        strip.text = ggtext::element_markdown(),
        legend.position = c(0.45,0.95),
        legend.background =  element_rect(fill='transparent'), 
        legend.spacing =  unit(0,"line"),
        legend.box.spacing =unit(0, "line"),
        legend.text = element_text(size=15),
        legend.direction = "horizontal",
        strip.background = element_blank())
dev.off()

# effects %>%
#   filter(predictor=="Treatment-P:GenotypeINV4") %>%
#   ggplot(aes(x=logFC, y=-log10(adj.P.Val), col= regulation, label=label))+
#   ggtitle("-P x inv4m") +
#   ylab("-log10(FDR)") +
#   xlab("log2(Fold Change)") +
#   geom_point(alpha=0.3) +
#   ggrepel::geom_text_repel(data=subset(effects %>%
#                                          filter(predictor=="Treatment-P:GenotypeINV4" & !is.na(label))),
#                            color ="black",
#                            fontface="italic",
#                            force=5,
#                            min.segment.length = 0.1) +
#   scale_color_manual( name="", values =c("#00AFBB", "grey","#bb0c00"), # to set the colours of our variable  
#                       labels = c("Downregulated","","Unregulated")) +
#   ggpubr::theme_classic2(base_size = 15) +
#   theme(legend.position = "top",
#         plot.title = element_text(face = "italic"),
#   )




  


write.csv(effects,file="~/Desktop/DEGs_FDR_for_table.csv")

effects %>% 
  group_by(predictor,regulation) %>%
  filter(!is.na(label)) %>%
  arrange(predictor,regulation,adj.P.Val) %>%
  mutate(logFC =  round(logFC,digits=2)) %>%
  mutate(neglogP =format(neglogP, digits = 2)) %>%
  select( predictor,regulation,locus_symbol,neglogP,logFC,locus_name) %>%
  kbl(caption="Table 1: Effects of diferent predictors on characterized maize genes",
      format= "latex",
      col.names = c("Predictor","Regulation","Locus","-log10(FDR)","logFC","Name"),
      align="r") %>%
  kable_classic(full_width = FALSE, html_font = "helvetica") %>%
  save_kable(file="~/Desktop/DEGs_named.tex",format="latex")

