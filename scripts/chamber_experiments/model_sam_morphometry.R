library(dplyr)
library(ggplot2)
library(ggfx)

microscopy <- read.csv("~/Desktop/MeristemMesurements_B73-INV_December2024.csv") %>% 
  dplyr::select(-starts_with("X")) %>%
  mutate(Ratio=Height/Width,
         radius=Width/2,
         arc.length = sqrt(radius^2+4*Height^2+(radius^2/2*Height)+asinh(2*Height/radius)),
         hr.ratio=Height/radius,
         shape.factor = 1000*Height/(radius^2),
         area.section= (4/3)*(Height*radius), 
         area.surface =  pi*radius*((radius^2+4*Height^2)^(3/2) -radius^3)/(6*Height^2),
         volume = pi*(Height*radius^2)/2)
microscopy$Genotype <- gsub("B73", "CTRL",microscopy$Genotype)
microscopy$Genotype <- gsub("INV", "Inv4m",microscopy$Genotype)
microscopy$Genotype <-factor(microscopy$Genotype, levels= c("CTRL","Inv4m"))




measured1 <-  c("Height","Width")
measured2 <-  c("Height","Width","SurfaceArea")
measured3 <-  c("Height","Width","Ratio")

base_cols <- c("Height","Width","Ratio","SurfaceArea")

paraboloid_cols <- c("Height","radius","hr.ratio", "shape.factor", "arc.length","area.section","area.surface","volume")

height_cols <-  c("Height","hr.ratio", "shape.factor")

selected_cols <- base_cols

selected_cols <- paraboloid_cols


microscopy %>%
  group_by(Genotype) %>%
  summarise_at(paraboloid_cols,.funs = median)

# quartz()
# microscopy %>%
#   ggplot(aes(x=Height, y=Width, color=Genotype)) +
#   geom_point()
to_plot <-  microscopy %>%
  tidyr::pivot_longer(
    cols=all_of(selected_cols),
    names_to = "pheno",
    values_to = "value") %>%
  mutate(pheno = factor(pheno, levels= selected_cols)) %>%
  group_by(pheno) %>%
  mutate(z_score = scale(value)) %>%
  mutate(Genotype_label = case_when(
    grepl("Inv4m",Genotype) ~ paste("<i>Inv4m</i>"), 
    grepl("CTRL",Genotype) ~ "CTRL", 
    .default = NA)
  ) %>%
  mutate(Genotype_label = factor(Genotype_label, levels=c("CTRL","<i>Inv4m</i>")))

table(to_plot$Genotype_label)

t_test <- to_plot %>%
  rstatix::t_test( z_score~Genotype_label, ref.group="CTRL",
                   detailed = TRUE,
                   alternative = "less"
                   )  %>%
  # rstatix::adjust_pvalue(method = "fdr") %>%
  rstatix::add_significance("p") %>%
  rstatix::add_y_position()

t_test$y.position <- 1.1*t_test$y.position

t_test %>% filter(p <0.05)
# t_test$y.position =c(165,145,1.25,20000 )

pal =c(CTRL="gold",`<i>Inv4m</i>`="purple4")

pheno_theme <- ggpubr::theme_classic2(base_size = 20) +
  ggplot2::theme(
    #plot.title = element_text(hjust = 0.5, face ="bold"),
    plot.title = ggtext::element_markdown(hjust = 0.5),
    axis.title.y = ggtext::element_markdown(),
    axis.title.x = element_blank(),
    axis.text.x = ggtext::element_markdown(color= "black", size=20),
    strip.background = element_blank(),
    strip.text = ggtext::element_markdown(size=20),
    legend.position = "none")


ylabels <-   c( Height = "Height<br><i>h</i>",
    radius ="Radius <br><i>r</i>",
    hr.ratio = " Ratio <br> <i>h/r</i>",
    shape.factor = "Shape factor<br> <i>h/r<sup>2</sup></i>", 
    arc.length = "Arc Length",
    area.section = "Section Area",
    area.surface = "Surface Area",
    volume= "Volume")

H0 <-"<i>H<sub>0</sub>: &mu;<sub>Inv4m</sub> &le; &mu;<sub>CTRL</sub></i>"
Ha <-"<i>H<sub>a</sub>: &mu;<sub>INV4m</sub> &gt; &mu;<sub>CTRL</sub></i>"
sam_title <- paste("SAM Morphology Paraboloid Model<br>",
                   "One tailed <i>t-test</i> (", H0,"|",Ha,")<br>",
                   "without multiple hypothesis correction")

quartz(width =11,height =8) 
to_plot %>%
  ggplot2::ggplot(aes( x = Genotype_label, y = z_score, col = Genotype_label )) +
  ggtitle(sam_title) +
  ylab("Z score") +
  ggplot2::geom_boxplot(width = 0.25, linewidth=1,alpha=0) %>% with_shadow(
    colour = "white",
    x_offset = 0,
    y_offset = 0,
    sigma = 1) +
  ggbeeswarm::geom_quasirandom(size=2) %>% with_shadow(
    colour = "white",
    x_offset = 0,
    y_offset = 0,
    sigma = 1) +
  ggpubr::stat_pvalue_manual(
    t_test, label = "p",
    size = 5,
    bracket.size = 0.8,
    hide.ns = TRUE) +
  scale_color_manual(values = pal) +
  scale_y_continuous(expand = expansion(mult=c(0.05, 0.2))) +
  facet_wrap(.~pheno, nrow=2, 
             scales = "free_y",
             labeller=as_labeller(ylabels)) +
  pheno_theme



t_test <-to_plot %>%
  rstatix::t_test( value~Genotype_label, ref.group="<i>Inv4m</i>",
                   detailed = TRUE,
                   alternative = "greater"
  )  %>%
  # rstatix::adjust_pvalue(method = "fdr") %>%
  rstatix::add_significance("p") %>%
  rstatix::add_y_position(scales = "free_y")

(30.7-28.2)/28.2

t_test$y.position <- 1.05*t_test$y.position


ylabels <-   c( Height = "Height",
                radius ="Radius",
                hr.ratio = "HR Ratio",
                shape.factor = "Shape Factor", 
                arc.length = "Arc Length",
                area.section = "Section Area",
                area.surface = "Surface Area",
                volume= "Volume")

yunit_labels <-   c( Height ="<br><i>h</i> [&mu; m]",
                radius ="<br><i>r</i> [&mu; m]",
                hr.ratio = "<br> <i>h/r</i>",
                shape.factor = "<br> <i>h/r<sup>2</sup></i> [nm<sup>-1</sup>]", 
                arc.length = "&mu; m",
                area.section = "&mu; m<sup>2</sup>",
                area.surface = "&mu; m<sup>2</sup>",
                volume= "&mu; m<sup>3</sup>")

plot_pheno <- function(x) {
  to_plot %>%
  filter(pheno==x) %>%
  ggplot2::ggplot(aes( x = Genotype_label, y = value, col = Genotype_label )) +
  ggtitle(ylabels[x]) +
  ylab(yunit_labels[x]) +
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
    t_test %>% filter(pheno==x), label = "p.signif",
    size = 10,
    bracket.size = 0.8,
    hide.ns = TRUE) +
  scale_color_manual(values = pal) +
  scale_y_continuous(expand = expansion(mult=c(0.05, 0.075))) +
  pheno_theme +
  theme( plot.margin = unit(c(1,-1,1,-1), 'lines'))
}

all_plots <- lapply(paraboloid_cols[1:4], FUN = plot_pheno)
class(all_plots)

quartz(width =12,height =4) 
ggpubr::ggarrange(plotlist = all_plots,
                  nrow = 1, ncol=4,
                  align = "v",
                  widths= c(1,1,1,1.1))






Y <- microscopy %>% dplyr::select(all_of(selected_cols)) %>% as.matrix() 

scaled_Y <- Y %>% scale()





Genotype  <- matrix(as.factor(microscopy$Genotype) %>% as.numeric())

mlm1 <- lm(Y ~ Genotype,)

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

effects$p.adjust <- p.adjust(effects$p.value, method = "fdr")

effects$CI.upper <- effects$Estimate  + qnorm(0.975)*effects$Std.Error
effects$CI.lower <- effects$Estimate  - qnorm(0.975)*effects$Std.Error


effects %>%
  filter(p.value< 0.05) %>%
  filter( predictor == "Genotype") 
result$`Response Height`

pt(coef(result[["Response Height"]])[2,3], mlm1$df, lower = FALSE)

#Cohen's d
effects[2,3]/sd(microscopy$Height)
#0.56

effects[2,3]/sd(microscopy$Height)



