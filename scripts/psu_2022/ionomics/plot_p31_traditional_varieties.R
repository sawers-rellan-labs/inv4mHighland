library(ggplot2)
library(dplyr)
library(forcats)
library(ggpubr)
library(ggfx)

P31 <- read.csv("/Users/fvrodriguez/Projects/CINVESTAV/p_data_landraces_b73.csv")
colnames(P31)
levels(P31$LANDRACE)


data <- P31 %>% 
  filter(LANDRACE %in% c("palomero","b73")) %>%
  mutate(LANDRACE = fct_recode(LANDRACE, "PT" = "palomero")) %>%
  mutate(LANDRACE = fct_recode(LANDRACE, "B73" = "b73" )) 

quartz()
ggplot(data, aes(LANDRACE, P31, color = LANDRACE ))+ geom_point( size = 4) +
  scale_color_manual(values=c("cornflowerblue","firebrick")) +
  ggpubr::theme_classic2(base_size = 20) +
  theme(legend.position = "none")


t.test <- 
  data %>%
    rstatix::t_test( P31 ~ LANDRACE, ref.group = "B73") %>%
  rstatix::adjust_pvalue(method = "BH") %>%
  rstatix::add_significance() %>% 
  rstatix::add_x_position(x = "LANDRACE") 

t.test$y.position <- 2750

pal <- c("gold","#4a0f82")


pheno_theme <- ggpubr::theme_classic2(base_size = 35) +
  ggplot2::theme(
    plot.title = element_text(hjust = 0.5, face ="bold"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(face="bold", color= "black", size=35),
    strip.background = element_blank(),
    legend.position = "none")

quartz(height=12, width=5)
data %>%
  droplevels() %>%
  ggplot2::ggplot(aes(x=LANDRACE, y=P31, color = LANDRACE )) +
  ggplot2::ggtitle("Leaf\nPhosphorus")  +
  ggplot2::ylab("mg/kg dry weight")  +
#  ylim(68,83)+
  ggplot2::geom_boxplot(width = 0.25,linewidth=1.5) %>% with_shadow(
    colour = "black",
    x_offset = 1,
    y_offset = 1,
    sigma = 1)+
  ggbeeswarm::geom_quasirandom(size= 5) %>% with_shadow(
    colour = "black",
    x_offset = 1,
    y_offset = 1,
    sigma = 1)+
  ggpubr::stat_pvalue_manual(
    t.test, 
    size = 10,
    bracket.size = 0.8,
    hide.ns = TRUE) +
  scale_color_manual(values = pal) +
  coord_cartesian(ylim = c(1000,3000)) +
  pheno_theme


summary(lm(data = data, P31~ LANDRACE))

data %>%
  group_by(LANDRACE) %>%
  summarize(mean = mean(P31), sd = sd(P31), cv = sd(P31)/mean(P31))
?summarise
