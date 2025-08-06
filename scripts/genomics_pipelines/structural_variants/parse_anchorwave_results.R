#Use R to draw a dotplot.
library(dplyr)
library(ggplot2)
library(svglite)

#Transform Coordinates using follow function.
changetoM <- function ( position ){
  position=position/1000000;
  paste(position, "M", sep="")}
#Read gene position, belong to which chromosome and so on
data =read.table("~/Desktop/PT_B73_cds.tab")
data$V1 <- gsub("chr","",data$V1)
data$V3 <- gsub("chr","",data$V3)
#Select all euchromosomes as factor.
data = data[which(data$V1 %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")),]
data = data[which(data$V3 %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")),]
data$V1 = factor(data$V1, levels=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"))
data$V3 = factor(data$V3, levels=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"))
#Using ggplot2 to plot a dotplot and beautify it.
figure1 <- data %>%
  filter(V4 > 170000000 & V4 < 190000000,V2 > 170000000 &  V2 < 190000000 ) %>%
  ggplot(aes(x=V4, y=V2)) +geom_point(size=0.5, aes(color=V5)) +
  facet_grid(V1 ~ V3, scales="free",space="free") +labs(x="PT", y="B73")+
  scale_x_continuous(labels=changetoM) + scale_y_continuous(labels=changetoM) +
  theme(axis.line = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
        axis.text.y = element_text( colour = "black"),
        legend.position='none',
        axis.text.x = element_text(angle=300, hjust=0, vjust=1, colour = "black"))
quartz()
figure1
# png("figure1.png")  
# figure1
# dev.off()
# pdf("figure1.pdf")
# figure1
# dev.off()
# svglite("figure1.svg")
# figure1
# dev.off()


anchors <- read.table("~/Desktop/anchors", header=TRUE)

# PT anchorwave
colnames(anchors)

anchors %>%
  filter(
    (referenceStart >= 172559094	& referenceStart <=  172884195) | 
    (referenceStart >= 188130243 & referenceStart <=  188231654)
)

break_pts <- anchors %>%
  filter(
    (referenceStart >= 172559094	& referenceStart <=  172884195) | 
      (referenceStart >= 188130243 & referenceStart <=  188231654)
  )

bounds <- c(break_pts$queryStart,break_pts$queryEnd) %>% sort

# Reference Upstream
bounds[6:7]

# Reference Downstream
bounds[2:3]



# TIL18 anchorwave

anchors <- read.table("~/Desktop/TIL18.anchors", header=TRUE)

break_pts <- anchors %>%
  filter(
    (referenceStart >= 172559094	& referenceStart <=  172884195) | 
      (referenceStart >= 188130243 & referenceStart <=  188231654)
  )

bounds <- c(break_pts$queryStart,break_pts$queryEnd) %>% sort

# Reference Upstream
bounds[6:7]

# Reference Downstream
bounds[2:3]

