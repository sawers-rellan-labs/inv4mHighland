library(dplyr)
library(ggplot2)
library(scales)

theme_unique_dark <- function (base_size = 12, base_family = "") {
  ret <- (theme_bw(base_size = base_size, base_family = base_family) +
            theme(text = element_text(colour = "white"),
                  title = element_text(color = "white"),
                  line = element_line(color = "white"),
                  rect = element_rect(fill = "black", color = "white"),
                  axis.ticks = element_line(color = "#969696"),
                  axis.title = element_text(color = "white"),
                  axis.text = element_text(color = "#eaeaea"),
                  axis.line = element_line(color = "#969696", linetype = 1),
                  legend.background = element_rect(fill = NULL, color = NULL),
                  legend.position = "bottom",
                  legend.key = element_rect(fill = NA, color = NA, linetype = 0),
                  strip.background = element_rect(fill=NA,colour=NA,size=NA,linetype = NULL),
                  strip.text = element_text(color="white",face="bold",vjust=.5,hjust=.5),
                  panel.background = element_rect(fill = "black", color = NULL),
                  panel.border = element_blank(),
                  panel.grid = element_line(color = "#252525"),
                  panel.grid.major = element_line(color = "#353535"),
                  panel.grid.minor = element_line(color = "#101010"),
                  plot.background = element_rect(fill = "black", colour = "black", linetype = 0)))
  ret
}

blastcols <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

knob <- read.table("~/Desktop/knob.blast", sep= "\t", header = FALSE) 

colnames(knob) <- blastcols 
knob$sstart <- as.integer(knob$sstart)
knob$send <- as.integer(knob$send)
scale_factor <- 333 # 333  max(knob_mex$bitscore)
knob$bitscaled <- scales::rescale(c(knob$bitscore,scale_factor),
                                     to=c(0,1)) %>% head(-1)

knob$bitscaled <- scales::rescale(knob$bitscore, to=c(0,1))
knob$queryg <- "B73"
knob$rep <- "knob180"

# TR-1_b73.blast      TR-1_hsu_2002.fasta TR-1_pt.blast       TR-1_til18.blast
tr1 <- read.table("~/Desktop/TR-1_b73.blast", sep= "\t", header = FALSE)
colnames(tr1) <- blastcols 
tr1$bitscaled <- scales::rescale(tr1$bitscore, to=c(0,1))
tr1$queryg <- "B73"
tr1$rep <- "TL-R1"

sum((tr1$qend -tr1$qstart) > 300)

to_plot <- rbind( knob,tr1)



inv4m_knobs <- knob %>%
  dplyr::filter(qstart > 150e6, qend < 200e6) %>%
  arrange(qstart) %>% 
  mutate(qlength= abs(qend-qstart))




library(ggplot2)
library(scales)



b73 <- knob %>%
  ggplot(aes(x=qstart, y=bitscaled)) +
  ggtitle("B73")+
  # annotate("rect", fill = "red", 
  #          xmin = 172561959, xmax = 172882309,
  #          ymin = -Inf, ymax = Inf) + 
  # annotate("rect", fill = "red", 
  #          xmin = 188131461, xmax = 188220418,
  #          ymin = -Inf, ymax = Inf) + 
  geom_vline(xintercept=c(172882309,188131461),
             #col = "grey75", 
             linewidth = 1.3)+
  geom_point(data =tr1, col = "gold") +
  geom_point(col="#3232a0") +
  coord_cartesian(xlim=c(170000000,250000000))+
  scale_x_continuous(labels = unit_format(unit = "", scale = 1e-6), breaks = c(175e6, 200e6,225e6, 250e6)) +
  ggpubr::theme_classic2(base_size = 20)+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(hjust=0.8),
        plot.title = element_text(hjust=0.5))  

quartz()
b73

  

knob_mex <- read.table("~/Desktop/knob_til18.blast", sep= "\t", header = FALSE) %>% arrange(V7)

colnames(knob_mex) <- blastcols
knob_mex$bitscaled <- scales::rescale(knob_mex$bitscore, to=c(0,1))
knob_mex$queryg <- "TIL18"
knob_mex$rep <- "TL-R1"

# TR-1_b73.blast      TR-1_hsu_2002.fasta TR-1_pt.blast       TR-1_til18.blast
tr1_mex <- read.table("~/Desktop/TR-1_til18.blast", sep= "\t", header = FALSE)
colnames(tr1_mex) <- blastcols 
tr1_mex$bitscaled <- scales::rescale(tr1_mex$bitscore, to=c(0,1))

tr1$queryg <- "TIL18"
tr1$rep <- "TL-R1"

sum((tr1$qend -tr1$qstart) > 300)

to_plot <- rbind(to_plot, knob_mex,tr1)



to_rle <- inv4m_knobs %>%
  filter(qstart > 172000000, qend <200000000) %>%
  arrange(qstart) %>%
  mutate(qmid = round((qstart+qend)/2)) %>% 
  mutate( d = qmid- lag(qmid,default = 0)) %>%
  mutate(is_linked= if_else(d<100000,1,0))

nrow(knob_mex)

knob_rle <- rle(to_rle$is_linked) 

knob_rle$lengths[2*(1:12)]
knob_rle$lengths[2*(1:12)-1]

knob_rle$lengths[knob_rle$values==1]

str(knob_rle)


to_rle$cluster <- rep(1:length(knob_rle$lengths),knob_rle$lengths)

table(to_rle$cluster)
to_rle[to_rle$cluster==2,] 
to_rle[to_rle$cluster==6,]
to_rle[to_rle$cluster==7,] %>% tail()

to_rle$qstart[to_rle$cluster==2] %>% range(na.rm = TRUE) %>% diff()
to_rle$qstart[to_rle$cluster==6] %>% range(na.rm = TRUE) %>% diff()

to_rle[to_rle$cluster==2,] %>% nrow()
to_rle$qlength[to_rle$cluster==2] %>% sum()
to_rle$qstart[to_rle$cluster==2] %>% range(na.rm = TRUE)
to_rle$qstart[to_rle$cluster==2] %>% range(na.rm = TRUE) %>% diff()
to_rle[to_rle$cluster==4,] %>% nrow()
to_rle$qlength[to_rle$cluster==4] %>% sum()
to_rle$qstart[to_rle$cluster==4] %>% range(na.rm = TRUE)
to_rle$qstart[to_rle$cluster==4] %>% range(na.rm = TRUE) %>% diff()
to_rle[to_rle$cluster==6,] %>% nrow()
to_rle$qlength[to_rle$cluster==6] %>% sum()
to_rle$qstart[to_rle$cluster==6] %>% range(na.rm = TRUE) 
to_rle$qstart[to_rle$cluster==6] %>% range(na.rm = TRUE) %>% diff()


library(rtracklayer)

library(GenomicRanges)

knob_ranges <- GenomicRanges::makeGRangesFromDataFrame(
 knob_pt%>%
  filter(qstart > 173000000, qend <200000000) %>%
  arrange(qstart) %>%
    mutate(strand=case_when(
      sign(send-sstart)>0 ~"+",
      sign(send-sstart)<0 ~"-",
      .default = "*"),
      source= "blast",
      type= "knob180_repeat"),
  seqnames.field = "qseqid",
  start.field = "qstart",
  end.field = "qend",
  strand.field = "strand",
  keep.extra.columns = TRUE,
  
)
export.gff(knob_ranges,"~/Desktop/pt_knobs.gff3", version="3")

mex_ranges <- GenomicRanges::makeGRangesFromDataFrame(
  knob_mex %>% 
    filter(qstart>180000000, qend < 195000000) %>%
    mutate(strand=case_when(
           sign(send-sstart)>0 ~"+",
           sign(send-sstart)<0 ~"-",
           .default = "*"),
           source= "blast",
           type= "knob180_repeat"),
  seqnames.field = "qseqid",
  start.field = "qstart",
  end.field = "qend",
  strand.field = "strand",
  keep.extra.columns = TRUE,
  
)
export.gff(mex_ranges,"~/Desktop/mex_knobs.gff3", version="3")



mex <- knob_mex %>%
  ggplot(aes(x=qstart, y=bitscaled)) +
  ggtitle("TIL18")+
  # annotate("rect", fill = "red", 
  #          xmin = 180269950, xmax = 180365316,
  #          ymin = -Inf, ymax = Inf) + 
  # annotate("rect", fill = "red", 
  #          xmin = 193570651, xmax = 193734606,
  #          ymin = -Inf, ymax = Inf) + 
  geom_vline(xintercept=c(180365316,193570651),
             #col = "grey75",
             linewidth = 1.3) +
  geom_point(data =tr1_mex, col= "gold") +
  geom_point(col= "#3232a0") +
  coord_cartesian(xlim=c(175000000,250000000))+
  scale_x_continuous(labels = unit_format(unit = "", scale = 1e-6), breaks = c(175e6, 200e6,225e6, 250e6)) +
  ggpubr::theme_classic2(base_size = 20)  +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(hjust=0.8),
        plot.title = element_text(hjust=0.5))   

quartz()
mex

knob_pt <- read.table("~/Desktop/knob_pt.blast", sep= "\t", header = FALSE)
colnames(knob_pt) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
head(knob_pt)
scale_factor <- max(knob_mex$bitscore) # 333
knob_pt$bitscaled <- scales::rescale(c(knob_pt$bitscore,scale_factor),
                                    to=c(0,1)) %>% head(-1)
knob_pt$gquery <- "PT"
knob_pt$rep <- "knob180"
tr1_pt <- read.table("~/Desktop/TR-1_pt.blast", sep= "\t", header = FALSE)
colnames(tr1_pt) <- blastcols 


scale_factor <- max(tr1_mex$bitscore) # 671
tr1_pt$bitscaled <- scales::rescale(c(tr1_pt$bitscore,scale_factor),
                                      to=c(0,1)) %>% head(-1)

tr1_pt$gquery <- "PT"
tr1_pt$rep <- "knob180"

to_plot <- rbind(to_plot,knob_pt,tr1_pt)



173369064-173486186
186925483-187092654

pt <- knob_pt %>%
  ggplot(aes(x=qstart, y=bitscaled)) +
  ggtitle("PT")+
  # annotate("rect", fill = "red", 
  #          xmin = 173369064, xmax = 173486186,
  #          ymin = -Inf, ymax = Inf) + 
  # annotate("rect", fill = "red", 
  #          xmin = 186925483, xmax = 187092654,
  #          ymin = -Inf, ymax = Inf) + 
  geom_vline(xintercept=c(173486186,186925483),
             #col = "grey75", 
             linewidth = 1.3) +
  geom_point(data =tr1_pt, col= "gold") +
  geom_point(color="#3232a0") +
  coord_cartesian(xlim=c(170000000,250000000), ylim=c(0,1))+
  scale_x_continuous(labels = unit_format(unit = "", scale = 1e-6), breaks = c(175e6, 200e6,225e6, 250e6)) +
  ggpubr::theme_classic2(base_size = 20) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(hjust=0.8),
        plot.title = element_text(hjust=0.5))  


quartz(height =7, width=7)
ggpubr::ggarrange(mex,pt, b73,  nrow=3, align = "hv")+
  ggtitle("Knob repeats delimit Inv4m")





