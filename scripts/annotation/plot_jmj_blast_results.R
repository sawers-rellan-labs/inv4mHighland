#plot blast hits
library(GenomicRanges)
library(ggbio)
library(dplyr)

jmj_annotation_table <- read.csv(file="~/Desktop/jmj_annotation_5_genes.csv")
exon_table <- jmj_annotation_table[jmj_annotation_table$type=="exon",]
exon_ranges<- with(exon_table,
               IRanges(start=start,end = end)
)

exons <- with(exon_table, 
              GRanges(
                seqnames=seqnames,
                ranges=exon_ranges,
                transcript = transcript_id,
                
                type="exon",
                track="jmj_annotation"
              )
)
exons$identity <-100


blast_cols <-c("qseqid",
               "sseqid",
               "identity",
               "length",
               "mismatch",
               "gapopen",
               "qstart",
               "qend",
               "sstart",
               "send",
               "evalue",
               "bitscore")


blast <- read.table("~/Desktop/B73_jmj_cluster_genomic_vs_B73_FLI-CDNA-jmj.blast", header=FALSE, sep="\t")
colnames(blast) <- blast_cols
blast$qseqid <- "4"
blast$qstart <- blast$qstart  + 177400000+1
blast$qend <- blast$qend  + 177400000+1

ranges <- with(blast,
IRanges(start=qstart,end = qend)
)


  
hits <- with(blast, 
  GRanges(
    type="blast_hit",
    seqnames=qseqid,
    ranges=ranges,
    transcript = sseqid,
    score = bitscore,
    identity = identity,
    track=sseqid
  )
)



blast_plot <- hits %>%
  subset(transcript !="EU945014.1") %>%
ggplot() +
  stat_stepping(aes(group = track, color= identity, fill = identity),
                facets=transcript ~.) +
  facet_wrap(transcript ~.,
             #scales="free_y",
             strip.position = "left",
             ncol=1)+
  scale_fill_continuous(name="% identity")+
  xlim(177400000,177590000)+
  theme(legend.position = "bottom",
        strip.text.y.left = element_text(angle =0),
        strip.background.y = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()
  )


gene_plot <- exons %>%
  ggplot() +
  stat_stepping(aes(group = track, color=transcript, fill = transcript),
                facets=transcript ~.)+
  facet_wrap(track ~ .,
             #scales="free_y",
             strip.position = "right",
             ncol=1)+
  scale_fill_brewer(palette = "Set1")+
  scale_color_brewer(palette = "Set1")+
  xlim(177400000,177590000)+
  theme(legend.position = "none",
        strip.text.y.right = element_text(angle =0),
        strip.background.y = element_blank()
  )



quartz()
ggpubr::ggarrange(blast_plot@ggplot,  gene_plot@ggplot,
                  ncol=1,common.legend = TRUE,
                  align = "hv",
                  heights = c(2,1))


#### cds candidates---


blast <- read.table("~/Desktop/B73_jmj_cluster_genomic_vs_candidates_v5.blast", header=FALSE, sep="\t")
colnames(blast) <- blast_cols
blast$qseqid <- "4"
blast$qstart <- blast$qstart  + 177400000+1
blast$qend <- blast$qend  + 177400000+1

ranges <- with(blast,
               IRanges(start=qstart,end = qend)
)



hits <- with(blast, 
             GRanges(
               type="blast_hit",
               seqnames=qseqid,
               ranges=ranges,
               transcript = sseqid,
               score = bitscore,
               identity = identity,
               track=sseqid
             )
)



blast_plot <- hits %>%
  subset(transcript !="EU945014.1") %>%
  ggplot() +
  stat_stepping(aes(group = track, color= identity, fill = identity),
                facets=transcript ~.) +
  xlim(177400000,177590000)+
  facet_wrap(transcript ~.,
             #scales="free_y",
             strip.position = "left",
             ncol=1)+
  scale_fill_continuous(name="% identity") +
  theme(legend.position = "none",
        strip.text.y.right = element_text(angle =0),
        strip.background.y = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()
  )

quartz()
ggpubr::ggarrange(blast_plot@ggplot,  gene_plot@ggplot,
                  ncol=1,common.legend = TRUE,
                  align = "hv")


# cDNA


blast <- read.table("~/Desktop/B73_jmj_cluster_genomic_vs_candidates_v5_cDNA.blast", header=FALSE, sep="\t")
colnames(blast) <- blast_cols
blast$qseqid <- "4"
blast$qstart <- blast$qstart  + 177400000+1
blast$qend <- blast$qend  + 177400000+1

ranges <- with(blast,
               IRanges(start=qstart,end = qend)
)



hits <- with(blast, 
             GRanges(
               type="blast_hit",
               seqnames=qseqid,
               ranges=ranges,
               transcript = sseqid,
               score = bitscore,
               identity = identity,
               track=sseqid
             )
)




transcript_order <- c("Zm00001eb191790_T001","Zm00001d051961_T002",
                      "Zm00001eb191790_T006",
                      "Zm00001eb191790_T013","Zm00001eb191820_T001")
hits$transcript <-factor(hits$transcript,levels=transcript_order)
blast_plot <- NULL
  
blast_plot <- hits %>%
  ggplot() +
  stat_stepping(aes(group = track, color= identity, fill = identity),
                facets=transcript ~.) +
  facet_wrap(transcript ~.,
             strip.position = "right",
             ncol=1)+
  xlim(177400000,177590000) +
  theme(
    legend.position = "left",
    strip.text.y.right = element_text(angle =0),
        strip.background.y = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()
  )

quartz(height =3,width=10)
ggpubr::ggarrange(blast_plot@ggplot,  gene_plot@ggplot,
                  ncol=1,
                  heights = c(2.5,1),
                  align = "hv")

