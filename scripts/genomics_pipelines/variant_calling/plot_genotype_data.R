#' Plot Genotype Data and Genomic Regions
#'
#' This script analyzes VCF genotype data and creates visualizations
#' showing variant distributions across genomic regions, with special
#' focus on the Inv4m inversion and introgression regions.
#'
#' @author Francisco Rodriguez
#' @date 2025-08-06

# Load required libraries -----------------------------------------------
library(vcfR)
library(dplyr)
library(ggplot2)
library(ggtext)
library(rtracklayer)
library(GenomicRanges)
library(IRanges)

# Source common configuration
source("../common_config.R")

# Helper functions -------------------------------------------------------

#' Load VCF file with error handling
#'
#' @param vcf_file Path to VCF file
#' @return vcfR object
load_vcf_data <- function(vcf_file) {
  if (!file.exists(vcf_file)) {
    stop("VCF file not found: ", vcf_file)
  }
  
  cat("Loading VCF file:", vcf_file, "\n")
  
  vcf_data <- tryCatch({
    read.vcfR(vcf_file, verbose = FALSE)
  }, error = function(e) {
    stop("Failed to load VCF file: ", e$message)
  })
  
  cat("Loaded", nrow(vcf_data), "variants\n")
  return(vcf_data)
}

#' Create genomic region visualization
#'
#' @param gene_annotation Gene annotation list from load_gene_annotation()
#' @param vcf_data VCF data (optional)
#' @return List of plots and data
create_genomic_region_plot <- function(gene_annotation, vcf_data = NULL) {
  
  # Extract key information
  genes <- gene_annotation$genes
  regions <- gene_annotation$regions
  gene_sets <- gene_annotation$gene_sets
  
  cat("Genomic regions summary:\n")
  cat("  Inv4m genes:", length(gene_sets$inv4m), "\n")
  cat("  Introgression genes:", length(gene_sets$introgression), "\n")
  cat("  Flanking genes:", length(gene_sets$flanking), "\n")
  
  # Create region annotation data frame
  region_data <- data.frame(
    region = c("inv4m", "introgression"),
    chr = c("4", "4"),
    start = c(
      start(gene_annotation$regions$inv4m),
      start(gene_annotation$regions$introgression)
    ),
    end = c(
      end(gene_annotation$regions$inv4m),
      end(gene_annotation$regions$introgression)
    ),
    gene_count = c(
      length(gene_sets$inv4m),
      length(gene_sets$introgression)
    )
  )
  
  return(list(
    regions = region_data,
    gene_sets = gene_sets,
    genes = genes
  ))
}

#' Analyze variant distribution across genomic regions
#'
#' @param vcf_data VCF data object
#' @param gene_annotation Gene annotation data
#' @return Summary statistics data frame
analyze_variant_distribution <- function(vcf_data, gene_annotation) {
  
  # Extract variant positions
  variant_pos <- data.frame(
    chr = getCHROM(vcf_data),
    pos = getPOS(vcf_data),
    stringsAsFactors = FALSE
  )
  
  # Filter to chromosome 4
  chr4_variants <- variant_pos[variant_pos$chr == "4", ]
  
  # Classify variants by region
  chr4_variants$region <- "outside"
  
  # Inv4m region
  inv4m_start <- start(gene_annotation$regions$inv4m)
  inv4m_end <- end(gene_annotation$regions$inv4m)
  
  inv4m_variants <- chr4_variants$pos >= inv4m_start & 
                   chr4_variants$pos <= inv4m_end
  chr4_variants$region[inv4m_variants] <- "inv4m"
  
  # Introgression region (excluding inv4m)
  introg_start <- start(gene_annotation$regions$introgression)
  introg_end <- end(gene_annotation$regions$introgression)
  
  flanking_variants <- (chr4_variants$pos >= introg_start & 
                       chr4_variants$pos <= introg_end) & 
                      !inv4m_variants
  chr4_variants$region[flanking_variants] <- "flanking"
  
  # Generate summary
  variant_summary <- chr4_variants %>%
    group_by(region) %>%
    summarise(
      variant_count = n(),
      .groups = "drop"
    )
  
  return(list(
    variants = chr4_variants,
    summary = variant_summary
  ))
}

# Main analysis pipeline -------------------------------------------------

#' Main function to plot genotype data
#'
#' @param vcf_file Path to VCF file (optional)
#' @param gff_file Path to GFF file (uses default if NULL)
#' @param output_prefix Prefix for output files
#' @export
plot_genotype_data <- function(vcf_file = NULL, gff_file = NULL, 
                              output_prefix = "genotype_analysis") {
  
  cat("Starting genotype data analysis...\n")
  
  # Load gene annotation
  cat("Loading gene annotation...\n")
  if (is.null(gff_file)) {
    gff_file <- REFERENCE_FILES$gff
  }
  
  gene_annotation <- load_gene_annotation(gff_file)
  
  # Create genomic region visualization
  cat("Creating genomic region analysis...\n")
  region_analysis <- create_genomic_region_plot(gene_annotation)
  
  # Save region information
  write.csv(
    region_analysis$regions,
    file = file.path(OUTPUT_DIR, paste0(output_prefix, "_regions.csv")),
    row.names = FALSE
  )
  
  # Save gene sets
  gene_sets_df <- data.frame(
    gene_id = c(
      gene_annotation$gene_sets$inv4m,
      gene_annotation$gene_sets$flanking
    ),
    region = c(
      rep("inv4m", length(gene_annotation$gene_sets$inv4m)),
      rep("flanking", length(gene_annotation$gene_sets$flanking))
    ),
    stringsAsFactors = FALSE
  )
  
  write.csv(
    gene_sets_df,
    file = file.path(OUTPUT_DIR, paste0(output_prefix, "_gene_sets.csv")),
    row.names = FALSE
  )
  
  # Analyze VCF data if provided
  if (!is.null(vcf_file)) {
    cat("Analyzing VCF data...\n")
    vcf_data <- load_vcf_data(vcf_file)
    
    variant_analysis <- analyze_variant_distribution(vcf_data, gene_annotation)
    
    # Save variant analysis
    write.csv(
      variant_analysis$summary,
      file = file.path(OUTPUT_DIR, paste0(output_prefix, "_variant_summary.csv")),
      row.names = FALSE
    )
    
    cat("Variant distribution summary:\n")
    print(variant_analysis$summary)
  } else {
    cat("No VCF file provided, skipping variant analysis\n")
    variant_analysis <- NULL
  }
  
  cat("Analysis complete!\n")
  cat("Results saved with prefix:", output_prefix, "\n")
  
  return(list(
    gene_annotation = gene_annotation,
    region_analysis = region_analysis,
    variant_analysis = variant_analysis
  ))
}

# Execute analysis if script is run directly
if (!interactive()) {
  results <- plot_genotype_data()
}


flanking_introgression_gene_ids <- shared_introgression_gene_ids[!(shared_introgression_gene_ids %in% inv4m_gene_ids)]

length(flanking_introgression_gene_ids)




inversion_limit <-round(c(inv4m_start,inv4m_end)/1000000,2)
inversion_limit[2] - inversion_limit[1]
# 15.2
introgression_limit <- round(c(157012149,195900523)/1000000,2) # checked manually

round(introgression_limit[2]-introgression_limit[1],1)
# 38.9


# Read VCF -----

vcf <- read.vcfR("~/Desktop/inv4m_PSU_imputed.vcf")
# Number of taxa: 13
# Number of sites: 19861
# Imputed genotypes by KNN imputation

gt <- extract.gt(vcf, element = c('GT'))

# gt[is.na(gt)] <-"0/0"
gt[gt=="0/0"] <-0
gt[gt=="0/1"] <-1
gt[gt=="1/0"] <-1
gt[gt=="1/1"] <-2
#gt[is.na(gt)] <-0



ngt <- matrix(as.integer(gt), ncol= ncol(gt))
dim(ngt)


colnames(ngt) <- colnames(gt)
rownames(ngt) <- rownames(gt)
# by_size <- c("3056","4083","3095","3113","4113","3083",
#              "4107","4122","3092","3080","3059","3131",
#              "4125","4080", "4131")
total <- ncol(ngt) * nrow(ngt)
table(ngt)/total
sum(is.na(ngt))/total
by_size <- c("3095","3113","4113","3083","4107","4122",
             "3092","3080","3059","3131","4125","4080", "4131")


inv4m <- c(rep("CTRL",6),rep("INV4",7))

# new_order <- data.frame(line= by_size, group = inv4m) %>%
#   mutate(group = factor(group)) %>%
#   arrange(group,line)
# by_rowid <-new_order$line

new_order <- c(sort(c("3095","3113","4113","3083","4107","4122")),
               sort(c("3092","3080","3059","3131","4125","4080", "4131")))

by_rowid <-new_order

ngt <- ngt[,by_rowid] %>% t()

colnames(ngt)[grep("S4_181",rownames(gt))]

S4_181558683_cor <- cor(ngt[,"S4_181558683" ,drop = FALSE] ,
    ngt, use= "pairwise.complete.obs") %>% t()
class(ngt)

apply(ngt,2,FUN=sum)

cor_test_results <- apply(ngt,2, 
  FUN=function(x){
  cor.test(x, y = ngt[,"S4_181558683"] ,
           method = "spearman")$p.value
  })
adjusted <- p.adjust(cor_test_results,method="fdr")
hist(adjusted)
significant <- names(asjusted)[which(ajusted <0.05)]   

length(significant)

inv4m_cor <-cbind(
  data.frame(marker = rownames(S4_181558683_cor)),
  S4_181558683_cor
) %>%
  rename(S4_181558683="r") %>%
  tidyr::separate(marker,into = c("CHR","BP"),sep = "_") %>%
  mutate(BP=as.numeric(BP)) %>%
  mutate(CHR=gsub("S","",CHR) %>% as.factor) %>%
  mutate(CHR=factor(CHR, levels=1:10 %>% as.character())) %>%
  mutate(p=cor_test_results) %>%
  mutate(FDR=p.adjust(p,method="fdr")) %>%
  tibble::rownames_to_column("ID")

inv4m_cor$R2 <-inv4m_cor$r^2


inv4m_cor$is_significant <- inv4m_cor$FDR < 0.005
sum(inv4m_cor$is_significant, na.rm = TRUE)

with(inv4m_cor,
     table(CHR)
)

8892/19861
1219/19861


with(inv4m_cor,
 table(CHR, is_significant)
)




# inv4m_cor <- inv4m_cor[!is.na(inv4m_cor$r),]

with(inv4m_cor, 
     hist(r))

correlated_SNP <- with(inv4m_cor %>% filter(is_significant),
 GRanges(seqnames=CHR,
  ID=ID,
  ranges=IRanges(start=BP, end=BP))
)

length(correlated_SNP)

7547/7683

7683-361

1-361/7683

outside_correlated <- inv4m_cor %>% 
 filter(is_significant) %>% 
 filter(!(CHR==4 & BP >= introgression_start & BP <= introgression_end))

nrow(outside_correlated)

anticorrelated_start <- 154283637 
anticorrelated_end <- 156985577  

anticorrelated <- correlated_SNP %>% 
  subset( seqnames==4 & start >= anticorrelated_start & end <= anticorrelated_end)

inside_correlated <- correlated_SNP %>% 
  subset( seqnames==4 & start >= introgression_start & end <= inv4m_end)

       




361-207

154283637 - 156985577

myGFF <- "/Users/fvrodriguez/ref/zea/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.59.chr.gff3"
B73<- rtracklayer::import(myGFF)  %>%
  subset(type=="mRNA"  & seqnames %in% 1:10)


transcripts <- as.data.frame(B73) 
colnames(transcripts)
transcripts$gene_id
transcripts$gene_id <- gsub("transcript:","",transcripts$ID)
transcripts$gene_id <- gsub("_.*","",transcripts$gene_id, perl =TRUE)

olap <-findOverlaps(outside_correlated,B73, ignore.strand=TRUE)
queryHits(olap) %>% sort() %>% unique() %>% length()
correlated_genes <- transcripts$gene_id[subjectHits(olap)] %>% sort() %>% unique()
length(correlated_genes)
new_label[correlated_genes]

olap <-findOverlaps(anticorrelated,B73, ignore.strand=TRUE)
queryHits(olap) %>% sort() %>% unique() %>% length()
correlated_genes <- transcripts$gene_id[subjectHits(olap)] %>% sort() %>% unique()
length(correlated_genes)
new_label[correlated_genes]


olap <-findOverlaps(inside_correlated ,B73, ignore.strand=TRUE)
queryHits(olap) %>% sort() %>% unique() %>% length()
correlated_genes <- transcripts$gene_id[subjectHits(olap)] %>% sort() %>% unique()
length(correlated_genes)
#new_label[correlated_genes]


# check with effects table
effects[effects$gene %in%correlated_genes,] %>% 
  select(gene:logFC,adj.P.Val,) %>%
  arrange(adj.P.Val) %>% print(n=100)

nrow(inv4m_cor)
quartz()
with(inv4m_cor, 
     hist(R2, main="R2 with inv4m tagging marker.\n 19861 markers in total", ylim=c(0,8000)))

with(inv4m_cor, 
     table(CHR))
     
with(inv4m_cor, 
     table(CHR) %>% prop.table() %>% sort() %>% round(3))

# 98% of perfectly correlated markers to inv4m are in CHR04
with(inv4m_cor %>% filter(R2==1), 
     table(CHR)%>% sort())

inv4m_cor$unlinked <-  abs(inv4m_cor$R2) <0.1
inv4m_cor$linked <- abs(inv4m_cor$R2) >0.5
inv4m_cor$linkage <- case_when(
  inv4m_cor$unlinked ~"unlinked",
  inv4m_cor$linked ~"linked",
  .default = "in_between"
)

with(inv4m_cor, 
     table(CHR,linkage))

with(inv4m_cor %>% filter(R2==1), 
     table(CHR))
     
inv4_distro <- data.frame(
  CHR= factor(1:10,levels= 1:10),
  all =with(inv4m_cor, 
      table(CHR)) %>% as.numeric(),
  inv4m = with(inv4m_cor %>% filter(is_significant), 
    table(CHR))  %>% as.numeric()
) %>%
  mutate(
    all_pct = 100*all/sum(all),
    inv4m_pct = 100*inv4m/sum(inv4m),
  )



p0 <-inv4_distro %>%
  tidyr::pivot_longer(cols = c("all_pct","inv4m_pct"),names_to = "marker_type",values_to = "pct") %>%
  ggplot(aes(y=pct,x=CHR, group=marker_type, fill=marker_type) ) +
  ylab("%") +
  #ggtitle("Distribution of Genotyped SNPs") +
  geom_bar(stat = "identity",position = "dodge")+
  scale_fill_manual(labels = c("all SNPs, n = 19861", "inv4m correlated, n = 7683"), values= c("gold","purple4")) +
  # scales::label_percent() +
  ggpubr::theme_classic2(base_size = 20)+
  theme(legend.position= c(0.7,0.8),
        legend.title = element_blank(),
        axis.title.x = element_blank())

# plot FDR
manhattan <- inv4m_cor %>% 
  
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(inv4m_cor, ., by=c("CHR"="CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot)

axisdf = manhattan %>%
  group_by(CHR) %>%
  summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

manhattan$FDR[manhattan$FDR==0] <- 1e-8

p1 <- manhattan %>%
ggplot(aes(x=BPcum, y=-log10(FDR))) +
  # Show all points
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("gold", "purple4"), 22 )) +
  geom_hline(yintercept= -log10(0.05), col= "red")+
 
  # custom X axis:
  scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
  coord_cartesian(ylim = c(0,10))+
  # scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  
  # Custom the theme:
  theme_bw(base_size = 20) +
  theme( 
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

p2 <- manhattan %>%
  ggplot(aes(x=BPcum, y=r)) +
  # Show all points
  geom_point( aes(color=as.factor(CHR)), size=1.3) +
  geom_point(data = manhattan %>% filter(FDR >=0.05),
             aes(color=as.factor(CHR)),
             fill= "white", shape=21, size=1.3) +
  scale_color_manual(values = rep(c("gold", "purple4"), 22 )) +
  
  # custom X axis:
  xlab("Chromosome") +
  scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
  coord_cartesian(ylim = c(-1.1,1.1))+
  # scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  
  # Custom the theme:
  theme_bw(base_size = 20) +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

# outside target




outside_cor <- inv4m_cor  %>%
  filter(!(CHR==4 & BP >= introgression_start & BP <=introgression_end))

quartz()
outside_cor %>%
  filter(CHR==4) %>%
  ggplot(aes(x=BP,y=r)) + geom_point()

outside_cor %>%
  filter(CHR==4 & R2 ==1) 


outside_distro <- data.frame(
  CHR= factor(1:10,levels= 1:10),
  all =with(outside_cor, 
            table(CHR) %>% prop.table() %>% as.numeric()%>% round(3)),
  #table(CHR) %>% as.numeric()),
  inv4m = with(outside_cor %>% filter(R2==1), 
               table(CHR) %>% prop.table()%>% as.numeric() %>% round(3))
  #table(CHR) %>% as.numeric())
)
outside_cor %>%
  filter(CHR==4) %>% pull(linkage) %>% table()
quartz()
outside_distro %>%
  tidyr::pivot_longer(cols = c("all","inv4m"),names_to = "marker_type",values_to = "fraction") %>%
  ggplot(aes(y=fraction,x=CHR, group=marker_type, fill=marker_type) )+
  geom_bar(stat = "identity",position = "dodge")+
  scale_fill_manual(values= c("gold","purple4")) +
  # scales::label_percent() +
  ggpubr::theme_classic2(base_size = 25)+
  theme(legend.position = "top")



quartz()
barplot(inv4_distro)


m  <- t(ngt)

get_genotype_runs <- function(x){
  result <-data.frame()
  
  for(ind_idx in by_rowid){
  result_ind <-data.frame()
  chr_levels  <-paste0("S",1:10)
  
    for(chr_idx in 1:10){
    chr_long <-NULL 
    chr_prefix <- paste0("S",chr_idx,"_")
    chr_name <- paste0("S",chr_idx)
    chr <- x[grepl(chr_prefix, rownames(x), perl = TRUE),]
    runs <- rle(chr[,ind_idx])
    chr_long <- data.frame(
      ind = ind_idx,
      chr=chr_name,
      pos= as.integer(gsub(chr_prefix,"",names(runs$values))),
      gt = runs$values) %>%
      group_by(ind, chr) %>%
      mutate(lag= lag(pos,default = 0)) %>%
      mutate(span= (pos-lag)/1e6) %>%
      mutate(chr= factor(chr, levels=chr_levels)) %>%
      mutate(ind=factor(ind, levels= by_rowid)) %>%
      mutate(missing = is.na(gt)) %>%
      mutate(missing_length = missing*span) %>%
      mutate(B73 = as.numeric(gt == 0)) %>%
      mutate(B73_length = B73*span) %>%
      mutate(Het = as.numeric(gt == 1)) %>%
      mutate(Het_length = Het*span) %>%    
      mutate(Mi21 = as.numeric(gt == 2)) %>%
      mutate(Mi21_length = Mi21*span)
    #chr_long$Highland[is.na(chr_long$Highland)] <- 0
    
    # sm <- gam(Highland~s(pos),family="binomial",data=chr_long)
    # chr_long$Highland_est <- predict(sm,newdata=chr_long,type="response")
    result_ind <- rbind(result_ind,chr_long)
    }
  
  result <- rbind(result,result_ind)

  }
  result
}

gt_runs <- get_genotype_runs(m)

with(gt_runs,
table(missing)
)

gt_runs$group <- case_when(gt_runs$ind %in% c("3092","3080","3059","3131","4125","4080", "4131") ~ "Inv4m",
                                  .default = "CTRL")

chr_length <- gt_runs %>%
  group_by(chr) %>%
  summarise(chr_length=max(pos)/1e6)
genome_size <- sum(chr_length$chr_length)

by_length <- gt_runs %>%
  inner_join(chr_length) %>%
  ungroup() %>%
  group_by(ind) %>%
  summarise(B73=sum(B73_length, na.rm = TRUE),
            Het=sum(Het_length, na.rm = TRUE),
            Mi21=sum(Mi21_length, na.rm = TRUE),
            missing=sum(missing_length, na.rm = TRUE),
            group = group[1]) %>%
  mutate(group = factor(group, levels= c("CTRL","Inv4m"))) %>%
  tidyr::pivot_longer(cols= c("B73","Het","Mi21", "missing"),names_to = "Genotype", values_to = "length") %>%
  mutate(pct = 100*length/genome_size)

bg_mean <- by_length %>%
  ungroup() %>%
  group_by(Genotype) %>%
  summarise_at(.vars = c("length","pct"), .funs = mean)

 by_length %>%
  ungroup() %>%
  group_by(group,Genotype) %>%
  summarise_at(.vars = c("length","pct"), .funs = mean)


quartz()
p3 <- by_length  %>%
  ggplot(aes(x= group, y=pct))+
  geom_hline(data = bg_mean, aes(yintercept = pct) , linetype = 2) +
  geom_boxplot(outlier.colour = "white", width=0.25) +
  ggbeeswarm::geom_quasirandom(size=2) +
  facet_wrap(.~Genotype, scales = "free_y", ncol=4) +
  xlab("Introgression group") +
  ylab("% of genome length") +
  ggpubr::theme_classic2(base_size = 20) +
  theme(strip.background = element_blank())

quartz()
by_length  %>%
  ggplot(aes(x= group, y=length))+
  geom_boxplot(outlier.colour = "white", width=0.25) +
  ggbeeswarm::geom_quasirandom(size=2) +
  facet_wrap(.~Genotype, scales = "free_y", ncol=4) +
  xlab("Introgression group") +
  ylab("Length [Mb]") +
  ggpubr::theme_classic2(base_size = 20) +
  theme(strip.background = element_blank())

out <- ggpubr::ggarrange(p0,p1,p2,p3,nrow=4, align = "hv", heights = c(2,2,3,3))
quartz(height=12, width=10)
out
ggsave(out, file="~/Desktop/SNP_distribution.svg", height=12, width=10)


# Mendelian segregation calculations
# AA Aa aa
#AA
#Aa
#aa
AA <- c(1, 1/2, 0,
        0, 1/2, 1,
        0,   0, 0) %>% 
  matrix(nrow = 3, byrow = TRUE)

aa <- c(0,   0, 0,
        1, 1/2, 0,
        0, 1/2, 1) %>% 
  matrix(nrow = 3, byrow = TRUE)

Aa <- c(1/2, 1/4,   0,
        1/2, 1/2, 1/2,
        0,   1/4, 1/2) %>% 
  matrix(nrow = 3, byrow = TRUE)

S  <- c(1, 1/4, 0,
        0, 1/2, 0,
        0, 1/4, 1) %>% 
  matrix(nrow = 3, byrow = TRUE)

f1 <- c(0,1,0)

bc2  <- AA %*% AA %*% f1

bc2s3 <- (S %*% S %*% S %*% bc2)[,1]

bc6s2 <- (S %*% S %*% S %*% S %*% S %*% S %*% bc2)[,1]



# For counting runs of fixed Mi21
# m[which(colSums(ngt) !=26),] <- 0

outside_cor

to_run <-m

# to_run[which(colSums(ngt) !=26),] <- 0

# get outside correlated run sizes
outside_correlated %>%
  filter(r<0) %>% tibble() %>% print(n=300)
#2.7 Mb
154283637 - 156985570
to_run[!rownames(to_run) %in% outside_correlated$ID,] <- 0
to_run[rownames(to_run) %in%  outide_correlated$ID,] <- 2

# get inside correalted run sizes
inside_correlated <- inv4m_cor %>% 
  filter(is_significant) %>% 
  filter(CHR==4 & BP >= introgression_start & BP <= introgression_end)
inside_correlated %>% tail()

to_run[!rownames(to_run) %in% inside_correlated$ID,] <- 0
to_run[rownames(to_run) %in% inside_correlated$ID,] <- 2

gt_runs <- get_genotype_runs(to_run)

gt_runs$group <- case_when(gt_runs$ind %in% c("3092","3080","3059","3131","4125","4080", "4131") ~ "CTRL",
                                       .default = "Inv4m")

chr_length <- gt_runs %>%
  group_by(chr) %>%
  summarise(chr_length=max(pos)/1e6)
genome_size <- sum(chr_length$chr_length)

by_length <- gt_runs %>%
  inner_join(chr_length) %>%
  ungroup() %>%
  group_by(ind) %>%
  summarise(B73=sum(B73_length, na.rm = TRUE),
            Het=sum(Het_length, na.rm = TRUE),
            Mi21=sum(Mi21_length, na.rm = TRUE),
            missing=sum(missing_length, na.rm = TRUE),
            group = group[1]) %>%
  mutate(group = factor(group, levels= c("CTRL","Inv4m"))) %>%
  tidyr::pivot_longer(cols= c("B73","Het","Mi21", "missing"),names_to = "Genotype", values_to = "length") %>%
  mutate(pct = 100*length/genome_size)

by_length %>%
  ungroup() %>%
  group_by(Genotype) %>%
  summarise_at(.vars = c("length","pct"), .funs = mean)

by_length %>%
  ungroup() %>%
  group_by(group,Genotype) %>%
  summarise_at(.vars = c("length","pct"), .funs = mean)


# Plotting Genotype runs #######################################################

#                            inv4m.  inv4m  zcn8
# dummy <- data.frame(chr = c("chr4", "chr4","chr8"), 
#                     y = c(173.9, 187.9,126.7))


dummy <- data.frame(chr = c("4", "4"), 
                    x = inversion_limit)

# ann_text <-data.frame(chr = c("4", "4"),  
#                        ind= c("3113","3059"),span = c(-5, -5))

# triangle <-data.frame( pos = c(0,NA),
#                        span =c(181.6,NA),
#                        ind=c(new_order[length(new_order)],NA),
#                        chr = c("4","4"),
#                        label = factor(c("Selection<br>marker","<i style='color:#43137D'>*inv4m*</i>"), 
#                                       levels= c("Selection<br>marker","<i style='color:#43137D'>*inv4m*</i>")))


triangle <-data.frame( pos = 0,
                       span =181.6,
                       ind= new_order[length(new_order)],
                       chr = "4",
                       label = factor("Selection marker<br>group"))


breaks <- c( 0, 1.5,2.5)
names(breaks ) <- c("B73/B73", "B73/Mi21","Mi21/Mi21")

selected3 <- c("3131","4131","4125")

bracket1 <- data.frame(x= c(inversion_limit[1],inversion_limit[1],inversion_limit[2],inversion_limit[2]),
                       y = c(0,-0.3,-0.3,0),
                       color="#0D0887FF")

bracket2 <- data.frame(x= c(introgression_limit[1],introgression_limit[1],inversion_limit[1],inversion_limit[1]),
                       y = c(-0.4,-0.7,-0.7,-0.4),
                       color="#CC4678FF")
bracket3 <- data.frame(x= c(inversion_limit[2],inversion_limit[2],introgression_limit[2],introgression_limit[2]),
                       y = c(-0.4,-0.7,-0.7,-0.4),
                       color="#CC4678FF")



pal

levels(factor(result$gt))

gt3_plot <- result  %>%
  filter(chr=="S4") %>%
  filter(ind %in% selected3 ) %>%
  droplevels() %>%
  ungroup() %>%
  mutate(ind =factor(ind, levels= selected3)) %>%
  mutate(chr=gsub("S","",chr)) %>%
  group_by(chr,ind) %>%
  ggplot( aes(x = span, y = ind, fill = gt)) + 
  ylab("Line")+
  geom_col(width = 1) +
  scale_y_discrete(expand=c(0,0))+
  scale_fill_viridis_b(breaks = breaks,
                       limits = c(-0.5, 2.5),
                       show.limits = TRUE,
                       na.value = "gray85", direction = -1)   +
  guides(fill = guide_legend(title= "Genotype",ncol = 3)) +
  coord_cartesian(xlim = c(150, 200))+
  annotate('point', x = 181.6 , y = 4,
           shape = 25, size = 7, color = "#0D0887FF", fill = "#0D0887FF") +
  annotate('text', x = c(inversion_limit,introgression_limit) , y = 4,
          label = rep("|",4), size = 7) +
  annotate('text', x = c(inversion_limit,introgression_limit) , y = 0,
           label = rep("|",4), size = 7) +
  # annotate('text', x = c(inversion_limit,introgression_limit) , y = -1.5,
  #          label = rep("|",4), size = 10, color="white")+
  ggpubr::theme_classic2(base_size = 25) +
  theme(legend.position = "top",
        #legend.title = element_blank(),
        legend.byrow = TRUE,
        legend.text = element_text(size = 20),
        legend.box.spacing = unit(0, "pt"),# The spacing between the plotting area and the legend box (unit)
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(color ="white"),
        axis.title.y = element_text(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        plot.background = element_rect(fill='transparent', color=NA)
  ) 

quartz()
gt3_plot 

save(gt3_plot, file = "~/Desktop/gt3_plot.rda")


p0 <-result  %>%
  filter(chr=="S4") %>%
  #filter(pos > 150e6 & pos <200e6) %>%
  ungroup() %>%
  group_by(chr,ind) %>%
  mutate(chr=gsub("S","",chr)) %>%
  ggplot( aes(x = span, y = ind, fill = gt)) + 
  ggtitle("Genotype",) +
  xlab("Position [Mb]")+
  ylab("Introgression Line")+
  geom_col(width = 1) +
  scale_fill_viridis_b(breaks = breaks,
                       limits = c(-0.5, 2.5),
                       show.limits = TRUE,
                       na.value = "gray85", direction = -1)   +
  labs(fill="Genotype")+
  facet_wrap(factor(chr, 1:10)~., 
             strip.position="left",
             ncol=1) +
 # geom_vline(data = dummy, aes(xinterceHighland = x), col = "#43137D", linetype="dashed", size=1) +
  #annotate('point', x = 181.6 , y = 14.7,
  #         shape = 25, size = 10, color = '#43137D', fill = '#43137D') +
  ggpubr::theme_classic2(base_size = 4)+
  coord_cartesian(xlim = c(150, 200))+
  theme(legend.position = "none",
        strip.background = element_blank(),
        # strip.text.y.left = element_text(size = 30, angle =0, color ="white",hjust =-1, vjust=-1),
        plot.title = element_text(size=30, hjust=1),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=30, color= "white"),
        axis.line.y = element_blank(),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 25, hjust= 0)
  )


quartz()


gt_plot <- result  %>%
  ungroup() %>%
  group_by(chr,ind) %>%
  mutate(chr=gsub("S","",chr)) %>%
  ggplot( aes(y = ind, x = span, fill = gt)) + 
  xlab("Position [Mb]")+
  ylab("Chromosome")+
  geom_col(width = 1) +
  scale_fill_viridis_b(breaks = breaks,
                       limits = c(-0.5, 2.5),
                       show.limits = TRUE,
                       na.value = "gray85", direction = -1)   +
  geom_point(data=triangle,
             mapping = aes(shape=label),
             position = position_nudge(y=2),
             size = 1, color = '#43137D', fill = '#43137D') +
  scale_shape_manual(values=25)+
  labs(fill="Genotype")+
  guides( ncol=2, 
    fill = guide_colorsteps( order = 1,
                             label.position = "left",
                             label.hjust = 1,
                             theme = theme(
                               legend.title=element_text(size=8, hjust= 0),
                               legend.text=element_text(size=8, vjust= 1.5),
                               legend.key.width  = unit(1, "lines"),
                               legend.key.height = unit(2.7, "lines")
                             )),
    shape = guide_legend( order = 2,
                          override.aes = list(size = 4),
                          label.position = "left",
                          label.hjust = 1,
                          theme = theme(
                            # legend.title=element_text(size=8, hjust= 0.88, col="white"),
                            legend.title=element_blank(),
                            legend.text=element_markdown(size=8),
                            legend.key.width  = unit(1, "lines"),
                          ))
  ) +
  facet_wrap(factor(chr, 1:10)~., 
             strip.position="left",
             ncol=1) +
  # geom_vline(data = dummy, aes(xinterceHighland = x), color= "#43137D", linetype="dashed") +
  coord_cartesian(xlim=c(0,450), expand = FALSE) +
  #coord_flip(expand = FALSE) +
  ggpubr::theme_classic2(base_size = 4)+
  theme(# legend.position = c(0.92,0.5), # Chr4
    legend.position = c(0.9,0.85), # Chr3
    #legend.direction = "vertical", legend.box = "horizontal",
    # legend.position = c(0.92,0.15), # inside bottom left
    legend.spacing.y = unit(0.1, "lines"),
    legend.box.just = "right",
    strip.background = element_blank(),
    strip.text.y.left = element_text(size = 12, angle =0),
    axis.text = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    # axis.title.y = element_blank(),
    axis.title.x = element_text(size = 14),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank()
  )


inlay <- p0 +
  ggpubr::theme_classic2(base_size = 8) +
  xlab("Chromosome 4") +
  annotate('point', x = 181.6 , y = 14.4,
          shape = 25, size = 4, color = '#43137D', fill = '#43137D') +
  annotate('text', x = c(inversion_limit,introgression_limit) , y = 14.4,
            label = rep("|",4), size = 4, color = '#43137D') +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_blank(),
        # axis.title.y = element_blank(),
        #axis.title. = element_text(size = 14),
        #axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        plot.background = element_rect(fill='transparent', color=NA)
        )


library(cowplot)
plot.with.inset <-
  ggdraw() +
  draw_plot(gt_plot) +
  draw_plot(inlay ,
            x = 0.45, y = 0.03, width = .6, height = 0.6,scale = 0.8)

png(file= "~/Desktop/Inv4mPSU_genotye.png",
     width = 7, height = 4,units = "in",res = 300)
plot.with.inset 
dev.off()







