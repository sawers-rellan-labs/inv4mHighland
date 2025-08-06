
urls <- c(
  B73 = file.path(genomeRepo,"Zm-B73-REFERENCE-NAM-5.0.fa"),
  parviglumis = file.path(genomeRepo,"Zv-TIL01-REFERENCE-PanAnd-1.0.fa"),
  mexicana = file.path(genomeRepo,"Zx-TIL18-REFERENCE-PanAnd-1.0.fa")
) 



translatedCDS  <- c(
  B73 = file.path(genomeRepo,"B73","Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.protein.fa"),
  parviglumis = file.path(genomeRepo,"parviglumis","Zv-TIL01-REFERENCE-PanAnd-1.0_Zv00001aa.1.protein.fa"),
  mexicana = file.path(genomeRepo,"mexicana","Zx-TIL18-REFERENCE-PanAnd-1.0_Zx00002aa.1.protein.fa")
) 

geneGff <- c(
  B73 = file.path(genomeRepo,"B73","Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3"),
  parviglumis = file.path(genomeRepo,"parviglumis","Zv-TIL01-REFERENCE-PanAnd-1.0_Zv00001aa.1.gff3"),
  mexicana = file.path(genomeRepo,"mexicana","Zx-TIL18-REFERENCE-PanAnd-1.0_Zx00002aa.1.gff3")
) 




library(GENESPACE)

genomeRepo <- "/rsstu/users/r/rrellan/DOE_CAREER/inv4m_synteny/gff"

wd <- "/rsstu/users/r/rrellan/DOE_CAREER/inv4m_synteny"

path2mcscanx <- "/usr/local/usrapps/maize/MCScanX/"

setwd(wd)

genomes2run<-c("B73","parviglumis","mexicana")

writeDirs <-  file.path(wd,genomes2run)

parsedPaths <- parse_annotations(
  rawGenomeRepo = genomeRepo,
  genomeDirs = genomes2run,
  genomeIDs = genomes2run,
  genespaceWd = wd,
  gffIdColumn = "ID",
  headerEntryIndex = 1,
  headerSep = " ",
  gffString = "gff3",
  faString = "fa")

path2mcscanx <- "/usr/local/usrapps/maize/MCScanX/" 
gpar <- init_genespace(
  wd = wd,
  ploidy = 1, 
  path2mcscanx = path2mcscanx)

out <- run_genespace(gpar, overwrite = T)



ripd <- plot_riparian(
  gsParam = out,
  refGenome = "B73", 
  useRegions = FALSE,
  pdfFile = file.path(wd, "riparian_plot.pdf"))

hits <- read_allBlast(
  filepath = file.path(out$paths$syntenicHits, 
                       "B73_vs_mexicana.allBlast.txt.gz"))
ggdotplot(hits = hits, 
          type = "all",
          outDir = wd,
          inversionColor = "green",
          verbose = FALSE)


roi <- data.frame(
  genome = c("B73","mexicana"),
  chr = c("chr4", "chr4"),
  start = c(170e6, 180e6),
  end = c(190e6, 200e6))

qreturn <- query_hits(gsParam = out, bed = roi, synOnly = TRUE)

xvx <- subset(qreturn[["B73, chr4: 1.7e+08-1.9e+08"]], genome2 == "mexicana")
gghits(xvx, useOrder = F)
dev.copy2pdf(file = file.path(wd, "inv4m_zoom.pdf"))

ripDat <- plot_riparian(
  gsParam = out, 
  highlightBed = roi,
  refGenome = "B73",
  inversionColor = "green",
  verbose = FALSE
)
dev.copy2pdf(file = file.path(wd, "invertions.pdf"))

ripDat <- plot_riparian(
  gsParam = out, 
  highlightBed = roi, 
  backgroundColor = NULL, 
  genomeIDs = c("B73", "mexicana"),
  inversionColor = "green",
  refGenome = "B73")
dev.copy2pdf(file = file.path(wd, "riparian_inv4m.pdf"))


roi <- data.frame(
  genome = c("B73","mexicana"),
  chr = c("chr4", "chr4"),
  start = c(0, 0),
  end = c(Inf, Inf))

ripDat <- plot_riparian(
  gsParam = out, 
  highlightBed = roi, 
  backgroundColor = NULL, 
  genomeIDs = c("B73", "mexicana"),
  inversionColor = "green",
  refGenome = "B73")
dev.copy2pdf(file = file.path(wd, "riparian_chr4.pdf"))
