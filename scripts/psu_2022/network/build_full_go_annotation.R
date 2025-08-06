# Methods description
#
# Once we had sets of differentially expressed genes for the three predictors (leaf, -P, \invfour) and two types of gene expression response (upregulated and downregulated), we proceeded to annotate them with gene ontology terms and KEGG pathways using \texttt{ClusterProfiler} \cite{yu2012, zicola2024}.  
# We started with the B73 NAM v5 gene ontology annotation from \cite{fattel2024} and added GO terms for each intermediate node in the gene ontology tree using the \texttt{ClusterProfiler} function \texttt{buildGOmap}. 
# Then we conducted gene over-representation analysis with the function \texttt{compareCluster}, using as universe/background the set of 24011 genes detected in at least one good quality leaf RNAseq library. 
# This function calculates the hypergeometric test for overrepresented ontology terms in the specified gene set and returns raw, and FDR-adjusted p-values.

library(clusterProfiler )
# from https://github.com/johanzi/GOMAP_maize_B73_NAM5
goa <- readRDS(file="GOMAP_maize_B73_NAM5/TERM2GENE_Fattel_2024.rds")
# Send this over to the server
gomap <-buildGOmap(goa)

saveRDS(gomap, file="TERM2GENE_full.rds")
