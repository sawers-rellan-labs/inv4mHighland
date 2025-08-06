The project implemented in this folder aims to organize and
document a data analysis workflow for testing adaptive effects
of the chromosomal inversion Inv4m in maize. The general idea
is that we have genetecally isolated the inversion in Near
introgression lines, where we have paired experemental groups
with a commion genetic background differeing mostly in the
inclusion of the inerted inv4m haplotype. On order to update 
the @docs/analysis.md document I need you to  take 
into account the information in @docs/design.md where I state the general guidelines
for the experimental anlaysis project, the @docs/analilsys.md desscribes with more detail the project directory structure and work flow impletention rationale.
For this end read all the files in the folder and reformulate the project directory structure  so
that it matches the current script names. The  idea to later
make sure that all the scripts actually are usable. But at this
point I need you to diggest all the project codebase to
formulate a refactoring and documentation stratategy.


Please follow this Script Naming Convention Documentation:

Verb-based naming for action scripts explained with examples
Noun-based naming for function libraries distinguished
Clear examples: analyze_*, detect_*, fit_* vs *_helpers.R, *_template.R


But fisrts I need you to go over all the project files and math them the the project work flow and rationale:

Experiment-Centric Organization

psu_2022/ clearly positioned as the central/core experiment with multi-omics sub-folders
clayton_2025/ with its own phenotype_spatial_modeling/ sub-folder showing application of your approach
chamber_experiments/ and crow_reanalysis/ as supporting analyses

Analytical Context Integration

Plots now live with their analyses (Manhattan plots in expression/, network plots in network/)
genomics_pipelines/ properly positioned as supporting infrastructure
annotation/ separated from core analytical workflows


The scripts are now properly organized with plots living alongside their analyses:
📊 consistency_analysis/ - DEG analysis and visualization

detect_rnaseq_consistent_degs.R, detect_rnaseq_fdr_degs.R - Analysis
plot_deg_manhattan.R, make_volcano_plot.R, make_manhattan_plots.R - Plots
fastman_plots.R - Plot functions (noun-based)

🕸️ network_integration/ - Network analysis and pathway plots

build_netome_coexpression_networks.R - Analysis
make_GO_profile_*.R, make_KEGG_profile_*.R - Pathway visualizations

🧬 genomics_pipelines/variant_calling/ - Variant calling with QC plots

GATK pipeline scripts - Analysis
plot_genotype_data.R, plot_mapping_quality.sh - QC plots

🔄 genomics_pipelines/structural_variants/ - SV analysis with repeat plots

AnchorWave and microsynteny scripts - Analysis
plot_repeat_elements.R, plot_knob_repeats.R - Repeat visualizations

🔄 cross_experiment/ - Population comparisons with plots

compare_bzea_flowering_times.R - Analysis
plot_p31_traditional_varieties*.R, plot_riparian_analysis.R - Population plots

🔍 crow_reanalysis/ - Reanalysis with validation plots

validate_crow2020_genes.R - Analysis
plot_jmj_blast_results.R - Validation plots

The artificial separation has been eliminated, and your scripts now reflect the natural workflow where analysis and visualization are integrated parts of the same analytical process!





Buit fisrt I need you to go through all the files in the project folders and summarize it.

match the currenbt file structure to the rationale as decribed below

Experiment-Centric Organization

psu_2022/ clearly positioned as the central/core experiment with multi-omics sub-folders
clayton_2025/ with its own phenotype_spatial_modeling/ sub-folder showing application of your approach
chamber_experiments/ and crow_reanalysis/ as supporting analyses

Analytical Context Integration

Plots now live with their analyses (Manhattan plots in expression/, network plots in network/)
genomics_pipelines/ properly positioned as supporting infrastructure
annotation/ separated from core analytical workflows

🎯 Framework Emphasis:
Your Methodological Innovation Highlighted

Spatial modeling prominently featured as core methodological advance
Consistency analysis vs mashr clearly distinguished
Cross-experiment integration framework outlined

Script Naming Convention Documentation

Verb-based naming for action scripts explained with examples
Noun-based naming for function libraries distinguished
Clear examples: analyze_*, detect_*, fit_* vs *_helpers.R, *_template.R

📋 Key Additions:
Workflow Descriptions

Each analytical pipeline clearly described with script references
Shows progression from PSU → Clayton → cross-experiment comparison
Links scripts to biological hypotheses

Implementation Strategy

Phase 1: Core framework using PSU 2022 scripts
Phase 2: Cross-experiment using Clayton 2025 spatial modeling
Phase 3: Crow reanalysis validation



