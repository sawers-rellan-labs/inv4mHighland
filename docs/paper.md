---
bibliography:
- Inv4mPhosphorus.bib
---

::: flushleft
**Introgression of a Mexican highland chromosomal inversion into
temperate maize accelerates flowering, promotes growth, and modulates a
cell proliferation gene network.**

Fausto Rodríguez-Zapata^1,2,\*^, Nirwan Tandukar^1,2^, Allison
Barnes^1,2^, Ruthie Stokes^1,2^, Alejandro Aragón-Raygoza ^1,3^, Meng
Li^3^, Sergio Pérez-Limón^3^, Melanie Perryman^3^, Miguel A. Piñeros^4^,
Josh Strable ^1,3^, Daniel Runcie^5^, Ruairidh Sawers ^3^ and Rubén
Rellán-Álvarez ^1,2,\*^\
**1** Department of Molecular and Structural Biochemistry, North
Carolina State University, Raleigh, NC, USA\
**2** Genetics and Genomics Program, North Carolina State University,
Raleigh, NC, USA\
**3** Department of Genetics, Development and Cell Biology, Iowa State
University, Ames, IA, USA\
**4** Department of Plant Science, Pennsylvania State University,
University Park, PA, USA\
**5** Robert W. Holley Center for Agriculture and Health, USDA-ARS,
Ithaca, NY, USA\
**6** Department of Plant Sciences, University of California, Davis, CA,
USA\

\*rrellan@ncsu.edu
:::

# Abstract {#abstract .unnumbered}

*Inv4m* is a chromosomal inversion prevalent in traditional maize
varieties adapted to the cold and often phosphorus-deficient Mexican
highlands. Field trials throughout Mexico have shown that, when grown at
high elevations, plants carrying the inversion flower faster and have
greater yield than plants without it. Although growth chamber
experiments indicate that *Inv4m* regulates the expression of
photosynthesis-related genes in response to cold, we have yet to know
the genes responsible for the adaptive effects of *Inv4m* in the field.
To identify *Inv4m*-regulated genes that underlie enhanced development
in the field, we bred B73-based Near Isogenic Lines (NILs) with either
the inversion or the standard karyotype. Additionally, we grew these
NILs in phosphorus-sufficient and deficient soils to test whether
*Inv4m* contributes to local adaptation through enhanced phosphorus
stress response. We measured plant reproductive and vegetative traits,
phosphorus, lipids, and gene expression in the leaves. Plants showed
classical responses to phosphorus starvation, including decreased
phosphorus and biomass accumulation, delayed flowering, and a switch
from phospholipid to glycolipid production. Notably, *Inv4m* plants
flowered earlier and grew taller regardless of phosphorus availability.
While increased leaf age and phosphorus deficiency resulted in
genome-wide expression changes, *Inv4m*'s effects were predominantly
confined to genes within the inversion. Our analyses suggest that
*Inv4m* introgression modulates a trans-coexpression network enriched in
cell proliferation and flower development genes, which includes DNA
replication fork proteins (*pcna2*, *mcm5*), histone demethylases
(*jmj2*, *jmj21*), and the FT florigen homolog *zcn26*. By
cross-referencing with a list of candidates from the literature, we
found other *Inv4m*-regulated genes associated with flowering time and
plant height. In a complementary growth chamber experiment, *Inv4m*
plants showed longer shoot apical meristems than controls, supporting
its effect on organ development. These findings provide insights into
*Inv4m*'s role in highland adaptation through the coordinated expression
of a developmental gene network.

# Author summary {#author-summary .unnumbered}

# Introduction {#introduction .unnumbered}

Maize might have preadapted to temperate zones through introgression
from its highland wild relative, *Zea mays ssp. mexicana* (shorthand
*mexicana*) [@yang2023]. Maize was originally domesticated in the
tropical lowlands of Mexico. Before its expansion into temperate
regions, maize was introduced to the Mexican highlands and the
Southwestern United States, where sympatry with highland teosinte likely
facilitated the introgression of adaptive alleles from *mexicana*.

However, not all highland-adaptive loci are present in temperate maize.
Highland-associated chromosomal inversions, such as *Inv4m* and *Inv9f*,
are prevalent in highland teosinte populations [@pyhajarvi2013] and
traditional Mexican maize varieties (TVs)
[@crow2020; @gonzalez-segovia2019-jy] but are rare in temperate maize.
Chromosomal inversions can contribute to local adaptation by preserving
locally adapted alleles across multiple loci and reducing recombination
within the inversion [@kirkpatrick2006b]. Genotyping of teosinte
populations using the Maize 50K chip revealed that *Inv4m* spans 13 Mb
and is predominantly found in *mexicana* populations [@pyhajarvi2013].
In Mexican TVs, variation in *Inv4m* is associated with elevation and
flowering time [@romero_navarro2017-cn]. Additionally, *Inv4m* shows
reduced genetic diversity, a clinal relationship with elevation, and is
nearly fixed in locations above 2500 m.a.s.l. [@crow2020]. The inversion
exhibits suppressed recombination, as confirmed in a biparental cross
[@gonzalez-segovia2019-jy]. *Inv4m* demonstrates classic patterns of
gene-by-environment interactions indicative of local adaptation. Plants
carrying the *Inv4m*-highland allele exhibit delayed flowering at low
elevations and earlier flowering at high elevations [@crow2020]. The
highland haplotype of *Inv4m* was introgressed from Zea mays ssp.
mexicana [@pyhajarvi2013; @calfee2021-mr; @hufford2013-gs], a wild maize
relative native to the Mexican highlands.

Despite strong evidence linking *Inv4m* to local adaptation, the
physiological processes and environmental factors underlying its
adaptive role remain unclear. Furthermore, the specific genes within
*Inv4m* that drive local adaptation are largely unidentified. Previous
research has shown that *Inv4m*-highland upregulates photosynthesis
genes in response to cold at the seedling stage [@crow2020] and is
associated with earlier flowering in the Mexican highlands, which likely
enhances fitness in environments with limited growth-degree accumulation
throughout the year [@romero_navarro2017-cn]. However, cold is not the
only limiting factor for plant growth in the highlands. Volcanic soils
(Andosols), which dominate the Mexican highlands, present an additional
constraint. Approximately 95% of natural Andosol profiles in Mexico are
found above 2000 m.a.s.l. [@paz-pellat2018; @inegi2013]. These soils are
characterized by high phosphorus retention [@krasilnikov2013], which
leads to low phosphorus availability for plant uptake
[@galvan-tejada2014]. MICH21, one of the Mexican highland maize
accessions analyzed by [@crow2020], originates from the Purépecha
Plateau, where Andosols and phosphorus-efficient TVs are common
[@paz-pellat2018; @galvan-tejada2014; @bayuelo-jimenez2011; @bayuelo-jimenez2014].
*Inv4m* may contribute to adaptation in the highlands by carrying
alleles that enhance the phosphorus starvation response (PSR). For
example, the phosphate transporter gene *ZmPho1;2a*, located within
*Inv4m*, is a strong candidate for adaptation to low phosphorus
availability [@salazar-vidal2016-rl].

In this study, we aimed to understand the physiological and molecular
effects of *Inv4m* and to identify candidate genes within the inversion
that could elucidate its adaptive role. Specifically, we tested whether
*Inv4m*-highland contributes to adaptation to low phosphorus
availability. To achieve this, we backcrossed MICH21, a Mexican highland
TV carrying *Inv4m*, into the B73 genetic background for eight
generations, generating Near Isogenic Lines (NILs) [@crow2020]. These
NILs were grown under temperate field conditions with two phosphorus
treatments to evaluate flowering time, height, and transcriptomic
responses.

We observed that *Inv4m* significantly reduces flowering time and
increases plant height, independent of phosphorus levels. We identified
a cluster of JUMONJI methyltransferases, which have higher copy numbers
in modern maize compared to highland teosinte and TVs, as potential
contributors to flowering time differences. Additionally, we found a
group of cell cycle-related genes that may underlie height differences.
While classical PSR genes were identified within the inversion, *Inv4m*
had no detectable effect on the phosphorus starvation response. These
findings provide insights into the genetic mechanisms underlying
*Inv4m*'s contribution to local adaptation and highlight potential
candidate genes driving its effects.

# Results {#results .unnumbered}

## *Inv4m* is a 15 Mb inversion delimited by breakpoint segments overlapping chromosomal knob repeats {#inv4m-is-a-15-mb-inversion-delimited-by-breakpoint-segments-overlapping-chromosomal-knob-repeats .unnumbered}

Complete genome sequences of selected Mexican Highland traditional
varieties allowed us to delimit the recombination breakpoints of *Inv4m*
by using an Anchorwave [@song2022] alignment of chromosome 4. With this
aim, we aligned both PT, a close relative to the *Inv4m* donor parent
Mi21, as a representative genome for the Mexican highlands, and TIL18
representing the highland teosinte *mexicana*, to B73, the recurrent
parent used for *Inv4m* introgression
(Fig [1](#fig::design){reference-type="ref" reference="fig::design"}A,
Supporting
Table [\[tab:breakpoints\]](#tab:breakpoints){reference-type="ref"
reference="tab:breakpoints"}). In each pairwise comparison, we found
that the anchor-defined *Inv4m* alignment is between the reference B73
plus strand and the highland genome minus strand, as expected. This
alignment allowed us to narrow down *Inv4m* location between genes
*Zm00001eb190470* and *Zm00001eb194800*, spanning 15.2 Mb and 432
annotated genes in the B73 v5 genome. This region corresponded to a 13.4
Mb inversion in the PT genome, bounded by genes *Zm00109aa017629* and
*Zm00109aa018009*; and, also a 13.2 Mb inversion in TIL18, extending
from genes *Zx00002aa015554* to *Zx00002aa015905*. However, the *Inv4m*
section was bordered by unaligned segments rather than by reference
orientation (plus/plus) alignments. These unaligned segments contained
no annotated genes, except for *Zm00001eb190490* in the upstream
breakpoint segment of B73. To explore the nature of the sequences
spanning the *Inv4m* breakpoints, we made a local similarity search
using LAST [@kielbasa2011] inside these *Inv4m* bounding unaligned
segments. We detected typical tandem array patterns in the dotplots
(blue and red rectangles in Fig [1](#fig::design){reference-type="ref"
reference="fig::design"}B) and we took a closer look at the available
repeat annotation of B73. The MaizeGDB annotation of B73 shows that the
upstream *Inv4m* breakpoint segment contained 438 `knob180` and 29
`TR-1` repeats making up $69\%$ and $5\%$ of the 636 annotated repeats
in this 323 kb section. These chromosomal knob 180 bp and TR-1 repeat
annotations matched the observed tandem array patterns observed in the
dotplots Fig [1](#fig::design){reference-type="ref"
reference="fig::design"}B Then we searched the two classes of knob
repeats in TIL18 and PT, which have no public repetitive element
annotation to date. For comparison purposes, we did a BLAST+
[@camacho2009] similarity search of the three chromosomes 4 using a
selected sequence of each repeat class as subjects (see methods). This
allowed us not only to detect both kinds of knob repeats in the
breakpoint segments but also inside the inversion and the rest of the
chromosome. We found in the long arm of chromosome 4 repeat arrays
corresponding to previously reported heterochromatic knobs, the ones
interspersed in the *Inv4m* breakpoints, and an internal knob repeat
array inside *Inv4m* Fig [1](#fig::design){reference-type="ref"
reference="fig::design"}A. Chromosomal knob repeats were the most
abundant kind of repeats within the *Inv4m* breakpoint sections Table
[\[tab:breakpoints\]](#tab:breakpoints){reference-type="ref"
reference="tab:breakpoints"} and the *Inv4m* proper. There is a notable
difference in the heterochromatic knobs between the long arms of the
analyzed chromosomes. TIL18 has a 30 Mb region of mixed repeats likely
corresponding to the 4L knob reported in
*mexicana*[@bilinski2018; @albert2010], PT is knobless, and B73 carries
a 4L knob consisting of *TR-1* repeats [@ghaffari2013; @albert2010]
(Fig [1](#fig::design){reference-type="ref" reference="fig::design"}A).
Although the downstream breakpoint had a lower content of knob repeats
with 37 `knob180` and no annotated `TR-1`'s, those made up $46\%$ out of
the 81 annotated repeats in 80 kb. In each genome, the breakpoint
segments on each side of the inversion showed asymmetries in size, knob
repeat number, presence of knob inverted repeat arrays, and closeness of
the knob repeats to the *Inv4m* anchor gene
(Fig [1](#fig::design){reference-type="ref" reference="fig::design"}B).
In particular, in B73 the upstream breakpoint segment is larger (323 kb
vs 89 kb), has a greater number of repeats (352 vs 26), presence of
inverted knob repeat arrays (blue and red rectangles in
Fig [1](#fig::design){reference-type="ref" reference="fig::design"} B
Upstream), and some knob repeats close to the inversion, while the
downstream breakpoint segment is smaller, has a lower number of repeats,
no inverted repeats of knob sequences (just red rectangles in
Fig [1](#fig::design){reference-type="ref" reference="fig::design"}B,
Downstream), and the few knob repeats present are closer to the
downstream external anchor gene rather than to *Inv4m*. We can observe
these four patterns in the opposite orientation in PT and TIL18 relative
to B73.

<figure>

<figcaption> <strong>(A)</strong> <em>Inv4m</em> breakpoints overlap
clusters of full 180 bp knob repeats with interspersed partial
<em>TR-1</em>s (normalized match score <span
class="math inline"> &lt; 50</span>). Each point represents a BLAST+
<span class="citation" data-cites="camacho2009"></span> hit with a
selected sequence of 180 bp knob or <em>TR-1</em> repeats. Around the
230 Mb mark regions corresponding to reported 4L cytogenetic knobs can
be observed in <em>mexicana</em> and B73 while absent in PT.
<strong>(B)</strong> Sequence self-similarity dotplots, using LAST <span
class="citation" data-cites="kielbasa2011"></span>, for the inversion
breakpoint segments as defined in Table <a href="#tab:breakpoints"
data-reference-type="ref"
data-reference="tab:breakpoints">[tab:breakpoints]</a>. Breakpoints are
upstream and downstream of <em>Inv4m</em> in each genome coordinate
system as depicted in (A). Forward matches in red, reverse in blue.
Rectangles indicate tandem repeats. <strong>(C)</strong> PT shows higher
(** <em>t-test</em> <span
class="math inline"><em>p</em> = 0.0022</span>) leaf phosphorus content
than B73 under phosphorus sufficiency (Puerto Vallarta, 2016).
<strong>(D)</strong> Near Isogenic Line population breeding scheme. The
<em>Inv4m</em> from Mi21, a highland traditional variety closely related
to PT, was introgressed into a temperate lowland genetic background
using B73 as recurrent parent. <strong>(E)</strong> Rock Springs, PA,
field station experimental setup at 5 ppm (-P), and 36 ppm (+P) soil
phosphorus [mg/kg Mehlich-3]. <strong>(F)</strong> NIL genotypes at
BC<span class="math inline"><sub>6</sub></span>S<span
class="math inline"><sub>2</sub></span>, of plants used in RNA and lipid
extractions, <span class="math inline"><em>n</em> = 13</span>, B73 NAM5
coordinates. <strong>(G)</strong> Zoom into <em>Inv4m</em>
introgression. Plants are sorted by genotype at the <em>Inv4m</em>
tagging SNP PZE04175660223 (181.6 Mb, downwards pointing triangle) then
by field row number. </figcaption>
</figure>

<figure id="fig::design">
<img src="figs/design.png" />
<figcaption><em><strong>Chromosomal Inversion <em>Inv4m</em>:
Delimitation, Introgression Breeding Design, Field Experimental Setup
and Population Genotype.</strong></em> </figcaption>
</figure>

## Divergent highland introgression is mostly confined to a 39 Mb segment spanning *Inv4m* {#divergent-highland-introgression-is-mostly-confined-to-a-39-mb-segment-spanning-inv4m .unnumbered}

To study the effects of *Inv4m* we repeatedly backcrossed Mi21, a
Mexican highland landrace containing *Inv4m* into B73 using a marker
inside the inversion. We then used RNA-Seq data to genotype the
resulting $\text{BC}_6\text{S}_2$ lines. The distribution of the 19861
QC-filtered SNPs is heavily biased towards chromosome 4. Almost half of
the variants, 8892 (44%), are on chromosome 4, with the other 9
chromosomes having an average of 1219 SNPs (6%) each (Supporting
Fig [8](#fig::SNPdistro){reference-type="ref"
reference="fig::SNPdistro"}A). The genotypes in the matrix of 19861 SNPs
$\times$ 13 individuals are 63.5% homozygous for B73, 8.96%
heterozygous, and 27.3% homozygous for Mi21, with 0.16% missing data. We
converted this genotype matrix to identical genotype run [@layer2016]
length in base pairs to estimate the proportion of the plant genome
introgressed from Mi21 (Fig [1](#fig::design){reference-type="ref"
reference="fig::design"}F and G). By genotype run length, the plants
are, on average, 85.6% homozygous for B73 (1821 Mb), 8.61% heterozygous
(183 Mb), and 5.5% homozygous for Mi21 (117 Mb), with 0.34 % missing
data (7.32 Mb) (Supporting Fig [8](#fig::SNPdistro){reference-type="ref"
reference="fig::SNPdistro"}D). As the genotypes come exclusively from
variant sites, the matrix has no markers fixed for the reference B73
allele. The NIL genomes are thus represented by a mosaic of 3 types of
sections. First, sections with fixed highland introgression, where the
alternate Mi21 allele is fixed. Second, sections of divergent highland
introgression consisting of markers significantly correlated with the
*Inv4m* tagging SNP PZE04175660223. And third, sections with random
introgressions uncorrelated with *Inv4m*. The fixed highland
introgression sections are made up of 126 SNPs marking 51 Mb (2.4% of
the genome) in genotype run length (average per plant). By manual
inspection we defined a region of divergent highland introgression
between positions 157012149 and 195900523 in B73 v5 coordinates,
bounding 39 Mb of shared introgression among *Inv4m* carrying lines.
With our selection for the highland allele of PZE04175660223 we have
produced plants with the introgression of the target *Inv4m* inversion
surrounded by 24 Mb of lingering unintended linkage drag
(Fig [1](#fig::design){reference-type="ref" reference="fig::design"}F
and G). Conversely, we have produced inversion-free lines in the CTRL
group by selecting for the B73 allele. In more detail, there are 7683
SNPs significantly correlated with *Inv4m*(Pearson correlation *t-test*,
$\textrm{\textit{FDR}} < 0.05$) constituting divergent highland
introgressions throughout the genome but mostly in *Inv4m* and flanking
regions (Supporting Fig [8](#fig::SNPdistro){reference-type="ref"
reference="fig::SNPdistro"}A). From these 7683, 7322 SNPs ($95.3\%$)
cover 35 Mb and 271 expressed genes inside the 39 Mb shared introgressed
region that includes *Inv4m*. The remaining 361 divergently introgressed
SNPs in the genomic background were dispersed throughout 10.4 Mb of
sequence covering 52 expressed genes. Notably, the majority of these
markers, 207 SNPs from 14 expressed genes, are negatively correlated
with PZE04175660223 and are located just upstream of the shared
introgression near *Inv4m* over a 2.7 Mb stretch, from position
154283637 to 15698557 (Fig [1](#fig::design){reference-type="ref"
reference="fig::design"}G and Supporting
Fig [8](#fig::SNPdistro){reference-type="ref"
reference="fig::SNPdistro"}C). In this segment, individuals selected for
*Inv4m* are homozygous for the B73 allele at all positions, while
individuals selected for the standard karyotype are frequently
heterozygous. In summary, each plant has an average of 117 Mb of
homozygous Mi21 genome per line, 51 Mb of fixed highland introgression,
35 Mb of divergent highland introgression colocalized and matching the
selection of *Inv4m*, and 10 Mb of other dispersed divergent highland
introgressions.

## Reproductive and vegetative traits show predominantly independent responses to phosphorus treatment and *Inv4m* introgression. {#reproductive-and-vegetative-traits-show-predominantly-independent-responses-to-phosphorus-treatment-and-inv4m-introgression. .unnumbered}

To evaluate the effects of *Inv4m* on plant growth, reproductive traits,
and its effect on plant phosphorus efficiency, we grew the NILs in two
levels of phosphorus in a field experiment in Rock Springs, PA. In
general, we found that the subject plants have a phenotypic response to
phosphorus treatment that is independent of the *Inv4m* introgression
(Fig [2](#fig::effects){reference-type="ref" reference="fig::effects"}A,
$H_0:\beta_i = 0$, *t-test*, $\textrm{\textit{FDR}} < 0.05$). Plants
grown in low phosphorus accumulate less vegetative biomass
($-20.4 \pm 5.63$ g stover dry weight at harvest, effect mean $\pm$
S.E.), reach maturity later ($3.63 \pm 0.54$ DTA, $3.3 \pm 0.46$ DTS),
and have reduced kernel weight ($-1.85 \pm 0.62$ g for a 50 kernel
sample). The negative effect of phosphorus deficiency on the
accumulation of shoot vegetative biomass is significant at all the
measured time points (Fig [2](#fig::effects){reference-type="ref"
reference="fig::effects"}A, Fig [9](#fig::growth){reference-type="ref"
reference="fig::growth"} A). In contrast, phosphorus deficiency did not
result in a significant difference in plant height between experimental
groups. An analysis of the fitted growth curves per row for vegetative
biomass, Fig [9](#fig::growth){reference-type="ref"
reference="fig::growth"} B-F, showed a delay in time to mid-mass in
phosphorus deficiency ($3.53 \pm 1.01$ days), in addition to the
decrease in estimated maximum dry stover weight ($-22.51 \pm 3.88$ g).
Overall phenotypic effects of phosphorus starvation were similar in both
*Inv4m* and CTRL plants (Fig [2](#fig::effects){reference-type="ref"
reference="fig::effects"}F-I), with no evidence of *Inv4m*$\times$ P
interaction effect except for cob diameter. In this case *Inv4m* ears
are more responsive to phosphorus starvation, where phosphorus
deficiency induced a reduction of $2.34 \pm 0.79$ cm in cob diameter
while the control plants showed no significant difference between
phosphorus treatments (Fig [2](#fig::effects){reference-type="ref"
reference="fig::effects"}E). Plants carrying *Inv4m* flowered earlier
($-1.27 \pm 0.50$ DTA,$-1.11 \pm 0.42$ DTS), and grew taller
($5.32  \pm 1.97$ cm) than plants in the control group. In addition to
this genotypic difference in height, we detected a reduced stover dry
weight at 40 DAP ($-3.14 \pm 1.21$ g) in *Inv4m* plants, and these later
three phenotypic responses show null dependency on soil phosphorus
treatment.

## Phosphorus deficiency modifies the plant mineral composition independently of *Inv4m* genotype {#phosphorus-deficiency-modifies-the-plant-mineral-composition-independently-of-inv4m-genotype .unnumbered}

Soil phosphorus depletion significantly reduced tissue phosphorus
content compared to phosphorus-sufficient conditions
($H_0:\beta_i \neq 0$, *p-value* $< 1 \times 10^{-9}$,
Fig [2](#fig::effects){reference-type="ref" reference="fig::effects"} A,
bottom) with a stronger effect in the stover ($-1569 \pm 119$ ppm, mean
$\pm$ S.E. as mg/kg dry weight) than in seed tissue ($-691 \pm 121$ ppm,
Fig [2](#fig::effects){reference-type="ref" reference="fig::effects"} H
vs I). This differential response resulted in a $2.08  \pm 0.22$
increase in the seed to stover phosphorus ratio
(Fig [10](#fig::ionome){reference-type="ref" reference="fig::ionome"} A
and B). Phosphorus deficiency also altered other mineral concentrations:
the stover showed increased sulfur ($111 \pm 32$ ppm) and zinc
($6.82 \pm 1.69$ ppm) content, while the seed to stover zinc ratio
decreased ($-0.32 \pm 0.07$). In seeds, magnesium content decreased
($-101 \pm 39.5$ ppm ) while calcium increased ($19.9 \pm 7.9$)
(Fig [2](#fig::effects){reference-type="ref" reference="fig::effects"}
A, and Fig [10](#fig::ionome){reference-type="ref"
reference="fig::ionome"} C-G). All ionomic responses to phosphorus
deficiency occurred independently of the plant *Inv4m* genotype.

<figure id="fig::effects">
<img src="figs/effects.png" />
<figcaption><em><strong><em>Inv4m</em> plants flower earlier and grow
higher independent of phosphorus treatment</strong></em>
<strong>(A)</strong> Significant effects of phosphorus treatment,
<em>Inv4m</em> genotype and interaction, in different plant phenotypes,
<span class="math inline"><em>n</em> = 64</span>. Standardized effect
and <span class="math inline">95%</span> confidence interval
(<em>t-test</em> for regression coefficients, <span
class="math inline"><em>H</em><sub>0</sub> : <em>β</em><sub><em>i</em></sub> = 0</span>,
<span class="math inline">$\textit{\textit{FDR}} &lt;0.05$</span>). The
phenotypes are sorted by maximum effect magnitude in each category the
the left. <strong>(B-G)</strong> Data distribution for significant
effects in original units. Yellow: reference group, Purple: effect
group. <strong>(B-D)</strong> Effects of <em>Inv4m</em>.
<strong>(E)</strong> Effect of <em>Inv4m</em><span
class="math inline">×</span> -P interaction. <strong>(F-I)</strong>
Effects of phosphorus starvation.<br />
pairwise <em>t-test</em> raw <em>p-value</em> * <span
class="math inline"> &lt; 0.05</span>, ** <span
class="math inline"> &lt; 0.01</span> *** <span
class="math inline"> &lt; 0.001</span>, *** <span
class="math inline"> &lt; 0.001</span>, **** <span
class="math inline"> &lt; 0.0001</span>.</figcaption>
</figure>

## Leaf stage and phosphorus starvation promote global changes in gene expression, while *Inv4m* introgression drives local transcriptional effects. {#leaf-stage-and-phosphorus-starvation-promote-global-changes-in-gene-expression-while-inv4m-introgression-drives-local-transcriptional-effects. .unnumbered}

Our analysis of leaf transcriptome data revealed distinct patterns of
gene expression driven by developmental stage, phosphorus availability,
and *Inv4m* introgression. A multidimensional scaling (MDS) of
$log_2[\text{CPM}]$ values captured 46% of expression variance in the
first three dimensions, each corresponding to specific experimental
factors (Fig [3](#fig::RNAseq){reference-type="ref"
reference="fig::RNAseq"} B and C). The first dimension explained 26% of
variance and is correlated to phosphorus treatment (Pearson $r=0.50$,
*t-test* $\textit{FDR} = 6.15 \times 10^{-4}$). The second dimension
explained 12% of observed expression variance and is inversely
correlated with leaf stage ($r=-0.62$,
$\textit{FDR} = 2.26 \times 10^{-5}$). And the third dimension with 8%
of explained variance is associated with *Inv4m*, separating samples by
genotype ($r=0.92$, $\textit{FDR} = 2.06\times 10^{-17}$). Both
phosphorus treatment and leaf stage had a genome-wide effect on gene
expression. Out of 24011 genes detected in at least one leaf sample, we
found 14330 differential expressed genes responding to leaf stage, and
10606 DEGs responding to soil phosphorus deficiency
(Fig [3](#fig::RNAseq){reference-type="ref" reference="fig::RNAseq"} D
and E). In contrast, significant transcriptional responses to the
*Inv4m* introgression were predominantly localized to the inversion and
its flanking regions (Fig [3](#fig::RNAseq){reference-type="ref"
reference="fig::RNAseq"} F). To characterize the distribution of *Inv4m*
responsive genes, we classified differentially expressed genes (DEGs)
according to their physical and genetic linkage to the inversion. We
defined *trans* regulated genes as those outside the inversion and its
surrounding linkage drag, while *cis* regulated genes were those within
the 39 Mb shared introgression spanning the inversion. The distribution
of DEGs responsive to *Inv4m* introgression in the genome by location
relative to the *Inv4m* locus is presented in
Table [\[tab:DEGs_distro\]](#tab:DEGs_distro){reference-type="ref"
reference="tab:DEGs_distro"}. Among the 970 DEGs responding to *Inv4m*
introgression, 646 were *trans* regulated while 324 were regulated in
*cis*. The *cis* regulated genes were proportionally distributed
(Fisher's exact test $p=0.28$) between the flanking regions (183 genes
in 24 Mb) and *Inv4m* proper (141 genes in 15 Mb). DEGs responding to
*Inv4m* were 40 times more likely to originate from the shared
introgression (Fisher's exact test $p=2.38 \times 10^{-300}$) and 45
times more likely to come from the inversion proper compared to *trans*
regulated genes (Fisher's exact test $p=2.73 \times 10^{-142}$). This
locus enrichment was even more pronounced among DEGs with high magnitude
response, $log_2[\text{FC}]>2$, hereafter named top DEGs. Of 139 top
DEGs, 96 were located in *cis* (61 in flanking regions, 35 in *Inv4m*
proper) and 43 in *trans*, representing 70% and 35% of top DEGs in the
shared introgression and inversion proper, respectively. These top DEGs
were 102 times more likely to originate from the shared introgression,
112 times more likely from the flanking segments, and 87 times more
likely from the *Inv4m* proper compared to *trans* regulated genes
(Fisher's exact test $p < 1 \times 10^{-48}$). The peak of significance
for differential expression due to *Inv4m* introgression is located
inside *Inv4m* proper and corresponds to the cluster of *JUMONJI* DNA
demethylases that includes *jmj2*, *jmj6* and
*jmj4*(Fig [3](#fig::RNAseq){reference-type="ref"
reference="fig::RNAseq"} F to I). However, the distributions of
*p-values* for differential expression due to *Inv4m* do not differ
between *Inv4m* proper and flanking regions, while both regions have
excess significance relative to the rest of the genome
(Kolmogorov--Smirnov test $p=0.29$ and $p\leq 2.26 \times 10^{-16}$
respectively, Supplementary Figure S ). Although a majority of the DEGs
in *Inv4m* proper are downregulated (82 genes, 58%), the proportion of
repressed genes is not different to that in the flanking the regions
(102 genes, 56%) (Fisher's exact test $p=0.73$) with a similar result
for top DEGs (Fisher's exact test $p=0.39$).

<figure id="fig::RNAseq">
<img src="figs/RNAseq.png" />
<figcaption><em><strong>Genome View of the Effect of Leaf Stage,
Phosphorus Treatment, and <em>Inv4m</em> in Gene
Expression</strong></em> <strong>(A)</strong> Experimental design.
Increasing number in leaf stage corresponds to older leaves. Younger
leaves are at a more apical position. <strong>(B-C)</strong> Gene
Expression Multidimensional Scaling. Samples cluster by leaf stage,
phosphorus treatment, and <em>Inv4m</em> genotype.
<strong>(D-F)</strong> Manhattan plots showing statistical significance
of differential gene expression for each experimental predictor. The
number of differentially expressed genes (<span
class="math inline">$\text{\textit{FDR}}&lt;0.05$</span>, red horizontal
line) is indicated after the colon. (D) Effect of increasing leaf stage
between two consecutive sampled leaves. (E) Effect of -P treatment and
(F) effect of <em><em>Inv4m</em></em> genotype relative to control
plants. <strong>(G-I)</strong> Zoom around the <em>Inv4m</em>
introgression detailing NIL genotype, gene expression significance, and
size of the <em>Inv4m</em> effect depicted in (F).(G) Example genotypes
of NILs homozygous for the highland allele of the <em>Inv4m</em> tagging
SNP PZE04175660223 (blue downward pointing triangle). (H) The most
significant DEGs due to <em>Inv4m</em> are <em>jmj2</em>, <em>jmj6</em>,
and <em>jmj9</em>, which belong to a gene cluster well within the
boundaries of the inversion proper. (I) Effects as <span
class="math inline"><em>l</em><em>o</em><em>g</em><sub>2</sub>(<em>F</em><em>C</em>)</span>
of <em>Inv4m</em> vs control plants, grey band <span
class="math inline">|<em>l</em><em>o</em><em>g</em><sub>2</sub>(<em>F</em><em>C</em>)| &lt; 2</span>,
and rolling window average trend lines (10 genes per window, green <span
class="math inline"><em>l</em><em>o</em><em>g</em><sub>2</sub>(<em>F</em><em>C</em>) &gt; 0</span>,
darkred <span
class="math inline"><em>l</em><em>o</em><em>g</em><sub>2</sub>(<em>F</em><em>C</em>) &lt; 0</span>).</figcaption>
</figure>

## Despite strong transcriptional effects, *Inv4m* differential expression lacks functional enrichment. {#despite-strong-transcriptional-effects-inv4m-differential-expression-lacks-functional-enrichment. .unnumbered}

The leaf gene expression and lipid profile reveal distinct molecular and
metabolic responses to leaf stage, phosphorus availability, and *Inv4m*
genotype (Fig [4](#fig::volcano){reference-type="ref"
reference="fig::volcano"}). Increased leaf stage, which correlates with
aging and senescence, leads to the downregulation of genes involved in
floral meristem determinacy (GO biological process enrichment, Fisher's
exact test *FDR* $=1.28\times 10 ^{-5}$). This set was constituted by
the MADS box transcription factors *zmm4*, *mads45*, *zmm15*, *mads67*;
*zap1* (a maize *APETALA1* homolog); and the TNF receptor-associated
factor, *traf7*. Conversely, samples from later leaf stages show
upregulation in a set of genes enriched in aging (Fisher's exact test,
*FDR* $=0.002$) that includes the diaminopimelate aminotransferase
*dapat3* and the *SAG12* cysteine protease homolog *ccp16*. The KEGG
pathway annotation shows that genes encoding photosynthetic antenna
proteins (Fisher's exact test, *FDR* $=0.015$) are significantly
downregulated in phosphorus deficiency (*lhcb10*) and increased leaf
stage (*lhcb1*, *lhcb10*) Fig [4](#fig::volcano){reference-type="ref"
reference="fig::volcano"} B. Downregulated genes under phosphorus
starvation do not show enrichment in specific biological processes.
However, the most significant PSR gene was *peamt2*, a
phosphoethanolamine N-methyltransferase, which is involved in
phospholipid biosynthesis, congruent with the reduction in phospholipid
synthesis under phosphorus-limiting conditions. On the other hand,
plants growing under phosphorus deficiency showed an evident
upregulation of classical phosphorus starvation response genes. Although
not included in the GO annotation of protein-coding genes, the most
prominent response is that of *pilncr1*, a long noncoding RNA spanning
*mir397*, a known master regulator of PSR. The upregulated
protein-coding genes showed enrichment in cellular response to phosphate
starvation (Fisher's exact test, *FDR* $=9.07 \times 10^{-11}$),
including *SPX* family transcription factors, and phosphate transporters
such as *phos1*, which facilitate phosphate uptake and redistribution
*pht1* *pht7*, and purple acid phosphatases *pap1* *pap14* that increase
phosphorus remobilization. Over-representation analysis of KEGG
metabolic pathways showed that phosphorus starvation upregulates genes
involved in glycerophospholipid metabolism (glycerophosphodiester
phosphodiesterases *gpx1*, *gpx3*, monogalactosyldiacylglycerol synthase
*mgd2*), phenylpropanoid and flavonoid biosynthesis, chitin degradation,
ABC transporters, and brassinosteroid biosynthesis (*brc2*). Conversely,
genes involved in galactolipid biosynthesis, which are crucial for
maintaining photosynthetic membranes, are significantly downregulated.

## Leaf lipids respond consistently to phosphorus starvation and leaf stage while response to *Inv4m* is dependent on leaf stage. {#leaf-lipids-respond-consistently-to-phosphorus-starvation-and-leaf-stage-while-response-to-inv4m-is-dependent-on-leaf-stage. .unnumbered}

Lipid profiling shows shifts associated with both leaf aging and
phosphorus starvation (Fig [4](#fig::volcano){reference-type="ref"
reference="fig::volcano"} B) . Increased leaf stage is linked to the
accumulation of digalactosyldiacylglycerol (DGGA), particularly DGGA36:3
and DGGA42:1, as well as a significant increase in triacylglycerols
(TGAs). The accumulation of TGAs suggests enhanced lipid storage,
possibly as an energy reservoir or stress adaptation mechanism during
senescence.Phosphorus starvation induces a well-characterized membrane
lipid remodeling response, shifting from phosphoglycerolipids (e.g.,
phosphatidylcholines \[PCs\], phosphatidylethanolamines \[PEs\],
lysophosphatidylethanolamines \[LPEs\], lysophosphatidylcholines
\[LPCs\], and phosphatidylglycerols \[PGs\]) to sugar-based glycolipids
such as DGGA. This shift is a widely observed adaptive mechanism in
plants under phosphorus limitation, reducing reliance on phosphorus-rich
membrane lipids while maintaining membrane integrity and function.
Additionally, phosphorus starvation also leads to an accumulation of
triacylglycerols (TGAs), which may serve as an alternative storage form
of fatty acids under stress conditions. *Inv4m* shows no apparent effect
on the differential gene expression or lipid production at this level.
Based on this nuanced analysis, *Inv4m* does indeed have significant
effects on lipid profiles when leaves are categorized by age
(bottom/older vs top/younger). The *Inv4m* genotype shows consistent
alterations across both leaf age groups, with notable decreases in
phosphatidylethanolamines (PE36:4 and PE36:6), decreased galactolipids
(DGDG36:2 and MGDG36:3), and a widespread increase in nearly all
detected triacylglycerols except TG56:2. Most intriguingly, the data
reveals a complex interaction between Inv4m genotype, phosphorus
availability, and leaf developmental stage. Under phosphorus starvation,
*Inv4m* older leaves specifically show enhanced accumulation of
monogalactosyldiacylglycerols (MGDG34:2 and MGDG34:3) and most
triacylglycerols (except TG56:2). However, younger leaves of the same
genotype exhibit an opposite response pattern, with decreased levels of
these same lipid species. This developmental stage-dependent response
suggests that Inv4m alters lipid remodeling mechanisms in a
tissue-specific manner, particularly influencing how plants manage
membrane composition and storage lipid accumulation during phosphorus
limitation. The consistent decrease in specific
phosphatidylethanolamines (PEs) across leaf stages in Inv4m plants,
coupled with altered galactolipid profiles, indicates that this genotype
fundamentally affects membrane lipid homeostasis, potentially modifying
the balance between phospholipids and non-phosphorus containing lipids.

<figure id="fig::volcano">
<img src="figs/volcano.png" />
<figcaption> <em><strong>Differential gene expression and lipid
metabolism</strong></em> <strong>(A)</strong> Gene expression effect,
see genes in Table <a href="#tab:gene_effects" data-reference-type="ref"
data-reference="tab:gene_effects">[tab:gene_effects]</a>. A larger
proportion of <em>Inv4m</em> DEGs have no characterized genes, even
though the <em>Inv4m</em> signal has greater statistical significance
than the other two conditions tested in the experiment.
<strong>(B)</strong> Lipid metabolite change and statistic al
significance. Classic responses to senescence and phosphorus starvation
are observed in the lipid profiles. DGGA has been previously linked to
leaf senescence in arabidopsis. In response to phosphorus starvation
plants use a greater amount of glycoglycerolipids, and lower the
synthesis of phosphoglycerolipids. No consistent lipid response to
<em>Inv4m</em> was detected in leaves of sampled stages.
<strong>(C)</strong> Metabolic pathway annotation of top DEGs.
<strong>(D)</strong> Gene Ontology annotation of top DEGs.</figcaption>
</figure>

<figure id="fig::    jmj">
<img src="figs/jmjcluster.png" />
<figcaption><em><strong>Jmj cluster and coexpression
network</strong></em> <strong>(A)</strong> <strong>(B)</strong>
</figcaption>
</figure>

## *Inv4m* plants show perturbations of a cell proliferation gene expression network in leaves and have elongated shoot apical meristems. {#inv4m-plants-show-perturbations-of-a-cell-proliferation-gene-expression-network-in-leaves-and-have-elongated-shoot-apical-meristems. .unnumbered}

A gene co-expression network analysis (Panel A) was performed to
investigate the trans-regulatory effects associated with the Inv4m
chromosomal inversion in maize leaves. This network comprised
approximately 24 genes and 28 co-expression edges, filtered to represent
interactions between genes within the 15 Mb Inv4m inversion on
chromosome 4 and genes located elsewhere in the genome. Of these, 12
genes were located within the inversion, and 12 genes were located
outside the inversion. We observed that 10 genes within the network were
upregulated and 5 were downregulated in Inv4m plants. Gene ontology (GO)
enrichment analysis (Panel B) of this Inv4m trans-coexpression network
revealed an over-representation of terms related to cell population
proliferation and regulation of flower development (Fisher's exact test
p-values here). Within the enriched group for cell population
proliferation, we observed that the expression of mrlk (meristematic
receptor-like kinase), PCNA2 (encoding a proliferative cell nuclear
antigen 2 protein, part of a clamp complex increasing processivity of
DNA polymerase delta), and MCM5 (encoding a minichromosome maintenance 5
protein, a helicase) were all perturbed in Inv4m plants. Specifically,
mrlk and PCNA2 were \[Upregulated/Downregulated - choose based on your
data\] in Inv4m plants, while MCM5 was \[Upregulated/Downregulated -
choose based on your data\]. We specifically noted the perturbation of
zcn26, an FT homolog (a gene known for its role in flowering time)
located outside the inversion, which was \[Upregulated or
Downregulated - specify from your data\] in Inv4m plants. Furthermore,
we identified disruptions in genes involved in protein methylation,
specifically jmj2 and jmj21 (JMJ domain-containing proteins, shown in
the network), suggesting a potential role for epigenetic regulation. To
determine whether these observed gene expression changes had phenotypic
consequences, we measured the shape and height of the vegetative
meristems in control and Inv4m individuals. We found that Inv4m plants
have shoots with a significant ($p<0.05$, one-tailed *t-test*,
$H_a: \mu_{\textit{Inv4}} > \mu_{\textit{CTRL}}$) 9% increase in SAM
height (9.65 $\mu$, 0.37 $\mu\text{m}$ lower bound with 95% confidence
), a 9% increase in shape coefficient (2.49 $\text{nm}^{-1}$, 0.59
$\text{nm}^{-1}$), indicating a elongated meristem (Panel C).

<figure id="fig::DEGnetwork">
<img src="figs/DEG_network.png" />
<figcaption><em><strong>Trans coexpression network of <em>Inv4m</em>
DEGs.</strong></em> <strong>(A)</strong> Edge polarity reflects causal
precedence of the genetic perturbation by <em>Inv4m</em> introgression.
<strong>(B)</strong> </figcaption>
</figure>

<figure id="fig::model">
<img src="figs/model.png" />
<figcaption> <em><strong>B73 phenotypic and gene expression response to
<em>Inv4m</em> introgression</strong></em> </figcaption>
</figure>

However, the variation at *Inv4m* seems to have no major effect in the
gene response to phosphorus. The *Inv4m* x P interaction t-test (*FDR*
$<0.005$) shows only one gene with differential response to phosphorus
depending on the *Inv4m* karyotype. This gene is found in near *Inv4m*,
but outside previously reported limits
(Fig [3](#fig::RNAseq){reference-type="ref" reference="fig::RNAseq"} C).
A single gene *aldh2*, 3 Mb upstream of the reported *Inv4m* limit,
shows a different response to phosphorus depending on the *Inv4m*
genotype. In phosphorus defficiency *aldh2* overexpresses in the CTRL
plants, but it is not responsive in the *Inv4m* plant
(Fig [3](#fig::RNAseq){reference-type="ref" reference="fig::RNAseq"} D).

# Discussion

Chromosomal knob repeats have been previously associated with megabase
sized inversions in maize and *Arabidopsis*. For instance, the maize
abnormal chromosome Ab10 shows, relative to N10, two inversions, 12 Mb
and 15 Mb in size, each closely preceded by one of three TR-1 repeat
arrays corresponding to cytogenetic chromomeres
[@liu2020c; @mroczek2006]. Also, in the maize line Zi330, *in situ*
hybridization with knob probes revealed two bands on the long arm of
chromosome 2: one corresponding to a TR-1 and the other to a 180 bp knob
repeat array [@yang2012]. Here the researchers found that, despite being
an inbred line, Zi330 had homologous chromosomes with inverted knob
repeat band order, resulting in heterozygosity at this locus. In
*Arabidopsis thaliana* Col-0, the 750 kb *hks4* chromosomal knob,
consisting of pericentromeric heterochromatin repeat arrays and other
interspersed repeats [@schmidt2020; @fransz2000], overlaps a 1.2 Mb
paracentric inversion with respect to Ler and *A. lyrata* [@zapata2016].

The effect of *Inv4m* in flowering time, previously reported as an
association study hit
[@romero_navarro2017-cn; @barnes2022; @gates2019-xu], was confirmed by
the experimental introgression of the isolated inverted haplotype into a
single genetic background (B73). Additionally, we observed that in this
particular background *Inv4m*, increases plant height at flowering.
These two effects are additive with the corresponding effects of
phosphorous treatment. The response of the experimental plants to
phosphorus was observed in more phenotypes and was of greater magnitude
than the effect of *Inv4m*(Fig [2](#fig::effects){reference-type="ref"
reference="fig::effects"}). In phosphorus deficiency the plants showed
delayed growth as both a later time to half maximum stover dry weight
and longer times to male and female flowering. Plants grew smaller in
the absence of phosphorus, accumulated less vegetative biomass as stover
dry weight at harvest, and had smaller kernels (KW50)
(Fig [2](#fig::effects){reference-type="ref" reference="fig::effects"}).
No significant *Inv4m*$\times$ P interaction was detected ($\beta$
t-test, $\textrm{\textit{FDR}} > 0.05$ for all phenotypes). Regarding
leaf transcriptome, the master regulator of the phosphorus starvation
response *PILNCR1-mir399* was over-expressed in phosphorus deficiency
independent of *Inv4m* genotype
(Fig [3](#fig::RNAseq){reference-type="ref" reference="fig::RNAseq"} D).
*PILNCR1-mir399* has been reported to be overexpressed in low
phosphorus, and in particular it is more responsive in P-inefficient
than in P-efficient genotypes [@du2018]. B73 has been reported to be
phosphorus inefficient with respect other inbreds like Mo17[@zhu2005],
and has been selected for yield under high input agriculture compared
with traditional varieties such as Mi21. However, *PILNCR1-mir399*
response to phosphorus in the experimental plants shows no dependency on
the variation at *Inv4m*. This suggests that the master regulator
phosphorus starvation response is not regulated by *Inv4m* in the
experimental conditions. An alternative explanation is that the
experimental P limitation was not hard enough to show a
genotype-dependent difference in response. No observable difference in
leaf anthocyanin accumulation or lower leaf senescence in response to
phosphorus deficiency was evident in the field.

We observed a clear over-representation of nitrogen-associated
biological processes and metabolic pathways in the 7373 genes responsive
to phosphorus treatment (Fisher exact test GO over-representation
adjusted *FDR* \< 1e-6). However, when we reduce the list to the top
1000 most significant, the over-represented gene processes and terms are
more related to phosphorus (Fisher exact test GO over-representation
adjusted *FDR* \< 1e-6). This suggests that the genes that show the
greater and/or more consistent effect have phosphorus-associated
functions while nitrogen-related genes make significant but smaller
adjustments in expression, maybe reflecting cross-talk between
phosphorus and nitrogen homeostasis[@torres-rodriguez2021]. I could not
find signals of photosynthesis or carbohydrate over representation in
the *Inv4m* responsive to inversion karyotype, as was previously found
in response to cold [@crow2020].

Overall, there is no evidence from this transcriptome experiment that
*Inv4m* contributes to local adaptation to phosphorus, in the sense that
its adaptive contribution does not depend on phosphorus amount in the
the soil. However, it is clear that *Inv4m* plants flower earlier and
are taller irrespective of the phosphorus status.

# Materials and methods {#materials-and-methods .unnumbered}

## *Inv4m* delimitation, breakpoints, and knob repeat analysis {#inv4m-delimitation-breakpoints-and-knob-repeat-analysis .unnumbered}

The genome sequences of TIL18 (Zx-TIL18-REFERENCE-PanAnd-1.0), PT
(Zm-PT-REFERENCE-HiLo-1.0), B73 (Zm-B73-REFERENCE-NAM-5.0), and
corresponding gff annotations of genes and repeats were obtained from
<https://download.maizegdb.org/> [@portwood2019]. We used Anchorwave
[@song2022] to align chromosome 4 of the inverted karyotypes, TIL18 and
PT, to chromosome 4 of the reference standard karyotype B73. We ran the
`genoAli` command with the `-IV` flag to select the optimized algorithm
for aligning inversions and translocations. The resulting `MAF` file
partitioned the chromosome into anchor delimited alignments, but
sections lacking alignment remained. Using LAST [@kielbasa2011] with
default parameters, we aligned each breakpoint segment sequence with
itself and then made dotplots of the similarity results with the
`last-dotplot` script for the six breakpoints in the three genomes. We
ran the `megablast` command, default parameters, with the breakpoint
sequences as queries against GenBank Accessions M32524.1 [@dennis1984]
representing the 180 bp knob and AY083942.1 (366 bp) [@hsu2003]
representing TR-1. For plotting results we normalized high scoring pair
bitscores ($x$) as: $$\begin{aligned}
\label{eq:repeatscore}
\text{Normalized match score} = 100\times\frac{x - x_{min}}{x_{max} - x_{min}}
\end{aligned}$$

## *Inv4m* Near Introgressed Lines, growth conditions, experimental design, and phenotype measurements {#inv4m-near-introgressed-lines-growth-conditions-experimental-design-and-phenotype-measurements .unnumbered}

To measure the effects of the *Inv4m* in plant field phenotypes and
their phosphorus starvation response transcriptome, we used a highland
traditional variety carrying the Highland haplotype of *Inv4m*
corresponding to the inverted karyotype. The accession Michoacán 21
(referred to as Mi21), from the Mexican Cónico group, was obtained from
the International Maize and Wheat Improvement Center (CIMMYT). In
contrast, the reference genome of the temperate inbred B73, the
recurrent parent for introgression, carries the lowland haplotype
corresponding to the standard non-inverted karyotype at *Inv4m*. From
the cross of Mi21 with B73 one F1 individual was backcrossed to B73 for
six generations. We selected lines carrying *Inv4m* with a diagnostic
SNP during each cycle using a cleaved amplified polymorphic sequence
(CAPS) marker. The marker SNP is PZE04175660223 located at position
4:181637780 in the NAM B73v5 *Zea mays* genome assembly. Amplification
of the polymorphic site was done with the following primer pair:
`CTGAGCAGGAGATGATGGCCACTC` and `GGAAAGGACATAAAAGAAAGGTGCA`, and
subsequently cleaved by *HinfI*. Plants were genotyped using the CASP
marker for selecting heterozygous plants at BC6S2 after selfing seeds of
*Inv4m* and CTRL homozygous individuals were selected for the field
trial.

Plants were planted on May 26 2022 at the Russell E. Larson Agricultural
Research Farm in Rock Springs, Pennsylvania (40°42'36\" N 77°57'0\" W,
366 m.a.s.l.) in soil classified as a Hagerstown silt loam (fine, mixed,
semiactive, mesic Typic Hapludalf). Experimental conditions were similar
to previously described [@strock2018]. The experiment had a complete
block design with two phosphorus (P) levels. Low-P fields (5 ppm
Melich-3 Phosphorus) and high-P fields (36 ppm Melich-3 Phosphorus) were
divided into smaller blocks. Three rows per block were planted with a
mean stand count of 8 plants per plot, and the plants from the center
row were selected for measurements to avoid border effects. Fields
received fertilization based on treatment requirements. Drip irrigation
was provided during dry periods. Each genotype was replicated four times
within its P treatment.

## Inductively Coupled Plasma Mass Spectrometry (ICP-MS) {#inductively-coupled-plasma-mass-spectrometry-icp-ms .unnumbered}

We followed the protocol described in [@baxter2014] for ionome analysis
of the flag leaf by Inductively Coupled Plasma--Mass Spectrometry.
Briefly, leaf tissue was grounded and then transferred to 110 mm
borosilicate glass test tube for digestion with 2.5 mL concentrated
nitric acid overnight at 105 °C. After two serial dilutions of this
digestion in ultra pure water (18.2M $\Omega$ Milli-Q system, Millipore)
1:4, and 0.9:4.1, a final aliquot of 1.2 mL was loaded into 96 well
autosampler trays. The samples were processed in an Agilent 7500 ICP-MS
system, and the signal was adjusted to machine drift. Collision data was
collected for Al, B, Ca, Fe, K, Mg, Mn, Mo, Ni, P, S, Zn, Na, and Cu.
The final concentrations of minerals were calculated as mg of element
per kg dry seed weight. Measurements beyond three times the
interquartile range were deemed as outliers and discarded for the
following analysis [@baxter2013]. Additionally, quantifications for Al,
B, Mo, Ni, Na, Cu, were discarded because the concentrations were too
close to the detection limit of the mass spectrometer as previously
reported [@baxter2013].

## Phenotype analysis {#phenotype-analysis .unnumbered}

For stover mass growth curves, a different plant at each time point 40,
50, 60, and harvest, 121 days after planting (DAP), was collected,
dried, and weighed for the same row. Stover dry mass data was fitted to
a logistic growth model using the R package `Growthcurver`
[@sprouffske2016]. Maximum Stover dry weight was estimated to be the
maximum over the four-time points and not dry weight at harvest. Ear
measurements were taken for one ear per row at harvest. We modeled the
individual phenotypes with the `R nlme` function as the response
variable in a mixed effects model with spatial structure. For each
phenotype $y$ we have:

$$\begin{aligned}
\label{eq:pheno_model}
y_{ijkr} = \beta_{0} + \beta_{1}\text{P}_i + \beta_{2} \textit{\textit{Inv4m}}_j + \beta_{3}[\textit{\textit{Inv4m}} \times \text{P}]_{ij} + u_k + \varepsilon_{ijkr}
\end{aligned}$$

Where the phenotype observation $y_{ijkr}$ corresponds to the plant $r$
in phosphorus treatment $i$ with genotype $j$ in block $k$. The fixed
effects coefficients $\beta_{0}$ for the overall mean, $\beta_{1}$ for
the effect of phosphorus treatment $i$, $\beta_{2}$ for the fixed effect
of genotype $j$. One random effect of the block $k$
$(u_{k}) \sim N(0, \sigma_u^2)$, and the residuals
$(\varepsilon_{ijkr}) \sim N(\mathbf{0}, \sigma_\varepsilon^2 \mathbf{\Sigma})$.
We added a correlation structure of the residuals $\mathbf{\Sigma}$
given by a spherical model [@pinheiro2000]:

$$\begin{aligned}
\label{eq:sp_model}
\Sigma_{pq} = \begin{cases}
1 & \text{if } d_{pq} = 0 \\
1 - \frac{3d_{pq}}{2\rho} + \frac{d_{pq}^3}{2\rho^3} & \text{if } 0 < d_{pq} < \rho \\
0 & \text{if } d_{pq} \geq \rho
\end{cases}
\end{aligned}$$

Where $d_{pq}$ is the Euclidean distance between plots $p$ and $q$ in
the field (based on row and column coordinates), and $\rho$ is the range
parameter of the spherical model. We corrected for multiple hypotheses
($H_0: \beta = 0$) by reporting *t-tests* for the fixed effects below
*FDR* = 0.05.

## Tissue sampling, RNA extraction, and sequencing {#tissue-sampling-rna-extraction-and-sequencing .unnumbered}

We sampled the plants at 63 DAP when we estimated them to be between v10
to v12 developmental stages. We took tissue from he first leaf with a
fully developed collar, or first leaf before the flag leaf, and every
other leaf below for a total of four sampled leaves per plant. These
leaves were numbered sequentially from 1 (most apical) to 4 (most
basal). We used four replicate plants per combination of P treatment and
*Inv4m* genotype for a total of 64 tissue samples. We took ten disc
samples from the leaf tips with a tissue puncher and immediately froze
the tissue in 1.5 mL tubes with two steel beads precooled with liquid
nitrogen and kept in dry ice until stored at -80°C. We extracted total
RNA with the QIAGEN RNAeasy Plant Mini Kit RNA extraction kit following
manufacturer procedures (QIAGEN 74904), and RNA samples were quantified
in nanodrop and sent to the NCSU Core Genomics Laboratory for
sequencing. Following QC in Bioanalyzer, Illumina libraries were
prepared and sequenced in a lane of Novaseq according to manufacturer
recommendations.

## Plant genotyping {#plant-genotyping .unnumbered}

We followed [@brouard2022] for GATK-based RNAseq genotyping of 15 plant
samples represented by 60 leaf libraries. Briefly, Illumina short reads
were mapped to the NAM5 Zea mays B73 genome [@hufford2021] using `STAR`
[@dobin2013], then we marked duplicates in the resulting BAM alignments,
split reads at intron-exon junctions and recalibrated sequence quality
per leaf library. At this point, we used HaploytypeCallerfor for
generating gvcfs per plant identified by field row id ($\sim 4$
libraries per plant). We did joint sample genotyping afterward with
`genotypeGVCFs`. Then we filtered for variant quality (
`window 35, cluster, QD < 2.0, FS > 30.0, SOR > 3.0, MQ < 40.0`) for the
genotypes and $50\%$ marker completion for individuals. This resulted in
200000 markers with $85\%$ complete data for 13 plants. Finally, we used
TASSEL5 K Nearest Neighbour imputing, producing a matrix of 19668
markers at $99.84\%$ completion. Shell scripts are available at the [
github repository](https://github.com/sawers-rellan-labs/\invfourRNA)

## Lipid extraction, identification, and quantification by HPCL-MS {#lipid-extraction-identification-and-quantification-by-hpcl-ms .unnumbered}

We used the lipid extraction by methyl-tert-butyl ether [@matyash2008]
and HPLC-MS methods as described in [@barnes2022]. First we ground the
frozen tissue samples using a SPEX Geno/Grinder (Metuchen, NJ, USA).
Then we added cold methanol (MeOH) to each sample, vortexed and kept
them on ice. Next, we added methyl tert-butyl ether (MTBE) and vortexed
again. After shaking the samples at 4°C, we added LC/MS grade water at
room temperature and vortexed. We then centrifuged the samples and
collected the supernatant from the upper organic phase, splitting it
into two aliquots for lipid profiling and pooling. We dried the samples
using a speed vacuum and resuspended them in a MeOH-Toluene mixture with
an internal standard. After vortexing and sonication, we transferred
aliquots into amber glass vials. We used Agilent UHPLC-QTOF MS/MS for
ionization and detection. We conditioned the columns by running \"no
sample injections\" followed by samples, pools, and blanks. We injected
the samples into the UHPLC-QTOF MS/MS in positive electrospray
ionization mode, maintaining specific source parameters. We corrected
retention times using Agilent MassHunter and Excel by extracting ion
chromatograms of internal standards, fitting a polynomial regression,
and calculating new retention times using a lipid library. For lipid
identification, we used MSDIAL, **ALLISON MS DIAL PROTOCOL GOES
HERE!!!!** MS/MS settings, and a post-identification file for accurate
m/z and retention times. We converted raw data to .abf format using the
Reifycs Abf converter, and filtered the alignment results based on
compound intensity. We curated data for duplicates, isotopes, and
ion-adducts using MS-FLO, and the curated data using the sum of all
known metabolite signals (mTIC). Finally, we used the lipid peak area
for further analysis.

## Differential gene expression and differential lipid analysis {#differential-gene-expression-and-differential-lipid-analysis .unnumbered}

We aligned reads to the maize Zm-B73-REFERENCE-NAM-5.0 genome using
`kallisto` [@bray2016]. The alternative transcript alignment was turned
into counts per gene per MB. We used `voom` to calculate variance
according to gene expression levels and counts were converted to
$log_2(\text{CPM})$. Lipid analysis followed a similar workflow, where
we calculated variance weights with `voom` for each lipid MS-spectra
peak area and transformed it to $log_2$ scale. We made a multivariate
multiple regression for gene expression and lipid MS signal separately
using `limma` [@ritchie2015]. For the log transformed expression/signal
$Y_{ijrs}$, from leaf $s$, in plant $r$, under phosphorus treatment $i$,
with genotype $j$, we have:

$$\begin{aligned}
% \label{eq:expression_model}
\begin{aligned}
Y_{ijrs} = {}& \beta_0 + \beta_{1}\text{Row}_l + \beta_{2}\text{Column}_{m} + \beta_3 \text{Leaf}_{s} +\beta_4 \text{P}_{i} + \beta_5 \textit{\textit{Inv4m}}_{j}+ \beta_6 [\text{P} \times \textit{\textit{Inv4m}}]_{ij} \\
& + \varepsilon_{ijrs}
\end{aligned}
\end{aligned}$$ with residuals: $$\begin{aligned}
\varepsilon_{ijrs} \sim \mathcal{N} (0,\phi\sigma^2)
\end{aligned}$$ We used the leaf stage ($\text{Leaf}_{s}$) as a
numerical variable with $s \in \{1,2,3,4\}$ instead of categorical. This
implies that $\beta_3$ represents the rate of change of expression with
increasing leaf stage number, while the rest of the coefficients were
defined as categorical in the same way as in equation
[\[eq:pheno_model\]](#eq:pheno_model){reference-type="eqref"
reference="eq:pheno_model"} and
[\[eq:sp_model\]](#eq:sp_model){reference-type="eqref"
reference="eq:sp_model"}. We adjusted the p-values for the t-tests of
the linear model coefficients as false discovery rates and genes whose
effect had a $FDR <0.05$ were deemed to be differentially expressed. For
phosphorus treatment ($\beta_4$) and *Inv4m* genotype ($\beta_5$) we
considered genes with an effect of $|log_2(\text{Fold Change})| >2$ as
top DEGs. In the case of the leaf effect, a gene was considered top DEG
if $|log_2(\text{Fold Change})|>0.7$, i.e. $>2.1$ between leaf stage 1
and leaf stage 4. R scripts and expression data are available at the [
github repository](https://github.com/sawers-rellan-labs/\invfourRNA).

## Gene Ontology and KEGG overrepresenattion analysis {#gene-ontology-and-kegg-overrepresenattion-analysis .unnumbered}

Once we had sets of differentially expressed genes for the three
predictors (leaf, -P, *Inv4m*) and two types of gene expression response
(upregulated and downregulated), we proceeded to annotate them with gene
ontology terms and KEGG pathways using `ClusterProfiler`
[@yu2012; @zicola2024]. We started with the B73 NAM v5 gene ontology
annotation from [@fattel2024] and added GO terms for each intermediate
node in the gene ontology tree using the `ClusterProfiler` function
`buildGOmap`. Then we conducted gene over-representation analysis with
the function `compareCluster`, using as universe/background the set of
24011 genes detected in at least one good quality leaf RNAseq library.
This function calculates the hypergeometric test for overrepresented
ontology terms in the specified gene set and returns raw, and
FDR-adjusted p-values. We then manually reviewed the combined 1700
significant (*FDR* $<0.05$) overrepresented GO term associations for the
6 predictor/regulation combinations, and we selected for illustration an
*ad hoc* subset with low semantic redundancy. Similarly, We tested for
KEGG pathways over representation using the `enrichKEGG` function from
`compareCluster`, which makes the same hypothesis tests on the NCBI
REFseq annotation of the B73 NAM assembly. Both types of
overrepresentation analysis were plotted with the package `DOSE`
[@yu2015].

## *Inv4mTrans* coexpression network filtering {#inv4mtrans-coexpression-network-filtering .unnumbered}

We then looked into further details of the biological functions of the
genes that are differentially expressed in the *Inv4m* plants with
respect to the controls. To isolate the effect of *Inv4m* in the gene
coexpression network of B73 we had to restrict our analysis to the
pairwise relationships between genes inside *Inv4m* and genes outside
the shared introgression while discounting the effects of the genes in
the flanking introgression. With this purpose, we started by defining a
set of coexpressed pairs of *Inv4m* DEGs where one was located in the
shared introgression and the other outside, and where the coexpressed
pairs showed a statistically significant correlation (*FDR*$<0,05$,
Pearson correlation *t-test*). Then we filtered from this set the
coexpression pairs that have not been previously reported for B73. For
this step, we used as the validated reported dataset the gene
coexpression network for B73 published in MaizeNetome [@feng2023]. This
network is based on gene expression quantified in 31 different tissues
and stages and contains significant pairwise gene expression
correlations detected with WGCNA [@langfelder2008]. As this network was
calculated with the `Zm-B73-REFERENCE-GRAMENE-4.0` genome we had to
uplift the v5 introgression coordinates (4:157012149-195900523) into v4
(4:155195539-193981713) for retrieval using the
[MaizeNetome](http://minteractome.ncpgr.cn/elementfetch.php) website.
Furthermore, we had to split the introgression in three intervals
corresponding to the left flank (4:155195539-170747954),
*Inv4m*(4:170747955-186145605) and right flank (4:186145606-195900523),
so each partition had less than the 500 genes allowed download limit.
From a total of 14248262 reported edges in the database, we found 640599
edges linking the expression of 717 genes inside the shared
introgression limits with 23748 genes throughout the genome. Then we
intersected this reported set of edges for B73 with our experiment's set
and obtained a network of 77 edges pairing 29 DEGs in the shared
introgression with 18 DEGs outside. The subgraph of the trans
coexpression network lof the *Inv4m* DEGs consisted of 30 edges linking
12 DEGs located in the *Inv4m* proper with 20 genes outside the shared
introgression (Figure [6](#fig::DEGnetwork){reference-type="ref"
reference="fig::DEGnetwork"} A). The trans coexpression network for the
flanking introgression DEGs consisted of 47 edges linking 17 genes
located in the flanking introgression segments with 13 genes outside the
shared introgression. At this point we made a Gene Ontology
overrepresentation analysis of the DEGs present in each set and compared
the annotation to the annotation of DEGs present in the *trans*
coexpression network of the flanking introgression DEGs, and of the
shared introgression taken as a whole (Figure
[6](#fig::DEGnetwork){reference-type="ref" reference="fig::DEGnetwork"}
B) following the procedure in the previous section. A final filtered
subnetwork of 4 DEGs within the *Inv4m* limit and 5 genes outside the
shared introgression showed no influence of *Inv4m* flanking DEGs
(Figure [6](#fig::DEGnetwork){reference-type="ref"
reference="fig::DEGnetwork"} A, dashed lines)

## Filtering of *Inv4m* DEGs by phenotype association {#filtering-of-inv4m-degs-by-phenotype-association .unnumbered}

As our data showed evidence of *Inv4m* accelerating flowering time and
increasing plant height, we put together a list of candidate genes
associated with these two phenotypes to tease out which DEGs were likely
contributors to the observed *Inv4m* effect in these traits. For
flowering time, we started with the list of 991 genes compiled by
[@wang2021] and 62 genes from [@li2023a]. Then we downloaded the maize
data from the GWAS atlas [@liu2023]
(`gwas_association_result_for_maize.txt.gz`) and selected genes that
overlapped association SNPs for the [Plant Phenotype and Trait
Ontology](https://ngdc.cncb.ac.cn/gwas/browse/ontology) term "days to
flowering trait\" `PPTO:0000155`. For this and the following candidate
gene list, we considered that a gene overlapped an association SNP if
the SNP was located within the 5 kb extended range of the gene model,
i.e. as described in the gff gene annotation $\pm 5$ kb. The final
source of associations for flowering time was the phenotypic plasticity
study in [@tibbs-cortes2024] from which we used 281 genes with
significant GWAS SNPs in the columns `DTS_slope`, `DTS_intcp`,
`DTA_slope`, `DTA_intcp`. For plant height, 27 genes from [@liu2023],
1210 genes with GWAS Atlas associations for the term "plant height\"
`PPTO:0000126`; and 39 genes overlapping phenotypic plasticity
association SNPs for `PH_slope` and `PH_intcp` [@tibbs-cortes2024]. The
final nonredundant list consisted of a total of 2224 candidate genes for
flowering time and 1272 candidates for plant height.

## Meristem clearing and size quantification {#meristem-clearing-and-size-quantification .unnumbered}

For vegetative meristem size quantification, maize seedlings were grown
at the North Carolina State University Phytotron in a Percival Model
LT-105 growth chamber (conditions: 29.4°C day/23.9°C night, relative
humidity 50%, 16 hours light/8 hours dark, light intensity of 412 $\mu$
mol at plant height; soil type: 1:1 Sun Gro Propagation Growing Mix
\[Canadian Sphagnum peat moss 50-65%, vermiculite, dolomitic lime,
0.0001% silicon dioxide\] : cement sand) in 24-well trays. Seedlings
were watered daily in the morning daily and fertilized three times per
week. Two-weeks after planting, seedlings were cut at soil level and
again 1cm above the soil cut. This 1cm tissue cassette of the shoot apex
was cut longitudinally in the medial plane, in a midrib to margin
orientation, by hand with a razor blade. Tissue was fixed in ice cold
and fresh FAA \[50% EtOH, 35% milliQ water, 10% formaldehyde (35%), 5%
glacial acetic acid (v/v)\] with a vacuum for 15 min. FAA was replaced,
and tissues were place overnight at 4°C on a rocker. For clearing,
tissues were then removed from FAA and dehydrated through a graded
ethanol series at room temperature for an hour each with gentle shaking:
50%, 70%, 85%, 95% EtOH (v/v)\], followed by 1:1 95% EtOH and Methyl
Salicylate, and finally, 100% Methyl Salicylate. Once in 100% Methyl
Salicylate, samples were left to shake at room temperature overnight.
Each shoot apex tissue cassette was placed on a microscope slide and
covered with a coverslip. Cleared shoot apices were imaged with
differential internal contrast using a Leica DM4B microscope equipped
with a DMC6200 digital camera. Each shoot apex image was measured using
ImageJ v2.14.0. The scale was properly set for each image. Width was
measured as a straight line at 0° anchored from the edge of P0. Height
was measured from the highest point of the meristem tip to the width
line at 90°. Surface area was measured using the polygon selection tool
by taking into consideration the width line and marking the edge of the
meristem dome. Data were analyzed using ANOVA and plots were generated
in R v.4.3.2 with the packages rstatix v.0.7.2, readxl v.1.4.3, ggplot2
v.3.5.1 and ggpubr v.0.6.0.

# Supporting information {#supporting-information .unnumbered}

![***Chromosome Distribution of the SNPs in the RNAseq samples, their
correlation with *Inv4m*, and genotype run length as percentage of the
genome.*** **(A)** Chromosome 4 contains 44% of all genotyped SNPs and
98% of *Inv4m* correlated SNPs (n=13). **(B)** Manhattan plot for the
significance (*t-test*) of the Pearson correlation of 19861 genotyped
SNPs with PZE04175660223. Red line: *FDR* $=0.05$ threshold. **(C)**
Pearson correlation ($r$) of the each SNP with the *Inv4m* tagging SNP
PZE04175660223, open circles: non-significant SNPs. SNPs with
significant correlation with *Inv4m*(full circles) are divergent
highland introgressions. **(D)** NIL genotype run length as percentage
of the genome. Each point represents a NIL, *Inv4m* plants are tagged by
the highland allele of PZE04175660223, CTRL plants by the reference
(B73) allele, dashed line is the overall mean for the 13
NILs.](figs/SNP_distribution.png){#fig::SNPdistro width="\\linewidth"}

<figure id="fig::growth">
<img src="figs/growth.png" />
<figcaption><strong><em>Plant growth curves show effect of phosphorus
treatment but no <em><em>Inv4m</em></em> effect.</em></strong>
<strong>(A)</strong> Stover dry weight observations.
<strong>(B)</strong> Stover dry weight, logistic growth fit curves per
field row. <strong>(C-F)</strong> Effects of phosphorus treatments on
logistic growth curve parameters, selected for their significance (<span
class="math inline">$\textrm{\textit{FDR}} &lt; 0.05$</span>) in a
linear mixed effect model, as pairwise comparisons.<br />
Student’s <em>t-test p-value</em>:* <span
class="math inline"> &lt; 0.05</span>, ** <span
class="math inline"> &lt; 0.01</span>, *** <span
class="math inline"> &lt; 0.001</span>, **** <span
class="math inline"> &lt; 0.0001</span>. </figcaption>
</figure>

<figure id="fig::ionome">
<img src="figs/ionome.png" style="width:95.0%" />
<figcaption><em><strong>Effect of Phosphorus Treatment in the Maize
Ionome.</strong></em> <strong>(A)</strong> Significant effects (linear
mixed model, <span class="math inline">$\textrm{\textit{FDR}} &lt;
0.05$</span>) of soil phosphorus treatment on mineral content per
tissue, and seed/stover mineral ratio. <strong>(B-G)</strong> Pairwise
comparisons between phosphorus treatments corresponding to the effect in
(A).<br />
Student’s <em>t-test p-value</em>:* <span
class="math inline"> &lt; 0.05</span>, ** <span
class="math inline"> &lt; 0.01</span>, *** <span
class="math inline"> &lt; 0.001</span>, **** <span
class="math inline"> &lt; 0.0001</span>. </figcaption>
</figure>

<figure id="fig::lipid_all_leaves">
<img src="figs/lipid_all_leaves.png" style="width:95.0%" />
<figcaption><em><strong>Effect of <em>Inv4m</em>, Phosphorus Treatment
and Interaction in lipids, all leaves.</strong></em> </figcaption>
</figure>

<figure id="fig::lipid_old_young">
<img src="figs/lipid_old_young.png" style="width:95.0%" />
<figcaption><em><strong>Effect of <em>Inv4m</em>, Phosphorus Treatment
and Interaction in lipids, old leaves vs young leaves.</strong></em>
</figcaption>
</figure>

# Acknowledgments {#acknowledgments .unnumbered}
