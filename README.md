# Acetate utilization over temperature gradient in Potter Cove sediments
As supplemental to submitted manuscript.\
Slurry RNA stable isotope probing incubations were conducted with <sup>13</sup>C-acetate over a temperature gradient.\
In this repository, the code is stored which was used to analyse and plot data acquired during the incubation and especially the analysis of the amplicon sequencing data and subsequent analyses and plotting.
<br/><br/>
## Measurements during incubation period
### Geochemical measurements
Dissolved ferrous iron and sulfate were measured at start and end point of incubations.\
Script [`geochemistry`](geochemistry.R) is used to calculate the accumulation rate and plot.

### <sup>13</sup>C-CO<sub>2</sub> in headspace
Isotope proportion of <sup>12</sup>C- and <sup>13</sup>C-CO<sub>2</sub> were measured in the incubation headspace over incubation period.\
Script [`delta13CO2`](delta13CO2.R) is used for plotting.
<br/><br/>
## Stable isotope probing
### Fractionation profiles
After incubation, RNA was extracted from incubations and density separated by centrifugation. Density and RNA concentration was measured in the resulting fractions.\
Script [`fractionation_profiles`](fractionation_profiles.R) is used to plot this data and mark which of the fractions were used for amplicon sequencing.

### quantitative PCR
Samples from the SIP experiment also used for amplicon sequencing were used for quantitative PCR (qPCR). 
The functional markergene *dsrA* for dissimilatory sulfate reduction and the bacterial 16S rRNA gene for quantification of bacteria in general were amplified.\
Script [`qPCR_calculation_plot`](qPCR_calculation_plot.R) is used to calculate copy numbers per ng used cDNA and plot it.

### Amplicon sequence data analysis - workflow
Involves all scripts starting with `S_` and some scripts which are sourced during this analysis in folder [Scripts to source](Scripts_to_source). \
The workflow is based on the pipeline by [C. Hassenrueck](http://doi.io-warnemuende.de/10.12754/misc-2022-0002) following the [dada2 pipeline](https://benjjneb.github.io/dada2/index.html).

The scripts provided here contain the used parameters, so when the rawdata is downloaded from ENA it can be directly used to replicate the analysis process.\
Raw sequence data can be found [here](https://www.ebi.ac.uk/ena/browser/view/PRJEB82428)

1. bash script [`S_01_dada2_seqprep`](S_01_dada2_seqprep.bash) for quality control, demultiplexing and primer clipping \
Mapping files for sequence analysis can be found [here](small_data)
2. R script [`S_02_dada2`](S_02_dada2.R) for sequence trimming, error correction, dereplication, denoising,  read merging, chimera removal and taxonomic classification
3. R script [`S_03_metacoder_import_data`](S_03_metacoder_import_data.R) for importing ASV table created in previous script into a `taxmap` format, remove mitochondrial and chloroplast sequences and doubletons and singletons\
   preparing RData objects and tables for plotting and dbRDA analysis
4. R script [`S_04_rarefaction-curves`](S_04_rarefaction-curves.R) for calculating and plotting rarefaction curves and removing samples with too low coverage if necessary
5. R script [`S_05_barplots`](S_05_barplots.R) for reformatting data and plotting the taxonomy in overview barplots\
   creating data output to use for line plots of taxa enriched in heavy fractions
6. R script [`S_06_plot_enr.taxa`](S_06_plot_enr.taxa.R) for calculating which ASVs incorporated the 13C-label and plotting
7. R script [`S_07_dbRDA`](S_07_dbRDA.R) for calculating distance-based redundancy analyses on complete and subsetted data sets and plotting

Files with some metadata needed to run the scripts can be found [here](small_data) and in the supplementary of the manuscript.

