# This script imports taxmap objects created in script metacoder_import_data.R produced through dada2 pipeline

# packages ####
library(tidyverse)
library(taxa)
library(metacoder) # version 0.3.7
library(phyloseq)
library(hues)
library(cowplot)
library(devEMF)
library(patchwork)


# import data ####
setwd("WORKINGDIRECTORY")

## path to scripts to load
sloc <- "SCRIPTLOCATION/Scripts_to_source/"

## project name to name files
prj <- "PROJECT"

## taxmap object from metacoder_import_data R script, samples with insufficient coverage removed
load(paste0("taxmap_", prj, ".RData"))

## metadata
mdata <- read_tsv("MDATAFILE.txt")


# create phyloseq object ####
objphy <- as_phyloseq(obj, otu_table = "asv_table", otu_id_col = "asv_id", sample_data = mdata, sample_id_col = "Name.Mapping.File")

ntaxa(objphy)                  # how many taxa (including subtaxa)
nsamples(objphy)               # how many samples
sample_names(objphy)           # check sample names
sample_variables(objphy)       # sample variables from mdata
otu_table(objphy)              # ASV count table
tax_table(objphy)              # taxonomy for ASVs
taxa_names(objphy)             # ASV IDs

rank_names(objphy)             # check rank names (need to be Kingdom, Phylum etc.)

## calculate relative abundance
opr <- transform_sample_counts(objphy, function(x) x/sum(x))

## melt phyloseq object into dataframe
phylo_melt <- psmelt(opr)

### save object phylo_melt used to perform all following operations and remove not needed objects
save(phylo_melt, file = paste0(prj, "_phylomelt.RData"))
rm(obj, objphy, opr)
load(paste0(prj, "_phylomelt.RData"))

# rename taxa depending on max abundance ####
## use function 'sort_abundant_taxa'
source(paste0(sloc, "filter_taxa_above_threshold_filt_others.R"))

## taxa above 5% in at least one samples
# ta5 <- sort_abundant_taxa(table = phylo_melt, abundance_threshold = 5)
# 47 ASVs, 36 genera, 30 families, 15 orders, 11 classes, 8 phyla
## missing unclassified Geobacteraceae

## taxa above 2% in at least one samples
# ta2 <- sort_abundant_taxa(table = phylo_melt, abundance_threshold = 2)
# 84 ASVs, 64 genera, 54 families, 32 orders, 22 classes, 17 phyla

## taxa above 4% in at least one samples
ta4 <- sort_abundant_taxa(table = phylo_melt, abundance_threshold = 4)
# 52 ASVs, 39 genera, 33 families, 18 orders, 13 classes, 8 phyla

# use 4% threshold, need uncl. Geopsychrobacteraceae
rm(ta5, ta2)

# sort and make nice names for taxa ####
## other columns needed as factor for plotting
cols_factor <- c("Temp", "Isotope", "Treatment")

## decide on rank level to plot: Genus and sorted by Class
View(subset(ta4$ASV_table_taxa_abv_thr, ta4$ASV_table_taxa_abv_thr$Sample == "Bacteria.2degrees.13C.Acetate.Fraction.Ultra-light"))

## reorder taxa for appearance in plot
class_unique <- (unique(ta4$ASV_table_taxa_abv_thr$Class_abt)) # 15 taxa to sort by
class_ordererd <- moveMe(class_unique, "other_Bacteria_<4% last")

genus_unique <- (unique(ta4$ASV_table_taxa_abv_thr$Genus_abt)) # 54 taxa to plot
genus_reorder <-  moveMe(genus_unique, "other_Bacteria_<4% last; 
                         other_o_Desulfobacterales_<4% after other_f_Desulfosarcinaceae_<4%;
                         Desulfuromonadales_f_unclassified after Sva1033_g_unclassified;
                         other_p_Desulfobacterota_<4% after Syntrophotalea;
                         other_c_Clostridia_<4% after Peptostreptococcales-Tissierellales_f_unclassified;
                         other_p_Planctomycetota_<4% after other_f_Pirellulaceae_<4%;
                         other_c_Gammaproteobacteria_<4% after Woeseia")
genus_renamed <- genus_reorder %>%  gsub("R76-B128", "R76-B128 (Kiritimatiellae)", .) %>% 
  gsub("other_o_Verrucomicrobiales_<4%", "o_Verrucomicrobiales", .) %>% 
  gsub("SG8-4_g_unclassified", "SG8-4_g_unclassified (Phycisphaerae)", .) %>% 
  gsub("SEEP-SRB4", "SEEP-SRB4 (Desulfocapsaceae)", .) %>% 
  gsub("Sva1033_g_unclassified", "f_Sva1033 (Desulfuromonadales)", .) %>% 
  gsub("Bacteroidetes BD2-2_g_unclassified", "f_Bacteroidetes BD2-2", .)
 

ta4n <- select(ta4$ASV_table_taxa_abv_thr, OTU, Sample, Class_abt, Genus_abt, Abundance, Treatment, Isotope, Fraction, Temp, Class,
               Genus_uncl, Family_uncl, Order_uncl, Class_uncl, Phylum_uncl) %>% 
  mutate(Class_abt = factor(Class_abt, levels = class_ordererd,
                          labels = class_ordererd),
         Genus_abt = factor(Genus_abt, levels = genus_reorder,
                          labels = genus_renamed), 
         Fraction = factor(Fraction, levels = c("Ultralight", "Light", "Midpoint", "Heavy", "Ultraheavy")),
         across(all_of(cols_factor), as.factor))

# check if for chosen second group level combinations of groups are more than for first group which is plotted
identical(length(unique(ta4n$Genus_abt)), 
          nrow(ta4n %>% group_by(Class_abt, Genus_abt) %>% select(Class_abt, Genus_abt) %>% distinct()))

## if false, check where additional combinations come from
taxa_comb <- ta4n %>% group_by(Class_abt, Genus_abt) %>% select(Class_abt, Genus_abt) %>% distinct()

## make modified Class_abt renaming some of the "other" entries in Genus_abt which cause the problem
ta4n_mod <- ta4n %>% 
  mutate(Class_abt_mod = if_else(Genus_abt == "other_Bacteria_<4%", "other_Bacteria_<4%",
                                  if_else(Genus_abt == "other_p_Desulfobacterota_<4%", "Desulfuromonadia",
                                          if_else(Genus_abt == "other_p_Planctomycetota_<4%", "Planctomycetes", Class_abt))))

## check again
identical(length(unique(ta4n_mod$Genus_abt)), 
          nrow(ta4n_mod %>% group_by(Class_abt_mod, Genus_abt) %>% select(Class_abt_mod, Genus_abt) %>% distinct()))
#is now fine

## need to name and order new column
class_mod_unique <- c(unique(ta4n_mod$Class_abt_mod))
class_mod_ordererd <- moveMe(class_mod_unique, "other_Bacteria_<4% last; Desulfuromonadia after Desulfobulbia")
ta4n_mod_f <- ta4n_mod %>% 
  mutate(Class_abt_mod = factor(Class_abt_mod, levels = class_mod_ordererd,
                                labels = class_mod_ordererd))

# sum ASV abundance by taxon ####
## need to decide which rank to plot and to group by which other rank
 # here plot genus, group by class
ta4s <- ta4n_mod_f %>% 
  group_by(Sample, Class_abt_mod, Genus_abt, Treatment, Isotope, Fraction, Temp) %>% 
  summarise(Abundance = sum(Abundance)) %>% 
  ungroup()

## check for max. abundance of plotted taxon, relevant for "other" taxa
plot_max <- ta4s %>% 
  group_by(Genus_abt) %>% 
  summarise(max = max(Abundance)) %>% 
  mutate(above_threshold = if_else(max < 0.04, FALSE, TRUE))


save(ta4n_mod_f, file = paste0(prj, ".rel.abd.ASV.long_ab4perc.RData"))
save(ta4s, file = paste0(prj, ".rel.abd.taxa.long_ab4perc.RData"))

write_tsv(ta4n_mod_f, "MS/raw_data_scripts/rawdata_ASV_plots.txt")

# plot taxa ####
## get automatic colors ####
# function to define multiple colour pallete based on subgroups from here
# https://stackoverflow.com/questions/49818271/stacked-barplot-with-colour-gradients-for-each-bar?rq=1
ColourPalleteMulti <- function(df, group, subgroup){
  
  # Find how many colour categories to create and the number of colours in each
  categories <- aggregate(as.formula(paste(subgroup, group, sep="~" )), df, function(x) length(unique(x)))
  category.start <- (scales::hue_pal(l = 40, direction = 1)(nrow(categories))) # Set the top of the colour pallete
  category.end  <- (scales::hue_pal(l = 100, direction = 1)(nrow(categories))) # set the bottom
  
  # Build Colour pallette
  colours <- unlist(lapply(1:nrow(categories),
                           function(i){
                             colorRampPalette(colors = c(category.start[i], category.end[i]))(categories[i,2])}))
  return(colours)
}

ColourPalleteMulti_cat <- function(df, group, subgroup){
  
  # Find how many colour categories to create and the number of colours in each
  categories <- aggregate(as.formula(paste(subgroup, group, sep="~" )), df, function(x) length(unique(x)))
  category.start <- (scales::hue_pal(l = 40, direction = 1)(nrow(categories))) # Set the top of the colour pallete
  category.end  <- (scales::hue_pal(l = 100, direction = 1)(nrow(categories))) # set the bottom
  
  # Build Colour pallette
  colours <- unlist(lapply(1:nrow(categories),
                           function(i){
                             colorRampPalette(colors = c(category.start[i], category.end[i]))(categories[i,2])}))
  return(categories)
}

colours_cat <- ColourPalleteMulti_cat(ta4s, "Class_abt_mod", "Genus_abt")


## create colour gradient with color pallet function
colours <- ColourPalleteMulti(ta4s, "Class_abt_mod", "Genus_abt")
### export generated color table
write(colours, paste0(prj, "_color-pallet-automatical.txt"))
### import hand edited color table
colours_edit <- scan(paste0(prj, "_color-pallet-edit.txt"), what = "")


## plot basic taxa ####
write_tsv(ta4s, "MS/raw_data_scripts/rawdata_barplots.txt")
load(paste0(prj, ".rel.abd.taxa.long_ab4perc.RData"))

plot_basic <- ggplot(ta4s, aes(x = Fraction, y = Abundance*100, fill = Genus_abt
                               , color = Genus_abt
                               )) + 
  geom_bar(position = position_stack(), stat = "identity", width = 0.8, linewidth = 0.01) +  # Stacked 100% barplot
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual("", values = colours_edit) +
  scale_color_manual("", values = colours_edit) +
  guides(fill = guide_legend(ncol = 3), color = guide_legend(ncol = 8)) +
  labs(y = "Relative abundance (%)") +
  facet_wrap(~Treatment*Isotope, nrow = 1) +
  theme_cowplot() +
  theme(text = element_text(size = 7), axis.text = element_text(size = 7),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
        axis.line = element_line(linewidth = 0), 
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1), 
        strip.background = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.text = element_text(size = 7), legend.position = "bottom",
        plot.title = element_text(size = 8, hjust = 0.5),
        legend.key.size = unit(0.3, "cm"), legend.key.spacing.y = unit(-0.02, "cm"))


## subset basic plot for different temperatures ####
plot_2deg <- plot_basic %+% subset(ta4s, ta4s$Temp == "2") +
  labs(title = "Incubation at 2 °C")

plot_5deg <- plot_basic %+% subset(ta4s, ta4s$Temp == "5") +
  labs(title = "Incubation at 5 °C")

plot_10deg <- plot_basic %+% subset(ta4s, ta4s$Temp == "10") +
  labs(title = "Incubation at 10 °C")

plot_15deg <- plot_basic %+% subset(ta4s, ta4s$Temp == "15") +
  labs(title = "Incubation at 15 °C")

plot_20deg <- plot_basic %+% subset(ta4s, ta4s$Temp == "20") +
  labs(title = "Incubation at 20 °C")

plot_25deg <- plot_basic %+% subset(ta4s, ta4s$Temp == "25") +
  labs(title = "Incubation at 25 °C")

plot_30deg <- plot_basic %+% subset(ta4s, ta4s$Temp == "30") +
  labs(title = "Incubation at 30 °C")


## combine plots for supplementary ####
comb1 <- plot_2deg + plot_5deg + guide_area() +
  plot_annotation(tag_levels = 'a') + 
  plot_layout(ncol = 1, guides = 'collect', heights = c(1,1,1.2)) &
  theme(legend.key.size = unit(0.3, "cm"), legend.box.margin = margin(l = -30, t = 5, b = 5, r = 2))

ggsave("MS/plots/supl/supl_FigS05_temp_SIP_barchart_comb1_2-5deg.png", plot = comb1,
       width = 16.5, height = 21, units = "cm")
ggsave("MS/plots/supl/supl_FigS05_temp_SIP_barchart_comb1_2-5deg.pdf", plot = comb1,
       width = 16.5, height = 21, units = "cm")
ggsave("MS/plots/supl/supl_FigS05_temp_SIP_barchart_comb1_2-5deg.emf", plot = comb1,
       width = 16.5, height = 21, units = "cm",
       device = {function(filename, ...) devEMF::emf(file = filename, ...)})

comb2 <- plot_10deg + plot_15deg + guide_area() +
  plot_annotation(tag_levels = 'a') + 
  plot_layout(ncol = 1, guides = 'collect', heights = c(1,1,1.2)) &
  theme(legend.key.size = unit(0.3, "cm"), legend.box.margin = margin(l = -30, t = 5, b = 5, r = 2))

ggsave("MS/plots/supl/supl_FigS06_temp_SIP_barchart_comb2_10-15deg.png", plot = comb2,
       width = 16.5, height = 21, units = "cm")
ggsave("MS/plots/supl/supl_FigS06_temp_SIP_barchart_comb2_10-15deg.pdf", plot = comb2,
       width = 16.5, height = 21, units = "cm")
ggsave("MS/plots/supl/supl_FigS06_temp_SIP_barchart_comb2_10-15deg.emf", plot = comb2,
       width = 16.5, height = 21, units = "cm",
       device = {function(filename, ...) devEMF::emf(file = filename, ...)})


comb3 <- plot_20deg + plot_25deg + guide_area() +
  plot_annotation(tag_levels = 'a') + 
  plot_layout(ncol = 1, guides = 'collect', heights = c(1,1,1.2)) &
  theme(legend.key.size = unit(0.3, "cm"), legend.box.margin = margin(l = -30, t = 5, b = 5, r = 2))

ggsave("MS/plots/supl/supl_FigS07_temp_SIP_barchart_comb3_20-25deg.png", plot = comb3,
       width = 16.5, height = 21, units = "cm")
ggsave("MS/plots/supl/supl_FigS07_temp_SIP_barchart_comb3_20-25deg.pdf", plot = comb3,
       width = 16.5, height = 21, units = "cm")
ggsave("MS/plots/supl/supl_FigS07_temp_SIP_barchart_comb3_20-25deg.emf", plot = comb3,
       width = 16.5, height = 21, units = "cm",
       device = {function(filename, ...) devEMF::emf(file = filename, ...)})


comb4 <- plot_30deg + guide_area() +
  plot_layout(ncol = 1, guides = 'collect', heights = c(1,1.2)) &
  theme(legend.key.size = unit(0.3, "cm"), legend.box.margin = margin(l = -30, t = 5, b = 5, r = 2))

ggsave("MS/plots/supl/supl_FigS08_temp_SIP_barchart_comb4_30deg.png", 
       plot = comb4, width = 17, height = 13, units = "cm")
ggsave("MS/plots/supl/supl_FigS08_temp_SIP_barchart_comb4_30deg.pdf", 
       plot = comb4, width = 17, height = 13, units = "cm")
ggsave("MS/plots/supl/supl_FigS08_temp_SIP_barchart_comb4_30deg.emf", 
       plot = comb4, width = 17, height = 13, units = "cm",
       device = {function(filename, ...) devEMF::emf(file = filename, ...)})

