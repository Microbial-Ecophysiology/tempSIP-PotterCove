# plot only taxa enriched in ultra-heavy and heavy 13C fractions
## need file created in script S_05_barplots.R

# load packages ####
library(tidyverse)
library(cowplot)
library(patchwork)
library(ggtext)

## import data ####
setwd("WORKINGDIRECTORY")

## project name to name files
prj <- "PROJECT"

## load data created in R script S_05_barplots.R
load(paste0(prj, ".rel.abd.ASV.long_ab4perc.RData"))
ta4n <- ta4n_mod_f


# calculations ####
## calculate mean in heavy and light fractions for 12C and 13C and differences ####
mean <- ta4n %>% 
  filter(Isotope == "13C" & Fraction %in% c("Heavy", "Ultraheavy")) %>% 
  group_by(OTU, Genus_abt, Class_abt_mod, Treatment, Temp, Genus_uncl, Family_uncl, Class_uncl, Order_uncl, Phylum_uncl) %>% 
  summarise(Abundance = mean(Abundance), .groups = "drop") %>%                    # calculate mean rel.abd per ASV of 13C heavy and ultra-heavy fractions
  left_join(
    (ta4n %>% 
       filter(Isotope == "13C" & Fraction %in% c("Light", "Ultralight")) %>% 
       group_by(OTU, Genus_abt, Class_abt_mod, Treatment, Temp, Genus_uncl, Family_uncl, Class_uncl, Order_uncl, Phylum_uncl) %>% 
       summarise(Abundance = mean(Abundance), .groups = "drop")),                    # calculate mean rel.abd per ASV of 13C light and ultra-light fractions
    by = c("OTU", "Genus_abt", "Class_abt_mod", "Treatment", "Temp", "Genus_uncl", "Family_uncl", "Class_uncl", "Order_uncl", "Phylum_uncl"),
    suffix = c(".13C.heavy", ".13C.light")
  ) %>% 
  left_join(
    left_join(
      (ta4n %>% 
         filter(Isotope == "12C" & Fraction %in% c("Heavy", "Ultraheavy")) %>% 
         group_by(OTU, Genus_abt, Class_abt_mod, Treatment, Temp, Genus_uncl, Family_uncl, Class_uncl, Order_uncl, Phylum_uncl) %>% 
         summarise(Abundance = mean(Abundance), .groups = "drop")),                  # calculate mean rel.abd per ASV of 13C light and ultra-light fractions
      (ta4n %>% 
         filter(Isotope == "12C" & Fraction %in% c("Light", "Ultralight")) %>% 
         group_by(OTU, Genus_abt, Class_abt_mod, Treatment, Temp, Genus_uncl, Family_uncl, Class_uncl, Order_uncl, Phylum_uncl) %>% 
         summarise(Abundance = mean(Abundance), .groups = "drop")),                   # calculate mean rel.abd per ASV of 13C light and ultra-light fractions
      by = c("OTU", "Genus_abt", "Class_abt_mod", "Treatment", "Temp", "Genus_uncl", "Family_uncl", "Class_uncl", "Order_uncl", "Phylum_uncl"),
      suffix = c(".12C.heavy", ".12C.light")
    ),
    by = c("OTU", "Genus_abt", "Class_abt_mod", "Treatment", "Temp", "Genus_uncl", "Family_uncl", "Class_uncl", "Order_uncl", "Phylum_uncl")
  ) %>% 
  left_join(
    ta4n %>% 
      # calculate max abundance per sample
      group_by(OTU, Treatment, Temp) %>% 
      summarise(max.abd = max(Abundance), .groups = "drop"),
    by = c("OTU", "Treatment", "Temp")
  ) %>% 
  mutate(
    # calculate differences
    dif.13C = Abundance.13C.heavy - Abundance.13C.light,
    dif.12C = Abundance.12C.heavy - Abundance.12C.light,
    dif.13C.12C = dif.13C - dif.12C)

## determine thresholds to use ####
# function rounding number up to next power of 10
roundUp <- function(x) 10^ceiling(log10(x))*10

# find for differences > 0 median and round to next power of 10
mean %>% filter(dif.13C > 0) %>% pull(dif.13C) %>% summary()
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0000000 0.0000041 0.0000087 0.0004474 0.0000208 0.4895800 
thr.dif.13C <- mean %>% filter(dif.13C > 0) %>% pull(dif.13C) %>% median() %>% roundUp()
# 1e-04

mean %>% filter(dif.13C.12C > 0) %>% pull(dif.13C.12C) %>% summary()
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
# 0.0000000 0.0000045 0.0000101 0.0004645 0.0000283 0.6541712
thr.dif.13C.12C <-  mean %>% filter(dif.13C.12C > 0) %>% pull(dif.13C.12C) %>% median() %>% roundUp()
# 0.001


## calculate labeled ASVs ####
# abundance threshold (not %) per sample below which labeling is considered false as not calculable reliably
min.abd.thr <- 0.005

ASVs_labeled_calc <- mean %>% 
  mutate(
    # calc orginal calc with threshold
    labeled_origcal_thr = if_else(max.abd > min.abd.thr & dif.13C > thr.dif.13C, T, F),
    # new calc considering 12C difference, implemented in revision process
    labeled_newcalc_thr = if_else(max.abd > min.abd.thr & dif.13C > 0 & dif.13C.12C > thr.dif.13C.12C, T, F),
    # modify Abundance column by applying new label calculation
    Abundance.13C.heavy_mod = if_else(labeled_newcalc_thr == T, Abundance.13C.heavy, 0)
  )

# look at how many ASVs in samples are labeled
sum_perc = list(same_labeling = ~sum(.x),
                perc_same_labeling = ~mean(.x))
summarise_labeled_methods <- ASVs_labeled_calc %>% 
  summarise(across(where(is.logical), sum_perc)) %>% 
  t() %>% as.data.frame() %>% rownames_to_column() %>% as_tibble()
print(summarise_labeled_methods)


## filter ASVs above threshold ####
list_asv_above_threshold_variable_calc <- function(input_table, rank, threshold, abundance_column, calc_col) {
  df1 <- input_table %>% 
    # need syntax .data[[]] here to work in function below where column name comes from vector as string
    filter(.data[[calc_col]]) %>% 
    group_by(across({{rank}})) %>% 
    summarise(max_abundance = max({{abundance_column}})) %>% 
    filter(max_abundance > threshold/100) %>% 
    arrange(desc(max_abundance))
  out_vector <- pull(df1, {{rank}})
  return(out_vector)
}

abundance_threshold <- 2
calc_used <- "labeled_newcalc_thr"
ASVs_ab_threshold <- list_asv_above_threshold_variable_calc(ASVs_labeled_calc, OTU, abundance_threshold, 
                                                            Abundance.13C.heavy_mod, calc_used)
# 31 ASVs
## export list of labeled ASVs above threshold
write(ASVs_ab_threshold, paste0("raw_data_scripts/", prj, "-list-labeled-ASVs-ab-", abundance_threshold, "perc.txt"))

## combine all calculated values in one table ####
df.comb <- ta4n %>% 
  # join ASV table with calculated mean values and label calculation results
  left_join(ASVs_labeled_calc, 
            by = c("OTU", "Genus_abt", "Class_abt_mod", "Treatment", "Temp", "Genus_uncl", "Family_uncl", "Class_uncl", "Order_uncl", "Phylum_uncl")) %>% 
  mutate(
    # include used threshold values in table
    used_threshold_origcalc = thr.dif.13C,
    used_threshold_newcalc = thr.dif.13C.12C,
    # calculation used
    used_calc = calc_used,
    # abundance threshold used
    used_abd.threshold = abundance_threshold,
    # ASV above used threshold
    ASV_is_ab_threshold = if_else(OTU %in% ASVs_ab_threshold, T, F),
    
    # include information if ASV was used in main plot or why not
    ASV_included_in_plot = if_else(OTU %in% c("sq1", "sq7", "sq10", "sq93", "sq12", "sq44",
                                              "sq36","sq19","sq18","sq5",
                                              "sq16", "sq3","sq26","sq6", 
                                              "sq21","sq4", "sq2", "sq29"), "main plot",
                                   if_else(OTU %in% c("sq77", "sq50"), "no clear label, only 13C in one sample", 
                                           # additional applied threshold: 9% per taxon, 5% per ASV
                                           if_else(ASV_is_ab_threshold, "supplementary plot, below main plot threshold", NA)
                                   ))
  )
# save table 
write_tsv(df.comb, "raw_data_scripts/rawdata_ASV_table_label_results.txt")
save(df.comb, file = "raw_data_scripts/rawdata_ASV_table_label_results.RData")



# plotting ####
## MS Fig 03 ####
### sort and order ASVs
data_asv <- df.comb %>% 
  # filter for labeled ASVs decided for main plot
  filter(ASV_included_in_plot == "main plot") %>% 
  # only keep columns for plotting
  select(OTU, Treatment, Temp, Abundance.13C.heavy_mod, Genus_uncl) %>% 
  distinct() %>% 
  rename(Abundance = Abundance.13C.heavy_mod) %>% 
  mutate(Temp = as.numeric(as.character(Temp)),
         OTU = factor(OTU, levels = c("sq1", "sq7", "sq10", "sq93", "sq12", "sq44",
                                      "sq36","sq19","sq18","sq5",
                                      "sq16", "sq3","sq26","sq6", 
                                      "sq21","sq4", "sq2", "sq29"),
                      labels = c("ASV 1", "ASV 7", "ASV 10", "ASV 93", "ASV 12", "ASV 44",
                                 "ASV 36","ASV 19","ASV 18","ASV 5",
                                 "ASV 16", "ASV 3","ASV 26","ASV 6", 
                                 "ASV 21","ASV 4", "ASV 2", "ASV 29"))) %>% 
  arrange(desc(Treatment))
### save txt file needed to reproduce the plot
write_tsv(data_asv, "data_plots/data_Fig03_rebuttal.txt", quote = "needed")


### plot with all ASVs for main plot
treat_label <- c("Acetate", "Acetate + iron oxide", "Acetate + sulfate")
treat_col <- c("#A5DB36FF", "#F8850FFF", "#39568CFF")

plot.asv.all <- 
  ggplot(data_asv, aes(x = Temp, y = Abundance*100, color = Treatment,
                       shape = Treatment, size = Treatment, linetype = Treatment, group = Treatment)) +
  geom_point() +
  scale_color_manual(values = treat_col, labels = treat_label) +
  scale_shape_discrete(labels = treat_label) +
  scale_size_manual(values = c(2,2.9,3.3), guide = "none") +
  scale_x_continuous(limits = c(0,31), expand = c(0,0)) +
  labs(x = "Temperature (°C)", 
       y = "Relative abundance in heavy fractions (%)") +
  guides(color = guide_legend("Substrate"), shape = guide_legend("Substrate"), 
         linetype = guide_legend("Substrate")) +
  facet_wrap(~OTU, nrow = 1) +
  theme_bw() +
  theme(panel.grid = element_blank(), panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        text = element_text(size = 7.8),
        legend.text = element_text(size = 7.8), legend.title = element_text(face = "bold"),
        axis.text = element_text(size = 6, colour = "black"),
        strip.background = element_rect(fill = NA, color = NA), strip.text = element_text(size = 7),
        plot.title = element_markdown(size = 7.8, hjust = 0.5, margin = margin(b = -2)),
        axis.ticks = element_line(linewidth = 0.5), axis.title = element_blank(),
        plot.margin = margin(b = 2, t = 3, r = 5.5, l = 5.5))

### subset plot for taxa
m.Sva <- plot.asv.all + subset(data_asv, data_asv$Genus_uncl == "Sva1033_g_unclassified") +
  labs(title = "**Sva1033**", tag = "a") +
  scale_y_continuous(limits = c(0,40), breaks = c(0,20,40))

m.Desumu <- plot.asv.all + subset(data_asv, data_asv$Genus_uncl == "Desulfuromusa") +
  labs(title = "***Desulfuromusa***", tag = "b") +
  scale_y_continuous(limits = c(0,56), breaks = c(0,25,50))

m.uDesulfu.16 <- plot.asv.all + subset(data_asv, data_asv$OTU == "ASV 16") +
  labs(title = "**unclassified**<br>***Desulfuromonadales***", tag = "c") +
  scale_y_continuous(limits = c(0,23), breaks = c(0,10,20))

m.Desumon <- plot.asv.all + subset(data_asv, data_asv$Genus_uncl == "Desulfuromonas") +
  labs(title = "***Desulfuromonas***", tag = "d") +
  scale_y_continuous(limits = c(0,23), breaks = c(0,10,20))

m.Sulfuri <- plot.asv.all + subset(data_asv, data_asv$Genus_uncl == "Sulfurimonas") +
  labs(title = "***Sulfurimonas***", tag = "e") +
  scale_y_continuous(limits = c(0,40), breaks = c(0,20,40))

m.Sulfuro <- plot.asv.all + subset(data_asv, data_asv$Genus_uncl == "Sulfurospirillum") +
  labs(title = "***Sulfurospirillum***", tag = "f") +
  scale_y_continuous(limits = c(0,40), breaks = c(0,20,40))


### combine plots
ms_design <- "
NAAAAAAAAAAAAAAAAAAAAAAA
NAAAAAAAAAAAAAAAAAAAAAAA
NAAAAAAAAAAAAAAAAAAAAAAA
NBBBBBBBBBBBBBBBBBBBBBBB
NBBBBBBBBBBBBBBBBBBBBBBB
NBBBBBBBBBBBBBBBBBBBBBBB
NDDDDDDEEEEEEEEEEEEEEEEE
NDDDDDDEEEEEEEEEEEEEEEEE
NDDDDDDEEEEEEEEEEEEEEEEE
NGGGGGGGGGGGHHHHHHHHHHHH
NGGGGGGGGGGGHHHHHHHHHHHH
NGGGGGGGGGGGHHHHHHHHHHHH
NOOOOOOOOOOOOOOOOOOOOOOO
"

ylab.ms <- ggplot(data.frame(l = plot.asv.all$labels$y, x = 1, y = 1)) +
  geom_text(aes(x, y, label = l), angle = 90, size = (7.8 / (14/5))) + 
  theme_void() +
  coord_cartesian(clip = "off")

xlab.ms <- ggplot(data.frame(l = plot.asv.all$labels$x, x = 1, y = 1)) +
  geom_text(aes(x, y, label = l), size = (7.8 / (14/5))) + 
  theme_void() +
  coord_cartesian(clip = "off")

MS.plot.asv <- m.Sva + m.Desumu + m.uDesulfu.16 + m.Desumon +
  m.Sulfuri + m.Sulfuro + 
  ylab.ms + xlab.ms +
  plot_layout(guides = "collect", design = ms_design) &
  theme(legend.position = "bottom", legend.margin = margin(t = -10))

ggsave("plots/MS_Fig03_temp_SIP_enr.species.png", plot = MS.plot.asv,
       width = 18, height = 17, units = "cm")
ggsave("plots/MS_Fig03_temp_SIP_enr.species.pdf", plot = MS.plot.asv,
       width = 18, height = 17, units = "cm")
ggsave("plots/MS_Fig03_temp_SIP_enr.species.emf", plot = MS.plot.asv,
       width = 18, height = 17, units = "cm",
       device = {function(filename, ...) devEMF::emf(file = filename, ...)})

## MS supl Figs labeled ASVs across fractions ####
### calculations ####
# reduce table to only entries of 31 labeled ASVs above threshold
df.sub <- df.comb %>% filter(ASV_is_ab_threshold) %>% 
  mutate(Fraction = factor(Fraction, levels = c("Ultralight", "Light", "Midpoint", "Heavy", "Ultraheavy"),
                           labels = c("UL", "L", "M", "H", "UH")),
         Treatment = factor(Treatment, labels = c("Acetate", "Acet. + iron ox.", "Acet. + sulfate")),
         Temp = factor(Temp, levels = c("2", "5", "10", "15", "20", "25", "30"), 
                       labels = c("2 °C", "5 °C", "10 °C", "15 °C", "20 °C", "25 °C", "30 °C")))

# df.sub in long format for calculation test results
df.sub.long <- df.sub %>% 
  select(OTU, Treatment, Temp, labeled_newcalc_thr) %>% 
  distinct() %>% 
  mutate(Labeled = factor(labeled_newcalc_thr, levels = c("TRUE", "FALSE"), labels = c("labeled", "not labeled")),
         x.calc = "UL",
         y.calc = Inf)

# plot relative abundance of all samples and ASVs
## variables to define later in function
text.size.axis = 6
text.size.rest = 7
text.size.facet = 7
point.size = 2
bar.line.size = 0.3
bar.size = 0.9


# function for abundance over fraction in treatment ~ temp table with labeling significant result
plot_asv_over_fraction <- function(OTU.filt, OTU.title) {
  plot_asv_frac <- ggplot() +
    geom_bar(data = (df.sub %>% filter(OTU == OTU.filt)),
             aes(x = Fraction, y = Abundance*100, fill = Isotope), stat = "identity", position = position_dodge(), color = "black", 
             linewidth = bar.line.size, width = bar.size) +
    scale_fill_manual(values = c("gray85", "black")) +
    geom_point(data = (df.sub.long %>% filter(OTU == OTU.filt)),
               aes(x = x.calc, 
                   y = max((df.sub %>% filter(OTU == OTU.filt) %>% pull(Abundance)))*98, 
                   shape = Labeled),
               size = point.size, fill = "white", color = "#E76BF3") +
    scale_shape_manual(values = c(19, 21), guide = guide_legend(title = "Calculation", nrow = 1), drop = F) +
    scale_y_continuous(expand = c(0,0), 
                       limits = c(0,
                                  max((df.sub %>% filter(OTU == OTU.filt) %>% pull(Abundance)))*110)) +
    labs(x = "Fractions", y = "Relative abundance (%)", title = OTU.title) +
    facet_grid(Treatment ~ Temp) +
    theme_bw() +
    theme(text = element_text(size = text.size.rest), axis.title = element_text(size = text.size.rest),
          axis.text = element_text(size = text.size.axis, colour = "black"),
          plot.title = element_text(size = text.size.rest, face = "bold", hjust = 0.5),
          legend.position = "bottom",
          legend.text = element_text(size = text.size.rest), 
          legend.title = element_text(size = text.size.rest, face = "bold"),
          legend.key.height = unit(3, "pt"), legend.key.width = unit(6, "pt"),
          strip.background = element_rect(fill = "white"),
          strip.text.x = element_text(size = text.size.facet),
          strip.text.y = element_text(size = text.size.axis))
  return(plot_asv_frac)
}

# unique ASVs
ASVs_ab_threshold <- df.sub %>% pull(OTU) %>% unique() %>% sub(pattern = "sq", replacement = "") %>% 
  as.numeric() %>% sort() %>% paste0("sq", .)
# unique ASVs with taxonomy
ASVs_ab_threshold_tax <- df.comb %>% filter(ASV_is_ab_threshold) %>% 
  mutate(ASV_tax = paste0(OTU, "_", Genus_uncl)) %>% pull(ASV_tax) %>% unique()

## plot and combine plots for all labeled ASVs
### plots of ASVs in main plot Fig03 ####
#  plot 1 Sva1033
plot_asv_Sva_1 <- plot_asv_over_fraction("sq1","ASV 1 Sva1033") +
  plot_asv_over_fraction("sq7","ASV 7 Sva1033") +
  plot_asv_over_fraction("sq10", "ASV 10 Sva1033") +
  plot_annotation(tag_levels = 'a') +
  plot_layout(guides = "collect", ncol = 1) &
  theme(legend.position = "bottom")
ggsave("plots/supl_FigSX_temp_SIP_lab.asv.over.fract.01.Sva.1.png", plot_asv_Sva_1, 
       height = 22, width = 16.5, units = "cm")
ggsave("plots/supl_FigSX_temp_SIP_lab.asv.over.fract.01.Sva.1.pdf", plot_asv_Sva_1, 
       height = 22, width = 16.5, units = "cm")
ggsave("plots/supl_FigSX_temp_SIP_lab.asv.over.fract.01.Sva.1.emf", plot_asv_Sva_1, 
       height = 22, width = 16.5, units = "cm",
       device = {function(filename, ...) devEMF::emf(file = filename, ...)})

#  plot 2 Sva1033
plot_asv_Sva_2 <- plot_asv_over_fraction("sq93", "ASV 93 Sva1033") +
  plot_asv_over_fraction("sq12","ASV 12 Sva1033") +
  plot_asv_over_fraction("sq44", "ASV 44 Sva1033") +
  plot_annotation(tag_levels = 'a') +
  plot_layout(guides = "collect", ncol = 1) &
  theme(legend.position = "bottom")
ggsave("plots/supl_FigSX_temp_SIP_lab.asv.over.fract.02.Sva.2.png", plot_asv_Sva_2, 
       height = 22, width = 16.5, units = "cm")
ggsave("plots/supl_FigSX_temp_SIP_lab.asv.over.fract.02.Sva.2.pdf", plot_asv_Sva_2, 
       height = 22, width = 16.5, units = "cm")
ggsave("plots/supl_FigSX_temp_SIP_lab.asv.over.fract.02.Sva.2.emf", plot_asv_Sva_2, 
       height = 22, width = 16.5, units = "cm",
       device = {function(filename, ...) devEMF::emf(file = filename, ...)})

# plot 3 Desulfuromusa
plot_asv_Desulfuromu_1 <- plot_asv_over_fraction("sq36", "ASV 36 Desulfuromusa") +
  plot_asv_over_fraction("sq19", "ASV 19 Desulfuromusa") +
  plot_asv_over_fraction("sq18","ASV 18 Desulfuromusa") +
  plot_asv_over_fraction("sq5","ASV 5 Desulfuromusa") +
  plot_annotation(tag_levels = 'a') +
  plot_layout(guides = "collect", ncol = 1) &
  theme(legend.position = "bottom", 
        # for having 4 subplots
        strip.text.y = element_text(size = 4), plot.margin = margin(t=3, r=5.5, l=5.5, b = 0))
ggsave("plots/supl_FigSX_temp_SIP_lab.asv.over.fract.03.Desulfuromusa.png", plot_asv_Desulfuromu_1, 
       height = 22, width = 16.5, units = "cm")
ggsave("plots/supl_FigSX_temp_SIP_lab.asv.over.fract.03.Desulfuromusa.pdf", plot_asv_Desulfuromu_1, 
       height = 22, width = 16.5, units = "cm")
ggsave("plots/supl_FigSX_temp_SIP_lab.asv.over.fract.03.Desulfuromusa.emf", plot_asv_Desulfuromu_1, 
       height = 22, width = 16.5, units = "cm",
       device = {function(filename, ...) devEMF::emf(file = filename, ...)})

# plot 4 uncl. Desulfuromonadales & Desulfuromonas
plot_asv_Desulfuromo <- plot_asv_over_fraction("sq16","ASV 16 uncl. Desulfuromonadales") +
  plot_asv_over_fraction("sq3","ASV 3 Desulfuromonas") +
  plot_asv_over_fraction("sq26", "ASV 26 Desulfuromonas") +
  plot_asv_over_fraction("sq6","ASV 6 Desulfuromonas") +
  plot_annotation(tag_levels = 'a') +
  plot_layout(guides = "collect", ncol = 1) &
  theme(legend.position = "bottom", 
        # for having 4 subplots
        strip.text.y = element_text(size = 4), plot.margin = margin(t=3, r=5.5, l=5.5, b = 0))
ggsave("plots/supl_FigSX_temp_SIP_lab.asv.over.fract.04.uncl.Desulfo.Desulfuromonas.png", plot_asv_Desulfuromo, 
       height = 22, width = 16.5, units = "cm")
ggsave("plots/supl_FigSX_temp_SIP_lab.asv.over.fract.04.uncl.Desulfo.Desulfuromonas.pdf", plot_asv_Desulfuromo, 
       height = 22, width = 16.5, units = "cm")
ggsave("plots/supl_FigSX_temp_SIP_lab.asv.over.fract.04.uncl.Desulfo.Desulfuromonas.emf", plot_asv_Desulfuromo, 
       height = 22, width = 16.5, units = "cm",
       device = {function(filename, ...) devEMF::emf(file = filename, ...)})

# plot 5 Sulfurimonas & Sulfurospirillum
plot_asv_Sulfuri.Sulfuro <- plot_asv_over_fraction("sq21","ASV 21 Sulfurimonas" ) +
  plot_asv_over_fraction("sq4","ASV 4 Sulfurimonas" ) +
  plot_asv_over_fraction("sq2", "ASV 2 Sulfurospirillum") +
  plot_asv_over_fraction("sq29", "sq29 Sulfurospirillum") +
  plot_annotation(tag_levels = 'a') +
  plot_layout(guides = "collect", ncol = 1) &
  theme(legend.position = "bottom",
        # for having 4 subplots
        strip.text.y = element_text(size = 4), plot.margin = margin(t=3, r=5.5, l=5.5, b = 0))
ggsave("plots/supl_FigSX_temp_SIP_lab.asv.over.fract.05.Sulfurimonas.Sulfurospir.png", plot_asv_Sulfuri.Sulfuro, 
       height = 22, width = 16.5, units = "cm")
ggsave("plots/supl_FigSX_temp_SIP_lab.asv.over.fract.05.Sulfurimonas.Sulfurospir.pdf", plot_asv_Sulfuri.Sulfuro, 
       height = 22, width = 16.5, units = "cm")
ggsave("plots/supl_FigSX_temp_SIP_lab.asv.over.fract.05.Sulfurimonas.Sulfurospir.emf", plot_asv_Sulfuri.Sulfuro, 
       height = 22, width = 16.5, units = "cm",
       device = {function(filename, ...) devEMF::emf(file = filename, ...)})

### plots of ASVs only in supplementary ####
# plot 6 other Sva1033 % uncl. Desulfuromonadales
plot_asv_Sva_3 <- plot_asv_over_fraction("sq43","ASV 43 Sva1033") +
  plot_asv_over_fraction("sq81","ASV 81 Sva1033") +
  plot_asv_over_fraction("sq57", "ASV 57 uncl. Desulfuromonadales") +
  plot_annotation(tag_levels = 'a') +
  plot_layout(guides = "collect", ncol = 1) &
  theme(legend.position = "bottom")
ggsave("plots/supl_FigSX_temp_SIP_lab.asv.over.fract.06.Sva.3.uncl.Desulf.png", plot_asv_Sva_3, 
       height = 22, width = 16.5, units = "cm")
ggsave("plots/supl_FigSX_temp_SIP_lab.asv.over.fract.06.Sva.3.uncl.Desulf.pdf", plot_asv_Sva_3, 
       height = 22, width = 16.5, units = "cm")
ggsave("plots/supl_FigSX_temp_SIP_lab.asv.over.fract.06.Sva.3.uncl.Desulf.emf", plot_asv_Sva_3, 
       height = 22, width = 16.5, units = "cm",
       device = {function(filename, ...) devEMF::emf(file = filename, ...)})

# plot 7 Trichloromoas & uncl. Geopsycho & Arcobacter
plot_asv_7 <- plot_asv_over_fraction("sq35","ASV 35 Trichloromonas") +
  plot_asv_over_fraction("sq56", "ASV 56 uncl. Geopsychrobacteraceae") +
  plot_asv_over_fraction("sq95","ASV 95 uncl. Arcobacteraceae") +
  plot_asv_over_fraction("sq134","ASV 134 uncl. Arcobacteraceae") +
  plot_annotation(tag_levels = 'a') +
  plot_layout(guides = "collect", ncol = 1) &
  theme(legend.position = "bottom",
        # for having 4 subplots
        strip.text.y = element_text(size = 4), plot.margin = margin(t=3, r=5.5, l=5.5, b = 0))
ggsave("plots/supl_FigSX_temp_SIP_lab.asv.over.fract.07.Trichl.uncl.Geops.Arco.png", plot_asv_7, 
       height = 22, width = 16.5, units = "cm")
ggsave("plots/supl_FigSX_temp_SIP_lab.asv.over.fract.07.Trichl.uncl.Geops.Arco.pdf", plot_asv_7, 
       height = 22, width = 16.5, units = "cm")
ggsave("plots/supl_FigSX_temp_SIP_lab.asv.over.fract.07.Trichl.uncl.Geops.Arco.emf", plot_asv_7, 
       height = 22, width = 16.5, units = "cm",
       device = {function(filename, ...) devEMF::emf(file = filename, ...)})

# plot 8 Fusibac & Maribac & Hoeflea & Colwellia
plot_asv_Fusi.Mari.Hof.Col <- plot_asv_over_fraction("sq17", "ASV 17 Fusibacter") +
  plot_asv_over_fraction("sq23","ASV 23 Maribacter") +
  plot_asv_over_fraction("sq40","ASV 40 Hoeflea") +
  plot_asv_over_fraction("sq105", "ASV 105 Colwellia") +
  plot_annotation(tag_levels = 'a') +
  plot_layout(guides = "collect", ncol = 1) &
  theme(legend.position = "bottom",
        # for having 4 subplots
        strip.text.y = element_text(size = 4), plot.margin = margin(t=3, r=5.5, l=5.5, b = 0))
ggsave("plots/supl_FigSX_temp_SIP_lab.asv.over.fract.08.Fusibac.Maribac.Hoflea.Colwel.png", plot_asv_Fusi.Mari.Hof.Col, 
       height = 22, width = 16.5, units = "cm")
ggsave("plots/supl_FigSX_temp_SIP_lab.asv.over.fract.08.Fusibac.Maribac.Hoflea.Colwel.pdf", plot_asv_Fusi.Mari.Hof.Col, 
       height = 22, width = 16.5, units = "cm")
ggsave("plots/supl_FigSX_temp_SIP_lab.asv.over.fract.08.Fusibac.Maribac.Hoflea.Colwel.emf", plot_asv_Fusi.Mari.Hof.Col, 
       height = 22, width = 16.5, units = "cm",
       device = {function(filename, ...) devEMF::emf(file = filename, ...)})

### plots of questionably labeled ASVs, only in supplementary ####
# plot 9 Marinobacer, Desulfobacter
plot_asv_Marino.Desulfob <- plot_asv_over_fraction("sq50","ASV 50 Marinobacter") +
  plot_asv_over_fraction("sq77","ASV 77 Desulfobacter") +
  plot_annotation(tag_levels = 'a') +
  plot_layout(guides = "collect", ncol = 1) &
  theme(legend.position = "bottom")
ggsave("plots/supl_FigSX_temp_SIP_lab.asv.over.fract.09.Marino.Desulfobac.png", plot_asv_Marino.Desulfob, 
       height = 15, width = 16.5, units = "cm")
ggsave("plots/supl_FigSX_temp_SIP_lab.asv.over.fract.09.Marino.Desulfobac.pdf", plot_asv_Marino.Desulfob, 
       height = 15, width = 16.5, units = "cm")
ggsave("plots/supl_FigSX_temp_SIP_lab.asv.over.fract.09.Marino.Desulfobac.emf", plot_asv_Marino.Desulfob, 
       height = 15, width = 16.5, units = "cm",
       device = {function(filename, ...) devEMF::emf(file = filename, ...)})

