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
## calculate average of heavy and ultra-heavy fraction for each labeled treatment
# with control if this value is higher than mean abundance in light fractions, if not set to 0
mean <- ta4n %>% 
  filter(Isotope == "13C" & Fraction %in% c("Heavy", "Ultraheavy")) %>% 
  group_by(OTU, Genus_abt, Class_abt_mod, Treatment, Temp, Genus_uncl, Family_uncl, Class_uncl, Order_uncl, Phylum_uncl) %>% 
  summarise(Abundance = mean(Abundance), .groups = "drop") %>%                    # calculate mean rel.abd per ASV of 13C heavy and ultra-heavy fractions
  mutate(Fraction = "Heavy-fractions") %>% 
  left_join(
    (ta4n %>% 
       filter(Isotope == "13C" & Fraction %in% c("Light", "Ultralight")) %>% 
       group_by(OTU, Genus_abt, Class_abt_mod, Treatment, Temp, Genus_uncl, Family_uncl, Class_uncl, Order_uncl, Phylum_uncl) %>% 
       summarise(Abundance = mean(Abundance), .groups = "drop") %>%                    # calculate mean rel.abd per ASV of 13C light and ultra-light fractions
       mutate(Fraction = "Light-fractions")),
    by = c("OTU", "Genus_abt", "Class_abt_mod", "Treatment", "Temp", "Genus_uncl", "Family_uncl", "Class_uncl", "Order_uncl", "Phylum_uncl")
  ) %>% 
  mutate(Abundance = if_else(Abundance.x < Abundance.y, 0, Abundance.x))        # check if ASVs are more abundant in heavy than in light fractions, if not set to 0



## find higher abundance ASVs
list_asv_above_threshold <- function(input_table, rank, threshold, abundance_column) {
  df1 <- input_table %>% 
    group_by(across({{rank}})) %>% 
    summarise(max_abundance = max({{abundance_column}})) %>% 
    filter(max_abundance > threshold/100) %>% 
    arrange(desc(max_abundance))
  out_vector <- pull(df1, {{rank}})
  return(out_vector)
}

abundance_threshold <- 2
ASV_ab_threshold <- list_asv_above_threshold(mean, OTU, abundance_threshold, Abundance)
# 31 ASVs above threshold

## from plots below, know that Fusibacter, Maribacter and Marinobacter ASVs are not really labeled
# remove them from ASV list
ASVs_not_labeled <- c("sq17", "sq23", "sq50")
ASV_ab_threshold_filt <- ASV_ab_threshold[!grepl(paste0(ASVs_not_labeled, collapse = "|"), ASV_ab_threshold)]

mean_filt <- mean %>% 
  mutate(ASV_mod = if_else(OTU %in% ASV_ab_threshold_filt, OTU, paste0("other ASVs < ", abundance_threshold, "%"))) %>% 
  group_by(ASV_mod, Genus_abt, Treatment, Temp, OTU) %>% 
  summarise(Abundance = sum(Abundance), .groups = "drop") %>% 
  left_join(unique(ta4n[,c("OTU", "Genus_uncl", "Family_uncl", "Class_uncl", "Order_uncl", "Phylum_uncl")]), by = c("OTU" = "OTU"))

unique((
  mean_filt %>% 
    subset(ASV_mod != "other ASVs < 2%") %>% 
    mutate(ASV_taxa = paste(Class_uncl, Order_uncl, Family_uncl, Genus_uncl, ASV_mod)) %>% 
    arrange(ASV_taxa)
)$ASV_taxa)

mean_filt_asvs <- mean_filt %>% 
  filter(ASV_mod != "other ASVs < 2%") %>% 
  mutate(ASV.taxa = paste(OTU, "\n", Genus_uncl),
         Temp = as.numeric(levels(Temp))[Temp],
         Treatment = factor(Treatment, levels = c("Acetate", "Acetate.Lepidocrocite", "Acetate.Sulfate"),
                            labels = c("Acetate", "Acetate + lepidocrocite", "Acetate + sulfate")),
         OTU.ASV = gsub("sq", "ASV ", OTU)) %>% 
  select(ASV.taxa, Treatment, Temp, Abundance, Genus_uncl, OTU, OTU.ASV) %>% 
  arrange(Genus_uncl)

asv_order <- c(unique(mean_filt_asvs$ASV.taxa))
mean_filt_asvs_sort <- mean_filt_asvs %>% 
  mutate(ASV.taxa = factor(ASV.taxa, levels = asv_order))

write_tsv(mean_filt_asvs_sort, "MS/raw_data_scripts/rawdata_enr_ASVs_13C_h.fraction.txt")

## filter for labeled ASVs > 2%
enr.asv <- ta4n %>% 
  filter(OTU %in% ASV_ab_threshold_filt) %>% 
  mutate(Treatment = factor(Treatment, levels = c("Acetate", "Acetate.Lepidocrocite", "Acetate.Sulfate"),
                            labels = c("Acetate", "Acetate + lepidocrocite", "Acetate + sulfate")),
         OTU.ASV = gsub("sq", "ASV ", OTU),
         OTU = factor(OTU), Fraction = factor(Fraction), OTU.ASV = factor(OTU.ASV),
         ASV.taxa = factor(paste(OTU, "\n", Genus_uncl)),
         Temp = as.numeric(as.character(Temp))) %>% 
  arrange(Genus_uncl) %>% 
  ungroup()

enr.asv.13C <- enr.asv %>% 
  filter(Isotope == "13C")

write_tsv(enr.asv.13C, "MS/raw_data_scripts/rawdata_enr_ASVs_across_fractions.txt")

# plotting ####
treat_label <- c("Acetate", "Acetate + lepidocrocite", "Acetate + sulfate")
treat_col <- c("#A5DB36FF", "#F8850FFF", "#39568CFF")

## MS supl Fig S09 plot with all enriched ASVs across temperature ####
all_asv <- ggplot(mean_filt_asvs_sort, aes(x = Temp, y = Abundance*100, color = Treatment,
                                           shape = Treatment, linetype = Treatment)) +
  geom_point(size = 1.6) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = treat_col, labels = treat_label) +
  scale_x_continuous(limits = c(0,31), expand = c(0,0)) +
  labs(x = "Temperature (°C)", y = "Relative abundance in heavy fractions (%)", 
       title = "ASVs > 2% relative abundance") +
  guides(color = guide_legend("Treatment"), 
         shape = guide_legend("Treatment"), linetype = guide_legend("Treatment")) +
  facet_wrap(~ASV.taxa, scale = "free", ncol = 4) +
  theme_bw() +
  theme(panel.grid = element_blank(), text = element_text(size = 7), axis.title = element_text(size = 7),
        legend.text = element_text(size = 7), axis.text = element_text(size = 6, colour = "black"),
        legend.position = "bottom", legend.margin = margin(t = -10),
        plot.title = element_text(size = 7, face = "bold", hjust = 0.5),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 6))

ggsave("MS/plots/supl/supl_FigS09_temp_SIP_all.enr.asv.abv2perc_acr.temp_filt.pdf", plot = all_asv,
       width = 16.5, height = 22, units = "cm")
ggsave("MS/plots/supl/supl_FigS09_temp_SIP_all.enr.asv.abv2perc_acr.temp_filt.png", plot = all_asv,
       width = 16.5, height = 22, units = "cm")
ggsave("MS/plots/supl/supl_FigS09_temp_SIP_all.enr.asv.abv2perc_acr.temp_filt.emf", plot = all_asv,
       width = 16.5, height = 22, units = "cm",
       device = {function(filename, ...) devEMF::emf(file = filename, ...)})


## MS supl plots enriched ASVs across fractions ####
## double check for each ASV if it actually was enriched
plot.asv.frac <- ggplot(enr.asv.13C, aes(x = Fraction, y = Abundance*100, group = Treatment,
                                         color = Treatment, shape = Treatment, linetype = Treatment)) +
  geom_point(size = 1.6) +
  geom_line(linewidth = 0.6) +
  scale_color_manual(values = treat_col, labels = treat_label) +
  scale_shape_discrete(labels = treat_label) +
  labs(x = "Fraction", 
       y = "Relative abundance (%)") +
  guides(color = guide_legend("Treatment"), shape = guide_legend("Treatment"),
         linetype = guide_legend("Treatment")) +
  facet_grid(OTU.ASV~Temp, scale = "free") +
  theme_bw() +
  theme(panel.grid = element_blank(), text = element_text(size = 7), axis.title = element_text(size = 7),
        legend.text = element_text(size = 7), axis.text = element_text(size = 6, colour = "black"),
        legend.position = "bottom", strip.background = element_rect(fill = "white"),
        legend.margin = margin(t = -10),
        strip.text = element_text(size = 6),
        plot.title = element_text(size = 8, hjust = 0.5, face = "bold"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


### subset ASVs over fractions for taxa ####
frac.Sva <- plot.asv.frac %+% subset(enr.asv.13C, Genus_uncl == "Sva1033_g_unclassified") +
  labs(title = "Sva1033")
ggsave("MS/plots/supl/supl_FigS11_temp_SIP_enr.asv.over.frac.Sva.png", plot = frac.Sva,
       width = 16.5, height = 17, units = "cm")
ggsave("MS/plots/supl/supl_FigS11_temp_SIP_enr.asv.over.frac.Sva.pdf", plot = frac.Sva,
       width = 16.5, height = 17, units = "cm")
ggsave("MS/plots/supl/supl_FigS11_temp_SIP_enr.asv.over.frac.Sva.emf", plot = frac.Sva,
       width = 16.5, height = 17, units = "cm",
       device = {function(filename, ...) devEMF::emf(file = filename, ...)})

frac.Desulfuromonas <- plot.asv.frac %+% subset(enr.asv.13C, Genus_uncl == "Desulfuromonas") +
  labs(title = expression(bolditalic("Desulfuromonas")))

frac.Trich <- plot.asv.frac %+% subset(enr.asv.13C, Genus_uncl == "Trichloromonas") +
  labs(title = expression(bolditalic("Trichloromonas")))

frac.Desulfuromusa <- plot.asv.frac %+% subset(enr.asv.13C, Genus_uncl == "Desulfuromusa") +
  labs(title = expression(bolditalic("Desulfuromusa")))

frac.Geopsy <- plot.asv.frac %+% subset(enr.asv.13C, Genus_uncl == "Geopsychrobacteraceae_g_unclassified") +
  labs(title = expression(bold("unclassified")~bolditalic("Geopsychrobacteraceae")))

frac.uncl.Des <- plot.asv.frac %+% subset(enr.asv.13C, Genus_uncl == "Desulfuromonadales_f_unclassified") +
  labs(title = expression(bold("unclassified")~bolditalic("Desulfuromonadales")))

frac.Sulfuri <- plot.asv.frac %+% subset(enr.asv.13C, Genus_uncl == "Sulfurimonas") +
  labs(title = expression(bolditalic("Sulfurimonas")))

frac.Sulfuro <- plot.asv.frac %+% subset(enr.asv.13C, Genus_uncl == "Sulfurospirillum") +
  labs(title = expression(bolditalic("Sulfurospirillum")))

frac.Arco <- plot.asv.frac %+% subset(enr.asv.13C, Genus_uncl == "Arcobacteraceae_g_unclassified") +
  labs(title = expression(bold("unclassified")~bolditalic("Arcobacteraceae")))

frac.Col <- plot.asv.frac %+% subset(enr.asv.13C, Genus_uncl == "Colwellia") +
  labs(title = expression(bolditalic("Colwellia")))

frac.Hof <- plot.asv.frac %+% subset(enr.asv.13C, Genus_uncl == "Hoeflea") +
  labs(title = expression(bolditalic("Hoeflea")))

frac.Desulfob <- plot.asv.frac %+% subset(enr.asv.13C, Genus_uncl == "Desulfobacter") +
  labs(title = expression(bolditalic("Desulfobacter")))

frac.Fusi <- plot.asv.frac %+% subset(enr.asv.13C, Genus_uncl == "Fusibacter") +
  labs(title = expression(bolditalic("Fusibacter")))
# Fusibacter not really labeled!

frac.Mari <- plot.asv.frac %+% subset(enr.asv.13C, Genus_uncl == "Maribacter") +
  labs(title = expression(bolditalic("Maribacter")))
# Maribacter not really labeled!

frac.Marino <- plot.asv.frac %+% subset(enr.asv.13C, Genus_uncl == "Marinobacter") +
  labs(title = expression(bolditalic("Marinobacter")))
# Marinobacter not labeled

## kick Fusibacter, Maribacter and Marinobacter out, re-run above code

### combine plots ####
frac.comb1 <- frac.Desulfuromonas + frac.Desulfuromusa + plot_annotation(tag_levels = 'a') +
  plot_layout(guides = "collect", ncol = 1, heights = c(3,4)) & 
  theme(legend.position = "bottom", legend.margin = margin(t = -5, b = -5))
ggsave("MS/plots/supl/supl_FigS12_temp_SIP_enr.asv.over.frac.Desulfuromonas.Desulfuromusa.png", plot = frac.comb1,
       width = 16.5, height = 19.5, units = "cm")
ggsave("MS/plots/supl/supl_FigS12_temp_SIP_enr.asv.over.frac.Desulfuromonas.Desulfuromusa.pdf", plot = frac.comb1,
       width = 16.5, height = 19.5, units = "cm")
ggsave("MS/plots/supl/supl_FigS12_temp_SIP_enr.asv.over.frac.Desulfuromonas.Desulfuromusa.emf", plot = frac.comb1,
       width = 16.5, height = 19.5, units = "cm",
       device = {function(filename, ...) devEMF::emf(file = filename, ...)})

frac.comb2.2 <- frac.Sulfuri + frac.Sulfuro + plot_annotation(tag_levels = 'a') +
  plot_layout(guides = "collect", ncol = 1, heights = c(2,2)) & 
  theme(legend.position = "bottom", legend.margin = margin(t = -5, b = -5))
ggsave("MS/plots/supl/supl_FigS13_temp_SIP_enr.asv.over.frac.Sulfuri.Sulfuro.png", plot = frac.comb2.2,
       width = 16.5, height = 14.5, units = "cm")
ggsave("MS/plots/supl/supl_FigS13_temp_SIP_enr.asv.over.frac.Sulfuri.Sulfuro.pdf", plot = frac.comb2.2,
       width = 16.5, height = 14.5, units = "cm")
ggsave("MS/plots/supl/supl_FigS13_temp_SIP_enr.asv.over.frac.Sulfuri.Sulfuro.emf", plot = frac.comb2.2,
       width = 16.5, height = 14.5, units = "cm",
       device = {function(filename, ...) devEMF::emf(file = filename, ...)})


frac.comb3 <- frac.Col + frac.Desulfob + frac.Geopsy + frac.Trich + plot_annotation(tag_levels = 'a') +
  plot_layout(guides = "collect", ncol = 1) & 
  theme(legend.position = "bottom", legend.margin = margin(t = -5, b = -5),
        plot.margin = margin(b = 0, t = 5.5, r = 5.5, l = 5.5))
ggsave("MS/plots/supl/supl_FigS14_temp_SIP_enr.asv.over.frac.Col.Desulfob.Geopsy.Tric.png", plot = frac.comb3,
       width = 16.5, height = 19, units = "cm")
ggsave("MS/plots/supl/supl_FigS14_temp_SIP_enr.asv.over.frac.Col.Desulfob.Geopsy.Tric.pdf", plot = frac.comb3,
       width = 16.5, height = 19, units = "cm")
ggsave("MS/plots/supl/supl_FigS14_temp_SIP_enr.asv.over.frac.Col.Desulfob.Geopsy.Tric.emf", plot = frac.comb3,
       width = 16.5, height = 19, units = "cm",
       device = {function(filename, ...) devEMF::emf(file = filename, ...)})


frac.comb4 <- frac.Arco + frac.uncl.Des + frac.Hof + plot_annotation(tag_levels = 'a') +
  plot_layout(guides = "collect", ncol = 1, heights = c(2,2,1)) & 
  theme(legend.position = "bottom", legend.margin = margin(t = -5, b = -5),
        plot.margin = margin(b = 0, t = 5.5, r = 5.5, l = 5.5))
ggsave("MS/plots/supl/supl_FigS15_temp_SIP_enr.asv.over.frac.Arco.unclDes.Hof.png", plot = frac.comb4,
       width = 16.5, height = 18, units = "cm")
ggsave("MS/plots/supl/supl_FigS15_temp_SIP_enr.asv.over.frac.Arco.unclDes.Hof.pdf", plot = frac.comb4,
       width = 16.5, height = 18, units = "cm")
ggsave("MS/plots/supl/supl_FigS15_temp_SIP_enr.asv.over.frac.Arco.unclDes.Hof.emf", plot = frac.comb4,
       width = 16.5, height = 18, units = "cm",
       device = {function(filename, ...) devEMF::emf(file = filename, ...)})


## MS Fig03 ####
### sort and order ASVs
data_asv <- mean_filt_asvs_sort %>% 
  filter(
    # for main plot remove low abundant Sva1033 ASV
    !OTU %in% c("sq81")) %>% 
  mutate(OTU = factor(OTU, levels = c("sq1", "sq43", "sq7", "sq10", "sq93", "sq12", "sq44",
                                      "sq36","sq19","sq18","sq5", "sq56",
                                      "sq16", "sq3","sq26","sq6", "sq35",
                                      "sq21","sq4", "sq2", "sq29", "sq77",
                                      "sq95","sq134", "sq105", "sq57", "sq40"),
                      labels = c("ASV 1", "ASV 43", "ASV 7", "ASV 10", "ASV 93", "ASV 12", "ASV 44",
                                 "ASV 36","ASV 19","ASV 18","ASV 5", "ASV 56",
                                 "ASV 16", "ASV 3","ASV 26","ASV 6", "ASV 35",
                                 "ASV 21","ASV 4", "ASV 2", "ASV 29", "ASV 77",
                                 "ASV 95","ASV 134", "ASV 105", "ASV 57", "ASV 40")))
### save txt file needed to reproduce the plot
write_tsv(data_asv, "MS/data_plots/data_Fig03.txt", quote = "needed")

### plot with all
plot.asv.all <- ggplot(data_asv, aes(x = Temp, y = Abundance*100, color = Treatment,
                                     shape = Treatment, linetype = Treatment, group = Treatment)) +
  geom_point(size = 1.6) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = treat_col, labels = treat_label) +
  scale_shape_discrete(labels = treat_label) +
  scale_x_continuous(limits = c(0,31), expand = c(0,0)) +
  labs(x = "Temperature (°C)", 
       y = "Relative abundance in heavy fractions (%)") +
  guides(color = guide_legend("Treatment"), shape = guide_legend("Treatment"), 
         linetype = guide_legend("Treatment")) +
  facet_wrap(~OTU, nrow = 1) +
  theme_bw() +
  theme(panel.grid = element_blank(), panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        text = element_text(size = 7.8),
        legend.text = element_text(size = 7.8), axis.text = element_text(size = 6, colour = "black"),
        strip.background = element_rect(fill = NA, color = NA), strip.text = element_text(size = 7),
        plot.title = element_markdown(size = 7.8, hjust = 0.5, margin = margin(b = -2)),
        axis.ticks = element_line(linewidth = 0.5), axis.title = element_blank(),
        plot.margin = margin(b = 2, t = 3, r = 5.5, l = 5.5))

### subset plot for taxa
m.Sva <- plot.asv.all %+% subset(data_asv, data_asv$Genus_uncl == "Sva1033_g_unclassified") +
  labs(title = "**Sva1033**", tag = "a") +
  scale_y_continuous(limits = c(0,40), breaks = c(0,20,40))

m.Desumu <- plot.asv.all %+% subset(data_asv, data_asv$Genus_uncl == "Desulfuromusa") +
  labs(title = "***Desulfuromusa***", tag = "b") +
  scale_y_continuous(limits = c(0,56), breaks = c(0,25,50))

m.uGeo <- plot.asv.all %+% subset(data_asv, data_asv$Genus_uncl == "Geopsychrobacteraceae_g_unclassified") +
  labs(title = "**unclassified**<br>***Geopsychrobacteraceae***", tag = "c") +
  scale_y_continuous(limits = c(0,5.6), breaks = c(0,2.5,5))

m.uDesulfu.16 <- plot.asv.all %+% subset(data_asv, data_asv$OTU == "ASV 16") +
  labs(title = "**unclassified**<br>***Desulfuromonadales***", tag = "d") +
  scale_y_continuous(limits = c(0,23), breaks = c(0,10,20))

m.Desumon <- plot.asv.all %+% subset(data_asv, data_asv$Genus_uncl == "Desulfuromonas") +
  labs(title = "***Desulfuromonas***", tag = "e") +
  scale_y_continuous(limits = c(0,23), breaks = c(0,10,20))

m.Tri <- plot.asv.all %+% subset(data_asv, data_asv$Genus_uncl == "Trichloromonas") +
  labs(title = "***Trichloromonas***", tag = "f") +
  scale_y_continuous(limits = c(0,5.6), breaks = c(0,2.5,5))

m.Sulfuri <- plot.asv.all %+% subset(data_asv, data_asv$Genus_uncl == "Sulfurimonas") +
  labs(title = "***Sulfurimonas***", tag = "g") +
  scale_y_continuous(limits = c(0,40), breaks = c(0,20,40))

m.Sulfuro <- plot.asv.all %+% subset(data_asv, data_asv$Genus_uncl == "Sulfurospirillum") +
  labs(title = "***Sulfurospirillum***", tag = "h") +
  scale_y_continuous(limits = c(0,40), breaks = c(0,20,40))

m.Desuba <- plot.asv.all %+% subset(data_asv, data_asv$Genus_uncl == "Desulfobacter") +
  labs(title = "***Desulfobacter***", tag = "i") +
  scale_y_continuous(limits = c(0,23), breaks = c(0,10,20))

m.Arc <- plot.asv.all %+% subset(data_asv, data_asv$Genus_uncl == "Arcobacteraceae_g_unclassified") +
  labs(title = "**unclassified**<br>***Arcobacteraceae***", tag = "j") +
  scale_y_continuous(limits = c(0,5.6), breaks = c(0,2.5,5))

m.Col <- plot.asv.all %+% subset(data_asv, data_asv$Genus_uncl == "Colwellia") +
  labs(title = "***Colwellia***", tag = "k") +
  scale_y_continuous(limits = c(0,5.6), breaks = c(0,2.5,5))

m.uDesulfu.57 <- plot.asv.all %+% subset(data_asv, data_asv$OTU == "ASV 57") +
  labs(title = "**unclassified**<br>***Desulfuromonadales***", tag = "l") +
  scale_y_continuous(limits = c(0,10.1), breaks = c(0,5,10))

m.Hof <- plot.asv.all %+% subset(data_asv, data_asv$Genus_uncl == "Hoeflea") +
  labs(title = "***Hoeflea***", tag = "m") +
  scale_y_continuous(limits = c(0,10.1), breaks = c(0,5,10))


ms_design <- "
NAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
NAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
NAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
NBBBBBBBBBBBBBBBBBBBBBBBCCCCCC
NBBBBBBBBBBBBBBBBBBBBBBBCCCCCC
NBBBBBBBBBBBBBBBBBBBBBBBCCCCCC
NDDDDDDEEEEEEEEEEEEEEEEEFFFFFF
NDDDDDDEEEEEEEEEEEEEEEEEFFFFFF
NDDDDDDEEEEEEEEEEEEEEEEEFFFFFF
NGGGGGGGGGGGHHHHHHHHHHHHIIIIII
NGGGGGGGGGGGHHHHHHHHHHHHIIIIII
NGGGGGGGGGGGHHHHHHHHHHHHIIIIII
NJJJJJJJJJJJKKKKKKLLLLLLMMMMMM
NJJJJJJJJJJJKKKKKKLLLLLLMMMMMM
NJJJJJJJJJJJKKKKKKLLLLLLMMMMMM
NOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
"

ylab.ms <- ggplot(data.frame(l = plot.asv.all$labels$y, x = 1, y = 1)) +
  geom_text(aes(x, y, label = l), angle = 90, size = (7.8 / (14/5))) + 
  theme_void() +
  coord_cartesian(clip = "off")

xlab.ms <- ggplot(data.frame(l = plot.asv.all$labels$x, x = 1, y = 1)) +
  geom_text(aes(x, y, label = l), size = (7.8 / (14/5))) + 
  theme_void() +
  coord_cartesian(clip = "off")

MS.plot.asv.2 <- m.Sva + m.Desumu + m.uGeo + m.uDesulfu.16 + m.Desumon + m.Tri + 
  m.Sulfuri + m.Sulfuro + m.Desuba + m.Arc + m.Col + m.uDesulfu.57 + m.Hof +
  ylab.ms + xlab.ms +
  plot_layout(guides = "collect", design = ms_design) &
  theme(legend.position = "bottom", legend.margin = margin(t = -10))

ggsave("MS/plots/MS_Fig03_temp_SIP_enr.species.png", plot = MS.plot.asv.2,
       width = 18, height = 21, units = "cm")
ggsave("MS/plots/MS_Fig03_temp_SIP_enr.species.pdf", plot = MS.plot.asv.2,
       width = 18, height = 21, units = "cm")
ggsave("MS/plots/MS_Fig03_temp_SIP_enr.species.emf", plot = MS.plot.asv.2,
       width = 18, height = 21, units = "cm",
       device = {function(filename, ...) devEMF::emf(file = filename, ...)})

## MS supl FigSX Fig03 same scale ####
### plot with all
plot.asv.all.supl <- ggplot(data_asv, aes(x = Temp, y = Abundance*100, color = Treatment,
                                          shape = Treatment, linetype = Treatment, group = Treatment)) +
  geom_point(size = 1.6) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = treat_col, labels = treat_label) +
  scale_shape_discrete(labels = treat_label) +
  scale_x_continuous(limits = c(0,31), expand = c(0,0)) +
  scale_y_continuous(limits = c(0,56), breaks = c(0,25,50)) +
  labs(x = "Temperature (°C)", 
       y = "Relative abundance in heavy fractions (%)") +
  guides(color = guide_legend("Treatment"), shape = guide_legend("Treatment"), 
         linetype = guide_legend("Treatment")) +
  facet_wrap(~OTU, nrow = 1) +
  theme_bw() +
  theme(panel.grid = element_blank(), panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        text = element_text(size = 7.8),
        legend.text = element_text(size = 7.8), axis.text = element_text(size = 6, colour = "black"),
        strip.background = element_rect(fill = NA, color = NA), strip.text = element_text(size = 7),
        plot.title = element_markdown(size = 7.8, hjust = 0.5, margin = margin(b = -2)),
        axis.ticks = element_line(linewidth = 0.5), axis.title = element_blank(),
        plot.margin = margin(b = 2, t = 3, r = 5.5, l = 5.5))

### subset plot for taxa
s.m.Sva <- plot.asv.all %+% subset(data_asv, data_asv$Genus_uncl == "Sva1033_g_unclassified") +
  labs(title = "**Sva1033**", tag = "a") +
  scale_y_continuous(limits = c(0,56), breaks = c(0,25,50))

s.m.Desumu <- plot.asv.all %+% subset(data_asv, data_asv$Genus_uncl == "Desulfuromusa") +
  labs(title = "***Desulfuromusa***", tag = "b") +
  scale_y_continuous(limits = c(0,56), breaks = c(0,25,50))

s.m.uGeo <- plot.asv.all %+% subset(data_asv, data_asv$Genus_uncl == "Geopsychrobacteraceae_g_unclassified") +
  labs(title = "**unclassified**<br>***Geopsychrobacteraceae***", tag = "c") +
  scale_y_continuous(limits = c(0,56), breaks = c(0,25,50))

s.m.uDesulfu.16 <- plot.asv.all %+% subset(data_asv, data_asv$OTU == "ASV 16") +
  labs(title = "**unclassified**<br>***Desulfuromonadales***", tag = "d") +
  scale_y_continuous(limits = c(0,56), breaks = c(0,25,50))

s.m.Desumon <- plot.asv.all %+% subset(data_asv, data_asv$Genus_uncl == "Desulfuromonas") +
  labs(title = "***Desulfuromonas***", tag = "e") +
  scale_y_continuous(limits = c(0,56), breaks = c(0,25,50))

s.m.Tri <- plot.asv.all %+% subset(data_asv, data_asv$Genus_uncl == "Trichloromonas") +
  labs(title = "***Trichloromonas***", tag = "f") +
  scale_y_continuous(limits = c(0,56), breaks = c(0,25,50))

s.m.Sulfuri <- plot.asv.all %+% subset(data_asv, data_asv$Genus_uncl == "Sulfurimonas") +
  labs(title = "***Sulfurimonas***", tag = "g") +
  scale_y_continuous(limits = c(0,56), breaks = c(0,25,50))

s.m.Sulfuro <- plot.asv.all %+% subset(data_asv, data_asv$Genus_uncl == "Sulfurospirillum") +
  labs(title = "***Sulfurospirillum***", tag = "h") +
  scale_y_continuous(limits = c(0,56), breaks = c(0,25,50))

s.m.Desuba <- plot.asv.all %+% subset(data_asv, data_asv$Genus_uncl == "Desulfobacter") +
  labs(title = "***Desulfobacter***", tag = "i") +
  scale_y_continuous(limits = c(0,56), breaks = c(0,25,50))

s.m.Arc <- plot.asv.all %+% subset(data_asv, data_asv$Genus_uncl == "Arcobacteraceae_g_unclassified") +
  labs(title = "**unclassified**<br>***Arcobacteraceae***", tag = "j") +
  scale_y_continuous(limits = c(0,56), breaks = c(0,25,50))

s.m.Col <- plot.asv.all %+% subset(data_asv, data_asv$Genus_uncl == "Colwellia") +
  labs(title = "***Colwellia***", tag = "k") +
  scale_y_continuous(limits = c(0,56), breaks = c(0,25,50))

s.m.uDesulfu.57 <- plot.asv.all %+% subset(data_asv, data_asv$OTU == "ASV 57") +
  labs(title = "**unclassified**<br>***Desulfuromonadales***", tag = "l") +
  scale_y_continuous(limits = c(0,56), breaks = c(0,25,50))

s.m.Hof <- plot.asv.all %+% subset(data_asv, data_asv$Genus_uncl == "Hoeflea") +
  labs(title = "***Hoeflea***", tag = "m") +
  scale_y_continuous(limits = c(0,56), breaks = c(0,25,50))


ms_design <- "
NAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
NAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
NAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
NBBBBBBBBBBBBBBBBBBBBBBBCCCCCC
NBBBBBBBBBBBBBBBBBBBBBBBCCCCCC
NBBBBBBBBBBBBBBBBBBBBBBBCCCCCC
NDDDDDDEEEEEEEEEEEEEEEEEFFFFFF
NDDDDDDEEEEEEEEEEEEEEEEEFFFFFF
NDDDDDDEEEEEEEEEEEEEEEEEFFFFFF
NGGGGGGGGGGGHHHHHHHHHHHHIIIIII
NGGGGGGGGGGGHHHHHHHHHHHHIIIIII
NGGGGGGGGGGGHHHHHHHHHHHHIIIIII
NJJJJJJJJJJJKKKKKKLLLLLLMMMMMM
NJJJJJJJJJJJKKKKKKLLLLLLMMMMMM
NJJJJJJJJJJJKKKKKKLLLLLLMMMMMM
NOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
"

s.ylab.ms <- ggplot(data.frame(l = plot.asv.all.supl$labels$y, x = 1, y = 1)) +
  geom_text(aes(x, y, label = l), angle = 90, size = (7.8 / (14/5))) + 
  theme_void() +
  coord_cartesian(clip = "off")

s.xlab.ms <- ggplot(data.frame(l = plot.asv.all.supl$labels$x, x = 1, y = 1)) +
  geom_text(aes(x, y, label = l), size = (7.8 / (14/5))) + 
  theme_void() +
  coord_cartesian(clip = "off")

MS.plot.asv.same.scale <- s.m.Sva + s.m.Desumu + s.m.uGeo + s.m.uDesulfu.16 + s.m.Desumon + s.m.Tri + 
  s.m.Sulfuri + s.m.Sulfuro + s.m.Desuba + s.m.Arc + s.m.Col + s.m.uDesulfu.57 + s.m.Hof +
  s.ylab.ms + s.xlab.ms +
  plot_layout(guides = "collect", design = ms_design) &
  theme(legend.position = "bottom", legend.margin = margin(t = -10))

ggsave("MS/plots/supl/supl_FigS10_temp_SIP_enr.species.scale.png", plot = MS.plot.asv.same.scale,
       width = 18, height = 21, units = "cm")
ggsave("MS/plots/supl/supl_FigS10_temp_SIP_enr.species.scale.pdf", plot = MS.plot.asv.same.scale,
       width = 18, height = 21, units = "cm")
ggsave("MS/plots/supl/supl_FigS10_temp_SIP_enr.species.scale.emf", plot = MS.plot.asv.same.scale,
       width = 18, height = 21, units = "cm",
       device = {function(filename, ...) devEMF::emf(file = filename, ...)})

