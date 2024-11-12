# plot fractionation profiles
 # including which fractions are sequenced

# load packages ####
library(tidyverse)
library(patchwork)
library(ggpubr)

# import data ####
setwd("WORKINGDIRECTORY")

## project name to name files
prj <- "PROJECT"

## raw data
raw <- read_tsv(paste0(prj, "_fractionation_data.txt"))
mdata <- read_tsv(paste0(prj, "_id_translation_treatment_no.txt"))

df <- left_join(raw, mdata, by = "Incubation_ID") %>% 
  mutate(seq.fraction = factor(case_when(Fraction %in% c(3,4) ~ "Ultraheavy",
                                  Fraction %in% c(5,6) ~ "Heavy",
                                  Fraction %in% c(7,8) ~ "Midpoint",
                                  Fraction %in% c(9,10) ~ "Light",
                                  Fraction %in% c(11,12) ~ "Ultralight"),
         levels = c("Ultralight", "Light", "Midpoint", "Heavy", "Ultraheavy", "NA"),
         labels = c("Ultralight", "Light", "Midpoint", "Heavy", "Ultraheavy", "")),
         Treatment = factor(Treatment, levels = c("Acetate", "Acetate.Lepidocrocite", "Acetate.Sulfate"),
                            labels = c("Acetate", "Acetate + lepidocrocite", "Acetate + sulfate")))

# plot ####
## ng RNA per fraction
ggplot(df, aes(x = Density, y = RNA_recovered_fraction, color = Isotope, linetype = Isotope)) +
  geom_path() +
  geom_point(aes(shape = seq.fraction), size = 2) +
  geom_text(label = df$Fraction, color = "black", size = 6, size.unit = "pt") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(limits = c(0,NA), expand = c(0,0)) +
  labs(x = "Density (g/ml)", y = "Total RNA per fraction (ng)",
       color = "Isotope", linetype = "Isotope", shape = "Sequenced fraction") +
  guides(color = guide_legend(order = 1), linetype = guide_legend(order = 1),
         shape = guide_legend(order = 2, nrow = 1)) +
  facet_grid(Temp~Treatment) +
  theme_bw() +
  theme(text = element_text(size = 8, color = "black"), 
                axis.title =  element_text(size = 7),
                axis.text = element_text(size = 7, color = "black"),
                legend.text = element_text(size = 7),
                axis.ticks = element_line(linewidth = 0.5),
                panel.grid = element_blank(),
        legend.position = "bottom")

## perc RNA recovered per fraction
isotope_label <- c("unlabeled", expression({}^{13}*"C"))
frac_plot <- ggplot(df, aes(x = Density, y = perc_recovered_RNA_per_fraction, color = Isotope, linetype = Isotope)) +
  geom_path(linewidth = 1.1) +
  geom_point(aes(shape = seq.fraction), size = 3.5, stroke = 0.8, fill = "white", na.rm = T) +
  # geom_text(label = df$Fraction, color = "black", size = 6, size.unit = "pt") +
  geom_label(label = df$Fraction, 
             label.padding = unit(0.1, "lines"), label.r = unit(0.5, "lines"), label.size = 0,
             color = "black", fill = "#ffffffcc", size = 5, size.unit = "pt", na.rm = T) +
  # scale_x_continuous(limits = c(1.75, NA)) +
  scale_y_continuous(breaks = c(0,25,50)) +
  scale_color_manual(values = c("black", "red3"), labels = isotope_label) +
  scale_linetype_manual(values = c("41","solid"), labels = isotope_label) +
  scale_shape_manual(values = c(23, 22, 21, 24, 25, NA),
                     labels = c("ultra-light", "light", "midpoint", "heavy", "ultra-heavy", "")) +
  coord_cartesian(xlim = c(1.75,NA)) +
  labs(x = "Density (g/ml)", y = "RNA recovered per fraction (%)",
       color = "Isotope", linetype = "Isotope", shape = "Sequenced fraction") +
  guides(color = guide_legend(order = 1), linetype = guide_legend(order = 1),
         shape = guide_legend(order = 2, nrow = 1)) +
  facet_grid(Temp~Treatment) +
  theme_bw() +
  theme(text = element_text(size = 8, color = "black"), 
        axis.title =  element_text(size = 7.5),
        axis.text = element_text(size = 7, color = "black"),
        legend.text = element_text(size = 7), legend.title = element_text(size = 7.5, face = "bold"),
        axis.ticks = element_line(linewidth = 0.5),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "white"), strip.text = element_text(size = 8, face = "bold"),
        legend.position = "bottom", legend.box = "vertical", legend.key.width = unit(1.2, "cm"),
        legend.byrow = T, legend.margin = margin(t = -10, b = -2), legend.key.spacing.x = unit(1, "mm"),
        legend.box.spacing = unit(25, "pt"), legend.box.margin = margin(t = -10))
frac_plot

ggsave("MS/plots/supl/supl_FigS02_temp_SIP_fraction_profile.png", frac_plot, width = 18, height = 22, units = "cm")
ggsave("MS/plots/supl/supl_FigS02_temp_SIP_fraction_profile.pdf", frac_plot, width = 18, height = 22, units = "cm")
ggsave("MS/plots/supl/supl_FigS02_temp_SIP_fraction_profile.emf", frac_plot, width = 18, height = 22, units = "cm",
       device = {function(filename, ...) devEMF::emf(file = filename, ...)})


