## script to perform dbRDAs of bacterial community data
## data used created in script S_03_metacoder_import_data.R

# load packages ####
library(tidyverse)
library(vegan)
library(patchwork)
library(viridis)
library(cowplot)
library(devEMF)

# import data ####
setwd("WORKINGDIRECTORY")

## path to scripts to load
sloc <- "SCRIPTLOCATION/Scripts_to_source/"

## project name to name files
prj <- "PROJECT"

# metadata
mdata <- read_tsv("MDATAFILE.txt")
mdata$Isotope <- factor(mdata$Isotope)
mdata$Fraction <- factor(mdata$Fraction)
mdata$Treatment <- factor(mdata$Treatment)

load(paste0(prj, "_for_dbRDA_asv-table_final.RData"))


## define functions ####
source(paste0(sloc, "mydbRDA_data.R"))
# column of metadata needs to be same as rownames rel.abd.table or will give Error message

# plots for MS ####

## some labels
treat_label <- c("Acetate", "Acetate + lepidocrocite", "Acetate + sulfate")


## all fractions, all temp - done ####
dbrda_all <- mydbRDA_data(mdata, rel.abd_forsubset, var1 = "Temp", var2 = "Treatment",
                          id_column = "Name.Mapping.File")

# test results printed to console:
# F(3,205) = 20.407, p = 0.001
# temp: 53.3643  0.001
# treat: 3.9205  0.001
# dbrda1: 87%, dbrda2: 9.9%
write_tsv(dbrda_all$dbRDA_plot, "MS/data_plots/data_Fig02A.txt")

plot_all <- ggplot(dbrda_all$dbRDA_plot, 
                   aes(x = dbRDA1, y = dbRDA2, color = Temp, shape = Treatment)) +
  stat_ellipse(aes(linetype = Treatment), type = "t", show.legend = T, level = 0.6,
               linewidth = 1) +
  scale_linetype_discrete(labels = treat_label) +
  scale_colour_viridis(breaks = c(2,5,10,15,20,25,30), name = "Temperature (°C)") +
  scale_shape_discrete(labels = treat_label) +
  scale_y_continuous(breaks = seq(-4,2,1)) +
  scale_x_continuous(breaks = seq(-1,1.5,1)) +
  geom_point() +
  coord_fixed() +
  guides(color = guide_colorbar(direction = "horizontal", title.position = "top")) +
  labs(title = "All fractions from\nall treatments", 
       x = "dbRDA1 (variation 87%)",
       y = "dbRDA2 (variation 10%)") +
  annotate("text", x = 0.18, y = -3.5, 
           label = "Temperature\nF(1,205) = 53.4 ***\nTreatment\nF(2,205) = 3.9 ***",
           hjust = 0, size = 7 / (14/5)) +
  theme_cowplot() +
  theme(text = element_text(size = 8), axis.text = element_text(size = 7),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        axis.line = element_line(linewidth = 0), legend.text = element_text(size = 8),
        plot.title = element_text(size = 8, hjust = 0.5))


## 13C all fractions, all temp - done ####
sub_13C <- mysubset(mdata, rel.abd_forsubset, variable = "Isotope", levels_to_compare = "13C",
                    id_column = "Name.Mapping.File")

dbrda_13C <- mydbRDA_data(sub_13C$mdata, sub_13C$rel.abd.table, var1 = "Temp", var2 = "Treatment",
                          id_column = "Name.Mapping.File")

# test results printed to console:
# Residual: 100
#           F       p
# model: 9.6198  0.001
# temp: 24.0971  0.001
# treat: 2.3714  0.001
# dbrda1: 84%, dbrda2: 14%

write_tsv(dbrda_13C$dbRDA_plot, "MS/data_plots/data_Fig02B.txt")

plot_13C <- ggplot(dbrda_13C$dbRDA_plot, 
                   aes(x = dbRDA1, y = dbRDA2, color = Temp, shape = Treatment)) +
  stat_ellipse(aes(linetype = Treatment), type = "t", show.legend = T, level = 0.6,
               linewidth = 1) +
  scale_linetype_discrete(labels = treat_label) +
  scale_colour_viridis(breaks = c(2,5,10,15,20,25,30), name = "Temperature (°C)") + 
  scale_shape_discrete(labels = treat_label) +
  scale_x_continuous(breaks = seq(-1,2,1)) +
  scale_y_continuous(breaks = seq(-3,2,1)) +
  geom_point() +
  coord_fixed() +
  guides(color = guide_colorbar(direction = "horizontal", title.position = "top")) +
  labs(title = "All fractions from\nlabeled treatments",
       x = "dbRDA1 (variation 84%)",
       y = "dbRDA2 (variation 14%)") +
  annotate("text", x = 0.42, y = -3, 
           label = "Temperature\nF(1,100) = 24.1 ***\nTreatment\nF(2,100) = 2.4 ***",
           hjust = 0, size = 7 / (14/5)) +
  theme_cowplot() +
  theme(text = element_text(size = 8), axis.text = element_text(size = 7),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        axis.line = element_line(linewidth = 0), legend.text = element_text(size = 8),
        plot.title = element_text(size = 8, hjust = 0.5))


## 13C heavy + u-heavy fractions, all temp - done ####
sub_13C_h_uheavy <- mysubset(sub_13C$mdata, sub_13C$rel.abd.table, variable = "Fraction",
                             levels_to_compare = c("Heavy", "Ultraheavy"), id_column = "Name.Mapping.File")

dbrda_13C_h_uheavy <- mydbRDA_data(sub_13C_h_uheavy$mdata, sub_13C_h_uheavy$rel.abd.table, var1 = "Temp",
                                   var2 = "Treatment", id_column = "Name.Mapping.File")

# test results printed to console:
# Residual: 37
#           F       p
# model: 7.2409  0.001
# temp: 16.2171  0.001
# treat: 2.7298  0.002
# dbrda1: 76%, dbrda2: 22%

write_tsv(dbrda_13C_h_uheavy$dbRDA_plot, "MS/data_plots/data_Fig02C.txt")

plot_13_h_uheavy <- ggplot(dbrda_13C_h_uheavy$dbRDA_plot, 
                           aes(x = dbRDA1, y = dbRDA2, color = Temp, shape = Treatment)) +
  stat_ellipse(aes(linetype = Treatment), type = "t", show.legend = T, level = 0.6,
               linewidth = 1) +
  scale_linetype_discrete(labels = treat_label) +
  scale_colour_viridis(breaks = c(2,5,10,15,20,25,30), name = "Temperature (°C)") + 
  scale_shape_discrete(labels = treat_label) +
  geom_point() +
  coord_fixed() +
  scale_y_continuous(limits = c(-1.95,1.6), breaks = seq(-1.8,2,1)) +
  scale_x_continuous(breaks = seq(-1,2,1)) +
  guides(color = guide_colorbar(direction = "horizontal", title.position = "top")) +
  labs(title = "Heavy + ultra-heavy fractions\n from labeled treatments",
       x = "dbRDA1 (variation 76%)",
       y = "dbRDA2 (variation 22%)") +
  annotate("text", x = -1.42, y = -1.75, 
           label = "Temperature: F(1,37) = 16.2 ***\nTreatment: F(2,37) = 2.7 **",
           hjust = 0, size = 7 / (14/5)) +
  theme_cowplot() +
  theme(text = element_text(size = 8), axis.text = element_text(size = 7),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        axis.line = element_line(linewidth = 0), legend.text = element_text(size = 8),
        plot.title = element_text(size = 8, hjust = 0.5),
        legend.key.width = unit(0.75, "cm"), legend.key.height = unit(0.5, "cm"))

## combine plots for MS ####
comb_plot_dbRDA2 <-  plot_all + theme(legend.position = "none") +
  plot_13C + theme(legend.position = "none") +
  (plot_13_h_uheavy / guide_area() + 
     plot_layout(guides = 'collect', heights = c(1,0.5))) +
  plot_annotation(tag_levels = 'a') &
  theme(plot.title = element_text(hjust = 0.5),
        legend.box.margin = margin(t = 0.8, unit = "cm"), legend.title = element_text(face = "bold"),
        legend.key.spacing.y = unit(-0.1, "cm"))


ggsave("MS/plots/MS_Fig02_temp_SIP_dbRDA.png",
       plot = comb_plot_dbRDA2, height = 12, width = 18, units = "cm")
ggsave("MS/plots/MS_Fig02_temp_SIP_dbRDA.pdf",
       plot = comb_plot_dbRDA2, height = 12, width = 18, units = "cm",
       useDingbats = F)
ggsave("MS/plots/MS_Fig02_temp_SIP_dbRDA.emf",
       plot = comb_plot_dbRDA2, height = 12, width = 18, units = "cm",
       device = {function(filename, ...) devEMF::emf(file = filename, ...)})


## 12C all fractions, all temp, for supplementary ####
sub_12C <- mysubset(mdata, rel.abd_forsubset, variable = "Isotope", levels_to_compare = "12C",
                    id_column = "Name.Mapping.File")

dbrda_12C <- mydbRDA_data(sub_12C$mdata, sub_12C$rel.abd.table, var1 = "Temp", var2 = "Treatment",
                          id_column = "Name.Mapping.File")

# test results printed to console:
# Residual: 101
#           F       p
# model: 13.67    0.001
# temp:  35.0324  0.001
# treat: 2.9892   0.003
# dbrda1: 87%, dbrda2: 9%

write_tsv(dbrda_12C$dbRDA_plot, "MS/data_plots/data_supl_FigS04_dbrda_12C.txt")

plot_12C <- ggplot(dbrda_12C$dbRDA_plot, 
                   aes(x = dbRDA1, y = dbRDA2, color = Temp, shape = Treatment)) +
  stat_ellipse(aes(linetype = Treatment), type = "t", show.legend = T, level = 0.6,
               linewidth = 1) +
  scale_linetype_discrete(labels = treat_label) +
  scale_colour_viridis(breaks = c(2,5,10,15,20,25,30), name = "Temperature (°C)") + 
  scale_shape_discrete(labels = treat_label) +
  geom_point() +
  coord_fixed() +
  scale_y_continuous(limits = c(-1.25,2.2)) +
  scale_x_continuous(breaks = seq(-1,2,1)) +
  guides(color = guide_colorbar(direction = "horizontal", title.position = "top")) +
  labs(title = "All fractions from\nall unlabeled treatments",
       x = "dbRDA1 (variation 87%)",
       y = "dbRDA2 (variation 9%)") +
  annotate("text", x = 0.23, y = 1.9, 
           label = "Temperature\nF(1,101) = 35.0 ***\nTreatment\nF(2,101) = 3.0 **",
           hjust = 0, size = 7 / (14/5)) +
  theme_cowplot() +
  theme(text = element_text(size = 8), axis.text = element_text(size = 7),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        plot.background = element_rect(fill = "white"),
        axis.line = element_line(linewidth = 0), legend.text = element_text(size = 8),
        plot.title = element_text(size = 8, hjust = 0.5),
        legend.title = element_text(face = "bold"),
        legend.key.size = unit(0.65, "cm"),
        legend.key.spacing.y = unit(-0.1, "cm"))

ggsave("MS/plots/supl/supl_FigS04_temp_SIP_12C_dbrda.png", plot_12C, 
       width = 10, height = 8, units = "cm")
ggsave("MS/plots/supl/supl_FigS04_temp_SIP_12C_dbrda.pdf", plot_12C, 
       width = 10, height = 8, units = "cm")
ggsave("MS/plots/supl/supl_FigS04_temp_SIP_12C_dbrda.emf", plot_12C, 
       width = 10, height = 8, units = "cm",
       device = {function(filename, ...) devEMF::emf(file = filename, ...)})


# models for subset of temperatures ####

## all fractions different temperature subsets ####

### all fractions temp 2 - 20 °C ####
sub_temp2_20 <- mysubset(mdata, rel.abd_forsubset, variable = "Temp",
                         levels_to_compare = c(2,5,10,15,20), 
                         id_column = "Name.Mapping.File")

dbrda_temp2_20 <- mydbRDA_data(sub_temp2_20$mdata, sub_temp2_20$rel.abd.table, 
                               var1 = "Temp",var2 = "Treatment", 
                               id_column = "Name.Mapping.File")

# test results printed to console:
# Residual: 145
#           F       p
# model: 10.091    0.001
# temp: 18.7423    0.001
# treat: 5.7649    0.001
# dbrda1: 63%, dbrda2: 30%


### all fractions temp 2 - 15 °C ####
sub_temp2_15 <- mysubset(mdata, rel.abd_forsubset, variable = "Temp",
                         levels_to_compare = c(2,5,10,15), 
                         id_column = "Name.Mapping.File")

dbrda_temp2_15 <- mydbRDA_data(sub_temp2_15$mdata, sub_temp2_15$rel.abd.table, 
                               var1 = "Temp",var2 = "Treatment", 
                               id_column = "Name.Mapping.File")

# test results printed to console:
# Residual: 115
#           F       p
# model: 7.6182    0.001
# temp: 13.6043    0.001
# treat: 4.6238    0.001
# dbrda1: 60%, dbrda2: 33%


### all fractions temp 2 - 10 °C ####
sub_temp2_10 <- mysubset(mdata, rel.abd_forsubset, variable = "Temp",
                         levels_to_compare = c(2,5,10), 
                         id_column = "Name.Mapping.File")

dbrda_temp2_10 <- mydbRDA_data(sub_temp2_10$mdata, sub_temp2_10$rel.abd.table, 
                               var1 = "Temp",var2 = "Treatment", 
                               id_column = "Name.Mapping.File")

# test results printed to console:
# Residual: 85
#           F       p
# model: 5.1237    0.001
# temp:  8.3267    0.001
# treat: 3.5205    0.001
# dbrda1: 55%, dbrda2: 36%


### all fractions temp 2 - 5 °C ####
sub_temp2_5 <- mysubset(mdata, rel.abd_forsubset, variable = "Temp",
                         levels_to_compare = c(2,5), 
                         id_column = "Name.Mapping.File")

dbrda_temp2_5 <- mydbRDA_data(sub_temp2_5$mdata, sub_temp2_5$rel.abd.table, 
                               var1 = "Temp",var2 = "Treatment", 
                               id_column = "Name.Mapping.File")

# test results printed to console:
# Residual: 56
#           F       p
# model: 2.781    0.003
# temp:  3.3769    0.013
# treat: 2.4830    0.007
# dbrda1: 51%, dbrda2: 40%


### all fractions temp 15 - 30 °C ####
sub_temp15_30 <- mysubset(mdata, rel.abd_forsubset, variable = "Temp",
                          levels_to_compare = c(15,20,25,30), 
                          id_column = "Name.Mapping.File")

dbrda_temp15_30 <- mydbRDA_data(sub_temp15_30$mdata, sub_temp15_30$rel.abd.table, 
                                var1 = "Temp",var2 = "Treatment", 
                                id_column = "Name.Mapping.File")

# test results printed to console:
# Residual: 116
#           F       p
# model: 15.007     0.001
# temp:  39.8219    0.001
# treat: 2.5995     0.004
# dbrda1: 89%, dbrda2: 8%


### all fractions temp 20 - 30 °C ####
sub_temp20_30 <- mysubset(mdata, rel.abd_forsubset, variable = "Temp",
                          levels_to_compare = c(20,25,30), 
                          id_column = "Name.Mapping.File")

dbrda_temp20_30 <- mydbRDA_data(sub_temp20_30$mdata, sub_temp20_30$rel.abd.table, 
                                var1 = "Temp",var2 = "Treatment", 
                                id_column = "Name.Mapping.File")

# test results printed to console:
# Residual: 86
#           F       p
# model: 14.563     0.001
# temp:  38.4987     0.001
# treat: 2.5952      0.007
# dbrda1: 89%, dbrda2: 8%


### all fractions temp 25 - 30 °C ####
sub_temp25_30 <- mysubset(mdata, rel.abd_forsubset, variable = "Temp",
                          levels_to_compare = c(25,30), 
                          id_column = "Name.Mapping.File")

dbrda_temp25_30 <- mydbRDA_data(sub_temp25_30$mdata, sub_temp25_30$rel.abd.table, 
                                var1 = "Temp",var2 = "Treatment", 
                                id_column = "Name.Mapping.File")

# test results printed to console:
# Residual: 56
#           F       p
# model: 18.127     0.001
# temp:  48.2551    0.001
# treat: 3.0632     0.005
# dbrda1: 89%, dbrda2: 7%, dbrda3: 4%


### subsetted fractions of 2 - 15 °C ####

## 13C all fractions temp 2 - 15 °C
sub_temp2_15 <- mysubset(mdata, rel.abd_forsubset, variable = "Temp",
                         levels_to_compare = c(2,5,10,15), 
                         id_column = "Name.Mapping.File")

sub_13C_temp2_15 <- mysubset(sub_temp2_15$mdata, sub_temp2_15$rel.abd.table,
                             variable = "Isotope", levels_to_compare = "13C",
                             id_column = "Name.Mapping.File")

dbrda_13C_temp2_15 <- mydbRDA_data(sub_13C_temp2_15$mdata, 
                                   sub_13C_temp2_15$rel.abd.table, 
                                   var1 = "Temp",var2 = "Treatment", 
                                   id_column = "Name.Mapping.File")

# test results printed to console:
# Residual: 55
#           F       p
# model: 4.3216    0.001
# temp:  6.7766    0.001
# treat: 3.0955    0.003
# dbrda1: 54%, dbrda2: 42%, dbrda3: 4%


## 13C heavy + ultra-heavy fractions temp 2 - 15 °C
sub_temp2_15 <- mysubset(mdata, rel.abd_forsubset, variable = "Temp",
                         levels_to_compare = c(2,5,10,15), 
                         id_column = "Name.Mapping.File")

sub_13C_temp2_15 <- mysubset(sub_temp2_15$mdata, sub_temp2_15$rel.abd.table,
                             variable = "Isotope", levels_to_compare = "13C",
                             id_column = "Name.Mapping.File")

sub_13C_h_temp2_15 <- mysubset(sub_13C_temp2_15$mdata, sub_13C_temp2_15$rel.abd.table,
                             variable = "Fraction", 
                             levels_to_compare = c("Heavy", "Ultraheavy"),
                             id_column = "Name.Mapping.File")

dbrda_13C_h_temp2_15 <- mydbRDA_data(sub_13C_h_temp2_15$mdata, 
                                     sub_13C_h_temp2_15$rel.abd.table, 
                                     var1 = "Temp",var2 = "Treatment", 
                                     id_column = "Name.Mapping.File")

# test results printed to console:
# Residual: 19
#           F       p
# model: 7.532    0.001
# temp:  10.6190   0.001
# treat: 5.9884    0.001
# dbrda1: 67%, dbrda2: 30%



### subsetted fractions of 2 - 10 °C ####

## 13C all fractions temp 2 - 10 °C 
sub_temp2_10 <- mysubset(mdata, rel.abd_forsubset, variable = "Temp",
                         levels_to_compare = c(2,5,10), 
                         id_column = "Name.Mapping.File")

sub_13C_temp2_10 <- mysubset(sub_temp2_10$mdata, sub_temp2_10$rel.abd.table,
                             variable = "Isotope", levels_to_compare = "13C",
                             id_column = "Name.Mapping.File")

dbrda_13C_temp2_10 <- mydbRDA_data(sub_13C_temp2_10$mdata, 
                                   sub_13C_temp2_10$rel.abd.table, 
                                   var1 = "Temp",var2 = "Treatment", 
                                   id_column = "Name.Mapping.File")

# test results printed to console:
# Residual: 40
#           F       p
# model: 2.9818    0.002
# temp:  4.6958    0.009
# treat: 2.1354    0.051
# dbrda1: 57%, dbrda2: 37%


## 13C heavy + ultra-heavy fractions temp 2 - 10 °C
sub_temp2_10 <- mysubset(mdata, rel.abd_forsubset, variable = "Temp",
                         levels_to_compare = c(2,5,10), 
                         id_column = "Name.Mapping.File")

sub_13C_temp2_10 <- mysubset(sub_temp2_10$mdata, sub_temp2_10$rel.abd.table,
                             variable = "Isotope", levels_to_compare = "13C",
                             id_column = "Name.Mapping.File")

sub_13C_h_temp2_10 <- mysubset(sub_13C_temp2_10$mdata, sub_13C_temp2_10$rel.abd.table,
                               variable = "Fraction", 
                               levels_to_compare = c("Heavy", "Ultraheavy"),
                               id_column = "Name.Mapping.File")

dbrda_13C_h_temp2_10 <- mydbRDA_data(sub_13C_h_temp2_10$mdata, 
                                     sub_13C_h_temp2_10$rel.abd.table, 
                                     var1 = "Temp",var2 = "Treatment", 
                                     id_column = "Name.Mapping.File")

# test results printed to console:
# Residual: 13
#           F       p
# model: 5.2649    0.001
# temp:  6.1471    0.002
# treat: 4.7160    0.001
# dbrda1: 64%, dbrda2: 32% dbrda3: 4%



### subsetted fractions of 2 - 5 °C ####

## 13C all fractions temp 2 - 5 °C
sub_temp2_5 <- mysubset(mdata, rel.abd_forsubset, variable = "Temp",
                         levels_to_compare = c(2,5), 
                         id_column = "Name.Mapping.File")

sub_13C_temp2_5 <- mysubset(sub_temp2_5$mdata, sub_temp2_5$rel.abd.table,
                             variable = "Isotope", levels_to_compare = "13C",
                             id_column = "Name.Mapping.File")

dbrda_13C_temp2_5 <- mydbRDA_data(sub_13C_temp2_5$mdata, 
                                  sub_13C_temp2_5$rel.abd.table, 
                                  var1 = "Temp",var2 = "Treatment", 
                                  id_column = "Name.Mapping.File")

# test results printed to console:
# Residual: 26
#           F       p
# model: 1.6074    0.139
# temp:  1.3652    0.204
# treat: 1.7285    0.145
# dbrda1: 62%, dbrda2: 28%, dbrda3: 10%


## 13C heavy + ultra-heavy fractions temp 2 - 5 °C
sub_temp2_5 <- mysubset(mdata, rel.abd_forsubset, variable = "Temp",
                        levels_to_compare = c(2,5), 
                        id_column = "Name.Mapping.File")

sub_13C_temp2_5 <- mysubset(sub_temp2_5$mdata, sub_temp2_5$rel.abd.table,
                            variable = "Isotope", levels_to_compare = "13C",
                            id_column = "Name.Mapping.File")

sub_13C_h_temp2_5 <- mysubset(sub_13C_temp2_5$mdata, sub_13C_temp2_5$rel.abd.table,
                               variable = "Fraction", 
                               levels_to_compare = c("Heavy", "Ultraheavy"),
                               id_column = "Name.Mapping.File")

dbrda_13C_h_temp2_5 <- mydbRDA_data(sub_13C_h_temp2_5$mdata, 
                                    sub_13C_h_temp2_5$rel.abd.table, 
                                    var1 = "Temp",var2 = "Treatment", 
                                    id_column = "Name.Mapping.File")

# test results printed to console:
# Residual: 8
#           F       p
# model: 2.885    0.014
# temp:  2.1087    0.094 
# treat: 3.2732    0.012 
# dbrda1: 76%, dbrda2: 19%



### subsetted fractions of 20 - 30 °C ####

## 13C all fractions temp 20 - 30 °C 
sub_temp20_30 <- mysubset(mdata, rel.abd_forsubset, variable = "Temp",
                          levels_to_compare = c(20,25,30), 
                          id_column = "Name.Mapping.File")

sub_13C_temp20_30 <- mysubset(sub_temp20_30$mdata, sub_temp20_30$rel.abd.table,
                              variable = "Isotope", levels_to_compare = "13C",
                              id_column = "Name.Mapping.File")

dbrda_13C_temp20_30 <- mydbRDA_data(sub_13C_temp20_30$mdata, 
                                    sub_13C_temp20_30$rel.abd.table, 
                                    var1 = "Temp",var2 = "Treatment", 
                                    id_column = "Name.Mapping.File")

# test results printed to console:
# Residual: 41
#           F       p
# model: 7.2148    0.001
# temp: 19.0103    0.001
# treat: 1.3171    0.196    
# dbrda1: 88%, dbrda2: 8%


## 13C heavy + ultra-heavy fractions temp 20 - 30 °C
sub_temp20_30 <- mysubset(mdata, rel.abd_forsubset, variable = "Temp",
                          levels_to_compare = c(20,25,30), 
                          id_column = "Name.Mapping.File")

sub_13C_temp20_30 <- mysubset(sub_temp20_30$mdata, sub_temp20_30$rel.abd.table,
                              variable = "Isotope", levels_to_compare = "13C",
                              id_column = "Name.Mapping.File")

sub_13C_h_temp20_30 <- mysubset(sub_13C_temp20_30$mdata, sub_13C_temp20_30$rel.abd.table,
                                variable = "Fraction", 
                                levels_to_compare = c("Heavy", "Ultraheavy"),
                                id_column = "Name.Mapping.File")
 
dbrda_13C_h_temp20_30 <- mydbRDA_data(sub_13C_h_temp20_30$mdata, 
                                      sub_13C_h_temp20_30$rel.abd.table, 
                                      var1 = "Temp",var2 = "Treatment", 
                                      id_column = "Name.Mapping.File")

# test results printed to console:
# Residual: 14
#           F       p
# model: 5.6436    0.002
# temp:  13.7792   0.001
# treat: 1.5758    0.183    
# dbrda1: 82%, dbrda2: 15%



### subsetted fractions of 25 - 30 °C ####

## 13C all fractions temp 25 - 30 °C
sub_temp25_30 <- mysubset(mdata, rel.abd_forsubset, variable = "Temp",
                          levels_to_compare = c(25,30), 
                          id_column = "Name.Mapping.File")

sub_13C_temp25_30 <- mysubset(sub_temp25_30$mdata, sub_temp25_30$rel.abd.table,
                              variable = "Isotope", levels_to_compare = "13C",
                              id_column = "Name.Mapping.File")

dbrda_13C_temp25_30 <- mydbRDA_data(sub_13C_temp25_30$mdata, 
                                    sub_13C_temp25_30$rel.abd.table, 
                                    var1 = "Temp",var2 = "Treatment", 
                                    id_column = "Name.Mapping.File")

# test results printed to console:
# Residual: 26
#           F       p
# model: 7.9428    0.001
# temp: 21.1318    0.001
# treat: 1.3483    0.220
# dbrda1: 89%, dbrda2: 6%


## 13C heavy + ultra-heavy fractions temp 25 - 30 °C
sub_temp25_30 <- mysubset(mdata, rel.abd_forsubset, variable = "Temp",
                          levels_to_compare = c(25,30), 
                          id_column = "Name.Mapping.File")

sub_13C_temp25_30 <- mysubset(sub_temp25_30$mdata, sub_temp25_30$rel.abd.table,
                              variable = "Isotope", levels_to_compare = "13C",
                              id_column = "Name.Mapping.File")

sub_13C_h_temp25_30 <- mysubset(sub_13C_temp25_30$mdata, sub_13C_temp25_30$rel.abd.table,
                                variable = "Fraction", 
                                levels_to_compare = c("Heavy", "Ultraheavy"),
                                id_column = "Name.Mapping.File")

dbrda_13C_h_temp25_30 <- mydbRDA_data(sub_13C_h_temp25_30$mdata, 
                                      sub_13C_h_temp25_30$rel.abd.table, 
                                      var1 = "Temp",var2 = "Treatment", 
                                      id_column = "Name.Mapping.File")

# test results printed to console:
# Residual: 8
#           F       p
# model: 7.7303    0.001
# temp:  19.2414   0.001
# treat: 1.9748    0.136    
# dbrda1: 83%, dbrda2: 12%

