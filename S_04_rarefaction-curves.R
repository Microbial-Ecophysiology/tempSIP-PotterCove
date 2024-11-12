# This script imports taxmap objects created in script metacoder_import_data.R produced through dada2 pipeline

# explore data and create rarefaction curves 

# server part ####

## packages ####
library(tidyverse)
library(taxa)
library(metacoder)
library(iNEXT)
library(cowplot)

## import data ####
setwd("WORKINGDIRECTORY")

## project name to name files
prj <- "PROJECT"

load(paste0("taxmap_", prj, ".RData"))

mdata <- read_tsv("MDATAFILE.txt")

sample.names <- mdata$Name.Mapping.File

asv_counts_clean <- as.data.frame(obj$data$asv_table[, sample.names])


## rarefaction curves ####
# calculate variety of statistics with iNEXT
# can vary no. of knots calculated. the more the longer it will take
stat_all <- iNEXT(asv_counts_clean, q = c(0,2), datatype = "abundance", knots = 100,
                  endpoint = max(colSums(asv_counts_clean)))

save(stat_all, file = paste0("iNext_", prj, ".RData"))


# PC part ####
## packages ####
library(tidyverse)
library(taxa)
library(metacoder)
library(iNEXT)
library(ggrepel)
library(cowplot)
library(RColorBrewer)
library(patchwork)
library(devEMF)


## import data ####
setwd("WORKINGDIRECTORY")

## project name to name files
prj <- "PROJECT"

## taxmap object
load(paste0("taxmap_", prj, ".RData"))

## metadata
mdata <- read_tsv("MDATAFILE.txt")
sample.names <- mdata$Name.Mapping.File

asv_counts_clean <- as.data.frame(obj$data$asv_table[, sample.names])


## basic data overview - counts per sample ####
hist(colSums(obj$data$asv_table[, sample.names]))


## total reads and ASVs per sample
mdata_add <- mdata %>% 
  left_join((
    colSums(obj$data$asv_table[, sample.names]) %>% 
      enframe(name = "Name.Mapping.File", value = "totReads")
  )) %>% 
  left_join((
    apply(obj$data$asv_table[, sample.names] , 2, function(x) sum(x > 0)) %>% 
      enframe(name = "Name.Mapping.File", value = "totASVs")
  ))

write.table(mdata_add, "MDATAFILE_totReads-totOTUS.txt", sep = "\t", 
            quote = F, row.names = FALSE)


sum(mdata_add$totReads >= 1000)
# 201
sum(mdata_add$totReads >= 2000)
# 187
sum(mdata_add$totReads < 1000)
# 8
sum(mdata_add$totReads < 2000)
# 22
min(mdata_add$totReads)
# 527

hist(mdata_add$totReads, breaks = 100)

subset(mdata_add, mdata_add$totReads == 527)$Name.Mapping.File
# "Bacteria.20degrees.12C.Acetate.Sulfate.Fraction.Light"



## rarefaction curves ####
### dataframe for plotting ####
load(paste0("iNext_", prj, ".RData"))

stat_all$DataInfo
stat_all$AsyEst

# export iNEXT data in separate data frame to plot manually
 # add mdata to make nicer plot
df <- fortify(stat_all, type = 1)

cols_factor <- c("Method", "Temp", "Treatment")
dfm <- left_join(df, mdata, by = c("Assemblage" = "Name.Mapping.File")) %>% 
  mutate(frac_iso = paste(Isotope, Fraction, sep = "_"),
         across(all_of(cols_factor), as.factor),
         Order.q = factor(Order.q, levels = c("0", "2"),
                        labels = c("Species richness", "Linearized\nsimpson diversity")),
         frac_iso = factor(frac_iso, levels = c("12C_Ultralight", "13C_Ultralight", "12C_Light", "13C_Light",
                                                "12C_Midpoint", "13C_Midpoint", "12C_Heavy", "13C_Heavy",
                                                "12C_Ultraheavy", "13C_Ultraheavy"),
                           labels = c("12C Ultra-light", "13C Ultra-light", "12C Light", "13C Light",
                                      "12C Midpoint", "13C Midpoint", "12C Heavy", "13C Heavy",
                                      "12C Ultra-heavy", "13C Ultra-heavy")))


### Rarefaction curves for MS ####
plot_o <- ggplot(dfm, aes(x = x/10000, y = y, color = Fraction.label, shape = Isotope_mod)) +
  geom_line(data = subset(dfm, dfm$Method == "Rarefaction"), 
            aes(group = Assemblage), linetype = "solid", linewidth = 0.6) +
  geom_point(data = subset(dfm, dfm$Method == "Observed"), 
             aes(), size = 1) +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Number of individuals x10000", y = "Species diversity", color = "Fraction", shape = "Treatment type") +
  facet_grid(Order.q ~ Substrate_long, scales = "free_y")
  
## function for plots
plot_subset <- function(plot_data, subset_var, subset_val) {
  plot <- ggplot(subset(plot_data, plot_data[[subset_var]] %in% subset_val), 
                 aes(x = x, y = y, color = frac_iso)) +
    geom_path(data = subset(plot_data, plot_data[["Method"]] == "Rarefaction" & 
                              plot_data[[subset_var]] %in% subset_val), 
              aes(group = Assemblage), linetype = "solid", linewidth = 0.3) +
    geom_point(data = subset(plot_data, plot_data[["Method"]] == "Observed" & 
                               plot_data[[subset_var]] %in% subset_val), 
               aes(), size = 0.6) +
    labs(x = "Number of individuals", y = "Species diversity", color = "Fraction") +
    scale_color_brewer(palette = "Paired") +
    guides(color = guide_legend(nrow = 2)) +
    facet_grid(Order.q~Treatment, scales = "free") + 
    theme_bw(base_line_size = 0.8, base_rect_size = 0.8) +
    theme(panel.grid.minor = element_blank(), 
          strip.background = element_rect(colour = "black", fill = "white"),
          text = element_text(size = 7, color = "black"), axis.text = element_text(size = 6, color = "black"),
          strip.text = element_text(size = 7), legend.text = element_text(size = 7),
          legend.position = "bottom", panel.background = element_blank(),
          legend.direction = "vertical", legend.title.position = "left",
          plot.title = element_text(face = "bold", hjust = 0.5, size = 8),
          legend.margin = margin(0,0,0,0), legend.box.margin = margin(-10,0,0,0),
          legend.key.spacing.y = unit(-0.1, "cm"))
  return(plot)
}

## make subset plots
plot_2 <- plot_subset(plot_data = dfm, subset_var = "Temp", subset_val = c("2")) +
  coord_cartesian(xlim = c(0, 50000)) +
  labs(title = "Incubation at 2 °C")

plot_5 <- plot_subset(plot_data = dfm, subset_var = "Temp", subset_val = c("5")) +
  coord_cartesian(xlim = c(0, 50000)) +
  labs(title = "Incubation at 5 °C")

plot_10 <- plot_subset(plot_data = dfm, subset_var = "Temp", subset_val = c("10")) +
  coord_cartesian(xlim = c(0, 50000)) +
  labs(title = "Incubation at 10 °C")

plot_15 <- plot_subset(plot_data = dfm, subset_var = "Temp", subset_val = c("15")) +
  coord_cartesian(xlim = c(0, 50000)) +
  labs(title = "Incubation at 15 °C")

plot_20 <- plot_subset(plot_data = dfm, subset_var = "Temp", subset_val = c("20")) +
  coord_cartesian(xlim = c(0, 50000)) +
  labs(title = "Incubation at 20 °C")

plot_25 <- plot_subset(plot_data = dfm, subset_var = "Temp", subset_val = c("25")) +
  coord_cartesian(xlim = c(0, 50000)) +
  labs(title = "Incubation at 25 °C")

plot_30 <- plot_subset(plot_data = dfm, subset_var = "Temp", subset_val = c("30")) +
  coord_cartesian(xlim = c(0, 50000)) +
  labs(title = "Incubation at 30 °C")



## rarefaction plots with 20000 as max limit ####
plot_2 <- plot_subset(plot_data = dfm, subset_var = "Temp", subset_val = c("2")) +
  coord_cartesian(xlim = c(0, 20000)) +
  labs(title = "Incubation at 2 °C")

plot_5 <- plot_subset(plot_data = dfm, subset_var = "Temp", subset_val = c("5")) +
  coord_cartesian(xlim = c(0, 20000)) +
  labs(title = "Incubation at 5 °C")

plot_10 <- plot_subset(plot_data = dfm, subset_var = "Temp", subset_val = c("10")) +
  coord_cartesian(xlim = c(0, 20000)) +
  labs(title = "Incubation at 10 °C")

plot_15 <- plot_subset(plot_data = dfm, subset_var = "Temp", subset_val = c("15")) +
  coord_cartesian(xlim = c(0, 20000)) +
  labs(title = "Incubation at 15 °C")

plot_20 <- plot_subset(plot_data = dfm, subset_var = "Temp", subset_val = c("20")) +
  coord_cartesian(xlim = c(0, 20000)) +
  labs(title = "Incubation at 20 °C")

plot_25 <- plot_subset(plot_data = dfm, subset_var = "Temp", subset_val = c("25")) +
  coord_cartesian(xlim = c(0, 20000)) +
  labs(title = "Incubation at 25 °C")

plot_30 <- plot_subset(plot_data = dfm, subset_var = "Temp", subset_val = c("30")) +
  coord_cartesian(xlim = c(0, 20000)) +
  labs(title = "Incubation at 30 °C")


## combine plots ####
comb1 <- plot_2 + plot_5 + plot_10 + guide_area() + plot_annotation(tag_levels ='a') +
  plot_layout(guides = "collect", ncol = 1, heights = c(1,1,1,0.2))
ggsave("MS/plots/supl/supl_FigS17_temp_SIP_rarefaction-curves_comb1_max20000.pdf", plot = comb1,
       width = 16.5, height = 22, units = "cm")
ggsave("MS/plots/supl/supl_FigS17_temp_SIP_rarefaction-curves_comb1_max20000.png", plot = comb1,
       width = 16.5, height = 22, units = "cm")
ggsave("MS/plots/supl/supl_FigS17_temp_SIP_rarefaction-curves_comb1_max20000.emf", plot = comb1,
       width = 16.5, height = 22, units = "cm",
       device = {function(filename, ...) devEMF::emf(file = filename, ...)})

comb2 <- plot_15 + plot_20 + plot_25 + guide_area() + plot_annotation(tag_levels ='a') +
  plot_layout(guides = "collect", ncol = 1, heights = c(1,1,1,0.2))
ggsave("MS/plots/supl/supl_FigS18_temp_SIP_rarefaction-curves_comb2_max20000.pdf", plot = comb2,
       width = 16.5, height = 22, units = "cm")
ggsave("MS/plots/supl/supl_FigS18_temp_SIP_rarefaction-curves_comb2_max20000.png", plot = comb2,
       width = 16.5, height = 22, units = "cm")
ggsave("MS/plots/supl/supl_FigS18_temp_SIP_rarefaction-curves_comb2_max20000.emf", plot = comb2,
       width = 16.5, height = 22, units = "cm",
       device = {function(filename, ...) devEMF::emf(file = filename, ...)})


plot30_2 <- plot_30 + theme(plot.margin = margin(10,10,10,10))
ggsave("MS/plots/supl/supl_FigS19_temp_SIP_rarefaction-curves_30deg_max20000.pdf", plot = plot30_2,
       width = 16.5, height = 8, units = "cm")
ggsave("MS/plots/supl/supl_FigS19_temp_SIP_rarefaction-curves_30deg_max20000.png", plot = plot30_2,
       width = 16.5, height = 8, units = "cm")
ggsave("MS/plots/supl/supl_FigS19_temp_SIP_rarefaction-curves_30deg_max20000.emf", plot = plot30_2,
       width = 16.5, height = 8, units = "cm",
       device = {function(filename, ...) devEMF::emf(file = filename, ...)})

