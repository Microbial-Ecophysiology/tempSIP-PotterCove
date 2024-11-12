# estimated iron and sulfate reduction rates by calculating difference 
  # between start and end point concentration

## load packages ####
library(tidyverse)
library(devEMF)
library(patchwork)

# import data ####
setwd("WORKINGDIRECTORY")

## project name to name files
prj <- "PROJECT"

## raw data
Fe <- read_tsv(paste0(prj, "_disFe_indiv..txt"))
SO4 <- read_tsv(paste0(prj, "_sulfate_indiv..txt"))
## label translation
mdata <- read_tsv(paste0(prj, "_id_translation_geochemistry.txt"))


## calculations ####
### Fe ####
Fe_calc <- Fe %>% 
  # join with metadata
  left_join(mdata, by = c("Treatment" = "short_label")) %>% 
  rename(TreatmentID = Treatment,
         Treatment = Treatment.y) %>% 
  # calculate rate by difference per treatment between time points
  group_by(TreatmentID) %>% 
  mutate(t1 = min(day), t2 = max(day)) %>% 
  mutate(measured_timepoint = case_when(day == t1 ~ "t1.m",
                                         day == t2 ~ "t2.m"),
         incubation_time = t2 - t1) %>% 
  ## spread data for this
  select(-day) %>% 
  pivot_wider(names_from = measured_timepoint, values_from = Fe_uM) %>% 
  ## calculate difference
  mutate(acc.Fe = t2.m - t1.m,
         Fe.rate = acc.Fe/incubation_time)


## calculate mean and sd
Fe_calc_m <- Fe_calc %>% 
  group_by(Treatment, Temp) %>% 
  summarise(Iron_red_rates_mean = mean(Fe.rate),
            Iron_red_rates_SD = sd(Fe.rate))

### Sulfate ####
SO4_calc <- SO4 %>% 
  # join with metadata
  left_join(mdata, by = c("Treatment" = "short_label")) %>% 
  rename(TreatmentID = Treatment,
         Treatment = Treatment.y) %>% 
  # calculate rate by difference per treatment between timepoints
  group_by(TreatmentID) %>% 
  mutate(t1 = min(day), t2 = max(day)) %>% 
  mutate(SO4_uM = SO4_mM*1000,
         measured_timepoint = case_when(day == t1 ~ "t1.m",
                                        day == t2 ~ "t2.m"),
         incubation_time = t2 - t1) %>% 
  ## spread data for this
  select(-day, -SO4_mM) %>% 
  pivot_wider(names_from = measured_timepoint, values_from = SO4_uM) %>% 
  ## calculate difference
  mutate(cons.SO4 = t1.m - t2.m,
         SO4.rate = cons.SO4/incubation_time)

## calculate mean and sd
SO4_calc_m <- SO4_calc %>% 
  group_by(Treatment, Temp) %>% 
  summarise(SRR_mean = mean(SO4.rate),
            SRR_SD = sd(SO4.rate))



## plot sulfate reduction rates ####
### save data for plot
write_tsv(SO4_calc_m, "MS/data_plots/data_supl_FigS03a_SRR.txt")

### plot
treat_label <- c("Acetate", "Acetate + lepidocrocite", "Acetate + sulfate")
treat_col <- c("#A5DB36FF", "#F8850FFF", "#39568CFF")

SRR_plot <- ggplot(SO4_calc_m, aes(x = Temp, y = SRR_mean, color = Treatment, 
                                   shape = Treatment)) +
  geom_hline(yintercept = 0, color = "grey", linewidth = 0.3) +
  geom_point(size = 1.6) +
  geom_path(aes(linetype = Treatment), linewidth = 0.6) +
  geom_errorbar(aes(ymin = SRR_mean -SRR_SD, ymax = SRR_mean + SRR_SD), 
                width = 0.6, linewidth = 0.5) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), breaks = c(0, 25, 50)) +
  scale_color_manual(values = treat_col, labels = treat_label) +
  scale_shape_discrete(labels = treat_label) +
  scale_linetype_discrete(labels = treat_label) +
  labs(x = "Temperature (°C)", y = "Net sulfate consumption (µM/day)") +
  coord_cartesian(ylim = c(-5,51), xlim = c(0,31)) +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        text = element_text(size = 8, color = "black"), 
        axis.title = element_text(size = 8), axis.text = element_text(size = 7, color = "black"),
        legend.text = element_text(size = 8), legend.title = element_text(face = "bold"),
        legend.position = "bottom", 
        legend.margin = margin(0,0,0,0))


## plot iron reduction rates ####
### save data for plot
write_tsv(Fe_calc_m, "MS/data_plots/data_supl_FigS03b_IRR.txt")

### plot
IRR_plot <- ggplot(Fe_calc_m, aes(x = Temp, y = Iron_red_rates_mean, 
                           color = Treatment, shape = Treatment)) +
  geom_hline(yintercept = 0, color = "grey", linewidth = 0.3) +
  geom_point(size = 1.6) +
  geom_path(aes(linetype = Treatment), linewidth = 0.6) +
  geom_errorbar(aes(ymin = Iron_red_rates_mean - Iron_red_rates_SD, 
                    ymax = Iron_red_rates_mean + Iron_red_rates_SD),
                width = 0.6, linewidth = 0.5) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), breaks = c(0, 5, 10)) +
  scale_color_manual(values = treat_col, labels = treat_label) +
  scale_shape_discrete(labels = treat_label) +
  scale_linetype_discrete(labels = treat_label) +
  labs(x = "Temperature (°C)", y = expression("Dissolved Fe"^{"2+"}*"accumulation (µM/day)")) +
  coord_cartesian(ylim = c(-1,10.5), xlim = c(0,31)) +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        text = element_text(size = 8, color = "black"), 
        axis.title = element_text(size = 8), axis.text = element_text(size = 7, color = "black"),
        legend.text = element_text(size = 8), legend.title = element_text(face = "bold"),
        legend.position = "bottom", 
        legend.margin = margin(0,0,0,0))
IRR_plot


## combine plots ####
comb <- (SRR_plot|IRR_plot)/guide_area() +
  plot_annotation(tag_levels = 'a') + plot_layout(guides = 'collect', heights = c(9,1)) 

ggsave("MS/plots/supl/supl_FigS03_temp_SIP_SRR-IRR.pdf", plot = comb, 
       height = 8, width = 18, units = "cm", useDingbats = F)
ggsave("MS/plots/supl/supl_FigS03_temp_SIP_SRR-IRR.png", plot = comb, 
       height = 8, width = 18, units = "cm")
ggsave("MS/plots/supl/supl_FigS03_temp_SIP_SRR-IRR.emf", plot = comb, 
       height = 8, width = 18, units = "cm",
       device = {function(filename, ...) devEMF::emf(file = filename, ...)})
