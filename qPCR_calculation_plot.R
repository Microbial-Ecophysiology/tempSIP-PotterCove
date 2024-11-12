# code calculating copies per ng cDNA and plotting the data
 #qPCR was performed at cDNA samples, which were also sequenced
 #dsrA qPCR from all samples
 #because some samples were below detection limit, for these samples 
 #the cDNA was controlled by performing a bacterial qPCR in addition

# load packages ####
library(tidyverse)
library(cowplot)

# load data ####
setwd("WORKINGDIRECTORY")

df1 <- read_tsv("rawdata_qPCR.txt")

# Calculation ####
## needed factors
vol_used_qPCR_ul <- 2
factor_to_dalton <- 6.0221E23
size_qPCR_standard_dalton <- 1199619.49

## copies per ng cDNA
 # (qPCR quantity (g) * factor_to_dalton)/size_qPCR_standard_dalton) 
 # / amount of cDNA in qPCR (ng)

df2 <- df1 %>% 
  mutate(
    no_cDNA_dsrA = 
      ((qPCR_dsrA_Quantity_ng * factor_to_dalton) / size_qPCR_standard_dalton) /
        (cDNA_conc*vol_used_qPCR_ul),
    no_cDNA_bac = 
      ((qPCR_bac_Quantity_ng * factor_to_dalton) / size_qPCR_standard_dalton) /
      (cDNA_conc*vol_used_qPCR_ul)
  )

# calculate mean and SD
df3 <- df2 %>%
  group_by(Treatment, pooled_Fraction, Incubation_ID, Temp) %>%
  mutate(no_cDNA_dsrA_mean = mean(no_cDNA_dsrA), no_cDNA_dsrA_sd = sd(no_cDNA_dsrA),
         no_cDNA_bac_mean = mean(no_cDNA_bac), no_cDNA_bac_sd = sd(no_cDNA_bac)) %>% 
  distinct(
    Incubation_ID, Treatment, Temp, pooled_Fraction,
    no_cDNA_dsrA_mean, no_cDNA_dsrA_sd,
    no_cDNA_bac_mean, no_cDNA_bac_sd) %>% 
  ungroup() %>% 
  mutate(pooled_Fraction = factor(pooled_Fraction, 
                                  levels = c("Ultralight", "Light", "Midpoint", "Heavy", "Ultraheavy")))


# supl MS plot ####
x <- c("Acetate", "Acetate +\nLepidocrocite", "Acetate +\nSulfate")
names(x) <- c("Acetate", "Acetate.Lepidocrocite", "Acetate.Sulfate")

## save data plot
df3 %>% 
  select(Incubation_ID, Treatment, Temp, pooled_Fraction, no_cDNA_dsrA_mean, no_cDNA_dsrA_sd,
         no_cDNA_bac_mean, no_cDNA_bac_sd) %>% 
  write_tsv("MS/data_plots/data_supl_FigS16_qPCR.txt")


## plot
plot_no_cDNA <- ggplot(df3, aes(x = pooled_Fraction, y = no_cDNA_dsrA_mean)) +
  geom_point(aes(color = no_cDNA_dsrA_mean == 0), size = 0.8) +
  geom_errorbar(aes(ymin = no_cDNA_dsrA_mean - no_cDNA_dsrA_sd, 
                    ymax = no_cDNA_dsrA_mean + no_cDNA_dsrA_sd,
                    color = no_cDNA_dsrA_mean == 0), linewidth = 0.3,
                width = 0.4) +
  scale_color_manual(values = setNames(c("red", "black"), c(T,F))) +
  geom_point(aes(x = pooled_Fraction, y = no_cDNA_bac_mean/800), color = "blue", size = 0.5) +
  geom_errorbar(aes(ymin = no_cDNA_bac_mean/800 - no_cDNA_bac_sd/800, 
                    ymax = no_cDNA_bac_mean/800 + no_cDNA_bac_sd/800), linewidth = 0.3,
                width = 0.4, color = "blue") +
  scale_y_continuous(sec.axis = sec_axis(~.*800, name = "bacterial 16S rRNA transcripts per ng cDNA")) +
  
  labs(x = "Temperature (°C)", y = expression(italic("dsrA")~"transcripts per ng cDNA")) +
  facet_grid(Treatment ~ Temp, labeller = labeller(Treatment = x)) +
  theme_bw() +
  theme(legend.position = "none", text = element_text(size = 8),
        strip.background = element_rect(color = "black", fill = "white"),
        axis.text.y = element_text(color = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, color = "black"),
        axis.title.y.right = element_text(color = "blue"),
        panel.grid = element_blank())

ggsave("MS/plots/supl/supl_FigS16_temp_SIP_dsrA_qPCR.png", 
       plot_no_cDNA, height = 8.5, width = 16.5, units = "cm")
ggsave("MS/plots/supl/supl_FigS16_temp_SIP_dsrA_qPCR.pdf", 
       plot_no_cDNA, height = 8.5, width = 16.5, units = "cm")
ggsave("MS/plots/supl/supl_FigS16_temp_SIP_dsrA_qPCR.emf", 
       plot_no_cDNA, height = 8.5, width = 16.5, units = "cm",
       device = {function(filename, ...) devEMF::emf(file = filename, ...)})
