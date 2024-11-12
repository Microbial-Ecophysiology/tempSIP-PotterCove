# plot delta 13C CO2 profiles of SIP incubations

# load packages ####
library(tidyverse)
library(devEMF)
library(patchwork)


# import data ####
setwd("WORKINGDIRECTORY")

## project name to name files
prj <- "PROJECT"

raw <- read_tsv(paste0(prj, "_CO2_isotope_data.txt"))
mdata <- read_tsv(paste0(prj, "_id_translation_geochemistry.txt"))


# calculate mean and sd ####
df <- full_join(raw, mdata, by = c("Treatment" = "short_label")) %>% 
  rename(TreatmentID = Treatment,
         Treatment = Treatment.y)

df_calc <- df %>% 
  group_by(Treatment, Temp, Isotope, day) %>% 
  summarise(mean = mean(d13C.12C.permil.vs.VPDB, na.rm = T),
            SD = sd(d13C.12C.permil.vs.VPDB, na.rm = T)) %>% 
  mutate(treat_iso = factor(case_when(Treatment == "Acetate" & Isotope == "12C" ~ "12C-acetate",
                               Treatment == "Acetate" & Isotope == "13C" ~ "13C-acetate",
                               Treatment == "Acetate.Lepidocrocite" & Isotope == "12C" ~ "12C-acetate + lepidocrocite",
                               Treatment == "Acetate.Lepidocrocite" & Isotope == "13C" ~ "13C-acetate + lepidocrocite",
                               Treatment == "Acetate.Sulfate" & Isotope == "12C" ~ "12C-acetate + sulfate",
                               Treatment == "Acetate.Sulfate" & Isotope == "13C" ~ "13C-acetate + sulfate"),
                            levels = c("13C-acetate", "12C-acetate", "13C-acetate + lepidocrocite", "12C-acetate + lepidocrocite",
                                       "13C-acetate + sulfate", "12C-acetate + sulfate")),
         Temp = factor(Temp, levels = c(2, 5, 10, 15, 20, 25, 30),
                       labels = c("2 °C", "5 °C", "10 °C", "15 °C", "20 °C", "25 °C", "30 °C")))


# plot ####
treat_label <- c( expression({}^{13}*"C-acetate"), expression({}^{12}*"C-acetate"),
                  expression({}^{13}*"C-acetate + lepidocrocite"), expression({}^{12}*"C-acetate + lepidocrocite"),
                  expression({}^{13}*"C-acetate + sulfate"), expression({}^{12}*"C-acetate + sulfate"))
treat_col <- c("#A5DB36FF","#b6ce85ff", "#F8850FFF","#da9f60ff", "#39568CFF", "#7893c8ff")

CO2_plot <- ggplot(df_calc, aes(x = day, y = mean, color = treat_iso, 
                                   shape = treat_iso)) +
  geom_point(size = 1.6, fill = "white", stroke = 1,) +
  geom_path(aes(linetype = treat_iso), linewidth = 0.6) +
  geom_errorbar(aes(ymin = mean - SD, ymax = mean + SD), 
                width = 0.6, linewidth = 0.5) +
  scale_color_manual(values = treat_col, labels = treat_label) +
  scale_shape_manual(values = c(19, 21, 17, 24, 15, 22), labels = treat_label) +
  scale_linetype_manual(values = c("solid","solid", "dashed","dashed", "longdash", "longdash"), labels = treat_label) +
  labs(x = "Incubation time (days)", y = expression(delta^{13}*"C (\u{2030} VPDB)"),
       color = "Treatment", shape = "Treatment", linetype = "Treatment") +
  facet_wrap(~Temp, nrow = 2, axes = "all") +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        text = element_text(size = 8, color = "black"), 
        axis.title = element_text(size = 8), axis.text = element_text(size = 7, color = "black"),
        strip.background = element_rect(fill = "white"), strip.text = element_text(size = 8, face = "bold"),
        legend.key.size = unit(0.3, "cm"), legend.key.spacing.y = unit(0.1, "cm"),
        legend.text = element_text(size = 7), legend.title = element_text(size = 7, face = "bold"),
        legend.position = "inside", legend.position.inside = c(1, 0),
        legend.justification = c("right", "bottom"))

ggsave("MS/plots/supl/supl_FigS01_temp_SIP_CO2_isotope.png", CO2_plot, width = 16, height = 8, units = "cm")
ggsave("MS/plots/supl/supl_FigS01_temp_SIP_CO2_isotope.pdf", CO2_plot, width = 16, height = 8, units = "cm")
ggsave("MS/plots/supl/supl_FigS01_temp_SIP_CO2_isotope.emf", CO2_plot, width = 16, height = 8, units = "cm",
       device = {function(filename, ...) devEMF::emf(file = filename, ...)})

