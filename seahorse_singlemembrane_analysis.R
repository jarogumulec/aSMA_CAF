# Seahorse analyser

#dependencies

library(readxl)
library(dplyr) 
library(tidyr) #na wider table


setwd("D:/OneDrive - MUNI/Experimenty/24-01-23 - Honza aSMA CAF/23-08-01 - Seahorse HGF CAF")



normalized <- read_excel("20230620_Honza_HGf_CAF_filtr_dle_signalu_pozadi_ex.xlsx", 
                                                                  sheet = "Normalized Rate")

unselected<- c("A02", "A06", "B04", "C02", "C05")


normalized.filtered <- normalized[!normalized$Well %in% unselected
                                  & normalized$Group != "Background"
                                  , ]

# Define the groups based on Measurement
normalized.filtered <- normalized.filtered %>%
  mutate(injection = case_when(
    Measurement %in% c(1, 2, 3) ~ "basal",
    Measurement %in% c(4, 5, 6) ~ "postOM",
    Measurement %in% c(7, 8, 9) ~ "post2Rot_AA",
    TRUE ~ NA_character_
  )) %>%
  group_by(Group, Well, injection)


normalised.filtered.means <- normalized.filtered %>%
  summarize(
    Mean_OCR = mean(OCR, na.rm = TRUE),
    Mean_ECAR = mean(ECAR, na.rm = TRUE),
    Mean_PER = mean(PER, na.rm = TRUE)
  )



#otocit tabulku aby byla jamka na radek - tidyr
normalised.filtered.means_wide <- normalised.filtered.means %>%
  pivot_wider(names_from = injection, values_from = c(Mean_OCR, Mean_ECAR, Mean_PER), names_sep = "_")

colnames(normalised.filtered.means_wide)


calculated_table <- normalised.filtered.means_wide %>%
  mutate(
    calc.basal.OCR = Mean_OCR_basal - Mean_OCR_post2Rot_AA,
    calc.leak.OCR = Mean_OCR_postOM - Mean_OCR_post2Rot_AA,
    post_OM.ECAR = Mean_ECAR_postOM - Mean_ECAR_basal
  )



