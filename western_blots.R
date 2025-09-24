
# WB celkova analyza


# dependent ----------

library(readxl)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(stringr)



# data load --------------

setwd("C:/Users/Jaro-work/OneDrive - MUNI/Experimenty/24-01-23 - Honza aSMA CAF/23-03-22 - Western blot ASMA PDPN CAFs")



## approach 1: from data Klara and kveta calculated -----
# this approach s unsuitable as it was normalized do membrane intensity average, not the HGFB

# # experiments 2023
# aSMA_PDPN.densitometry_230918 <- read_excel("23-09-18 - PDPN WB denzitometrie/densitometry.xlsx")
# aSMA_PDPN.densitometry_230918 <- aSMA_PDPN.densitometry_230918[c(1:11, 13:23, 25:35),c(1:4, 19:20)]
# 
# 
# 
# # experiments 2024
# aSMA_PDPN.densitometry_241007_Kveta <- read_excel("24-10-07 - Klara westerny leto 24/densitometry_Klára_léto 2024.xlsx", 
#                                            sheet = "Květa")
# aSMA_PDPN.densitometry_241007_Kveta <- aSMA_PDPN.densitometry_241007_Kveta[c(1:8, 11:18, 21:28),c(1:4, 19:20)]
# 
# 
# aSMA_PDPN.densitometry_241007_Klara <- read_excel("24-10-07 - Klara westerny leto 24/densitometry_Klára_léto 2024.xlsx", 
#                                                   sheet = "Klára")
# aSMA_PDPN.densitometry_241007_Klara <- aSMA_PDPN.densitometry_241007_Klara[c(1:6, 11:16, 21:26),c(1:4, 19:20)]


## approach 2: to HGFB of every membrane - non-normalised to housekeep -----------


# # downside: The HGFB will be 1.0 without variability at all. statistic needed would be One sample test
# 
# densitometry_HGFBnormalised <- read_excel("densitometry_vse_dohromady.xlsx") 
# 
# # removes empty cases
# densitometry_HGFBnormalised <- densitometry_HGFBnormalised[densitometry_HGFBnormalised$id != "" & !is.na(densitometry_HGFBnormalised$id), ]
# #removes unrelevant columns
# #densitometry_HGFBnormalised <- densitometry_HGFBnormalised[,c(1:3, 5, 22:24)] # pre 9/2025
# densitometry_HGFBnormalised <- densitometry_HGFBnormalised[,c(1:3, 5, 33:35)]
# 
# 
# colnames(densitometry_HGFBnormalised) <- c("id", "run", "replicate", "patient", 
#                         "PDPN.normalised.int.raw.vals", 
#                         "ASMA.normalised.int.raw.vals", 
#                         "TUB.normalised.int.raw.vals")
# 
# colnames(densitometry_HGFBnormalised)


## approach 3: to HGFB of every membrane - normalised to housekeep (actb+tub avg) 25-09-24 -----------


# downside: The HGFB will be 1.0 without variability at all. statistic needed would be One sample test
# but now allready with normalisation to housekeep

densitometry_HGFBnormalised <- read_excel("densitometry_vse_dohromady.xlsx", 
                                         sheet = "dodelane_25-09-18")



# removes empty cases
densitometry_HGFBnormalised <- densitometry_HGFBnormalised[densitometry_HGFBnormalised$id != "" & !is.na(densitometry_HGFBnormalised$id), ]
#removes unrelevant columns
densitometry_HGFBnormalised <- densitometry_HGFBnormalised[,c(1:3, 5, 40:41)]


colnames(densitometry_HGFBnormalised) <- c("id", "run", "replicate", "patient", 
                                           "ASMA.normalised.int.raw.vals",
                                           "PDPN.normalised.int.raw.vals")

colnames(densitometry_HGFBnormalised)













# rename patients for consistency
densitometry_HGFBnormalised <- densitometry_HGFBnormalised %>%
  mutate(
    patient = case_when(
      # Rename 3-digit patients to 'CAF' prefix
      str_detect(patient, "^\\d{3}$") ~ paste0("CAF ", patient),
      # Rename 2-digit patients to 'CAF 0XX' format
      str_detect(patient, "^\\d{2}$") ~ paste0("CAF 0", patient),
      # Keep HGFB unchanged
      patient == "HGFB" ~ "hGF",
      # Default to original value (just in case)
      TRUE ~ patient
    )
  )

write.csv(densitometry_HGFBnormalised, "../densitometry_HGFBnormalised_housekeepnormalised_25-09.csv", row.names = FALSE)
#


# plot aSMA--------------------


## stats to calculate first because of one sample adjujsted test ggpubr cannot do alone

# Step 1: Perform one-sample t-tests manually for each patient group, excluding the constant group (HGFB)
test_results <- densitometry_HGFBnormalised %>%
  group_by(patient) %>%
  summarize(
    p.val = if(all(ASMA.normalised.int.raw.vals == 1)) 1 else t.test(ASMA.normalised.int.raw.vals, mu = 1)$p.value
  )

# Step 2: Apply Benjamini-Hochberg (BH) correction (optional step to create p.adj table)
test_results <- test_results %>%
  mutate(p.adj = p.adjust(p.val, method = "BH"))

# Step 3: Merge the p-values with the original data for plotting
densitometry_HGFBnormalised <- densitometry_HGFBnormalised %>%
  left_join(test_results, by = "patient")

# Step 4: Remove rows with NA values (in case there are any)
densitometry_HGFBnormalised <- densitometry_HGFBnormalised %>%
  filter(!is.na(ASMA.normalised.int.raw.vals))

# Step 5: Create a distinct dataset for patient-level p-values
p_value_data <- test_results %>%
  filter(!is.na(p.val)) %>%
  distinct(patient, p.val)

# Step 6: Define a y-limit for the text to be placed slightly above the maximum values in each group
y_max <- max(densitometry_HGFBnormalised$ASMA.normalised.int.raw.vals, na.rm = TRUE) * 1.1



p_value_data <- p_value_data %>%
  mutate(sig_label = case_when(
    p.val < 0.0001 ~ "****",
    p.val < 0.001  ~ "***",
    p.val < 0.01   ~ "**",
    p.val < 0.05   ~ "*",
    TRUE           ~ ""
  ))

# 
# Step 7: Create the plot
ggplot(densitometry_HGFBnormalised, aes(x = patient, y = ASMA.normalised.int.raw.vals)) +
  geom_point(position = position_jitter(seed = 1, width = 0.3), size = 0.05, aes(colour = factor(patient))) +
  geom_boxplot(outlier.shape = NA, fill = "0", lwd = 0.5) +
  theme_bw() +
  theme(
    axis.text.x = element_text(color = "black", size = 6.3, angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(color = "black", size = 6.3),
    axis.title  = element_text(size = 8)
  ) +
  xlab(NULL) +
  ylab("aSMA (A.U.)") +
  theme(legend.position = "none") +
  scale_y_continuous(limits = c(0, NA)) +
  scale_y_log10() +
  geom_text(
    data = p_value_data,
    aes(x = patient, y = y_max, label = sig_label),
    vjust = -1,
    size = 2
  )



ggsave("asma.svg", plot = last_plot(),
       width = 7.5, height = 3.5, units = "cm", dpi = 300, scale = 1, limitsize = TRUE)


# Optional: Check the p-value table
print(test_results)
write.csv(test_results, "../densitometry_asma_test_results_housekeepnormalised.csv", row.names = FALSE)


# pplot PDPN -----

# Step 1: Perform one-sample t-tests manually for each patient group, excluding the constant group (hGF)
test_results_PDPN <- densitometry_HGFBnormalised %>%
  group_by(patient) %>%
  summarize(
    p.val = if(all(PDPN.normalised.int.raw.vals == 1)) 1 else t.test(PDPN.normalised.int.raw.vals, mu = 1)$p.value
  )

# Step 2: Apply Benjamini-Hochberg (BH) correction (optional step to create p.adj table)
test_results_PDPN <- test_results_PDPN %>%
  mutate(p.adj = p.adjust(p.val, method = "BH"))


# Step 3: Merge the p-values with the original data for plotting
densitometry_HGFBnormalised_PDPN <- densitometry_HGFBnormalised %>%
  left_join(test_results_PDPN, by = "patient")

# Step 4: Remove rows with NA values (in case there are any)
densitometry_HGFBnormalised_PDPN <- densitometry_HGFBnormalised_PDPN %>%
  filter(!is.na(PDPN.normalised.int.raw.vals))

# Step 5: Create a distinct dataset for patient-level p-values
p_value_data_PDPN <- test_results_PDPN %>%
  filter(!is.na(p.val)) %>%
  distinct(patient, p.val)

p_value_data_PDPN <- p_value_data_PDPN %>%
  mutate(sig_label = case_when(
    p.val < 0.0001 ~ "****",
    p.val < 0.001  ~ "***",
    p.val < 0.01   ~ "**",
    p.val < 0.05   ~ "*",
    TRUE           ~ ""
  ))



# Step 6: Define a y-limit for the text to be placed slightly above the maximum values in each group
y_max_PDPN <- max(densitometry_HGFBnormalised_PDPN$PDPN.normalised.int.raw.vals, na.rm = TRUE) * 1.1




# Step 7: Create the plot for PDPN
ggplot(densitometry_HGFBnormalised_PDPN, aes(x = patient, y = PDPN.normalised.int.raw.vals)) +
  geom_point(position = position_jitter(seed = 1, width = 0.3), size = 0.05, aes(colour = factor(patient))) +
  geom_boxplot(outlier.shape = NA, fill = "0", lwd = 0.5) +
  theme_bw() +
  theme(
    axis.text.x = element_text(color = "black", size = 6.3, angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(color = "black", size = 6.3),
    axis.title = element_text(size = 8)
  ) +
  xlab(NULL) +
  ylab("PDPN (A.U.)") +
  theme(legend.position = "none") +
  scale_y_continuous(limits = c(0, NA)) +
  scale_y_log10() +
  # Step 8: Manually annotate with the unadjusted p-values (one per patient)
  geom_text(
    data = p_value_data_PDPN,
    aes(x = patient, y = y_max_PDPN, label = sig_label),
    vjust = -1,
    size = 2
  )

# Save the plot as SVG
ggsave("pdpn.svg", plot = last_plot(),
       width = 7.5, height = 3.5, units = "cm", dpi = 300, scale = 1, limitsize = TRUE)

# Optional: Check the p-value table
print(test_results_PDPN)
write.csv(test_results_PDPN, "../densitometry_pdpn_test_results_housekeepnormalised.csv", row.names = FALSE)









