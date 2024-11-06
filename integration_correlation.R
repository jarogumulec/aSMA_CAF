

library(readxl)
library(tidyr)
library(dplyr)
library(stringr)
library(reshape2)
library(ggplot2)
library(purrr) 
library(pheatmap)
library(Hmisc) # stat p vals for ccorr
library(grid) # hvezdicky do korelacni matice

# nacteni vsech dat co byly udelany -------------

setwd("c:/Users/Jaro-work/OneDrive - MUNI/Experimenty/24-01-23 - Honza aSMA CAF")
# vse ulozeno do rootu

## cytokiny -----------
cytokine <- read.csv("normalised_cytokine.csv")
cytokine <- cytokine[, c("treatment", "cytokine", "normalised.expression.avg")]

# averaging
cytokine_avg <- cytokine %>%
  group_by(treatment, cytokine) %>%
  summarize(mean_normalised_expression = mean(normalised.expression.avg, na.rm = TRUE))

# remove Pos Neg Blank
cytokine_avg <- cytokine_avg %>%
  filter(!cytokine %in% c("Pos", "Neg", "BLANK"))

# reshaped
cytokine_avg <- dcast(cytokine_avg, treatment ~ cytokine, value.var = "mean_normalised_expression")

rm(cytokine)
# Rename Group column to patient
cytokine_avg <- cytokine_avg %>%
  rename(patient = treatment)




## seahorse -------
seahorse <- read.csv("Seahorse_complete_calculated.csv")
# only the needed ones
seahorse <- seahorse[, c("Group", "calc.basal.OCR", "calc.leak.OCR", "calc.post_OM.ECAR", "Mean_ECAR_basal", "calc.OCR_to_ecar", "calc.OCR_ATP_linked")]
# average
seahorse_avg <- seahorse %>%
  group_by(Group) %>%
  summarize(across(where(is.numeric), mean, na.rm = TRUE))

# Rename Group column to patient
seahorse_avg <- seahorse_avg %>%
  rename(patient = Group)
rm(seahorse)
colnames(seahorse_avg) <- c("patient", "basal OCR", "leak OCR", "post-OM ECAR",   "basal ECAR", "OCR to ECAR", "ATP-linked  OCR")
colnames(seahorse_avg) 

## western ----------
western <- read.csv("densitometry_HGFBnormalised.csv")
# only needed ones
western <- western[, c("patient", "PDPN.normalised.int.raw.vals", "ASMA.normalised.int.raw.vals")]
# average
western_avg <- western %>%
  group_by(patient) %>%
  summarize(across(where(is.numeric), mean, na.rm = TRUE))
rm(western)

colnames(western_avg) <- c("patient", "PDPN", "ASMA")



## AFM -------

AFM <- read.csv("23-09-18 - AFM CAF/AFM_filtered.csv") # 156BR2 filtered out. done in AFM.R
#AFM <- read_excel("23-09-18 - AFM CAF/AFM_caf_2023-2024.xlsx") # this would load unfiltered dataset





# only needed ones
AFM <- AFM[, c("patient", "mean")]
# average
AFM_avg <- AFM %>%
  group_by(patient) %>%
  summarize(across(where(is.numeric), mean, na.rm = TRUE))
rm(AFM)

## Rename Group column to patient
#AFM_avg <- AFM_avg %>%
#  rename(Youngmodulus = mean)

colnames(AFM_avg) <- c("patient", "Young modulus")

## merging ---------


tables_to_merge <- list(AFM_avg, cytokine_avg, seahorse_avg, western_avg)
# Merge all tables by the 'patient' column using full join
merged_data <- reduce(tables_to_merge, full_join, by = "patient")




# correlation table ---------

# Select only numeric columns from merged_data
numeric_data <- merged_data %>%
  select(where(is.numeric))

# Calculate the correlation matrix
correlation_matrix <- cor(numeric_data, use = "pairwise.complete.obs")

# View the correlation matrix
print(correlation_matrix)


# correl matrix --------




# Create a correlation matrix of numeric columns
numeric_data <- merged_data %>%
  select(where(is.numeric))
correlation_matrix <- cor(numeric_data, use = "pairwise.complete.obs")

# Plot the heatmap with clustering
color_palette <- colorRampPalette(c("blue", "white", "red"))(50)

# Plot the heatmap with the centered color scale at 0
pheatmap(correlation_matrix,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         color = color_palette,
         breaks = seq(-1, 1, length.out = 51),  # Ensures that 0 is centered in the color scale
         main = "Clustered Correlation Matrix",
         fontsize_row = 6,
         fontsize_col = 6)




# corr matrix signif only spravny clustering -------

#todo








# correl matrix signif only -------------
# blby ze clustruje jinak ale je to prehledny 
# nesignifikantni to totiz udela jako r=0, coz ovlivni clustering

# Step 1: Calculate correlation matrix with p-values
cor_results <- Hmisc::rcorr(as.matrix(numeric_data))  # Gives both correlation and p-values
correlation_matrix <- cor_results$r
p_value_matrix <- cor_results$P

# Step 2: Set non-significant correlations (e.g., p > 0.05) to NA
significance_level <- 0.05
correlation_matrix[p_value_matrix > significance_level] <- NA

# Step 3: Plot the heatmap with only significant correlations
color_palette <- colorRampPalette(c("blue", "white", "red"))(50)

# Replace NA values with 0 for clustering
correlation_matrix_clust <- correlation_matrix
correlation_matrix_clust[is.na(correlation_matrix_clust)] <- 0

pheatmap(correlation_matrix_clust,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         color = color_palette,
         breaks = seq(-1, 1, length.out = 51),
         main = "Significant Correlations Only",
         fontsize_row = 6,
         fontsize_col = 6,
         na_col = "grey")  # Color for non-significant correlations



