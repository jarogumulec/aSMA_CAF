

library(readxl)
library(tidyr)
library(purrr) 
library(stringr)
library(reshape2)
library(ggplot2)
library(pheatmap)
library(Hmisc) # stat p vals for ccorr
library(grid) # hvezdicky do korelacni matice
library(dplyr)  #this is last because summarise woud be masked grom Hmisc
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

write.csv(merged_data, "merged_data_for_correlation.csv", row.names = FALSE)



# correlation table ---------

# Select only numeric columns from merged_data
numeric_data <- merged_data %>%
  select(where(is.numeric))

# Calculate the correlation matrix
correlation_matrix <- cor(numeric_data, use = "pairwise.complete.obs")

# View the correlation matrix
print(correlation_matrix)


# correl matrix --------


## V1 all data shown even insignif-------------


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




## v2 show R in signif, spravny clustering -------

# Modified cor_pvalues function to handle columns with many NAs
# Adjusted function to handle NA values more robustly
cor_pvalues <- function(mat) {
  n <- ncol(mat)
  p_matrix <- matrix(NA, n, n)
  colnames(p_matrix) <- colnames(mat)
  rownames(p_matrix) <- colnames(mat)
  
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      # Extract columns and filter out NA values
      x <- mat[[i]]
      y <- mat[[j]]
      valid_idx <- complete.cases(x, y)
      
      # Use only complete cases for x and y
      x <- x[valid_idx]
      y <- y[valid_idx]
      
      # Check if there are enough observations to perform the test
      if (length(x) > 2) {  # Minimum of 3 data points for a correlation test
        test <- cor.test(x, y)
        p_matrix[i, j] <- test$p.value
        p_matrix[j, i] <- test$p.value
      }
    }
  }
  return(p_matrix)
}

# Calculate the p-values matrix
p_matrix <- cor_pvalues(numeric_data)

# Mask non-significant correlations in the correlation matrix
alpha <- 0.05
masked_correlation_matrix <- correlation_matrix
masked_correlation_matrix[p_matrix > alpha] <- NA  # Mask non-significant correlations

# Plot the heatmap with masked non-significant correlations
pheatmap(correlation_matrix,
         display_numbers = ifelse(!is.na(masked_correlation_matrix), round(masked_correlation_matrix, 2), ""),  # Display significant values only
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         color = color_palette,
         breaks = seq(-1, 1, length.out = 51),  # Ensures that 0 is centered in the color scale
         main = "Clustered Correlation Matrix (Significant Correlations Only)",
         fontsize_row = 6,
         fontsize_col = 6)



## v3 all + * for signif-----------

# Adjust the masked correlation matrix to show * for significant values, excluding the diagonal
display_matrix <- ifelse(!is.na(masked_correlation_matrix), "•", "")
diag(display_matrix) <- ""  # Remove asterisk or value for the diagonal (self-correlations)

# Plotting the square correlation matrix
p <- pheatmap(correlation_matrix,
              display_numbers = display_matrix,  # Show only * for significant correlations
              clustering_distance_rows = "euclidean",
              clustering_distance_cols = "euclidean",
              clustering_method = "complete",
              color = color_palette,
              breaks = seq(-1, 1, length.out = 51),  # Ensures that 0 is centered in the color scale
              main = "Clustered Correlation Matrix (Significant Correlations Only)",
              fontsize_row = 6,
              fontsize_col = 6)

# Use grid to enforce the aspect ratio to 1:1
grid.newpage()
pushViewport(viewport(width = unit(1, "snpc"), height = unit(1, "snpc"))) # square normalized device coordinates
print(p, newpage = FALSE)
popViewport()



## V4 just 1 half of corrmatrix + * -----------
# this is handled by corrplot which is cool, 
# + it but problem is


library(corrplot)

# Define the significance level
alpha <- 0.01

# Plot with hierarchical clustering, showing all correlations but only marking significant ones with an asterisk

corrplot(correlation_matrix, 
         method = "color", 
         type = "upper",                   # Display only the upper triangle
         order = "hclust",                 # Perform hierarchical clustering
         hclust.method = "complete",       # Specify the clustering method
         p.mat = p_matrix,                 # Pass in the p-value matrix
         sig.level = alpha,                # Set significance level
         insig = "label_sig",              # Show asterisks for significant correlations only
         pch = "*",                        # Use asterisk symbol
         pch.cex = 1.2,                    # Adjust size of asterisk
         pch.col = "black",                # Asterisk color
         addCoef.col = NULL,               # Remove correlation coefficient text
         col = colorRampPalette(c("blue", "white", "red"))(200),  # Color scale centered at 0
         tl.col = "black",                 # Text label color
         tl.srt = 45,                      # Rotate text labels for readability
         tl.cex = 0.63)                    # Set text size for labels (approximately 6.3)



## v4 signif only wrng clustering  -------------
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





