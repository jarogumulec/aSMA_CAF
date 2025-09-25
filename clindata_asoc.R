# clinicaldata


# dependencies-------------
library(googlesheets4)
library(dplyr)

# Authenticate with Google
gs4_auth()

#  tim autentifikuji tidyverse pro pristup do mych dokumentu!


# 0 data load -------------
sheet_url <- "https://docs.google.com/spreadsheets/d/10m3U4kt0p-Ru3_-IJJdfAAiGcGOgfCxFrod-U2nlTcY/edit#gid=620669745"
clindata <- read_sheet(sheet_url)
rm(sheet_url)


# list of patients which Mary analysed on mRNA level
# here potentially exclude those with duplex onco and others


# load list from correlation table
setwd("c:/Users/Jaro-work/OneDrive - MUNI/Experimenty/24-01-23 - Honza aSMA CAF")
merged_data <- read.csv("merged_data_for_correlation.csv")
CAF_patients <- levels(as.factor(merged_data$patient))


# Convert CAF labels to M labels and remove leading zeros

patientlist <- gsub("CAF 0*", "M", CAF_patients[CAF_patients != "hGF"])



## 0.1 preprocessing table --------------

# rename columns
{
colnames(clindata)[1] <- "id"  # rename some columns
colnames(clindata)[2] <- "rodnecislo.cast"
colnames(clindata)[3] <- "age.at.dg"
colnames(clindata)[4] <- "gender"
colnames(clindata)[5] <- "date.surgery"
colnames(clindata)[6] <- "sample.primary.site" # move to columns which are just controlones
colnames(clindata)[7] <- "sample.secondary.site" # move to columns which are just controlones
colnames(clindata)[8] <- "location"

colnames(clindata)[15] <- "therapy.strategy"
colnames(clindata)[16] <- "therapy.adjuv.radio"
colnames(clindata)[17] <- "therapy.adjuv.chemo"
colnames(clindata)[18] <- "therapy.finished"
colnames(clindata)[19] <- "therapy.success"

colnames(clindata)[20] <- "recurrence"
colnames(clindata)[21] <- "date.recurrence"
colnames(clindata)[22] <- "exitus"
colnames(clindata)[23] <- "exitus.cancerrelated"
colnames(clindata)[24] <- "date.exitus"
colnames(clindata)[25] <- "date.lastcontrol"

colnames(clindata)[27] <- "duplex.location"
colnames(clindata)[28] <- "smoker"

colnames(clindata)[29] <- "serum" # move to columns which are just controlones

colnames(clindata)[30] <- "date.serumsample" # rather to be used intead of surgery for survival?
colnames(clindata)[32] <- "RFS.days"
colnames(clindata)[33] <- "OS.days"
}
# correcting mixed type columns
# this is creepy correction but didnt found better way co correct it:

cols_to_convert <- c("pT", "pN", "pM", "G", "therapy.adjuv.radio", "duplex", "Comment", "smoker")

for (col in cols_to_convert) {
  clindata[[col]] <- as.character(clindata[[col]])
}


# filtering only those which Are done on correlation table
clindata.subset <- clindata[clindata$id %in% patientlist, ]


{
  rm(cols_to_convert) # removing temp vals
  rm(col) # removing temp vals
  rm(clindata)
  rm(patientlist)
  rm(CAF_patients)

}




# Convert "M" notation back to "CAF" notation with leading zeros for 2-digit numbers
# back to be consistent with merged_data
clindata.subset$id <- gsub("^M(\\d{2})$", "CAF 0\\1", clindata.subset$id) # For two-digit numbers
clindata.subset$id <- gsub("^M(\\d{3})$", "CAF \\1", clindata.subset$id)  # For three-digit numbers

# Rename the column to "patient"
colnames(clindata.subset)[colnames(clindata.subset) == "id"] <- "patient"






# Replace NULL with NA in the entire data frame
clindata.subset$duplex[clindata.subset$duplex == "NULL"] <- 0


## 0.2 checking table for integrity ------
checktable <- data.frame(id = clindata.subset$patient,
                         days_difference = as.numeric(difftime(clindata.subset$date.serumsample, clindata.subset$date.surgery, units = "days")),
                         clindata.subset$date.serumsample,
                         clindata.subset$date.surgery,
                         clindata.subset$sample.primary.site,
                         clindata.subset$sample.secondary.site,
                         clindata.subset$duplex,
                         clindata.subset$therapy.strategy
                         )



## 0.3 final filterout cancer duplex -------------
# not necesarry now

# clindata.subset <- clindata.subset[!(clindata.subset$id %in% c("M134", "M99")), ]



# 1 clinical characteristics screen and simplification -------------


# colnames(clindata.subset) # just for 
# now based on previous, simplify

table(clindata.subset$p16)
table(clindata.subset$location)

#clindata.subset$location.simple <- ifelse(clindata.subset$location %in% "oro", "oro", "other")
#table(clindata.subset$location.simple)
table(clindata.subset$location)


table(clindata.subset$gender)
table(clindata.subset$pT)

clindata.subset <- clindata.subset %>%
  mutate(pT.simple = case_when(
    pT %in% c("1", "2") ~ "1-2",
    pT %in% c("3", "4a") ~ "3-4",
    TRUE ~ NA_character_
  ))

table(clindata.subset$pT.simple)


table(clindata.subset$pM)
table(clindata.subset$pN)

# originally better 0-1 vs 2-3 but not optimal counts got
# clindata.subset <- clindata.subset %>%
#   mutate(pN.simple = case_when(
#     pN %in% c("0", "1") ~ "0-1",
#     pN %in% c("2", "3b") ~ "2-3",
#     TRUE ~ NA_character_
#   ))

clindata.subset <- clindata.subset %>%
  mutate(pN.simple = case_when(
    pN %in% c("0") ~ "0",
    pN %in% c("1", "2", "2a", "2b", "3b") ~ "1+",
    TRUE ~ NA_character_
  ))



#age
median(clindata.subset.selected$age.at.dg)
quantile(clindata.subset.selected$age.at.dg, 0.25, na.rm = TRUE)
quantile(clindata.subset.selected$age.at.dg, 0.75, na.rm = TRUE)



table(clindata.subset$pN.simple)

table(clindata.subset$therapy.strategy)
table(clindata.subset$smoker)
table(clindata.subset$Stage)

clindata.subset <- clindata.subset %>%
  mutate(Stage.simple = case_when(
    Stage %in% c("I", "II") ~ "I-II",
    Stage %in% c("III", "IVA", "IVB") ~ "III-IV",
    TRUE ~ NA_character_
  ))

table(clindata.subset$Stage.simple)


table(clindata.subset$recurrence)
table(clindata.subset$exitus)
table(clindata.subset$exitus.cancerrelated)

# numbers of patients with particular characteristic
sum(clindata.subset$therapy.adjuv.radio == 1) # pac s adjuv radio
sum(clindata.subset$therapy.adjuv.radio == 0 & clindata.subset$therapy.adjuv.chemo == 0) # bez adjuvance







# 2 filterout complicated vals not used for clustering -------

# create table where just those for statistics are used


colnames(clindata.subset)


# Select specified columns to create a new table
selected_columns <- c("patient", "age.at.dg", "gender", "location", "pM", "G", 
                      "p16", "therapy.adjuv.radio", "therapy.adjuv.chemo", 
                      "therapy.finished", "therapy.success", "recurrence", 
                      "smoker", "pT.simple", "pN.simple", "Stage.simple")

# Create the new table with selected columns
clindata.subset.selected <- clindata.subset %>%
  select(all_of(selected_columns))

write.csv(clindata.subset.selected, "clindata.csv", row.names = FALSE)

# selecting complete columns (not therapy chemo which has some NAs and will not be used in correlations)
selected_columns <- c("patient", "age.at.dg", "gender", "location", "pM", "G", 
                      "p16", "smoker", "pT.simple", "pN.simple", "Stage.simple")


clindata.subset.selected <- clindata.subset %>%
  select(all_of(selected_columns))
write.csv(clindata.subset.selected, "clindata_complete.csv", row.names = FALSE)

rm(checktable, selected_columns)





# 3 Analysis w measured data -  categorical factors  ----------

# Assuming merged_data and clindata_subset_selected are your two tables

merged_clin_data <- merge(merged_data, clindata.subset.selected, by = "patient")



rm(clindata.subset, merged_data, clindata.subset.selected)
write.csv(merged_clin_data, "merged_clindata_complete.csv", row.names = FALSE)


# MAY START JUST HERE - LOADING PROCESSED -------
merged_clin_data <- read.csv("merged_clindata_complete.csv")

## 3.1 counting sigificances -----


numeric_vars <- c("Young.modulus", "GCSF", "GM.CSF", "GRO.a", "GRO.a.b.y.", "IFN.y", "IL.10", "IL.13", "IL.15", 
                  "IL.1a", "IL.2", "IL.3", "IL.5", "IL.6", "IL.7", "IL.8", "MCP.1", "MCP.2", "MCP.3", "MIG", 
                  "RANTES", "TGF.b1", "TNF.a", "TNF.b", "basal.OCR", "leak.OCR", "post.OM.ECAR", "basal.ECAR", 
                  "OCR.to.ECAR", "ATP.linked..OCR", "PDPN", "ASMA")

categorical_vars <- c("location", "p16", "smoker", "pT.simple", "pN.simple", "Stage.simple")

#"gender", "pM", "G"

# alternativelly, for ASMA Status only

numeric_vars <- c("Young.modulus", "GCSF", "GM.CSF", "GRO.a", "GRO.a.b.y.", "IFN.y", "IL.10", "IL.13", "IL.15", 
                  "IL.1a", "IL.2", "IL.3", "IL.5", "IL.6", "IL.7", "IL.8", "MCP.1", "MCP.2", "MCP.3", "MIG", 
                  "RANTES", "TGF.b1", "TNF.a", "TNF.b", "basal.OCR", "leak.OCR", "post.OM.ECAR", "basal.ECAR", 
                  "OCR.to.ECAR", "ATP.linked..OCR", "PDPN")

categorical_vars <- c("location", "p16", "smoker", "pT.simple", "pN.simple", "Stage.simple", "ASMA_2group")

# excl missing 167 for ASMA Status - just for this 
merged_clin_data <- merged_clin_data %>%
  filter(patient != "CAF 167")




# Initialize a data frame to store the test results
results <- data.frame(
  Numeric_Variable = character(),
  Categorical_Variable = character(),
  t_statistic = numeric(),
  p_value = numeric(),
  N = numeric(),
  stringsAsFactors = FALSE
)

# Loop over each combination of numeric and categorical variable
for (num_var in numeric_vars) {
  for (cat_var in categorical_vars) {
    # Subset data for complete cases of the current pair of variables
    subset_data <- merged_clin_data %>% select(all_of(c(num_var, cat_var))) %>% na.omit()
    
    if (nrow(subset_data) > 2 && length(unique(subset_data[[cat_var]])) == 2) {  # Minimum 2 observations per group for t-test
      # Perform t-test
      t_test_result <- t.test(as.formula(paste(num_var, "~", cat_var)), data = subset_data)
      
      # Add the result with sample size (N)
      results <- rbind(results, data.frame(
        Numeric_Variable = num_var,
        Categorical_Variable = cat_var,
        t_statistic = t_test_result$statistic,
        p_value = t_test_result$p.value,
        N = nrow(subset_data)
      ))
    }
  }
}

# Adjust p-values using BH correction
results$p_adj <- p.adjust(results$p_value, method = "BH")

# Filter significant results (e.g., p_adj < 0.05)
#significant_results <- subset(results, p_adj < 0.05) # tvl, nic nevychazi signif po adjustment
significant_results <- subset(results, p_value < 0.05)

# Display significant results
print(significant_results)




## 3.2 heatmap of signif only -----------

# Re-cast significant results into wide format, ensuring correct names
heatmap_data <- dcast(significant_results, Numeric_Variable ~ Categorical_Variable, value.var = "p_adj")

# Melt the data back to long format for ggplot
melted_data <- melt(heatmap_data, id.vars = "Numeric_Variable", variable.name = "Categorical_Variable", value.name = "Adj_P_Value")

# Plot the heatmap
ggplot(melted_data, aes(x = Categorical_Variable, y = Numeric_Variable, fill = Adj_P_Value)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "blue", high = "red", na.value = "grey50", name = "Adj. p-value") +
  theme_minimal() +
  labs(title = "Significant Associations between Clinical and Measured Variables",
       x = "Clinical Variables", y = "Measured Variables") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

rm(melted_data, heatmap_data, significant_results, cat_var, num_var, subset_data)

## 3.3. effect sizes ------

library(ggplot2)
library(dplyr)


# Calculate Cohen's d and approximate CI for Cohen's d
# !!! pozor, davam -t_statistic / sqrt(N), to minus jsem tam dopsal, protoze to ukazovalo naopak,
# !!! asi to pocitalo 0 vs 1. ja ale vzdy chci 1 vs 0.

results <- results %>%
  mutate(
    # Calculate Cohen's d
    Effect_Size = -t_statistic / sqrt(N),  # !!! davam minus t_statistic - predtim bylo plus - to to snad resi
    # Estimate confidence intervals for Cohen's d (approximate)
    CI_low = Effect_Size - 1.96 * (1 / sqrt(N)),
    CI_high = Effect_Size + 1.96 * (1 / sqrt(N)),
    # Determine significance based on p-value
    Significance = ifelse(p_value < 0.05, "Significant", "Not Significant")
  )


forest_data <- results %>%
  arrange(Categorical_Variable, Effect_Size)

# Plot the forest plot



##3.4 forrest plot identical order  --------
# directions are now corrected by "-t_statistic"

comparison_labels <- c(
  "location" = "location (oro vs lar)",
  "p16" = "p16 (1 vs 0)",
  "smoker" = "smoker (1 vs 0)",
  "pT.simple" = "pT.simple (3-4 vs 1-2)",
  "pN.simple" = "pN.simple (1+ vs 0)",
  "Stage.simple" = "Stage.simple (III-IV vs I-II)",
  "ASMA_2group" = "ASMA status"
)

# Plot with custom facet labels
ggplot(forest_data, aes(x = reorder(Numeric_Variable, Effect_Size), y = Effect_Size, color = Significance)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "grey")) +
  facet_wrap(~ Categorical_Variable, scales = "free_y", nrow = 1, labeller = labeller(Categorical_Variable = comparison_labels)) +
  coord_flip() +
  theme_minimal() +
  labs(title = "Forest Plot of Effect Sizes by Clinical Variable",
       x = "Measured Variables", y = "Effect Size (Cohen's d)") +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        strip.text = element_text(size = 10),
        legend.position = "top")


## 3.5 Forrest plot ordered lo-hi------------

# Define reorder_within for facet-specific ordering
reorder_within <- function(x, by, within, fun = mean, sep = "___", ...) {
  new_x <- paste(x, within, sep = sep)
  stats::reorder(new_x, by, FUN = fun, ...)
}

# Define scale_x_reordered to handle the custom facet ordering
scale_x_reordered <- function(..., sep = "___") {
  ggplot2::scale_x_discrete(labels = function(x) gsub(paste0(sep, ".*$"), "", x), ...)
}


# Plot with unique sorting for each facet based on effect size
ggplot(forest_data, aes(x = reorder_within(Numeric_Variable, Effect_Size, Categorical_Variable), 
                        y = Effect_Size, color = Significance)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "grey")) +
  facet_wrap(~ Categorical_Variable, scales = "free_y", nrow = 1, 
             labeller = labeller(Categorical_Variable = comparison_labels)) +
  coord_flip() +
  theme_minimal() +
  labs(title = "Forest Plot of Effect Sizes by Clinical Variable",
       x = "Measured Variables", y = "Effect Size (Cohen's d)") +
  theme(
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    strip.text = element_text(size = 10),
    legend.position = "top"
  ) +
  scale_x_reordered()  # Ensures that each facet is sorted independently








## 3.6 check direction of effect -----
# this is a relevant checkup whether the direciton of the effect is correct


library(dplyr)

# Function to calculate and verify the means for each level of each categorical variable
check_effect_direction <- function(data, numeric_vars, cat_vars) {
  check_results <- list()
  
  for (num_var in numeric_vars) {
    for (cat_var in cat_vars) {
      # Calculate mean for each level of the categorical variable
      mean_summary <- data %>%
        group_by(!!sym(cat_var)) %>%
        summarize(mean_value = mean(!!sym(num_var), na.rm = TRUE))
      
      # Extract the mean values for each level
      levels_ordered <- mean_summary[[cat_var]]
      means <- mean_summary$mean_value
      
      # Calculate Cohen's d and determine direction
      cohen_d <- (means[1] - means[2]) / sd(data[[num_var]], na.rm = TRUE)
      
      # Save the results for this numeric-categorical variable pair
      check_results[[paste0(num_var, "_", cat_var)]] <- data.frame(
        Numeric_Variable = num_var,
        Categorical_Variable = cat_var,
        Level_1 = levels_ordered[1],
        Mean_Level_1 = means[1],
        Level_2 = levels_ordered[2],
        Mean_Level_2 = means[2],
        Computed_Cohen_d = cohen_d,
        Effect_Direction = ifelse(cohen_d > 0, paste0(levels_ordered[1], " > ", levels_ordered[2]), paste0(levels_ordered[2], " > ", levels_ordered[1]))
      )
    }
  }
  
  # Combine all results into one data frame
  do.call(rbind, check_results)
}

# Define your numeric and categorical variables
numeric_vars <- c("Young.modulus", "GCSF", "GM.CSF", "GRO.a", "GRO.a.b.y.", "IFN.y", "IL.10", "IL.13", "IL.15", "IL.1a", 
                  "IL.2", "IL.3", "IL.5", "IL.6", "IL.7", "IL.8", "MCP.1", "MCP.2", "MCP.3", "MIG", "RANTES", "TGF.b1", 
                  "TNF.a", "TNF.b", "basal.OCR", "leak.OCR", "post.OM.ECAR", "basal.ECAR", "OCR.to.ECAR", "ATP.linked..OCR", 
                  "PDPN", "ASMA")
cat_vars <- c("location", "p16", "smoker", "pT.simple", "pN.simple", "Stage.simple")

# Run the check
control_table <- check_effect_direction(merged_clin_data, numeric_vars, cat_vars)

# View or print the control table to verify results
print(control_table)

rm(forest_data, results, t_test_result, categorical_vars, numeric_vars, comparison_labels, reorder_within, scale_x_reordered)

# 4 Analysis w measured data -continuous  ----------

# that is age. display as forrest plot as previosuly.

# Initialize an empty data frame to store correlation results
correlation_results <- data.frame(Variable = character(), 
                                  Correlation = numeric(), 
                                  p_value = numeric(),
                                  CI_low = numeric(),
                                  CI_high = numeric(),
                                  Significance = character(),
                                  stringsAsFactors = FALSE)

# Calculate correlation for each variable against age.at.dg
for (var in numeric_vars) {
  if (var %in% names(merged_clin_data)) {
    # Pearson correlation
    cor_test <- cor.test(merged_clin_data[[var]], merged_clin_data$age.at.dg, 
                         use = "pairwise.complete.obs", method = "pearson")
    
    # Store results with 95% confidence intervals for correlation
    correlation_results <- correlation_results %>%
      add_row(Variable = var, 
              Correlation = cor_test$estimate, 
              p_value = cor_test$p.value,
              CI_low = cor_test$conf.int[1],
              CI_high = cor_test$conf.int[2],
              Significance = ifelse(cor_test$p.value < 0.05, "Significant", "Not Significant"))
  }
}

# Order by correlation values for plot
correlation_results <- correlation_results %>%
  arrange(Correlation)

# Create the forest plot
ggplot(correlation_results, aes(x = reorder(Variable, Correlation), y = Correlation, color = Significance)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.2) + # Confidence intervals as "error bars"
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") + # Reference line at 0
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "grey")) +
  coord_flip() +
  theme_minimal() +
  labs(title = "Correlation of Age with Measured Parameters",
       x = "Measured Variables", 
       y = "Correlation Coefficient (r)") +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        legend.position = "top") +
  scale_y_continuous(limits = c(-1, 1))  # Set limits for correlation coefficients
