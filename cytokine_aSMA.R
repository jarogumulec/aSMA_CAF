# Cytokine Assay - new verison 
# reformatin done in R 23-05-15
# now more calculations are done in R
# now also blanks and neg. controls are used to substract also from pos controls (i think its better)
# 23-08-11 modified with triplicate membranes
# 23-09-08 - Kveta updated the BR-BR3, the previous calc from 8/2023 fas not correct (2x BR1 nebo neco takoveho)


library(readxl)
library(tidyr)
library(dplyr)
library(stringr)
library(reshape2)
library("gplots") # heatmap.2
library(ggplot2)
library(multcomp)


#library("heatmap.plus")
#library("RColorBrewer") 
library("heatmaply") # althouth i primary use gplots for heatmap.2, this for viridis and normalisation





# skript pro vypocet -----------------


# Define a function for data processing
process_replicate <- function(replicate_data) {
  # Clean column names
  colnames(replicate_data) <- gsub("\\.\\.\\..*", "", colnames(replicate_data))
  
  # Remove unneeded columns
  replicate_data <- replicate_data[, -c(2:33)]
  
  # Remove lines with empty "treatment"
  replicate_data <- replicate_data[!is.na(replicate_data$treatment), ]
  
  # Pivot to long format
  replicate_data_long <- replicate_data %>%
    pivot_longer(cols = -treatment,
                 names_to = "cytokine",
                 values_to = "normalised.expression")
  
  # Calculate average for technical replicates
  replicate_data_averaged <- replicate_data_long %>%
    group_by(cytokine, treatment) %>%
    summarize(normalised.expression.avg = mean(normalised.expression))
  
  # Filter out "BLANK", "Pos", "Neg"
  #replicate_data_averaged <- replicate_data_averaged %>%
  #  filter(cytokine != "BLANK" & treatment != "Pos" & treatment != "Neg")
  
  # Return the processed data
  return(replicate_data_averaged)
}

#


## data load - 2023/08 - not correct BR1-3 - OBSOLETE------------------
# setwd("c:/Users/admin/OneDrive - MUNI/Experimenty/20-08-11 - Autophagy FaDu Cotreatment/23-05-12 - cytokiny aSMA CAF pacienti na FaDu")
# replicate_1 <- read_excel("23-05-12 - první replika/v2_odec_blank_23-05-12 - CAF aSMA.xlsx", skip = 1)
# replicate_2 <- read_excel("23-07-26 - druha replika/cytokiny Kveta_proR.xlsx", skip = 1)
# replicate_3 <- read_excel("23-08-09 - třetí replika/treti_replika_calc.xlsx", skip = 1)

# data load - 2023/09 - BR1-3------------------
#setwd("c:/Users/admin/OneDrive - MUNI/Experimenty/20-08-11 - Autophagy FaDu Cotreatment/23-05-12 - cytokiny aSMA CAF pacienti na FaDu/23-09-06 - Kveta prehodnoceno BR1-3")

# presun do slozky CAFs - pryc z puvodni autofagie
setwd("c:/Users/Jaro-work/OneDrive - MUNI/Experimenty/24-01-23 - Honza aSMA CAF/23-05-12 - cytokiny aSMA CAF pacienti na FaDu/23-09-06 - Kveta prehodnoceno BR1-3")

replicate_1 <- read_excel("BR1/BR1.xlsx", skip = 1, sheet = "calc")
replicate_2 <- read_excel("BR2/BR2.xlsx", skip = 1, sheet = "calc")
replicate_3 <- read_excel("BR3/BR3.xlsx", skip = 1, sheet = "calc")

# 2024-10 jeste dalsi pacienti - ti ale nebyli delani v triplikatu ale jen 2X
setwd("c:/Users/Jaro-work/OneDrive - MUNI/Experimenty/24-01-23 - Honza aSMA CAF/23-05-12 - cytokiny aSMA CAF pacienti na FaDu/24-10-11 Kveta_dalsi_pacienti_jen_2BR")

group_2_replicate_1 <- read_excel("BR1 2group CAFs E13.xlsx", skip = 1, sheet = "calc")
group_2_replicate_2 <- read_excel("BR2 2group CAFs E14.xlsx", skip = 1, sheet = "calc")




# # testovaci 2023/09 var 2 - list pruimer cele desky a nikoli jen pozitivni kontroly - koncenzus toto finalne nebrat:
# replicate_1 <- read_excel("BR1/BR1.xlsx", skip = 1, sheet = "calc_avg_all")
# replicate_2 <- read_excel("BR2/BR2.xlsx", skip = 1, sheet = "calc_avg_all")
# replicate_3 <- read_excel("BR3/BR3.xlsx", skip = 1, sheet = "calc_avg_all")





# data process --------------

replicate_1_data <- process_replicate(replicate_1)
replicate_2_data <- process_replicate(replicate_2)
replicate_3_data <- process_replicate(replicate_3)
#2024 added
group_2_replicate_1_data <- process_replicate(group_2_replicate_1)
group_2_replicate_2_data <- process_replicate(group_2_replicate_2)
group_2_replicate_2_data$treatment <- as.character(group_2_replicate_2_data$treatment)



# alternative + normalizovat mezi membranami -----
# pokud nepouzito, preskocit na line # merge replicates to one file 

# 2024-10-14
# i kdyz prumer na membranu drive zatracen, je to asi to nelepsi co jde udelat aby byly membrany porovnatelne.

normalize_to_one <- function(data, column, treatment_column) {
  # Calculate the mean within each treatment
  treatment_means <- tapply(data[[column]], data[[treatment_column]], mean, na.rm = TRUE)
  
  # Normalize each value within its treatment group
  data[[column]] <- mapply(function(value, treatment) {
    value / treatment_means[treatment]
  }, data[[column]], data[[treatment_column]])
  
  return(data)
}


# Normalize each replicate so that each treatment has an average of 1.0
replicate_1_data <- normalize_to_one(replicate_1_data, "normalised.expression.avg", "treatment")
replicate_2_data <- normalize_to_one(replicate_2_data, "normalised.expression.avg", "treatment")
replicate_3_data <- normalize_to_one(replicate_3_data, "normalised.expression.avg", "treatment")
group_2_replicate_1_data <- normalize_to_one(group_2_replicate_1_data, "normalised.expression.avg", "treatment")
group_2_replicate_2_data <- normalize_to_one(group_2_replicate_2_data, "normalised.expression.avg", "treatment")
group_2_replicate_2_data$treatment <- as.character(group_2_replicate_2_data$treatment)



aggregate(normalised.expression.avg ~ treatment, data = replicate_1_data, mean)
aggregate(normalised.expression.avg ~ treatment, data = replicate_2_data, mean)
aggregate(normalised.expression.avg ~ treatment, data = replicate_3_data, mean)
aggregate(normalised.expression.avg ~ treatment, data = group_2_replicate_1_data, mean)
aggregate(normalised.expression.avg ~ treatment, data = group_2_replicate_2_data, mean)



# merge replicates to one file ---------------

# merge
# Create an empty data frame to hold the combined data
normalised_cytokine <- data.frame()

# Process and combine replicate 1 data
replicate_1_data$replicate <- "BR1"
normalised_cytokine <- rbind(normalised_cytokine, replicate_1_data)

# Process and combine replicate 2 data
replicate_2_data$replicate <- "BR2"
normalised_cytokine <- rbind(normalised_cytokine, replicate_2_data)

# Process and combine replicate 3 data
replicate_3_data$replicate <- "BR3"
normalised_cytokine <- rbind(normalised_cytokine, replicate_3_data)


# Tady budu davat dalsi sady 2024
# Process and combine group 2 replicate 1 data
group_2_replicate_1_data$replicate <- "G2BR1"
normalised_cytokine <- rbind(normalised_cytokine, group_2_replicate_1_data)

# Process and combine group 2 replicate 1 data
group_2_replicate_2_data$replicate <- "G2BR2"
#bylo to numberic a nutno to predelat na character jako ostatni - proste tam byly jen cisla
group_2_replicate_2_data$treatment <- as.character(group_2_replicate_2_data$treatment)
normalised_cytokine <- rbind(normalised_cytokine, group_2_replicate_2_data)




# sort
normalised_cytokine <- normalised_cytokine %>%
  arrange(treatment, cytokine, replicate)

rm(replicate_1, replicate_1_data, replicate_2, replicate_2_data, replicate_3, replicate_3_data)
rm(group_2_replicate_2_data, group_2_replicate_1_data, group_2_replicate_2, group_2_replicate_1)





# table ready for heatmap -----

# Create treatment column based on replicate and level
normalised_cytokine$heatmap_treatment <- paste(normalised_cytokine$treatment, normalised_cytokine$replicate, sep = "_")


# save this table:

# zapise tabulku pro dalsi porovnavani.
write.csv(normalised_cytokine, "../normalised_cytokine.csv", row.names = FALSE)




normalised_cytokine_hm <- normalised_cytokine[, c("heatmap_treatment", "cytokine", "normalised.expression.avg")]

# Reshape the data for the heatmap
normalised_cytokine_hm <- dcast(normalised_cytokine_hm, heatmap_treatment ~ cytokine, value.var = "normalised.expression.avg")

# Remove the replicate column from the data
#heatmap_data$replicate <- NULL

rownames(normalised_cytokine_hm) <- normalised_cytokine_hm[, 1]  # Set the first column as row names
normalised_cytokine_hm <- normalised_cytokine_hm[, -1]  # Remove the first column


# heatmaps with replicates -----------------


# heatmap of unscaled data


breaks <- seq(0, 5, length.out = 101)


normalised_cytokine_hm <- as.data.frame(normalised_cytokine_hm) # makes data frame
unscaled_data <- as.matrix(t(normalised_cytokine_hm))


heatmap.2(unscaled_data, trace="none", density="none",  col=viridis(100), breaks = breaks,#col=bluered(255),
          scale="none",   Colv = F, dendrogram = "row", keysize = 0.9)

rm(unscaled_data)



#scaling of the table

# remove of pos and neg ctrls
normalised_cytokine_hm <- normalised_cytokine_hm %>%
  select(-c("Pos", "Neg", "BLANK"))


scaled_data <- as.data.frame(scale(normalised_cytokine_hm, center = TRUE, scale = TRUE))

scaled_data <- as.matrix(scaled_data)
scaled_data <- t(scaled_data)


# heatmap of scaled data
heatmap.2(scaled_data, trace="none", density="none",  col=bluered(255),
          scale="none",   Colv = F, dendrogram = "row", keysize = 0.9)


rm(scaled_data)



# statistics on aSMA levels -----------------------

# averaging the replicartees 

# Group by cytokine and treatment, and calculate the average
normalised_cytokine_avg <- normalised_cytokine %>%
  group_by(cytokine, treatment) %>%
  summarize(normalised.expression.avg = mean(normalised.expression.avg))



#Create a new column based on treatment values
normalised_cytokine_avg_stat <- normalised_cytokine_avg %>%
  mutate(aSMA.group = case_when(
    treatment %in% c(98, 101, 104, 95) ~ "aSMA_high",
    treatment %in% c(89, 97, 105, "HGFB") ~ "aSMA_low",
    TRUE ~ "Other"
  ))

normalised_cytokine_avg_stat <- normalised_cytokine_avg_stat %>%
  filter(!(cytokine %in% c("BLANK", "Pos", "Neg")))


# List of unique cytokines
unique_cytokines <- unique(normalised_cytokine_avg_stat$cytokine)

# Initialize an empty list to store results
results_list <- list()

# Loop through unique cytokines and perform t-tests
for (cytokine in unique_cytokines) {
  subset_data <- normalised_cytokine_avg_stat[normalised_cytokine_avg_stat$cytokine == cytokine, ]
  t_test_result <- t.test(normalised.expression.avg ~ aSMA.group, data = subset_data)
  
  # Store results in the list
  results_list[[cytokine]] <- list(
    Cytokine = cytokine,
    T_Statistic = t_test_result$statistic,
    P_Value = t_test_result$p.value
  )
}

# Convert the list of results into a data frame
results_df <- do.call(rbind, results_list)






# heatmaps averaged replicates -----------------

#normalised_cytokine_avg_hm <- normalised_cytokine_avg[, c("heatmap_treatment", "cytokine", "normalised.expression.avg")]
# Reshape the data for the heatmap

normalised_cytokine_avg_hm <- dcast(normalised_cytokine_avg, treatment ~ cytokine, value.var = "normalised.expression.avg")

# Remove the replicate column from the data
#heatmap_data$replicate <- NULL

rownames(normalised_cytokine_avg_hm) <- normalised_cytokine_avg_hm[, 1]  # Set the first column as row names


# THIS STEP NECESARRY ALSO FOR aSMA CORRELATION
normalised_cytokine_avg_hm <- normalised_cytokine_avg_hm[, -1]  # Remove the first column

rm(normalised_cytokine_avg) # odstrani mezivypocet prumeru


# heatmap of unscaled data

normalised_cytokine_avg_hm <- as.data.frame(normalised_cytokine_avg_hm)
unscaled_data <- as.matrix(t(normalised_cytokine_avg_hm))


# aby nezacinaly scales v -5, ale v expresi 0
breaks <- seq(0, 5, length.out = 101)

heatmap.2(unscaled_data, trace="none", density="none",   col=viridis(100), breaks = breaks,
          scale="none",   Colv = T, dendrogram = "both", keysize = 0.9)

rm(unscaled_data)




#scaling of the table

# remove of pos and neg ctrls
normalised_cytokine_avg_hm <- normalised_cytokine_avg_hm %>%
  select(-c("Pos", "Neg", "BLANK"))

scaled_data <- as.data.frame(scale(normalised_cytokine_avg_hm, center = TRUE, scale = TRUE))

## alternativne - odstraneni pacienta 98
# normalised_cytokine_avg_hm_98 <- normalised_cytokine_avg_hm[!(rownames(normalised_cytokine_avg_hm) == "98"), ]
# scaled_data <- as.data.frame(scale(normalised_cytokine_avg_hm_98, center = TRUE, scale = TRUE))

scaled_data <- as.matrix(scaled_data)
scaled_data <- t(scaled_data)


# heatmap of scaled data
heatmap.2(scaled_data, trace="none", density="none",  col=bluered(255),
          scale="none",   Colv = T, dendrogram = "both", keysize = 0.9)


rm(scaled_data)


# statistika 101 x 97 vsech genu ---------------


subset_101_97 <- subset(normalised_cytokine, treatment %in% c("101", "97"))

subset_101_97 <- subset_101_97 %>%
  filter(!(cytokine %in% c("Pos", "BLANK", "Neg")))


# List to store the results
results_list <- list()

# Get unique cytokines
unique_cytokines <- unique(subset_101_97$cytokine)

# Loop through unique cytokines and perform t-tests
results_list <- lapply(unique_cytokines, function(cytokine) {
  subset_data <- subset_101_97[subset_101_97$cytokine == cytokine, ]
  t_test_result <- t.test(normalised.expression.avg ~ treatment, data = subset_data)
  return(list(Cytokine = cytokine, T_Statistic = t_test_result$statistic, P_Value = t_test_result$p.value))
})

# Combine the results into a data frame
t_test_results_df <- do.call(rbind, results_list)

# Print the results
print(t_test_results_df)

# to je blby - porovnani 97x101 nic nevychazi signifikantne: 3 skupiny blbý....





# Boxplot panel všechny ----------

# Group data by cytokine for faceting

grouped_data <- normalised_cytokine %>%
  group_by(cytokine)




# Create the facetted ggplot
facetted_plot <- ggplot(grouped_data, aes(x = treatment, y = normalised.expression.avg, fill = treatment, color = treatment)) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x=element_text(color = "black", size=11, angle=30, vjust=1, hjust=1)) +
  xlab("treatment") +
  ylab("absorbance (AU)") +
  theme(legend.position = "none") +
  facet_wrap(~ cytokine, scales = "free_y", ncol = 3)

# Print the facetted plot
print(facetted_plot)









