# Seahorse analyser

#dependencies

library(readxl)
library(dplyr) 
library(tidyr) #na wider table
library(stringr) # extract file date part

setwd("D:/OneDrive - MUNI/Experimenty/24-01-23 - Honza aSMA CAF/23-08-01 - Seahorse HGF CAF")



#load all in folder ---------------

# Get a list of all .xlsx files in the working directory
file_list <- list.files(pattern = "\\.xlsx$", full.names = TRUE)

# Initialize an empty data frame to store the combined data
combined_data <- data.frame()


for (file in file_list) {
  current_data <- read_excel(file, sheet = "Normalized Rate")
  
  # Extract the numeric portion from the filename (assuming it's always in the format YYYYMMDD)
  numeric_part <- str_extract(basename(file), "\\d{8}")
  
  current_data$filedate <- numeric_part  # Add a new column with the extracted numeric portion
  combined_data <- rbind(combined_data, current_data)
}

# removes temp files
rm(current_data, file, file_list, numeric_part)

# udela sloupecek Well a filedate spojene aby bylo mozne pouzivat radek jako jednotku experimentu
# a prejmenuje to vsechno aby se to jmenovalo well a udela zalohu Well puvodni
combined_data <- transform(combined_data,
                           Well_orig = Well,
                           Well = paste(filedate, Well, sep = "_"))


# filtering in combined data set based on unselected list ---------

# tohle se zatim musi definovat manualne
unselected_list <- list(
  "20230620" = c("A02", "A06", "B04", "C02", "C05"),
  "20230627" = c("A02", "A06", "B06", "C05", "D03"),
  "20230801" = c("A06", "B02", "B04", "C03")
)

# nevim proc ale nejak to tu bylo treva
unique_well_orig <- unique(combined_data$Well_orig)
unique_filedate <- unique(combined_data$filedate)

# jen kontrola
cat("Unique Well_orig values:", unique_well_orig, "\n")
cat("Unique filedate values:", unique_filedate, "\n")


# Initialize an empty data frame to store the filtered data
filtered_combined_data <- data.frame()

# Loop through each filedate, filter the data, and append to filtered_combined_data
for (filedate in unique_filedate) {
  unselected <- unselected_list[[filedate]]
  
  # Filter the data based on the conditions
  filtered_data <- combined_data[!(combined_data$Well_orig %in% unselected) &
                                   (combined_data$filedate == filedate) &
                                   (combined_data$Group != "Background"), ]
  
  filtered_combined_data <- rbind(filtered_combined_data, filtered_data)
}

# rename  filtered_combined_data to simply filtered and rm all dat unnecesarry shit
filtered <- filtered_combined_data
rm(combined_data, filtered_combined_data, filtered_data, unselected_list, filedate, unique_filedate, unique_well_orig, unselected)




## load 1 (original solution) --------------
# normalized <- read_excel("20230620_Honza_HGf_CAF_filtr_dle_signalu_pozadi_ex.xlsx", 
#                                                                   sheet = "Normalized Rate")
# unselected<- c("A02", "A06", "B04", "C02", "C05")
# normalized.filtered <- normalized[!normalized$Well %in% unselected
#                                   & normalized$Group != "Background"
#                                   , ]



# Define the groups based on Measurement
filtered <- filtered %>%
  mutate(injection = case_when(
    Measurement %in% c(1, 2, 3) ~ "basal",
    Measurement %in% c(4, 5, 6) ~ "postOM",
    Measurement %in% c(7, 8, 9) ~ "post2Rot_AA",
    TRUE ~ NA_character_
  )) %>%
  group_by(Group, Well, injection, filedate, Well_orig)


filtered.means <- filtered %>%
  summarize(
    Mean_OCR = mean(OCR, na.rm = TRUE),
    Mean_ECAR = mean(ECAR, na.rm = TRUE),
    Mean_PER = mean(PER, na.rm = TRUE)
  )



#otocit tabulku aby byla jamka na radek - tidyr
filtered.means_wide <- filtered.means %>%
  pivot_wider(names_from = injection, values_from = c(Mean_OCR, Mean_ECAR, Mean_PER), names_sep = "_")


rm(filtered.means) # uz netreba, budu pouzivat bud jen filtered nebo wide

filtered.means_wide <- filtered.means_wide %>%
  mutate(
    calc.basal.OCR = Mean_OCR_basal - Mean_OCR_post2Rot_AA,
    calc.leak.OCR = Mean_OCR_postOM - Mean_OCR_post2Rot_AA,
    post_OM.ECAR = Mean_ECAR_postOM - Mean_ECAR_basal,
    post_OM.ECAR.perc = Mean_ECAR_postOM / Mean_ECAR_basal
  )






