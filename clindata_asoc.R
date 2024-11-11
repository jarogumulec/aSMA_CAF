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
  rm(merged_data) 
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



# 1 clinical characteristics -------------




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

