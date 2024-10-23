
# WB celkova analyza


# dependent ----------

library(readxl)



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


## approach 2: to HGFB of every membrane -----------

# downside: The HGFB will be 1.0 without variability at all. statistic needed would be One sample test

densitometry_HGFBnormalised <- read_excel("densitometry_vse_dohromady.xlsx")

# removes empty cases
densitometry_HGFBnormalised <- densitometry_HGFBnormalised[densitometry_HGFBnormalised$id != "" & !is.na(densitometry_HGFBnormalised$id), ]
#removes unrelevant columns
densitometry_HGFBnormalised <- densitometry_HGFBnormalised[,c(1:5, 22:24)]

write.csv(densitometry_HGFBnormalised, "../densitometry_HGFBnormalised.csv", row.names = FALSE)
#



# # if needed name from id to row name
# aSMA_PDPN.densitometry <- as.data.frame(aSMA_PDPN.densitometry)
# rownames(aSMA_PDPN.densitometryf) <- aSMA_PDPN.densitometry$id
