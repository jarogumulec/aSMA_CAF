

library(readxl)

# nacteni vsech dat co byly udelany -------------

setwd("c:/Users/Jaro-work/OneDrive - MUNI/Experimenty/24-01-23 - Honza aSMA CAF")
# vse ulozeno do rootu

# cytokiny
cytokine <- read.csv("normalised_cytokine.csv")
# seahorse 
seahorse <- read.csv("Seahorse_complete_calculated.csv")
# western
western <- read.csv("densitometry_HGFBnormalised.csv")









# OLD - correlation a data load --------------



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











# a SMA korelace -------------
# not yet adapted for it



## load data and do matching ----



aSMA_PDPN.densitometry <- read_excel("../23-09-18 - PDPN WB denzitometrie/densitometry.xlsx")
aSMA_PDPN.densitometry <- aSMA_PDPN.densitometry[c(1:11, 13:23, 25:35),c(1:4, 19:20)]
rownames(aSMA_PDPN.densitometry) <- aSMA_PDPN.densitometry$id 





# set of desired patients - only those which are in cytokine array


#desired_patients <- c("HGFB", "M101", "M104", "F105","M89","M97", "M98")


# alternatively also 95

desired_patients <- c("101", "104", "105","89", "95", "97", "98", "HGFB")
#desired_patients <- c("M101", "M104", "F105","M89", "M95", "M97", "M98", "HGFB")


sequence <- c("101_BR1", "101_BR2", "101_BR3", "104_BR1", "104_BR2","104_BR3","105_BR1","105_BR2","105_BR3","89_BR1","89_BR2",
              "89_BR3","95_BR1","95_BR2","95_BR3","97_BR1","97_BR2","97_BR3", "98_BR2","HGFB_BR1","HGFB_BR2","HGFB_BR3")



# Subset the table based on the desired patient order

aSMA_PDPN.densitometry.subset <- aSMA_PDPN.densitometry[aSMA_PDPN.densitometry$sample.PDPN %in% desired_patients, ]

#asma <- aSMA_densitometryCAF[aSMA_densitometryCAF$patient %in% desired_patients, ]

# Reorder the subsetted table based on the desired patient order

aSMA_PDPN.densitometry.subset <- aSMA_PDPN.densitometry.subset[match(sequence, aSMA_PDPN.densitometry.subset$id), ]


# #name control 
# row.names(normalised_cytokine_avg_hm)
# colnames(normalised_cytokine_avg_hm)
# asma$patient


# # remove POS and NEG from 
# 
# normalised_cytokine_avg_hm_without_pos_c <- normalised_cytokine_avg_hm[, -c(1, 21, 22)]

# remove unrelated column from ASMA
asma <- aSMA_PDPN.densitometry.subset$`ASMA.normalised.int.1-0`
PDPN <- aSMA_PDPN.densitometry.subset$`PDPN.normalised.int.1-0`




## correlation ASMA-------
correlations <- cor(asma, normalised_cytokine_hm)
correlations <- as.data.frame(t(correlations))
correlations$cytokine <- rownames(correlations)


plt <- 
  ggplot(correlations, aes(x = reorder(cytokine, V1), y = V1)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(x = "Cytokine", y = "Correlation") +
  ggtitle("Correlation with aSMA, BR-pair-wise") +
  coord_flip() +
  theme_bw() +
  theme(axis.text.x=element_text(color = "black", size=6.3, vjust=0.5, hjust=0.5), # musi byt vjust 0.5 !!
        axis.text.y=element_text(color = "black", size=6.3),
        axis.title = element_text(size = 8))   


plt

ggsave("plot.svg", plot = plt,
       width = 4, height = 6, units = "cm", dpi = 300, scale = 1, limitsize = TRUE)




## correlation PDPN-------
correlations <- cor(PDPN, normalised_cytokine_hm)
correlations <- as.data.frame(t(correlations))
correlations$cytokine <- rownames(correlations)



plt <- 
  ggplot(correlations, aes(x = reorder(cytokine, V1), y = V1)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(x = "Cytokine", y = "Correlation") +
  ggtitle("Correlation with PDPN, BR-pair-wise") +
  coord_flip() +
  theme_bw() +
  theme(axis.text.x=element_text(color = "black", size=6.3, vjust=0.5, hjust=0.5), # musi byt vjust 0.5 !!
        axis.text.y=element_text(color = "black", size=6.3),
        axis.title = element_text(size = 8))   


plt

ggsave("plot.svg", plot = plt,
       width = 4, height = 6, units = "cm", dpi = 300, scale = 1, limitsize = TRUE)

