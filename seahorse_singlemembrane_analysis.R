# Seahorse analyser

#dependencies

library(readxl)
library(dplyr) 
library(tidyr) #na wider table
library(stringr) # extract file date part
library(ggplot2)
library(ggpubr)

setwd("D:/OneDrive - MUNI/Experimenty/24-01-23 - Honza aSMA CAF/23-08-01 - Seahorse HGF CAF")

# nove
setwd("C:/Users/Jaro-work/OneDrive - MUNI/Experimenty/24-01-23 - Honza aSMA CAF/23-08-01 - Seahorse HGF CAF")


#load all in folder ---------------

# Get a list of all .xlsx files in the working directory
file_list <- list.files(pattern = "\\.xlsx$", full.names = TRUE)

# Initialize an empty data frame to store the combined data
combined_data <- data.frame()


# in following alternatively put "Normalized Rate" OR "Rate"

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




# Filtration 1: whole plates ----------------
#new 2024-05-15
# here i filter whole plate if error

#define dates to filter
filter_dates <- c(20230704)

combined_data <- combined_data %>%
  filter(!filedate %in% filter_dates)



# Filtration 2: specific wells based on unselected list ---------

# tohle se zatim musi definovat manualne
unselected_list <- list(
  "20230620" = c("A02", "A06", "B04", "C02", "C05"),
  "20230627" = c("A02", "A06", "B06", "C05", "D03"),
  "20230801" = c("A06", "B02", "B04", "C03"),
  "20240426" = c("A02", "A06"),
  "20240429" = c("A02", "A05", "B03", "C04", "C05", "C06")
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


# add group set which shows both date and treatment
filtered <- transform(filtered, Group_exp = paste(Group, filedate, sep = "_"))


rm(combined_data, filtered_combined_data, filtered_data, unselected_list, filedate, unique_filedate, unique_well_orig, unselected)




## Filtration 3: outliers defined ex-post -------------------------

# use with care, i explicitly defined what to exclude because of outliers!


## backuphelper
#filtered_bckp_pre2ndfilter <- filtered
#filtered <- filtered_bckp_pre2ndfilter


filtered <- filtered[!(filtered$Group_exp %in% 
                         c("CAF 104_20230627", "CAF 105_20230627", "CAF 89_20230627", "HGFb_20230627")), ]






## load 1 (original solution) --------------
# normalized <- read_excel("20230620_Honza_HGf_CAF_filtr_dle_signalu_pozadi_ex.xlsx", 
#                                                                   sheet = "Normalized Rate")
# unselected<- c("A02", "A06", "B04", "C02", "C05")
# normalized.filtered <- normalized[!normalized$Well %in% unselected
#                                   & normalized$Group != "Background"
#                                   , ]



# Calculating subtracked OCR and ECAR ---------------

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
# do calculations
filtered.means_wide <- filtered.means_wide %>%
  mutate(
    calc.basal.OCR = Mean_OCR_basal - Mean_OCR_post2Rot_AA,
    calc.leak.OCR = Mean_OCR_postOM - Mean_OCR_post2Rot_AA,
    calc.post_OM.ECAR = Mean_ECAR_postOM - Mean_ECAR_basal,
    calc.post_OM.ECAR.perc = Mean_ECAR_postOM / Mean_ECAR_basal,
    calc.OCR_to_ecar = Mean_OCR_basal / Mean_ECAR_basal
  )

# and make category for averaging repetitions experiment-wise treatment-wise
filtered.means_wide <- transform(filtered.means_wide,
                           Group_exp = paste(Group, filedate, sep = "_"))


# now filtered is original for timelapse plotting and wide is for stats



# dopoc ATP-linked -----------------
#2024-04
# jeste vypocita ATP-linked po odecteni leaku protoze leak je vysoky
filtered.means_wide <- filtered.means_wide %>%
  mutate(
    calc.OCR_ATP_linked = calc.basal.OCR - calc.leak.OCR
  )



# plot replicate-wise--------------------


ggplot(filtered.means_wide, aes(x=Group_exp, y=calc.basal.OCR)) +
  geom_point(position = position_jitter(seed = 1, width = 0.3), size=0.05, aes(colour = factor(Group))) +
  #scale_colour_manual(values = mypalette) +
  geom_boxplot(outlier.shape = NA, fill="0", lwd = 0.5) + # , aes(colour = factor(treatments))
  theme_bw() +
  theme(axis.text.x=element_text(color = "black", size=6.3, angle=90, vjust=0.5, hjust=1), # musi byt vjust 0.5 !!
        axis.text.y=element_text(color = "black", size=6.3),
        axis.title = element_text(size = 8))  + 
  xlab(NULL) +
  ylab("Basal OCR") +
  theme(legend.position = "none") +
  scale_y_continuous(
  #  labels = scales::percent_format(scale = 100),
    limits = c(0, NA)) #+
  stat_compare_means(aes(label = ..p.signif..),
                     method = "t.test", ref.group = "control", size=2, vjust=1)

ggsave("bocr.svg", plot = last_plot(),
       width = 8, height = 4.27, units = "cm", dpi = 300, scale = 1, limitsize = TRUE)





ggplot(filtered.means_wide, aes(x=Group_exp, y=calc.leak.OCR)) +
  geom_point(position = position_jitter(seed = 1, width = 0.3), size=0.05, aes(colour = factor(Group))) +
  geom_boxplot(outlier.shape = NA, fill="0", lwd = 0.5) + # , aes(colour = factor(treatments))
  theme_bw() +
  theme(axis.text.x=element_text(color = "black", size=6.3, angle=90, vjust=0.5, hjust=1), # musi byt vjust 0.5 !!
        axis.text.y=element_text(color = "black", size=6.3),
        axis.title = element_text(size = 8))  + 
  xlab(NULL) +
  ylab("Proton leak") +
  theme(legend.position = "none") +
  scale_y_continuous(
    #  labels = scales::percent_format(scale = 100),
    limits = c(-1, NA)) #+
stat_compare_means(aes(label = ..p.signif..),
                   method = "t.test", ref.group = "control", size=2, vjust=1)

ggsave("leakocr.svg", plot = last_plot(),
       width = 5, height = 4.27, units = "cm", dpi = 300, scale = 1, limitsize = TRUE)


ggplot(filtered.means_wide, aes(x=Group_exp, y=Mean_ECAR_basal)) +
  geom_point(position = position_jitter(seed = 1, width = 0.3), size=0.05, aes(colour = factor(Group))) +
  geom_boxplot(outlier.shape = NA, fill="0", lwd = 0.5) + # , aes(colour = factor(treatments))
  theme_bw() +
  theme(axis.text.x=element_text(color = "black", size=6.3, angle=90, vjust=0.5, hjust=1), # musi byt vjust 0.5 !!
        axis.text.y=element_text(color = "black", size=6.3),
        axis.title = element_text(size = 8))  + 
  xlab(NULL) +
  ylab("Basal ECAR") +
  theme(legend.position = "none") +
  scale_y_continuous(
    #  labels = scales::percent_format(scale = 100),
    limits = c(0, NA)) #+
stat_compare_means(aes(label = ..p.signif..),
                   method = "t.test", ref.group = "control", size=2, vjust=1)

ggsave("becar.svg", plot = last_plot(),
       width = 5, height = 4.27, units = "cm", dpi = 300, scale = 1, limitsize = TRUE)


ggplot(filtered.means_wide, aes(x=Group_exp, y=calc.post_OM.ECAR.perc)) +
  geom_point(position = position_jitter(seed = 1, width = 0.3), size=0.05, aes(colour = factor(Group))) +
  geom_boxplot(outlier.shape = NA, fill="0", lwd = 0.5) + # , aes(colour = factor(treatments))
  theme_bw() +
  theme(axis.text.x=element_text(color = "black", size=6.3, angle=90, vjust=0.5, hjust=1), # musi byt vjust 0.5 !!
        axis.text.y=element_text(color = "black", size=6.3),
        axis.title = element_text(size = 8))  + 
  xlab(NULL) +
  ylab("post-OM ECAR") +
  theme(legend.position = "none") +
  scale_y_continuous(
      labels = scales::percent_format(scale = 100),
    limits = c(NA, NA)) #+
stat_compare_means(aes(label = ..p.signif..),
                   method = "t.test", ref.group = "control", size=2, vjust=1)

ggsave("OMecar.svg", plot = last_plot(),
       width = 5, height = 4.27, units = "cm", dpi = 300, scale = 1, limitsize = TRUE)

ggplot(filtered.means_wide, aes(x=Group_exp, y=calc.OCR_to_ecar)) +
  geom_point(position = position_jitter(seed = 1, width = 0.3), size=0.05, aes(colour = factor(Group))) +
  geom_boxplot(outlier.shape = NA, fill="0", lwd = 0.5) + # , aes(colour = factor(treatments))
  theme_bw() +
  theme(axis.text.x=element_text(color = "black", size=6.3, angle=90, vjust=0.5, hjust=1), # musi byt vjust 0.5 !!
        axis.text.y=element_text(color = "black", size=6.3),
        axis.title = element_text(size = 8))  + 
  xlab(NULL) +
  ylab("OCR/ECAR") +
  theme(legend.position = "none") +
  scale_y_continuous(
    #labels = scales::percent_format(scale = 100),
    limits = c(NA, NA)) #+
stat_compare_means(aes(label = ..p.signif..),
                   method = "t.test", ref.group = "control", size=2, vjust=1)

ggsave("ocrtoecar.svg", plot = last_plot(),
       width = 5, height = 4.27, units = "cm", dpi = 300, scale = 1, limitsize = TRUE)







ggplot(filtered.means_wide, aes(x=Group_exp, y=calc.OCR_ATP_linked)) +
  geom_point(position = position_jitter(seed = 1, width = 0.3), size=0.05, aes(colour = factor(Group))) +
  #scale_colour_manual(values = mypalette) +
  geom_boxplot(outlier.shape = NA, fill="0", lwd = 0.5) + # , aes(colour = factor(treatments))
  theme_bw() +
  theme(axis.text.x=element_text(color = "black", size=6.3, angle=90, vjust=0.5, hjust=1), # musi byt vjust 0.5 !!
        axis.text.y=element_text(color = "black", size=6.3),
        axis.title = element_text(size = 8))  + 
  xlab(NULL) +
  ylab("ATP-linked OCR") +
  theme(legend.position = "none") +
  scale_y_continuous(
    #  labels = scales::percent_format(scale = 100),
    limits = c(0, NA)) #+
stat_compare_means(aes(label = ..p.signif..),
                   method = "t.test", ref.group = "control", size=2, vjust=1)

ggsave("atplinkedocr.svg", plot = last_plot(),
       width = 5, height = 5.27, units = "cm", dpi = 300, scale = 1, limitsize = TRUE)







# plot stats from wide table (paper, var 2) -----------


# prejmenovani duplikatu 95
filtered.means_wide_forstatistics <- filtered.means_wide %>%
  mutate(Group = recode(Group, "CAF 95-1" = "CAF 95", "CAF 95-2" = "CAF 95", "HGfb" = "HGFb"))


ggplot(filtered.means_wide_forstatistics, aes(x=Group, y=calc.basal.OCR)) +
  geom_point(position = position_jitter(seed = 1, width = 0.3), size=0.05, aes(colour = factor(Group))) +
  geom_boxplot(outlier.shape = NA, fill="0", lwd = 0.5) + # , aes(colour = factor(treatments))
  theme_bw() +
  theme(axis.text.x=element_text(color = "black", size=6.3, angle=90, vjust=0.5, hjust=1), # musi byt vjust 0.5 !!
        axis.text.y=element_text(color = "black", size=6.3),
        axis.title = element_text(size = 8))  + 
  xlab(NULL) +
  ylab("Basal OCR") +
  theme(legend.position = "none") +
  scale_y_continuous(
    #  labels = scales::percent_format(scale = 100),
    limits = c(0, NA)) +
  stat_compare_means(aes(label = ..p.signif..),
                     method = "t.test", ref.group = "HGFb", size = 2, vjust = 1,
                     method.args = list(p.adjust.method = "BH"))

ggsave("bocr_grouped.svg", plot = last_plot(),
       width = 4, height = 3, units = "cm", dpi = 300, scale = 1, limitsize = TRUE)







ggplot(filtered.means_wide_forstatistics, aes(x=Group, y=calc.OCR_to_ecar)) +
  geom_point(position = position_jitter(seed = 1, width = 0.3), size=0.05, aes(colour = factor(Group))) +
  geom_boxplot(outlier.shape = NA, fill="0", lwd = 0.5) + # , aes(colour = factor(treatments))
  theme_bw() +
  theme(axis.text.x=element_text(color = "black", size=6.3, angle=90, vjust=0.5, hjust=1), # musi byt vjust 0.5 !!
        axis.text.y=element_text(color = "black", size=6.3),
        axis.title = element_text(size = 8))  + 
  xlab(NULL) +
  ylab("OCR to ECAR ratio") +
  theme(legend.position = "none") +
  scale_y_continuous(
    #    labels = scales::percent_format(scale = 100),
    limits = c(0, NA)) +
  stat_compare_means(aes(label = ..p.signif..),
                     method = "t.test", ref.group = "HGFb", size = 2, vjust = 1,
                     method.args = list(p.adjust.method = "BH"))

ggsave("ocrtoecar_grouped.svg", plot = last_plot(),
       width = 4, height = 3, units = "cm", dpi = 300, scale = 1, limitsize = TRUE)





ggplot(filtered.means_wide_forstatistics, aes(x=Group, y=calc.OCR_ATP_linked)) +
  geom_point(position = position_jitter(seed = 1, width = 0.3), size=0.05, aes(colour = factor(Group))) +
  geom_boxplot(outlier.shape = NA, fill="0", lwd = 0.5) + # , aes(colour = factor(treatments))
  theme_bw() +
  theme(axis.text.x=element_text(color = "black", size=6.3, angle=90, vjust=0.5, hjust=1), # musi byt vjust 0.5 !!
        axis.text.y=element_text(color = "black", size=6.3),
        axis.title = element_text(size = 8))  + 
  xlab(NULL) +
  ylab("OCR ATP-linked") +
  theme(legend.position = "none") +
  scale_y_continuous(
    #  labels = scales::percent_format(scale = 100),
    limits = c(0, NA)) +
  stat_compare_means(aes(label = ..p.signif..),
                     method = "t.test", ref.group = "HGFb", size = 2, vjust = 1,
                     method.args = list(p.adjust.method = "BH"))

ggsave("ocr_atplinked_grouped.svg", plot = last_plot(),
       width = 4, height = 3, units = "cm", dpi = 300, scale = 1, limitsize = TRUE)





ggplot(filtered.means_wide_forstatistics, aes(x=Group, y=calc.leak.OCR)) +
  geom_point(position = position_jitter(seed = 1, width = 0.3), size=0.05, aes(colour = factor(Group))) +
  geom_boxplot(outlier.shape = NA, fill="0", lwd = 0.5) + # , aes(colour = factor(treatments))
  theme_bw() +
  theme(axis.text.x=element_text(color = "black", size=6.3, angle=90, vjust=0.5, hjust=1), # musi byt vjust 0.5 !!
        axis.text.y=element_text(color = "black", size=6.3),
        axis.title = element_text(size = 8))  + 
  xlab(NULL) +
  ylab("Proton leak") +
  theme(legend.position = "none") +
  scale_y_continuous(
    #  labels = scales::percent_format(scale = 100),
    limits = c(-1, 30)) +
  stat_compare_means(aes(label = ..p.signif..),
                     method = "t.test", ref.group = "HGFb", size = 2, vjust = 1,
                     method.args = list(p.adjust.method = "BH"))

ggsave("leakocr_grouped.svg", plot = last_plot(),
       width = 4, height = 3, units = "cm", dpi = 300, scale = 1, limitsize = TRUE)


ggplot(filtered.means_wide_forstatistics, aes(x=Group, y=Mean_ECAR_basal)) +
  geom_point(position = position_jitter(seed = 1, width = 0.3), size=0.05, aes(colour = factor(Group))) +
  geom_boxplot(outlier.shape = NA, fill="0", lwd = 0.5) + # , aes(colour = factor(treatments))
  theme_bw() +
  theme(axis.text.x=element_text(color = "black", size=6.3, angle=90, vjust=0.5, hjust=1), # musi byt vjust 0.5 !!
        axis.text.y=element_text(color = "black", size=6.3),
        axis.title = element_text(size = 8))  + 
  xlab(NULL) +
  ylab("Basal ECAR") +
  theme(legend.position = "none") +
  scale_y_continuous(
    #  labels = scales::percent_format(scale = 100),
    limits = c(0, NA)) +
  stat_compare_means(aes(label = ..p.signif..),
                     method = "t.test", ref.group = "HGFb", size = 2, vjust = 1,
                     method.args = list(p.adjust.method = "BH"))

ggsave("becar_grouped.svg", plot = last_plot(),
       width = 4, height = 3, units = "cm", dpi = 300, scale = 1, limitsize = TRUE)




# plot average (paper) -------------
# average individual experiments (1 nr per biol replicate)

average_data <- filtered.means_wide %>%
  group_by(Group, Group_exp) %>%
  summarize(
    avg_calc_basal_OCR = mean(calc.basal.OCR),
    avg_calc_leak_OCR = mean(calc.leak.OCR),
    avg_basal_ECAR = mean(Mean_ECAR_basal),
    avg_post_OM_ECAR = mean(calc.post_OM.ECAR),
    avg_ratio_ECAR = mean(calc.post_OM.ECAR.perc),
    avg_ratio_OCR_to_ECAR = mean(calc.OCR_to_ecar),
    avg_calc.OCR_ATP_linked = mean(calc.OCR_ATP_linked)
  )

# prejmenovani duplikatu 95
average_data <- average_data %>%
  mutate(Group = recode(Group, "CAF 95-1" = "CAF 95", "CAF 95-2" = "CAF 95", "HGfb" = "HGFb"))





##

ggplot(average_data, aes(x=Group, y=avg_calc_basal_OCR)) +
  geom_point(position = position_jitter(seed = 1, width = 0.3), size=0.05, aes(colour = factor(Group))) +
  geom_boxplot(outlier.shape = NA, fill="0", lwd = 0.5) + # , aes(colour = factor(treatments))
  theme_bw() +
  theme(axis.text.x=element_text(color = "black", size=6.3, angle=90, vjust=0.5, hjust=1), # musi byt vjust 0.5 !!
        axis.text.y=element_text(color = "black", size=6.3),
        axis.title = element_text(size = 8))  + 
  xlab(NULL) +
  ylab("Basal OCR") +
  theme(legend.position = "none") +
  scale_y_continuous(
    #  labels = scales::percent_format(scale = 100),
    limits = c(0, NA)) +
  stat_compare_means(aes(label = ..p.signif..),
                     method = "t.test", ref.group = "HGFb", size = 2, vjust = 1,
                     method.args = list(p.adjust.method = "BH"))
  
  
  
   
ggsave("bocr_mean.svg", plot = last_plot(),
       width = 4, height = 3, units = "cm", dpi = 300, scale = 1, limitsize = TRUE)





ggplot(average_data, aes(x=Group, y=avg_calc.OCR_ATP_linked)) +
  geom_point(position = position_jitter(seed = 1, width = 0.3), size=0.05, aes(colour = factor(Group))) +
  geom_boxplot(outlier.shape = NA, fill="0", lwd = 0.5) + # , aes(colour = factor(treatments))
  theme_bw() +
  theme(axis.text.x=element_text(color = "black", size=6.3, angle=90, vjust=0.5, hjust=1), # musi byt vjust 0.5 !!
        axis.text.y=element_text(color = "black", size=6.3),
        axis.title = element_text(size = 8))  + 
  xlab(NULL) +
  ylab("OCR ATP-linked") +
  theme(legend.position = "none") +
  scale_y_continuous(
    #  labels = scales::percent_format(scale = 100),
    limits = c(0, NA)) +
stat_compare_means(aes(label = ..p.signif..),
                   method = "t.test", ref.group = "HGFb", size = 2, vjust = 1,
                   method.args = list(p.adjust.method = "BH"))

ggsave("ocr_atplinked.svg", plot = last_plot(),
       width = 4, height = 3, units = "cm", dpi = 300, scale = 1, limitsize = TRUE)





ggplot(average_data, aes(x=Group, y=avg_calc_leak_OCR)) +
  geom_point(position = position_jitter(seed = 1, width = 0.3), size=0.05, aes(colour = factor(Group))) +
  geom_boxplot(outlier.shape = NA, fill="0", lwd = 0.5) + # , aes(colour = factor(treatments))
  theme_bw() +
  theme(axis.text.x=element_text(color = "black", size=6.3, angle=90, vjust=0.5, hjust=1), # musi byt vjust 0.5 !!
        axis.text.y=element_text(color = "black", size=6.3),
        axis.title = element_text(size = 8))  + 
  xlab(NULL) +
  ylab("Proton leak") +
  theme(legend.position = "none") +
  scale_y_continuous(
    #  labels = scales::percent_format(scale = 100),
    limits = c(0, 30)) +
  stat_compare_means(aes(label = ..p.signif..),
                     method = "t.test", ref.group = "HGFb", size = 2, vjust = 1,
                     method.args = list(p.adjust.method = "BH"))

ggsave("leakocr_mean.svg", plot = last_plot(),
       width = 4, height = 3, units = "cm", dpi = 300, scale = 1, limitsize = TRUE)


ggplot(average_data, aes(x=Group, y=avg_basal_ECAR)) +
  geom_point(position = position_jitter(seed = 1, width = 0.3), size=0.05, aes(colour = factor(Group))) +
  geom_boxplot(outlier.shape = NA, fill="0", lwd = 0.5) + # , aes(colour = factor(treatments))
  theme_bw() +
  theme(axis.text.x=element_text(color = "black", size=6.3, angle=90, vjust=0.5, hjust=1), # musi byt vjust 0.5 !!
        axis.text.y=element_text(color = "black", size=6.3),
        axis.title = element_text(size = 8))  + 
  xlab(NULL) +
  ylab("Basal ECAR") +
  theme(legend.position = "none") +
  scale_y_continuous(
    #  labels = scales::percent_format(scale = 100),
    limits = c(0, NA)) +
stat_compare_means(aes(label = ..p.signif..),
                   method = "t.test", ref.group = "HGFb", size = 2, vjust = 1,
                   method.args = list(p.adjust.method = "BH"))

ggsave("becar_mean.svg", plot = last_plot(),
       width = 4, height = 3, units = "cm", dpi = 300, scale = 1, limitsize = TRUE)


ggplot(average_data, aes(x=Group, y=avg_ratio_ECAR)) +
  geom_point(position = position_jitter(seed = 1, width = 0.3), size=0.05, aes(colour = factor(Group))) +
  geom_boxplot(outlier.shape = NA, fill="0", lwd = 0.5) + # , aes(colour = factor(treatments))
  theme_bw() +
  theme(axis.text.x=element_text(color = "black", size=6.3, angle=90, vjust=0.5, hjust=1), # musi byt vjust 0.5 !!
        axis.text.y=element_text(color = "black", size=6.3),
        axis.title = element_text(size = 8))  + 
  xlab(NULL) +
  ylab("post-OM ECAR ratio") +
  theme(legend.position = "none") +
  scale_y_continuous(
      labels = scales::percent_format(scale = 100),
    limits = c(0.1, NA)) +
stat_compare_means(aes(label = ..p.signif..),
                   method = "t.test", ref.group = "HGFb", size = 2, vjust = 1,
                   method.args = list(p.adjust.method = "BH"))

ggsave("OMecar_mean.svg", plot = last_plot(),
       width = 4, height = 3, units = "cm", dpi = 300, scale = 1, limitsize = TRUE)

ggplot(average_data, aes(x=Group, y=avg_ratio_OCR_to_ECAR)) +
  geom_point(position = position_jitter(seed = 1, width = 0.3), size=0.05, aes(colour = factor(Group))) +
  geom_boxplot(outlier.shape = NA, fill="0", lwd = 0.5) + # , aes(colour = factor(treatments))
  theme_bw() +
  theme(axis.text.x=element_text(color = "black", size=6.3, angle=90, vjust=0.5, hjust=1), # musi byt vjust 0.5 !!
        axis.text.y=element_text(color = "black", size=6.3),
        axis.title = element_text(size = 8))  + 
  xlab(NULL) +
  ylab("OCR to ECAR ratio") +
  theme(legend.position = "none") +
  scale_y_continuous(
#    labels = scales::percent_format(scale = 100),
    limits = c(0, NA)) +
stat_compare_means(aes(label = ..p.signif..),
                   method = "t.test", ref.group = "HGFb", size = 2, vjust = 1,
                   method.args = list(p.adjust.method = "BH"))

ggsave("ocrtoecar_mean.svg", plot = last_plot(),
       width = 4, height = 3, units = "cm", dpi = 300, scale = 1, limitsize = TRUE)



# time curves for ECAR and OCR --------

#this is an analogue for output from Wave/ from excel, but here its integrated



##time - showing replicates-----
#vytvori group date identifikátor
filtered <- transform(filtered, Group_exp = paste(Group, filedate, sep = "_"))

summary_data <- filtered %>%
  group_by(Group_exp, Group, Measurement, filedate) %>%
  summarise(
    Mean_OCR = mean(OCR),
    SD_OCR = sd(OCR)
  )

# Create the line chart with whiskers
# rovnice prepocitava z cisla measurement 
ggplot(summary_data, aes(x = (8.552 * Measurement - 7.0683), y = Mean_OCR, group = Group_exp, color = Group)) +
  geom_line() +
    geom_point(position = position_jitter(seed = 1, width = 0.3), size = 0.05, aes(colour = factor(Group))) +
  geom_errorbar(aes(ymin = Mean_OCR - SD_OCR, ymax = Mean_OCR + SD_OCR), width = 5, position = position_dodge(width = 1)) +
  theme_bw() +
  labs(x = "Time (min)", y = "Mean OCR") +
  theme(
    axis.text.x = element_text(color = "black", size = 6.3),
    axis.text.y = element_text(color = "black", size = 6.3),
    axis.title = element_text(size = 8),
    legend.position = "top",
    legend.text = element_text(size = 6.3)
  ) +
  scale_y_continuous(limits = c(0, NA))

ggsave("ocr_time1.svg", plot = last_plot(),
       width = 6, height = 6, units = "cm", dpi = 300, scale = 1, limitsize = TRUE)




## graf prumeru podle replik-----------


# prejmenovani duplikatu 95
filtered.95 <- filtered %>%
  mutate(Group = recode(Group, "CAF 95-1" = "CAF 95", "CAF 95-2" = "CAF 95", "HGfb" = "HGFb"))


#prumer podle Group a ne experimentu jak driv
summary_data <- filtered.95 %>%
  group_by(Group, Measurement) %>%
  summarise(
    Mean_OCR = mean(OCR),
    SD_OCR = sd(OCR)
  )


# rovnice prepocitava z cisla measurement 
ggplot(summary_data, aes(x = (8.552 * Measurement - 7.0683), y = Mean_OCR, group = Group, color = Group)) +
  geom_line() +
  geom_point(position = position_jitter(seed = 1, width = 0.3), size = 0.05, aes(colour = factor(Group))) +
  geom_errorbar(aes(ymin = Mean_OCR - SD_OCR, ymax = Mean_OCR + SD_OCR), width = 5, position = position_dodge(width = 1)) +
  theme_bw() +
  labs(x = "Time (min)", y = "Mean OCR") +
  theme(
    axis.text.x = element_text(color = "black", size = 6.3),
    axis.text.y = element_text(color = "black", size = 6.3),
    axis.title = element_text(size = 8),
    legend.position = "top",
    legend.text = element_text(size = 6.3)
  ) +
  scale_y_continuous(limits = c(0, NA))



### testovaci - SE misto SD pro tento graf---------------
ggsave("ocr_time2.svg", plot = last_plot(),
       width = 6, height = 6, units = "cm", dpi = 300, scale = 1, limitsize = TRUE)

# verze se SE misto SD
summary_data <- filtered.95 %>%
  group_by(Group, Measurement) %>%
  summarise(
    Mean_OCR = mean(OCR),
    SD_OCR = mean_se(OCR)
  )


ggplot(summary_data, aes(x = (8.552 * Measurement - 7.0683), y = Mean_OCR, group = Group, color = Group)) +
  geom_line() +
  geom_point(position = position_jitter(seed = 1, width = 0.3), size = 0.05, aes(colour = factor(Group))) +
  geom_errorbar(aes(ymin = SD_OCR$ymin, ymax = SD_OCR$ymax), width = 5, position = position_dodge(width = 1)) +
  theme_bw() +
  labs(x = "Time (min)", y = "Mean OCR") +
  theme(
    axis.text.x = element_text(color = "black", size = 6.3),
    axis.text.y = element_text(color = "black", size = 6.3),
    axis.title = element_text(size = 8),
    legend.position = "top",
    legend.text = element_text(size = 6.3)
  ) +
  scale_y_continuous(limits = c(0, NA))