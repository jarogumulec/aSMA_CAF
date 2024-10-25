
# AFM celkova analyza


# dependent ----------

library(readxl)
library(ggplot2)
library(ggpubr)



# data load --------------

setwd("C:/Users/Jaro-work/OneDrive - MUNI/Experimenty/24-01-23 - Honza aSMA CAF/23-09-18 - AFM CAF")
AFM <- read_excel("AFM_caf_2023-2024.xlsx")


# plot AFM--------------------
ggplot(AFM, aes(x=patient, y=mean)) +
  geom_point(position = position_jitter(seed = 1, width = 0.3), size=0.05, aes(colour = factor(patient))) +
  geom_boxplot(outlier.shape = NA, fill="0", lwd = 0.5) + # , aes(colour = factor(treatments))
  theme_bw() +
  theme(axis.text.x=element_text(color = "black", size=6.3, angle=90, vjust=0.5, hjust=1), # musi byt vjust 0.5 !!
        axis.text.y=element_text(color = "black", size=6.3),
        axis.title = element_text(size = 8))  + 
  xlab(NULL) +
  ylab("Young modulus (kPa)") +
  theme(legend.position = "none") +
  scale_y_continuous(
    #  labels = scales::percent_format(scale = 100),
    limits = c(0, 50)) + # dve ustrelene u 156 tam jsou tak odrizle, i tak staci, aby bylo videt
  stat_compare_means(aes(label = ..p.signif..),
                     method = "t.test", ref.group = "hGF", size = 2, vjust = 1,
                     method.args = list(p.adjust.method = "BH"))




ggsave("AFM.svg", plot = last_plot(),
       width = 7.5, height = 3.5, units = "cm", dpi = 300, scale = 1, limitsize = TRUE)





