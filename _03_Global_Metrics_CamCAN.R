library(tidyverse)
library(janitor)
library(plyr)
library(readxl)
library(car)
library(Rmisc)
library(tidylog)
library(plotly)
library(jsonlite)

rm(list = ls())

# Import processed data---------------------------------------------------------
source("_01_DataManip_CamCAN.R")

################################################################################
listfile_components <- list.files(getwd(), pattern = "*.txt")[grep("components", list.files(getwd(), pattern = "*.txt"))]
n_subj <- 645
n_threshold <- 5

components <- ldply(listfile_components, read.table, header = T, sep = "\t") %>%
  mutate(threshold = rep(c(.1, .12, .15, .17, .2), each = n_subj)) %>%
  relocate(threshold, .after = (X)) %>%
  plyr::rename(c("X" = "Subj_ID")) %>%
  replace("Subj_ID", rep(seq_len(n_subj), times = n_threshold)) %>%
  pivot_longer(
    cols = !c("Subj_ID", "threshold"),
    names_to = "Region",
    values_to = "components"
  )

# # Get proportion of LLC per threshold across subjects
LLC <- components %>%
  group_by(Subj_ID) %>%
  count(threshold, components) %>%
  group_by(Subj_ID, threshold) %>%
  mutate(prop = prop.table(n)) %>%
  slice_max(prop, n = 1) %>%
  group_by(threshold) %>%
  summarise_at(vars(prop), mean) %>%
  plyr::rename(c("prop" = "Largest Connected Component"))

data_full_per_subject_LLC <- merge(data_full_per_subject, LLC, by = "threshold") %>%
  arrange(Subj_ID, threshold)

# Evolution of Global metrics - what is the optimal threshold?
evo <- data_full_per_subject_LLC %>%
  group_by(threshold) %>%
  summarise_at(vars(Eglob, Clustering_coeff_glob, Eloc, `Largest Connected Component`, Q_consensus), funs(mean)) %>%
  pivot_longer(
    cols = !c("threshold"),
    names_to = "Metrics"
  )

plotly::ggplotly(
  evo %>%
    ggplot(aes(threshold, value, color = Metrics)) +
    geom_line() +
    geom_point(size = 2) +
    geom_jitter(height = 0.05, alpha = 0.2) +
    ggtitle("Evolution of global metrics as a function of the threshold") +
    xlab("Threshold") +
    ylab("") +
    scale_x_continuous(breaks = c(0.1, 0.12, 0.15, 0.17, 0.2)) +
    scale_y_continuous(breaks = seq(0, 1, 0.05)) +
    coord_cartesian(ylim = c(0.25, 1)) +
    geom_rect(
      aes(
        xmin = 0.14,
        xmax = 0.16,
        ymin = 0.25,
        ymax = 1.01
      ),
      fill = "red", alpha = 0.2, color = "red", linewidth = 0.1
    ) +
    ggpubr::theme_pubclean()
)


# AMI for all threshold consensus clustering
source("_NMI&AMI_functions.R")
# Make sure this is indexed on 131 observations only
data_AMI <- data_full_per_region %>% subset(threshold == "0.15")

c <- AMI_func(factor(data_AMI$Consensus_vector_0.1), factor(data_AMI$Consensus_vector_0.15))
d <- AMI_func(factor(data_AMI$Consensus_vector_0.12), factor(data_AMI$Consensus_vector_0.15))
e <- AMI_func(factor(data_AMI$Consensus_vector_0.17), factor(data_AMI$Consensus_vector_0.15))
f <- AMI_func(factor(data_AMI$Consensus_vector_0.2), factor(data_AMI$Consensus_vector_0.15))
# Mean AMI with partitions at other thresholds
(c + d + e + f) / 5

# 75% of AMI between the affiliation vectors at different thresholds
