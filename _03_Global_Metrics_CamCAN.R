# Import processed data---------------------------------------------------------
source("_01_DataManip_CamCAN.R")

################################################################################
library(plotly)

data_full_per_subject_LLC <- merge(data_full_per_subject, LLC_threshold, by = "threshold") %>%
  arrange(Subj_ID, threshold)

# Evolution of Global metrics - what is the optimal threshold?
evo <- data_full_per_subject_LLC %>%
  group_by(threshold) %>%
  summarise_at(vars(Eglob, Clustering_coeff_glob, Eloc, `Largest Connected Component`), funs(mean)) %>%
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


# OMST overall FC correlation with age
setwd(paste0(getwd(), "/OMST"))
mean_cor <- as.data.frame(fromJSON("mean_cor.json")) %>%
  mutate(Subj_ID = rep(seq_len(628)))
setwd(str_replace(getwd(), "\\/OMST", ""))

mean_cor_data <- merge(data_full_per_subject %>% 
                         filter(threshold == "0.15") %>% 
                         dplyr::select(Subj_ID, Age, Age_group),
                       mean_cor,
                       by = "Subj_ID"
                       ) %>% 
  rename(mean_FC = `fromJSON("mean_cor.json")`)

mod <- lm(mean_FC~Age, data = mean_cor_data)
summary(mod)
performance::check_model(mod)
performance::check_outliers(mod)
effectsize::eta_squared(mod, partial = TRUE)

cor.test(mean_cor_data$Age, mean_cor_data$mean_FC)
# Balance Integration/Segregation
plot(mean_cor_data$Age, mean_cor_data$mean_FC, pch = 19, col = "darkblue")
# Regression line
abline(lm(mean_cor_data$mean_FC ~ mean_cor_data$Age), col = "red", lwd = 3)

# Young adults have a lower overall FC which mean that proportional thresholding select on average
# correlation which have a higher probability of being spurious correlations and artificially inflate the network differences
# (van den Heuvel, 2017)



# setwd(paste0(getwd(), "/OMST"))
# Eglob_OMST <- as.data.frame(fromJSON("Eglob.json")) %>%
#   mutate(Subj_ID = rep(seq_len(628)))
# setwd(str_replace(getwd(), "\\/OMST", ""))
# 
# Eglob_OMST_data <- merge(data_full_per_subject %>%  
#                          dplyr::select(Subj_ID, Age, Age_group),
#                        Eglob_OMST,
#                        by = "Subj_ID"
# ) %>% 
#   rename(Eglob_OMST = `fromJSON("Eglob.json")`)
# 
# mod <- lm(Eglob_OMST~Age, data = Eglob_OMST_data)
# summary(mod)

setwd(paste0(getwd(), "/OMST"))
cost_OMST <- as.data.frame(fromJSON("globalcosteff_max.json")) %>%
  mutate(Subj_ID = rep(seq_len(628))) %>% 
  rename(GE_Cost = `fromJSON("globalcosteff_max.json")`) %>% 
  merge(., data_full_per_subject %>%
          dplyr::select(Subj_ID, Age, Age_group), 
        by = "Subj_ID")

setwd(str_replace(getwd(), "\\/OMST", ""))

mod <- lm(GE_Cost~Age, data = cost_OMST)
summary(mod)
performance::check_outliers(cost_OMST[1])

cor.test(cost_OMST$Age, cost_OMST$GE_Cost)
# Balance Integration/Segregation
plot(cost_OMST$Age, cost_OMST$GE_Cost, pch = 19, col = "darkblue")
# Regression line
abline(lm(cost_OMST$GE_Cost ~ cost_OMST$Age), col = "red", lwd = 3)
text(paste(
  "Correlation between Age and max GE-Cost; t(626) = 5.9, p < .001:",
  round(cor.test(cost_OMST$Age, cost_OMST$GE_Cost)$estimate, 2)
), x = 65, y = 0.45)

