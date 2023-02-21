##########################################################################################
# Script for Brain-Cognition analyses

# Written by CG - 2023
##########################################################################################
source("_06_TFProfile_CamCAN.R")

################################################################################
# DESCRIPTIVES
################################################################################

participants <- read_excel("meta_data_628/participant_data_T1.xlsx")[, 1:2] %>% 
  dplyr::rename(Subj_ID = Subject) %>% 
  replace("Subj_ID", seq_len(628)) 

CAMCAN_cognitive_data <- read_excel("meta_data_628/CognitiveData_CamCAN_Apr2022.xlsx") %>% 
  filter(Observations %in% participants$Observations) %>%
  dplyr::select(-gender_code) %>%
  dplyr::rename(Age_CogData = Age) %>% 
  dplyr::select(c(
    Observations,
    Age_CogData,
    MMSE,
    # Cattell Fluid intelligence
    Cattell,
    # Proverb comprehension (abstraction & EF)
    Proverbs_Summary__Score,
    # Picture-picture priming (word production)
    Picture__Primming_Summary_ACC_baseline_all,
    # Tip-of-the-tongue
    TOT_Summary_ToT_ratio
  )) %>%
  dplyr::rename(Proverb = Proverbs_Summary__Score) %>% 
  dplyr::rename(Naming = Picture__Primming_Summary_ACC_baseline_all) %>%
  dplyr::rename(ToT_Ratio = TOT_Summary_ToT_ratio) %>% 
  mutate_at(vars(MMSE), funs(as.numeric(.)))

CAMCAN_cognitive_data_supp <- read_excel("meta_data_628/CognitiveData_CamCAN_Supplement.xlsx") %>% 
  filter(Observations %in% participants$Observations) %>%
  dplyr::select(c(
    Observations,
    Hotel_Task,
    # Sentence Comprehension
    Sentence_Comprehension_c,
    # Memory
    Story_Recall
  )) %>% mutate_at(vars(Hotel_Task, Sentence_Comprehension_c, Story_Recall), funs(as.numeric(.)))

All_data <- merge(participants, CAMCAN_cognitive_data, by = "Observations") %>% 
  merge(., CAMCAN_cognitive_data_supp, by = "Observations") %>% 
  na.omit() %>% 
  merge(., TFP_General_imputed_stats, by = "Subj_ID") 

# All_data <- merge(participants, CAMCAN_cognitive_data, by = "Observations") %>% 
#   merge(., CAMCAN_cognitive_data_supp, by = "Observations") %>% 
#   na.omit() %>% 
#   merge(., Gradient_stats %>% dplyr::select(Subj_ID, Age, Age_group, `1st_network`,
#                                             ilr_modular_inter, ilr_modular_intra, ilr_internodal), by = "Subj_ID")

################################################################################
# ILR TRANSFORMATION 
################################################################################
Cog_data_ILR <- All_data %>%  
  mutate(ToT_Ratio_inverse = ToT_Ratio*(-1)) %>% 
  mutate(Hotel_Task_inverse = Hotel_Task*(-1)) %>% 
  dplyr::select(-c(ToT_Ratio, Hotel_Task))

# Cog_data_ILR <- All_data %>%  
#   mutate(ToT_Ratio_inverse = ToT_Ratio*(-1)) %>% 
#   mutate(Hotel_Task_inverse = Hotel_Task*(-1)) %>% 
#   dplyr::select(-c(ToT_Ratio, Hotel_Task)) %>% 
#   pivot_longer(c(ilr_modular_inter, ilr_modular_intra, ilr_internodal),
#                names_to = "balances", values_to = "ilr") %>% 
#   unite(BalancexRSN, "1st_network", "balances", remove = FALSE) %>% 
#   dplyr::select(-c(balances, `1st_network`)) %>% 
#   spread(BalancexRSN, ilr)

# cor_data <- Cog_data_ILR %>% dplyr::select(ilr_modular_1, ilr_modular_2, ilr_internodal,
#                                       # Cattell Fluid intelligence
#                                       Cattell,
#                                       # Hotel taks for executive functions (EF)
#                                       Hotel_Task,
#                                       # Proverb comprehension (abstraction & EF)
#                                       Proverb,
#                                       # Picture-picture priming (word production)
#                                       Naming,
#                                       # Sentence comprehension
#                                       Sentence_Comprehension_c,
#                                       # Tip-of-the-tongue
#                                       ToT_Ratio,
#                                       Story_Recall,
#                                       Verbal_Fluency,
#                                       Age) %>% na.omit()
# 
# corrplot::corrplot(cor(cor_data))
#
# library(psych)
#
# test <- psych::corr.test(cor_data, adjust = "fdr", ci = T, method = "pearson")
# test
# print(test, short = FALSE)


################################################################################
# CANONICAL CORRELATION  ANALYSIS
################################################################################
library(CCA)
library(CCP)

Data_CCA <- Cog_data_ILR %>% na.omit()

cog_measures <- Data_CCA[,c(4:9, 26:27)] %>% scale(.) %>% as.data.frame()

rs_measures <- Data_CCA[,23:25] %>% scale(.) %>% as.data.frame()

# pca_loadings <- FactoMineR::PCA(rs_measures)
# plot(pca_loadings)
# pca_loadings$eig
# contrib <- pca_loadings$var$contrib %>% as.data.frame()
# rs_measures_pca <- pca_loadings$ind$coord[,1:4]

Data_CCA %>% pivot_longer(
  MMSE:Story_Recall,
  names_to = "Cognitive_assessment",
  values_to = "performance"
) %>% 
  group_by(Cognitive_assessment) %>% 
  mutate(performance = as.numeric(scale(performance))) %>%  
  ggplot(aes(Age, performance, color = Cognitive_assessment)) +
  geom_hline(yintercept = 0, color = "red") +
  geom_jitter(height = 0.05, alpha = 0.1) +
  geom_smooth(linewidth = 2, method = "gam", alpha = .3) +
  scale_x_continuous(breaks = seq(20, 90, 15)) +
  coord_cartesian(ylim = c(-1, 1)) +
  scale_color_brewer(palette = "Paired") +
  theme_pubclean(base_size = 18) +
  theme(plot.title.position = "plot") +
  labs(y = "Normalized score") +
  ggtitle("Cognitive performances across the adult lifespan") +
  facet_wrap(~Cognitive_assessment, scale = "free")

# Relationship between topological integrative mechanisms and Efficiency is fully mediated by Age
# med <- robmed::test_mediation(Balance_eff~m(Age) + ilr_modular_1, data = Data_CCA, robust = "MM")
# summary(med)

desc_cca <- matcor(rs_measures, cog_measures)
img.matcor(desc_cca, type = 2)

cc_results <- cancor(rs_measures, cog_measures)
cc_results$cor


Brain_Mode <- as.matrix(rs_measures) %*% cc_results$xcoef[, 1]
Behavioral_Mode <- as.matrix(cog_measures) %*% cc_results$ycoef[, 1]


cca_df <- Data_CCA %>% 
  mutate(Brain_Mode=Brain_Mode*(-1),
         Behavioral_Mode=Behavioral_Mode)

RSA <- cca_df %>% 
  mutate(Age_binned = Hmisc::cut2(Age, g = 5, levels.mean = TRUE)) %>% 
  group_by(Age_binned) %>% group_split() %>% 
  map_dfr(. %>% 
            mutate(cor_res = cor.test(Brain_Mode, Behavioral_Mode)$estimate)) %>% 
  dplyr::select(Age_binned, cor_res, Brain_Mode, Behavioral_Mode) %>% 
  mutate(Similarity = -1*cor_res) %>% 
  distinct(Age_binned, .keep_all = TRUE)

RSA %>% 
  ggplot(aes(as.numeric(Age_binned), Similarity)) + 
  geom_line(linewidth = 2) +
  scale_x_continuous(breaks = seq(20, 90, 10)) +
  scale_y_continuous(breaks = seq(0, 0.2, 0.05)) +
  ggtitle("Similarity between canonical variates across the adult lifespan") +
  labs(x = "Age") +
  theme_classic2(base_size = 18)
  


flexplot(Behavioral_Mode~Brain_Mode + Age, cca_df, method = "lm")
flexplot(Behavioral_Mode~Age, cca_df, method = "lm")
flexplot(Brain_Mode~Age, cca_df, method = "lm")

# Canonical loadings - correlation between variables and canonical variates
cc_loadings <- cc(rs_measures, cog_measures)
cc_results2 <- comput(rs_measures, cog_measures, cc_loadings)
cc_results2[3:6]

# tests of canonical dimensions
rho <- cc_loadings$cor
# Define number of observations, number of variables in first set, and number of variables in the second set.
n <- dim(rs_measures)[1]
p <- length(rs_measures)
q <- length(cog_measures)

# Calculate p-values using the F-approximations of different test statistics:
p.asym(rho, n, p, q, tstat = "Wilks")

cca_df %>% 
  ggplot(aes(x=Brain_Mode, Behavioral_Mode))+
  geom_hline(yintercept =0, color = "gray") +
  geom_vline(xintercept =0, color = "gray") +
  geom_point(aes(colour = Age), size = 4) + 
  geom_smooth(method = "lm", color = "black", alpha = .4, size = 2) +
  scale_color_viridis(option = "plasma") +
  theme_classic2(base_size = 18) +
  
  labs(x = "Brain Mode\n (Functional Segregation - Integration)", y = "Behavioral Mode\n (Worse - Better performance)",
       title = "Canonical correlation between brain/behavioral modes") + 
  annotate("text", fontface = "bold",
           x = -0.07, y = -0.13, label = "Correlation = .32\n Wilks'Lambda: 0.88, F(24, 1631) = 3.18, p < .001",
           color = "black", size = 4
  ) +
  theme(plot.title.position = "plot")



plot_loadings_x <- cc_loadings$scores$corr.X.yscores[,1] %>% as.data.frame() %>% 
  rownames_to_column("Balances") %>% 
  plyr::rename(c("." = "loading")) %>% 
  arrange(loading)

plot_loadings_x$Balances <- factor(plot_loadings_x$Balances) %>%
  fct_reorder(plot_loadings_x$loading, .desc = FALSE)

loadings_x <- ggplot(plot_loadings_x, aes(x = Balances, y = loading))+
  geom_col(aes(fill = loading), alpha = .8) +
  scale_fill_distiller(palette = "RdBu", direction = -1) +
  coord_flip() +
  labs(y = "Loadings", x = element_blank()) + 
  theme_pubr(base_size = 18)

plot_loadings_y <- cc_loadings$scores$corr.Y.yscores[,1] %>% as.data.frame() %>% 
  rownames_to_column("Cognitive_assessment") %>% 
  plyr::rename(c("." = "loading")) %>% 
  arrange(loading)

plot_loadings_y$Cognitive_assessment <- factor(plot_loadings_y$Cognitive_assessment) %>%
  fct_reorder(plot_loadings_y$loading, .desc = TRUE)

loadings_y <- ggplot(plot_loadings_y, aes(x = Cognitive_assessment, y = loading))+
  geom_col(aes(fill = loading), alpha = .8) +
  scale_fill_distiller(palette = "RdBu", direction = 1) +
  coord_flip() +
  labs(y = "Loadings", x = element_blank()) + 
  theme_pubr(base_size = 18)

gridExtra::grid.arrange(loadings_x, loadings_y, nrow = 1)
# 
# library(widyr)
# library(ggraph)
# library(igraph)
# 
# adjacency_to_2col <- function(data) {
#   crossdata <- lapply(rownames(data), function(x) sapply(colnames(data), function(y) list(x, y, data[x, y])))
#   crossdatatmp <- matrix(unlist(crossdata), nrow = 3)
#   crossdatamat <- t(crossdatatmp)
#   crossdatadf <- as.data.frame(crossdatamat, stringsAsFactors = F)
#   crossdatadf[, 3] <- as.numeric(crossdatadf[, 3])
#   return(crossdatadf %>% na.omit())
# }
# 
# correlations <- 
#   desc_cca$XYcor %>% as.data.frame() %>% 
#   adjacency_to_2col(.) %>% 
#   plyr::rename(c("V3" = "correlation")) %>% 
#   filter(correlation != 1) %>% 
#   mutate(corr_abs = abs(correlation)) %>% 
#   arrange(desc(corr_abs))
# 
# correlations %>% 
#   head(length(correlations$corr_abs)*.8) %>% 
#   graph_from_data_frame(directed = TRUE) %>% 
#   ggraph(layout = "auto") +
#   geom_edge_link(aes(edge_alpha = correlation)) +
#   geom_node_point() +
#   geom_node_text(aes(label = name), repel = TRUE) +
#   theme_pubclean()

################################################################################
# MEDIATION ANALYSIS
################################################################################

library(robustbase)
library(robmed)

# MEDIATION

# STEP 1: Direct effect DV ~ IV
mod1 <- lm(Behavioral_Mode~ Age, data = cca_df)
summary(mod1)
performance::check_model(mod1)
flexplot(Behavioral_Mode ~ Age, cca_df, method = "lm")

 # STEP 2: Indirect effect Mediator ~ IV
mod2 <- lm(Brain_Mode~ Age, data = cca_df)
summary(mod2)
performance::check_model(mod2)
flexplot(Brain_Mode~Age, data = cca_df, method = "lm")

# Step 3: DV ~ mediator
mod3 <- lm(Behavioral_Mode~ Brain_Mode, cca_df)
summary(mod3)
performance::check_model(mod3)
flexplot(Behavioral_Mode~ Brain_Mode, cca_df, method = "lm")

# Step 4: DV ~ mediator + controlling for IV
mod4 <- lm(Behavioral_Mode~ Brain_Mode + Age, cca_df)
summary(mod4)
performance::check_model(mod4)
# --> Partial mediation
flexplot(Behavioral_Mode~Age + Brain_Mode, cca_df, method = "lm")


cca_df_med <- cca_df %>% dplyr::select(Behavioral_Mode, Brain_Mode, Age) %>% scale(.) %>% as.data.frame()

psych::mediate(Behavioral_Mode ~(Brain_Mode) + Age, data = cca_df_med) %>% summary()
# med <- robmed::test_mediation(Behavioral_Mode~m(Brain_Mode) + Age, data = cca_df_med, robust = "MM")
# summary(med)

library(lavaan)
library(sem)
library(semPlot)


library(processR)

cca_df_med <- cca_df %>% dplyr::select(Behavioral_Mode, Brain_Mode, Age, G1) %>% 
  mutate(AgexG1 = Age*G1,
         Brain_ModexG1 = Brain_Mode * G1) %>% scale(.) %>% as.data.frame()

labels=list(X="Age",M="Brain_Mode",Y="Behavioral_Mode",W="G1")
par(mfrow = c(2,1), mar=c(0,0,0,0), oma=c(0,0,0,0))
pmacroModel(7,labels=labels)
statisticalDiagram(7, labels=labels)

quantile(cca_df_med$G1, probs = c(0.16, 0.5, 0.84))
# 16%        50%        84% 
# -1.2560748  0.3430827  1.0692292 

mod <- "
  Behavioral_Mode ~ cprime*Age + b*Brain_Mode
  Brain_Mode ~ a1*Age + a2*G1 + a3*AgexG1
  indirect := a3*b
  direct := cprime
  total := direct + indirect
"

MedMod <- '
Behavioral_Mode ~ b*Brain_Mode + cprime*Age
Brain_Mode ~ a1*Age + a2*G1 + a3*AgexG1

# Define simple slopes and conditional indirect effects using :=

IndMedMod:= a3*b

# simple slope of Brain_Mode on Age is a1+a3*G1 (Hayes et al., 2017)
aLow: = a1+a3*(-1.26)
aMedian: = a1+a3*0.34
aMean: = a1+a3*(0)
aHigh: = a1+a3*1.07


# conditional indirect effects is b*(a1+a3*G1)
abLow: = b*aLow
abMedian: = b*aMedian
abMean: = b*aMean
abHigh: = b*aHigh
'

set.seed(2000)
fit <- lavaan::sem(MedMod, data = cca_df_med, 
                   se = "bootstrap", bootstrap = 100, 
                   fixed.x=FALSE, meanstructure = TRUE)

summary(fit, standardized = FALSE, fit.measures = TRUE, rsquare = TRUE, ci = TRUE)
parameterEstimates(fit, level = 0.95, boot.ci.type = "bca.simple", standardized = FALSE)
