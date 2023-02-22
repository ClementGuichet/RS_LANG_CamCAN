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

################################################################################
# ILR TRANSFORMATION 
################################################################################
Cog_data_ILR <- All_data %>%  
  mutate(ToT_Ratio_inverse = ToT_Ratio*(-1)) %>% 
  mutate(Hotel_Task_inverse = Hotel_Task*(-1)) %>% 
  dplyr::select(-c(ToT_Ratio, Hotel_Task))

################################################################################
# CANONICAL CORRELATION  ANALYSIS
################################################################################
library(CCA)
library(CCP)

Data_CCA <- Cog_data_ILR %>% na.omit()

cog_measures <- Data_CCA[,c(4:9, 26:27)] %>% scale(.) %>% as.data.frame()

rs_measures <- Data_CCA[,23:25] %>% scale(.) %>% as.data.frame()

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

################################################################################
# FIGURES
################################################################################
library(flexplot)

Brain_Mode <- as.matrix(rs_measures) %*% cc_results$xcoef[, 1]
Behavioral_Mode <- as.matrix(cog_measures) %*% cc_results$ycoef[, 1]

cca_df <- Data_CCA %>% 
  mutate(Brain_Mode=Brain_Mode*(-1),
         Behavioral_Mode=Behavioral_Mode)

flexplot(Brain_Mode~Age, cca_df, method = "lm")
flexplot(Behavioral_Mode~Age, cca_df, method = "lm")

cca_df %>% 
  ggplot(aes(x=Brain_Mode, Behavioral_Mode))+
  geom_hline(yintercept =0, color = "gray") +
  geom_vline(xintercept =0, color = "gray") +
  geom_point(aes(colour = Age), size = 4) + 
  geom_smooth(method = "lm", color = "black", alpha = .4, size = 2) +

  labs(x = "Brain Mode", 
       y = "Behavioral Mode\n (Worse - Better performance)",
       title = "Canonical correlation between brain/behavioral modes") + 
  annotate("text", fontface = "bold",
           x = -0.15, y = -0.12, label = "Correlation = .29",
           color = "black", size = 4
  ) +
  theme(plot.title.position = "plot") +
  scale_color_viridis(option = "plasma") +
  theme_classic2(base_size = 18)


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


################################################################################
# MEDIATION ANALYSIS
################################################################################

library(robustbase)
library(robmed)

# MEDIATION

# STEP 1: Direct effect DV ~ IV
mod1 <- lm(Behavioral_Mode~ Age, data = cca_df)
summary(mod1)
flexplot(Behavioral_Mode ~ Age, cca_df, method = "lm")

 # STEP 2: Indirect effect Mediator ~ IV
mod2 <- lm(Brain_Mode~ Age, data = cca_df)
summary(mod2)
flexplot(Brain_Mode~Age, data = cca_df, method = "lm")

# Step 3: DV ~ mediator
mod3 <- lm(Behavioral_Mode~ Brain_Mode, cca_df)
summary(mod3)
flexplot(Behavioral_Mode~ Brain_Mode, cca_df, method = "lm")

# Step 4: DV ~ mediator + controlling for IV
mod4 <- lm(Behavioral_Mode~ Brain_Mode + Age, cca_df)
summary(mod4)
# --> Partial mediation
flexplot(Behavioral_Mode~Age + Brain_Mode, cca_df, method = "lm")

cca_df_med <- cca_df %>% dplyr::select(Behavioral_Mode, Brain_Mode, Age, Balance_eff, Eloc) %>% scale(.) %>% as.data.frame()


psych::mediate(Behavioral_Mode ~(Balance_eff) + Age, data = cca_df_med) %>% summary()
# med <- robmed::test_mediation(Behavioral_Mode~m(Brain_Mode) + Age, data = cca_df_med, robust = "MM")
# summary(med)
library(processR)

labels=list(X="Age",M="Brain_Mode",Y="Behavioral_Mode")
par(mfrow = c(1,1), mar=c(0,0,0,0), oma=c(0,0,0,0))
pmacroModel(4,labels=labels)
statisticalDiagram(4, labels=labels)

library(lavaan)
library(sem)
library(semPlot)

library(processR)

write.csv(cca_df_med, 'file_mediation.csv')
cca_df_med <- cca_df %>% dplyr::select(Behavioral_Mode, Brain_Mode, Age, G1, Eloc) %>% 
  mutate(AgexG1 = Age*G1,
         Brain_ModexG1 = Brain_Mode * G1) %>% scale(.) %>% as.data.frame()

labels=list(X="Age",M1="Brain_Mode", M2 = "Eloc", M3 = "Balance I/S", Y="Behavioral_Mode")
par(mfrow = c(1,1), mar=c(0,0,0,0), oma=c(0,0,0,0))
pmacroModel(6,labels=labels)
statisticalDiagram(6, labels=labels)


mod_bis <- "
# Mediation Effect
Behavioral_Mode ~ b1*Brain_Mode+b2*Eloc+c1*Age
Brain_Mode ~ a1*Age
Eloc ~ a2*Age+d1*Brain_Mode
ind1:=a1*b1
ind2:=a2*b2
secondInd1:=a1*d1*b2
total1:=c1+a1*b1+a2*b2+a1*d1*b2
"

set.seed(2023)
fit <- lavaan::sem(mod_bis, data = cca_df_med, 
                   se = "bootstrap", bootstrap = 1000, 
                   fixed.x=FALSE, meanstructure = TRUE)

summary(fit, standardized = FALSE, fit.measures = TRUE, rsquare = TRUE, ci = TRUE)
parameterEstimates(fit, level = 0.95, boot.ci.type = "bca.simple", standardized = FALSE)


# Representational similarity


All_data_RSN <- merge(participants, CAMCAN_cognitive_data, by = "Observations") %>%
  merge(., CAMCAN_cognitive_data_supp, by = "Observations") %>%
  na.omit() %>%
  merge(., Gradient_stats %>% dplyr::select(Subj_ID, Age, `1st_network`,
                                            Integration, Peripherisation, Polyvalent_interfaces), by = "Subj_ID")

Cog_data_ILR_RSN <- All_data_RSN %>% 
  na.omit() %>% 
  mutate(ToT_Ratio_inverse = ToT_Ratio*(-1)) %>%
  mutate(Hotel_Task_inverse = Hotel_Task*(-1)) %>%
  dplyr::select(-c(ToT_Ratio, Hotel_Task)) %>%
  pivot_longer(c(Integration, Peripherisation, Polyvalent_interfaces),
               names_to = "topological_balances", values_to = "value") %>%
  unite(BalancexRSN, "1st_network", "topological_balances", remove = FALSE) %>%
  dplyr::select(-c(topological_balances, `1st_network`)) %>%
  spread(BalancexRSN, value) 
  

cog_measures_RSN <- Cog_data_ILR_RSN[,c(4:9, 11:12)] %>% scale(.) %>% as.data.frame()

rs_measures_RSN <- Cog_data_ILR_RSN[,13:45] %>% scale(.) %>% as.data.frame()


# pca_loadings <- FactoMineR::PCA(rs_measures_RSN)
# plot(pca_loadings)
# pca_loadings$eig
# contrib <- pca_loadings$var$contrib %>% as.data.frame()
# rs_measures_pca <- pca_loadings$ind$coord[,1:4]


ncors_integration <- 1-cor(cog_measures_RSN, rs_measures_RSN %>% dplyr::select(ends_with("Integration")))
ncors_integration_plot <- ncors_integration >= quantile(ncors_integration, na.rm = TRUE)[3]

library(circlize)

circlize::chordDiagram(ncors_integration_plot)

ncors_peri <- 1-cor(cog_measures_RSN, rs_measures_RSN %>% dplyr::select(ends_with("Peripherisation")))
ncors_peri_plot <- ncors_peri >= quantile(ncors_peri, na.rm = TRUE)[3]

library(circlize)

circlize::chordDiagram(ncors_peri_plot)
