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
  dplyr::rename(Picture_Priming = Picture__Primming_Summary_ACC_baseline_all) %>%
  dplyr::rename(ToT_Ratio = TOT_Summary_ToT_ratio) %>% 
  mutate_at(vars(MMSE), funs(as.numeric(.)))

CAMCAN_cognitive_data_supp <- read_excel("meta_data_628/CognitiveData_CamCAN_Supplement.xlsx") %>% 
  filter(Observations %in% participants$Observations) %>%
  dplyr::select(c(
    Observations,
    Hotel_Task,
    # Sentence Comprehension
    Sentence_Comprehension_c,
    Verbal_Fluency,
    # Memory
    Story_Recall
  )) %>% mutate_at(vars(Hotel_Task, Sentence_Comprehension_c, Verbal_Fluency, Story_Recall), funs(as.numeric(.)))

All_data <- merge(participants, CAMCAN_cognitive_data, by = "Observations") %>% 
  merge(., CAMCAN_cognitive_data_supp, by = "Observations") %>% 
  merge(., TFP_General, by = "Subj_ID")

################################################################################
# ILR TRANSFORMATION 
################################################################################
library(zCompositions)
library(compositions)

data_coda_modular <- All_data %>%
  dplyr::select(Connector, Satellite, Provincial, Peripheral) %>% 
  acomp(.) %>% 
  # Preserves the ratios between non-zero components
  cmultRepl(., output = "prop")

data_coda_interareal <- All_data %>%
  dplyr::select(Global_Bridge, Local_Bridge, Super_Bridge, Not_a_Bridge) %>%
  acomp(.) %>% 
  cmultRepl(., output = "prop")

data_coda_all <- cbind(data_coda_modular, data_coda_interareal,
                       All_data %>% dplyr::select(-c(Connector, Satellite, Provincial, Peripheral,
                                                    Global_Bridge, Local_Bridge, Super_Bridge, Not_a_Bridge)))

Cog_data_ILR <- data_coda_all %>%  
  
  mutate(ToT_Ratio_inverse = ToT_Ratio*(-1)) %>% 
  mutate(Hotel_Task_inverse = Hotel_Task*(-1)) %>%
  
  mutate(ilr_modular_all_1 = (((4/4)^0.5)*log((Connector*Peripheral)^0.5/(Provincial*Satellite)^0.5))) %>% 
  mutate(ilr_modular_all_2 = (((4/4)^0.5)*log((Connector*Satellite)^0.5/(Peripheral*Provincial)^0.5))) %>% 
  
  # Inter-modular integration
  mutate(ilr_modular_1 = (((1/2)^0.5)*log(Connector/(Provincial)))) %>%

  # Intra-module integration
  mutate(ilr_modular_2 = (((2/3)^0.5)*log((Connector*Peripheral)^0.5/(Satellite)))) %>% 
  
  mutate(ilr_internodal_all_1 = (((3/4)^0.5)*log((Global_Bridge)/((Super_Bridge*Not_a_Bridge*Local_Bridge)^(1/3))))) %>%  
  mutate(ilr_internodal_all_2 = (((2/3)^0.5)*log((Not_a_Bridge)/((Super_Bridge*Local_Bridge)^(1/2))))) 

  
# cor_data <- Cog_data_ILR %>% dplyr::select(ilr_modular_1, ilr_modular_2, ilr_internodal,
#                                       # Cattell Fluid intelligence
#                                       Cattell,
#                                       # Hotel taks for executive functions (EF)
#                                       Hotel_Task,
#                                       # Proverb comprehension (abstraction & EF)
#                                       Proverb,
#                                       # Picture-picture priming (word production)
#                                       Picture_Priming,
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

# lm(ilr_internodal~Age, Cog_data_ILR) %>% summary()
# flexplot(ilr_internodal~Age, Cog_data_ILR, method = "lm")


################################################################################
# CANONICAL CORRELATION  ANALYSIS
################################################################################
library(CCA)
library(CCP)
library(flexplot)

Data_CCA <- Cog_data_ILR %>% na.omit() %>% 
  dplyr::select(MMSE, Cattell, ToT_Ratio_inverse, Hotel_Task_inverse,
                Sentence_Comprehension_c, Verbal_Fluency, Proverb, Story_Recall,
                Picture_Priming,
                
                ilr_modular_all_1, ilr_modular_all_2, ilr_internodal_all_1, ilr_internodal_all_2,
                Age
                )

cog_measures <- Data_CCA[,1:9] %>% scale(.) %>% as.data.frame()

rs_measures <- Data_CCA[,c(10, 12)] %>% scale(.) %>% as.data.frame()

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

flexplot(Behavioral_Mode~Brain_Mode, cca_df, method = "lm")
flexplot(Behavioral_Mode~Age, cca_df, method = "lm")
flexplot(Brain_Mode~Age, cca_df, method = "lm")


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
           x = -0.075, y = -0.13, label = "Correlation = .32\n Wilks'Lambda: 0.86, F(32, 2074) = 2.65, p < .001",
           color = "black", size = 4
  ) +
  theme(plot.title.position = "plot")

# Canonical loadings - correlation between variables and canonical variates
cc_loadings <- cc(rs_measures, cog_measures)
cc_results2 <- comput(rs_measures, cog_measures, cc_loadings)
cc_results2[3:6]


# tests of canonical dimensions
rho <- cc_loadings$cor
## Define number of observations, number of variables in first set, and number of variables in the second set.
n <- dim(rs_measures)[1]
p <- length(rs_measures)
q <- length(cog_measures)

## Calculate p-values using the F-approximations of different test statistics:
p.asym(rho, n, p, q, tstat = "Wilks")


library(widyr)
library(ggraph)
library(igraph)

adjacency_to_2col <- function(data) {
  crossdata <- lapply(rownames(data), function(x) sapply(colnames(data), function(y) list(x, y, data[x, y])))
  crossdatatmp <- matrix(unlist(crossdata), nrow = 3)
  crossdatamat <- t(crossdatatmp)
  crossdatadf <- as.data.frame(crossdatamat, stringsAsFactors = F)
  crossdatadf[, 3] <- as.numeric(crossdatadf[, 3])
  return(crossdatadf %>% na.omit())
}

correlations <- 
  desc_cca$XYcor %>% as.data.frame() %>% 
  adjacency_to_2col(.) %>% 
  plyr::rename(c("V3" = "correlation")) %>% 
  filter(correlation != 1) %>% 
  mutate(corr_abs = abs(correlation)) %>% 
  arrange(desc(corr_abs))

correlations %>% 
  head(length(correlations$corr_abs)*.8) %>% 
  graph_from_data_frame(directed = TRUE) %>% 
  ggraph() +
  geom_edge_link(aes(edge_alpha = correlation)) +
  geom_node_point() +
  geom_node_text(aes(label = name), repel = TRUE) +
  theme_pubclean()

################################################################################
# MEDIATION ANALYSIS
################################################################################

library(robustbase)
library(robmed)

# mod1 <- lmrob(scale(Picture_Priming)~ scale(ilr_modular_1)*scale(Age), data = Cog_data_ILR, method = "MM")
# summary(mod1)
# 
# plot_bin <- Cog_data_ILR %>% 
#   mutate(Ratio = Connector/Provincial) %>% 
#   mutate(
#   bin_modular = ifelse(ilr_modular_1 < quantile(.$ilr_modular_1)[3], "+Provincial","+Connector")
# )
# flexplot(Picture_Priming~ Age + bin_modular, plot_bin, method = "lm")
# 
# # Individuals who are older have significantly worse performances in the Picture priming task (i.e., word production). 
# # Moreover, there was a significant interaction with the modular balance such that having more Connector than Provincial hubs
# # past a certain ratio is associated with worse performances after 60 years old. In contrasts, older subjects who retain a normal
# # balance show a slower decline.
# # 
# # In short, functional disorganization of the language connectome is linked with a faster decline in word production across the lifespan
# 
# mod1 <- lmrob(scale(Picture_Priming)~ scale(ilr_modular_2)*scale(Age), data = Cog_data_ILR, method = "MM")
# summary(mod1)
# 
# flexplot(Picture_Priming~ Age + ilr_modular_2, Cog_data_ILR, method = "lm")
# 
# 
# # Same thing but point of inflexion is at 40 rather than at 60
# 
# 


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
med <- robmed::test_mediation(Behavioral_Mode~m(Brain_Mode) + Age, data = cca_df_med, robust = "MM")
summary(med)

library(lavaan)
library(sem)
library(semPlot)

mod <- "
  Brain_Mode ~ a*Age
  Behavioral_Mode ~ c1*Age + b*Brain_Mode
  indirect := a*b
  direct := c1
  total := direct + indirect
"

fit <- lavaan::sem(mod, data = cca_df_med)

summary(fit, standardized = FALSE, fit.measures = TRUE, rsquare = TRUE, ci = TRUE)
parameterEstimates(fit, level = 0.95, boot.ci.type = "bca.simple", standardized = FALSE)
semPaths(fit, "std", layout = "tree", edge.label.cex = 1.25, fade = TRUE)

mod_moderated <- "
  Behavioral_Mode ~ c1*Age + b*Brain_Mode  + w1*Age:Brain_Mode
  direct := c1
  total := direct + b + w1
"

fit_moderated <- lavaan::sem(mod_moderated, data = cca_df_med)

mod__mediated_moderated <- "
  Brain_Mode ~ a*Age
  Behavioral_Mode ~ c1*Age + b*Brain_Mode + w1*Age:Brain_Mode
  indirect := a*b
  indirect_moderated := indirect + w1
  direct := c1
  total := direct + indirect
  total_moderated := total + w1
"

fit_mediated_moderated <- lavaan::sem(mod__mediated_moderated, data = cca_df_med)

anova(fit, fit_moderated, fit_mediated_moderated)


