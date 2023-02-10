##########################################################################################
# Script for Brain-Cognition analyses

# Written by CG - 2023
##########################################################################################
source("_06_TFProfile_CamCAN.R")

# 2 RECONFIGURATION MECHANISM
# Integration within module and integration between modules

# MECHANISM WITHIN:  Satellite reconfigure into Connectors or Peripheral, meaning they either get integrated into a module or left on their own
# MECHANISM BETWEEN: Half of Provincial hubs reconfigure into Connector

geometric_all <- TFP_General %>% 
  summarize_at(vars(Connector:Super_Bridge), funs(geomMeanExtension(., epsilon = epsilon)))

TFP_General %>%  
  mutate(Connector = log(Connector / geometric_all$Connector)) %>%
  mutate(Provincial = log(Provincial / geometric_all$Provincial)) %>%
  mutate(Satellite = log(Satellite / geometric_all$Satellite)) %>%
  mutate(Peripheral = log(Peripheral / geometric_all$Peripheral)) %>%
  dplyr::select(Subj_ID, Age, Connector:Satellite) %>%
  pivot_longer(
    cols = !c("Subj_ID", "Age"),
    names_to = "Functional_role",
    values_to = "Score"
  ) %>%
  ggplot(aes(Age, Score, color = Functional_role)) +
  geom_hline(yintercept = 0, color = "red") +
  geom_jitter(height = 0.05, alpha = 0.1) +
  geom_smooth(linewidth = 2, method = "gam", alpha = .3) +
  scale_x_continuous(breaks = seq(20, 90, 5)) +
  scale_y_continuous(breaks = seq(-0.2, 0.2, 0.1)) +
  coord_cartesian(ylim = c(-0.2, 0.2)) +
  scale_color_brewer(palette = "PuOr") +
  geom_vline(xintercept = 52, color = "red", linewidth = 1.5, alpha = 1) +
  theme_pubr() +
  ggtitle("Evolution of modular functional roles across adult lifespan")

TFP_General %>%  
  mutate(Global_Bridge = log(Global_Bridge / geometric_all$Global_Bridge)) %>%
  mutate(Local_Bridge = log(Local_Bridge / geometric_all$Local_Bridge)) %>%
  mutate(Super_Bridge = log(Super_Bridge / geometric_all$Super_Bridge)) %>%
  mutate(Not_a_Bridge = log(Not_a_Bridge / geometric_all$Not_a_Bridge)) %>%
  dplyr::select(Subj_ID, Age, Global_Bridge:Super_Bridge) %>%
  pivot_longer(
    cols = !c("Subj_ID", "Age"),
    names_to = "Functional_role",
    values_to = "Score"
  ) %>%
  ggplot(aes(Age, Score, color = Functional_role)) +
  geom_hline(yintercept = 0, color = "red") +
  geom_jitter(height = 0.05, alpha = 0.1) +
  geom_smooth(linewidth = 2, method = "gam", alpha = .3) +
  scale_x_continuous(breaks = seq(20, 90, 5)) +
  scale_y_continuous(breaks = seq(-0.2, 0.2, 0.1)) +
  coord_cartesian(ylim = c(-0.2, 0.2)) +
  scale_color_brewer(palette = "PuOr") +
  geom_vline(xintercept = 52, color = "red", linewidth = 1.5, alpha = 1) +
  theme_pubr() +
  ggtitle("Evolution of internodal functional roles across adult lifespan") +
  facet_wrap(~Functional_role)


################################################################################
list_TFP_RSN <- TFP_RSN %>% 
  group_by(`1st_network`) %>% group_split()
  
list_tmp <- list()
for (i in 1:length(list_TFP_RSN)) {
  library(compositions)
  
  tmp_raw <- rbindlist(list_TFP_RSN[i]) %>% arrange(Subj_ID) 
  
  tmp_coda_modular <- tmp_raw %>%
    dplyr::select(Connector, Satellite, Provincial, Peripheral) %>% 
    acomp(.) %>% 
    # Preserves the ratios between non-zero components
    cmultRepl(., output = "prop")
  
  tmp_coda_interareal <- tmp_raw %>%
    dplyr::select(Global_Bridge, Local_Bridge, Super_Bridge, Not_a_Bridge) %>%
    acomp(.) %>%
    cmultRepl(., output = "prop")
  
  tmp_geometric_all <- tmp_raw_imputed %>% 
    summarize_at(vars(Connector:Not_a_Bridge), funs(geomMeanExtension(., epsilon = epsilon)))
  
  tmp_raw_imputed <- cbind(tmp_raw %>% dplyr::select(Subj_ID, `1st_network`, Age),
                               tmp_coda_modular, 
                               tmp_coda_interareal)
  
  tmp_final <- tmp_raw_imputed %>% 
    mutate(Connector = log(Connector / tmp_geometric_all$Connector)) %>%
    mutate(Provincial = log(Provincial / tmp_geometric_all$Provincial)) %>%
    mutate(Satellite = log(Satellite / tmp_geometric_all$Satellite)) %>%
    mutate(Peripheral = log(Peripheral / tmp_geometric_all$Peripheral)) %>% 
    mutate(Global_Bridge = log(Global_Bridge / tmp_geometric_all$Global_Bridge)) %>%
    mutate(Local_Bridge = log(Local_Bridge / tmp_geometric_all$Local_Bridge)) %>%
    mutate(Super_Bridge = log(Super_Bridge / tmp_geometric_all$Super_Bridge)) %>%
    mutate(Not_a_Bridge = log(Not_a_Bridge / tmp_geometric_all$Not_a_Bridge)) 
    
  list_tmp[[i]] <- tmp_final
}

plot_TFP_RSN <- rbindlist(list_tmp)


RSN <- "SMN"

plot_TFP_RSN %>% 
  filter(grepl(RSN, `1st_network`)) %>% 
  dplyr::select(Subj_ID, Age, Connector:Peripheral) %>%
  pivot_longer(
    cols = !c("Subj_ID", "Age"),
    names_to = "Functional_role",
    values_to = "Score"
  ) %>%
  ggplot(aes(Age, Score, color = Functional_role)) +
  geom_hline(yintercept = 0, color = "red") +
  geom_jitter(height = 0.05, alpha = 0.1) +
  geom_smooth(linewidth = 2, method = "gam", alpha = .3) +
  scale_x_continuous(breaks = seq(20, 90, 5)) +
  scale_y_continuous(breaks = seq(-1, 1, 0.1)) +
  coord_cartesian(ylim = c(-1, 1)) +
  scale_color_brewer(palette = "PuOr") +
  geom_vline(xintercept = 52, color = "red", linewidth = 1.5, alpha = 1) +
  theme_pubr() +
  ggtitle("Evolution of modular functional roles across adult lifespan")
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
  dplyr::rename(ToT_Ratio = TOT_Summary_ToT_ratio)

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
  merge(., data_TFP_analysis, by = "Subj_ID")

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

Data_CCA <- Cog_data_ILR %>% na.omit() %>% 
  dplyr::select(Cattell, ToT_Ratio_inverse, Hotel_Task_inverse,
                Sentence_Comprehension_c, Verbal_Fluency, Proverb, Story_Recall,
                Picture_Priming,
                
                ilr_modular_all_1, ilr_modular_all_2, ilr_internodal_all_1, ilr_internodal_all_2,
                Age
                )

cog_measures <- Data_CCA[,1:8] %>% scale(.) %>% as.data.frame()

rs_measures <- Data_CCA[,9:11] %>% scale(.) %>% as.data.frame()

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


p <- cca_df %>% 
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

p
# plotly::ggplotly(p)

flexplot(Behavioral_Mode~Brain_Mode + Age_group, cca_df, method = "lm")


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
  rename(c("V3" = "correlation")) %>% 
  filter(correlation != 1) %>% 
  arrange(correlation)

correlations %>% 
  head(30) %>% 
  graph_from_data_frame() %>% 
  ggraph() +
  geom_edge_link(aes(edge_alpha = correlation)) +
  geom_node_point() +
  geom_node_text(aes(label = name), repel = TRUE) +
  theme_pubclean()

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
flexplot(Behavioral_Mode~Brain_Mode + Age, cca_df, method = "lm")


cca_df_med <- cca_df %>% dplyr::select(Behavioral_Mode, Brain_Mode, Age) %>% scale(.) %>% as.data.frame()

psych::mediate(Behavioral_Mode ~(Brain_Mode) + Age, data = cca_df_med) %>% summary()
med <- robmed::test_mediation(Behavioral_Mode~m(Brain_Mode) + Age, data = cca_df_med, robust = "MM")
summary(med)

library(lavaan)
library(sem)
library(semPlot)

mod <- "
  Behavioral_Mode ~ dir*Age
  Brain_Mode ~ ind1*Age
  Behavioral_Mode ~ ind2*Brain_Mode
  mediation := ind1*ind2
  direct := dir
  total := (ind1*ind2) + dir
"

fit <- lavaan::sem(mod, data = cca_df_med)

summary(fit, standardized = FALSE, fit.measures = TRUE, rsquare = TRUE, ci = TRUE)
parameterEstimates(fit, level = 0.95, boot.ci.type = "bca.simple", standardized = FALSE)
# semPaths(fit, "std", layout = "tree2", edge.label.cex = 1.25, fade = TRUE)
