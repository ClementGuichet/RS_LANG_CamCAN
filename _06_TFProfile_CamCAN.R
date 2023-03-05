##########################################################################################
# Script for Age-related analyses and visualization for the topologico-functional profile

# Written by CG - 2023
##########################################################################################
source("_05_Hub_detection_CamCAN.R")
source("_radarplotting_function.R")
source("_geometricmeanCruz.R")

################################################################################
# Investigating the evolution of graph-based metrics ----
################################################################################

############################################################################
# Global Disruption
############################################################################

mod <- lm(Balance_eff ~ Age, TFP_General)
summary(mod)
effectsize::eta_squared(mod)

TFP_General %>% 
  pivot_longer(c(Eglob, Eloc, Balance_eff), names_to = "Efficiencies", values_to = "value") %>% 
  group_by(Efficiencies) %>% 
  mutate(value = as.numeric(scale(value))) %>% 
  ggplot(aes(Age, value, color = Efficiencies)) +
  geom_hline(yintercept = 0, color = "red") +
  geom_vline(xintercept = 54, color = "red", linewidth = 1.5, alpha = 1) +
  geom_jitter(height = 0.05, alpha = 0.05, size = 4) +
  geom_smooth(linewidth = 2, method = "lm", alpha = .1) +
  scale_x_continuous(breaks = seq(20, 90, 5)) +
  scale_y_continuous(breaks = seq(-0.5, 0.7, 0.2)) +
  coord_cartesian(ylim = c(-0.5, 0.7)) +
  scale_color_brewer(palette = "Dark2", 
                     labels = c("Balance I/S", "Global efficiency", "Local efficiency")) +
  theme_pubclean(base_size = 18) +
  theme(plot.title.position = "plot") +
  labs(y = "Normalized efficiencies") +
  ggtitle("Functional network dynamics across the adult lifespan", 
          subtitle = "Inflection point: 54 yo")


############################################################################
# Nodal Disruption
############################################################################

############################################################################
# Bayesian non-parametric multiplicative replacement
############################################################################
library(zCompositions)
library(compositions)

tmp_coda_modular <- TFP_General %>%
  dplyr::select(Connector, Satellite, Provincial, Peripheral) %>% 
  acomp(.)
  # Preserves the ratios between non-zero components
  # cmultRepl(., output = "prop")

tmp_coda_internodal <- TFP_General %>%
  dplyr::select(Global_Bridge, Local_Bridge, Super_Bridge, Not_a_Bridge) %>%
  acomp(.) %>%
  cmultRepl(., output = "prop")

TFP_General_imputed <- cbind(TFP_General %>% dplyr::select(-c(Connector, Satellite, Provincial, Peripheral,
                                                              Global_Bridge, Local_Bridge, Super_Bridge, Not_a_Bridge)),
                         tmp_coda_modular, 
                         tmp_coda_internodal)

############################################################################
# Standardized geometric mean across all individuals for each topological role
############################################################################

# For modular roles

geometric_all <- TFP_General_imputed %>% 
  summarize_at(vars(Connector:Not_a_Bridge), funs(geometricmean(.)))

TFP_General_imputed %>%  
  mutate(Connector = log(Connector / geometric_all$Connector)) %>%
  mutate(Provincial = log(Provincial / geometric_all$Provincial)) %>%
  mutate(Satellite = log(Satellite / geometric_all$Satellite)) %>%
  mutate(Peripheral = log(Peripheral / geometric_all$Peripheral)) %>%
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
  scale_y_continuous(breaks = seq(-0.15, 0.15, 0.1)) +
  coord_cartesian(ylim = c(-0.15, 0.15)) +
  scale_color_brewer(palette = "PuOr") +
  geom_vline(xintercept = 54, color = "red", linewidth = 1.5, alpha = 1) +
  theme_pubclean(base_size = 18) +
  theme(plot.title.position = "plot",
        legend.title = element_blank()) +
  labs(y = "Centered log-ratios") +
  
  ggtitle("Lifespan dynamics of modular topological roles")


# For internodal roles

TFP_General_imputed %>%  
  mutate(Global_Bridge = log(Global_Bridge / geometric_all$Global_Bridge)) %>%
  mutate(Local_Bridge = log(Local_Bridge / geometric_all$Local_Bridge)) %>%
  mutate(Super_Bridge = log(Super_Bridge / geometric_all$Super_Bridge)) %>%
  mutate(Not_a_Bridge = log(Not_a_Bridge / geometric_all$Not_a_Bridge)) %>%
  dplyr::select(Subj_ID, Age, Global_Bridge:Not_a_Bridge) %>%
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
  scale_y_continuous(breaks = seq(-0.15, 0.15, 0.1)) +
  coord_cartesian(ylim = c(-0.15, 0.15)) +
  scale_color_brewer(palette = "PuOr",
                     labels = c("Global Interface", "Local_Interface", "Sub-Interface", "Super-Interface")) +
  geom_vline(xintercept = 54, color = "red", linewidth = 1.5, alpha = 1) +
  theme_pubclean(base_size = 18) +
  theme(plot.title.position = "plot",
        legend.title = element_blank()) +
  labs(y = "Centered log-ratios") +
  
  ggtitle("Lifespan dynamics of interface topological roles")
# facet_wrap(~Functional_role)

############################################################################
# DEFINING ILRs (ISOMETRIC LOG-RATIOS)
############################################################################

TFP_General_imputed_stats <- TFP_General_imputed %>% 
  # Integration
  mutate(Integration = ((3/4)^0.5)*log(Connector/(Provincial*Satellite*Peripheral)^(1/3))) %>%
  # Segregation
  mutate(Segregation = ((3/4)^0.5)*log(Provincial/(Connector*Satellite*Peripheral)^(1/3))) %>% 
  # Peripheralization
  mutate(Peripheralization = ((3/4)^0.5)*log(Peripheral/(Provincial*Satellite*Connector)^(1/3))) %>% 
  # Internodal role reconfiguration
  mutate(Polyvalence = (((2/3)^0.5)*log((Super_Bridge)^(1)/((Global_Bridge*Local_Bridge)^(1/2))))) 

mgcv::gam(Integration~s(Age), data = TFP_General_imputed_stats) %>% 
  summary()

mgcv::gam(Segregation~s(Age), data = TFP_General_imputed_stats) %>% 
  summary()

mgcv::gam(Peripheralization~s(Age), data = TFP_General_imputed_stats) %>% 
  summary()

mgcv::gam(Polyvalence~s(Age), data = TFP_General_imputed_stats) %>% 
  summary()

TFP_General_imputed_stats %>%  
  pivot_longer(c(Integration, 
                 Segregation,
                 Peripheralization,
                 Polyvalence),
               names_to = "topological_balances",
               values_to = "balance") %>% 
  group_by(topological_balances) %>% mutate(balance = as.numeric(scale(balance))) %>% 
  ggplot(aes(Age, balance, color = topological_balances)) +
  geom_hline(yintercept = 0, color = "red") +
  geom_hline(yintercept = 0, color = "red") +
  
  geom_jitter(height = 0.05, alpha = 0.15, size = 4) +
  geom_smooth(linewidth = 2, method = "gam", alpha = .1) +
  scale_x_continuous(breaks = seq(20, 90, 10)) +
  scale_y_continuous(breaks = seq(-0.4, 0.4, 0.1)) +
  geom_vline(xintercept = 54, color = "red", linewidth = 1.5, alpha = 1) +
  coord_cartesian(ylim = c(-0.4, 0.4)) +
  theme_pubr(base_size = 18) +
  theme(plot.title.position = "plot",
        legend.title = element_blank()) +
  labs(y = "Isometric log-ratios \n (z-scored)") +
  scale_color_brewer(palette = "Dark2") +
  ggtitle("Lifespan dynamic of the topological landscape")


################################################################################
# Decomposition at the RSN-level
################################################################################
# 
# list_TFP_RSN <- TFP_RSN %>% 
#   group_by(`1st_network`) %>% group_split()
# 
# list_tmp <- list()
# list_raw_imputed <- list()
# for (i in 1:length(list_TFP_RSN)) {
#   library(zCompositions)
#   library(compositions)
#   
#   tmp_raw <- rbindlist(list_TFP_RSN[i]) %>% arrange(Subj_ID) 
#   
#   tmp_coda_modular <- tmp_raw %>%
#     dplyr::select(Connector, Satellite, Provincial, Peripheral)
#   
#   if (min(tmp_coda_modular) == 0) {
#     tmp_coda_modular_bis <- tmp_coda_modular %>% 
#       acomp(.) %>% 
#       cmultRepl(., output = "prop")
#   } else {
#     tmp_coda_modular_bis <- tmp_coda_modular
#   }
#   
#   
#   tmp_coda_internodal <- tmp_raw %>%
#     dplyr::select(Global_Bridge, Local_Bridge, Super_Bridge, Not_a_Bridge)
#   
#   if (min(tmp_coda_internodal) == 0) {
#     tmp_coda_internodal_bis <- tmp_coda_internodal %>% 
#       acomp(.) %>%
#       cmultRepl(., output = "prop")
#   } else {
#     tmp_coda_internodal_bis <- tmp_coda_internodal
#   }
#   
#   tmp_raw_imputed <- cbind(tmp_raw %>% dplyr::select(-c(Connector, Satellite, Provincial, Peripheral,
#                                                         Global_Bridge, Local_Bridge, Super_Bridge, Not_a_Bridge)),
#                            tmp_coda_modular_bis, 
#                            tmp_coda_internodal_bis)
#   
#   list_raw_imputed[[i]] <- tmp_raw_imputed
#   
#   
#   tmp_geometric_all <- tmp_raw_imputed %>% 
#     summarize_at(vars(Connector:Not_a_Bridge), funs(geometricmean(.)))
#   
#   tmp_final <- tmp_raw_imputed %>% 
#     mutate(Connector = log(Connector / tmp_geometric_all$Connector)) %>%
#     mutate(Provincial = log(Provincial / tmp_geometric_all$Provincial)) %>%
#     mutate(Satellite = log(Satellite / tmp_geometric_all$Satellite)) %>%
#     mutate(Peripheral = log(Peripheral / tmp_geometric_all$Peripheral)) %>% 
#     mutate(Global_Bridge = log(Global_Bridge / tmp_geometric_all$Global_Bridge)) %>%
#     mutate(Local_Bridge = log(Local_Bridge / tmp_geometric_all$Local_Bridge)) %>%
#     mutate(Super_Bridge = log(Super_Bridge / tmp_geometric_all$Super_Bridge)) %>%
#     mutate(Not_a_Bridge = log(Not_a_Bridge / tmp_geometric_all$Not_a_Bridge)) 
#   
#   list_tmp[[i]] <- tmp_final
# }
# 
# TFP_RSN_imputed <- rbindlist(list_raw_imputed)
# TFP_RSN_CLR <- rbindlist(list_tmp)


# RSN <- "DMN"
# 
# TFP_RSN_CLR %>%
#   filter(grepl(RSN, `1st_network`)) %>%
#   # dplyr::select(Subj_ID, Age, Connector:Peripheral) %>%
#   dplyr::select(Subj_ID, Age, Connector:Peripheral) %>%
#   pivot_longer(
#     cols = !c("Subj_ID", "Age"),
#     names_to = "Functional_role",
#     values_to = "Score"
#   ) %>%
#   ggplot(aes(Age, Score, color = Functional_role)) +
#   geom_hline(yintercept = 0, color = "red") +
#   geom_jitter(height = 0.05, alpha = 0.1) +
#   geom_smooth(linewidth = 2, method = "gam", alpha = .3) +
#   scale_x_continuous(breaks = seq(20, 90, 5)) +
#   scale_y_continuous(breaks = seq(-1, 1, 0.1)) +
#   coord_cartesian(ylim = c(-0.2, 0.2)) +
#   scale_color_brewer(palette = "PuOr") +
#   theme_pubclean(base_size = 18) +
#   theme(plot.title.position = "plot") +
#   labs(y = "Centered log-ratios") +
#   ggtitle("Evolution of modular functional roles across adult lifespan")


################################################################################
# GRADIENT ANALYSIS
################################################################################
# 
# Gradient_stats <- TFP_RSN_imputed %>% 
#   # Integration
#   mutate(Integration = ((3/4)^0.5)*log(Connector/(Provincial*Satellite*Peripheral)^(1/3))) %>% 
#   # Peripherisation
#   mutate(Peripherisation = ((3/4)^0.5)*log(Peripheral/(Provincial*Satellite*Connector)^(1/3))) %>% 
#   # Internodal role reconfiguration
#   mutate(Polyvalent_interfaces = (((2/3)^0.5)*log((Super_Bridge)^(1)/((Global_Bridge*Local_Bridge)^(1/2)))))
# 
# # Gradient loadings
# Gradient_stats$`1st_network` <- factor(Gradient_stats$`1st_network`) %>% 
#   fct_reorder(Gradient_stats$G1, .desc = FALSE)
# 
# ggplot(Gradient_stats %>% group_by(`1st_network`) %>% 
#          summarize_at(vars(G1), mean), aes(`1st_network`, G1)) +
#   geom_col(aes(fill = G1)) +
#   scale_fill_distiller(palette = "RdBu", direction = -1) +
#   scale_y_continuous(breaks = seq(-8, 8, 2)) +
#   coord_flip() +
#   geom_hline(yintercept = 0, color = "red") +
#   theme_classic2(base_size = 18)
# 
# 
# # Integration
# 
# mod_lin <- mgcv::gam(Integration~Age*G1,
#                      data = Gradient_stats, method = "REML")
# 
# mod_gam_sAge <- mgcv::gam(Integration~s(Age, by = G1),
#                      data = Gradient_stats, method = "REML")
# 
# mod_gam_sG1 <- mgcv::gam(Integration~s(G1, by = Age),
#                      data = Gradient_stats, method = "REML")
# 
# mod_gam_full <- mgcv::gam(Integration~s(Age) + s(G1, k = 20) + 
#                        ti(Age, G1),
#                      data = Gradient_stats, method = "REML")
# 
# 
# AIC(mod_lin, mod_gam_sAge, mod_gam_sG1, mod_gam_full)
# summary(mod_gam_full)
# 
# mgcv::vis.gam(mod_gam_full, view = c("Age", "G1"),
#         plot.type = "persp", theta = 55, phi = 15,
#         n.grid = 30, lwd = 0.4, 
#         color = "topo", ticktype = "detailed")
# 
# Gradient_stats %>% 
#   filter(!(grepl("VMM|PMM", `1st_network`))) %>% 
#   group_by(`1st_network`) %>% 
#   mutate(Integration_scaled = as.numeric(scale(Integration))) %>% 
#   ungroup() %>% 
#   ggplot(aes(Age, Integration_scaled, color = forcats::fct_rev(`1st_network`))) +
#   geom_hline(yintercept = 0, color = "red") +
#   geom_jitter(height = 0.05, alpha = 0.05, size = 4) +
#   geom_smooth(linewidth = 2, method = "gam", alpha = 0) +
#   scale_x_continuous(breaks = seq(20, 90, 10)) +
#   scale_y_continuous(breaks = seq(-0.2, 0.3, 0.1)) +
#   coord_cartesian(ylim = c(-0.2, 0.3)) +
#   scale_color_brewer(palette = "RdBu") +
#   theme_classic2(base_size = 18) +
#   theme(plot.title.position = "plot") +
#   labs(y = "Integraion\n (z-scored ILR)",
#        color = "RSN") +
#   ggtitle("Modulation of the Integration mechanism by RSN across the adult lifespan")
# 
# Gradient_stats %>% 
#   mutate(G1_binned = Hmisc::cut2(G1, g = 4, levels.mean = TRUE)) %>% 
#   filter(!(grepl("VMM|PMM", `1st_network`))) %>% 
#   group_by(G1_binned) %>% 
#   mutate(Integration_scaled = as.numeric(scale(Integration))) %>% 
#   ungroup() %>% 
#   ggplot(aes(Age, Integration_scaled, color = forcats::fct_rev(G1_binned))) +
#   geom_hline(yintercept = 0, color = "red") +
#   geom_jitter(height = 0.05, alpha = 0.05, size = 4) +
#   geom_smooth(linewidth = 2, method = "gam", alpha = 0) +
#   scale_x_continuous(breaks = seq(20, 90, 10)) +
#   scale_y_continuous(breaks = seq(-0.2, 0.3, 0.1)) +
#   coord_cartesian(ylim = c(-0.2, 0.3)) +
#   scale_color_brewer(palette = "RdBu") +
#   theme_classic2(base_size = 18) +
#   theme(plot.title.position = "plot") +
#   labs(y = "Integration\n (z-scored ILR)",
#        color = "G1") +
#   ggtitle("Modulation of the Integration mechanism by G1\n across the adult lifespan")
# 
# 
# # Peripherisation
# 
# mod_lin <- mgcv::gam(Peripherisation~Age*G1,
#                      data = Gradient_stats, method = "REML")
# 
# mod_gam_sAge <- mgcv::gam(Peripherisation~s(Age, by = G1),
#                           data = Gradient_stats, method = "REML")
# 
# mod_gam_sG1 <- mgcv::gam(Peripherisation~s(G1, by = Age),
#                          data = Gradient_stats, method = "REML")
# 
# mod_gam_full <- mgcv::gam(Peripherisation~s(Age) + s(G1, k = 20) + 
#                             ti(Age, G1),
#                           data = Gradient_stats, method = "REML")
# 
# 
# AIC(mod_lin, mod_gam_sAge, mod_gam_sG1, mod_gam_full)
# summary(mod_gam_full)
# 
# Gradient_stats %>% 
#   filter(!(grepl("VMM|PMM", `1st_network`))) %>% 
#   group_by(`1st_network`) %>% 
#   mutate(Peripherisation_scaled = as.numeric(scale(Peripherisation))) %>% 
#   ungroup() %>% 
#   ggplot(aes(Age, Peripherisation_scaled, color = forcats::fct_rev(`1st_network`))) +
#   geom_hline(yintercept = 0, color = "red") +
#   geom_jitter(height = 0.05, alpha = 0.05, size = 4) +
#   geom_smooth(linewidth = 2, method = "gam", alpha = 0) +
#   scale_x_continuous(breaks = seq(20, 90, 10)) +
#   scale_y_continuous(breaks = seq(-0.2, 0.2, 0.1)) +
#   coord_cartesian(ylim = c(-0.2, 0.2)) +
#   scale_color_brewer(palette = "RdBu") +
#   theme_classic2(base_size = 18) +
#   theme(plot.title.position = "plot") +
#   labs(y = "Peripherisation\n (z-scored ILR)",
#        color = "G1") +
#   ggtitle("Modulation of the Peripherisation mechanism by RSN across the adult lifespan")
# 
# 
# Gradient_stats %>% 
#   mutate(G1_binned = Hmisc::cut2(G1, g = 4, levels.mean = TRUE)) %>% 
#   group_by(G1_binned) %>% 
#   mutate(Peripherisation_scaled = as.numeric(scale(Peripherisation))) %>% 
#   ungroup() %>% 
#   ggplot(aes(Age, Peripherisation_scaled, color = forcats::fct_rev(G1_binned))) +
#   geom_hline(yintercept = 0, color = "red") +
#   geom_jitter(height = 0.05, alpha = 0.05, size = 4) +
#   geom_smooth(linewidth = 2, method = "gam", alpha = 0) +
#   scale_x_continuous(breaks = seq(20, 90, 10)) +
#   scale_y_continuous(breaks = seq(-0.2, 0.2, 0.1)) +
#   coord_cartesian(ylim = c(-0.2, 0.2)) +
#   scale_color_brewer(palette = "RdBu") +
#   theme_classic2(base_size = 18) +
#   theme(plot.title.position = "plot") +
#   labs(y = "Peripherisation\n (z-scored ILR)",
#        color = "G1") +
#   ggtitle("Modulation of the Peripherisation mechanism by G1 across the adult lifespan")
# 
# 
# 
# # Polyvalent_Interfaces 
# mod_lin <- mgcv::gam(Polyvalent_interfaces~Age*G1,
#                      data = Gradient_stats, method = "REML")
# 
# mod_gam_sAge <- mgcv::gam(Polyvalent_interfaces~s(Age, by = G1),
#                           data = Gradient_stats, method = "REML")
# 
# mod_gam_sG1 <- mgcv::gam(Polyvalent_interfaces~s(G1, by = Age),
#                          data = Gradient_stats, method = "REML")
# 
# mod_gam_full <- mgcv::gam(Polyvalent_interfaces~s(Age) + s(G1, k = 20) + 
#                             ti(Age, G1),
#                           data = Gradient_stats, method = "REML")
# 
# 
# AIC(mod_lin, mod_gam_sAge, mod_gam_sG1, mod_gam_full)
# summary(mod_gam_full)




################################################################################
# Decomposition at the Community-level
################################################################################

list_TFP_ComStruct <- TFP_ComStruct %>% 
  group_by(Consensus_vector_0.15) %>% group_split()

list_tmp <- list()
list_raw_imputed <- list()
for (i in 1:length(list_TFP_ComStruct)) {
  library(zCompositions)
  library(compositions)
  
  tmp_raw <- rbindlist(list_TFP_ComStruct[i]) %>% arrange(Subj_ID) 
  
  tmp_coda_modular <- tmp_raw %>%
    dplyr::select(Connector, Satellite, Provincial, Peripheral)
  
  if (min(tmp_coda_modular) == 0) {
    tmp_coda_modular_bis <- tmp_coda_modular %>% 
      acomp(.) %>% 
      cmultRepl(., output = "prop")
  } else {
    tmp_coda_modular_bis <- tmp_coda_modular
  }
  
  
  tmp_coda_internodal <- tmp_raw %>%
    dplyr::select(Global_Bridge, Local_Bridge, Super_Bridge, Not_a_Bridge)
  
  if (min(tmp_coda_internodal) == 0) {
    tmp_coda_internodal_bis <- tmp_coda_internodal %>% 
      acomp(.) %>%
      cmultRepl(., output = "prop")
  } else {
    tmp_coda_internodal_bis <- tmp_coda_internodal
  }
  
  tmp_raw_imputed <- cbind(tmp_raw %>% dplyr::select(-c(Connector, Satellite, Provincial, Peripheral,
                                                        Global_Bridge, Local_Bridge, Super_Bridge, Not_a_Bridge)),
                           tmp_coda_modular_bis, 
                           tmp_coda_internodal_bis)
  
  list_raw_imputed[[i]] <- tmp_raw_imputed
  
  
  tmp_geometric_all <- tmp_raw_imputed %>% 
    summarize_at(vars(Connector:Not_a_Bridge), funs(geometricmean(.)))
  
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

TFP_ComStruct_imputed <- rbindlist(list_raw_imputed)
TFP_ComStruct_CLR <- rbindlist(list_tmp)

TFP_ComStruct_imputed %>%
  dplyr::select(Subj_ID, Age, Connector:Peripheral, Consensus_vector_0.15) %>%
  pivot_longer(
    cols = !c("Subj_ID", "Age", "Consensus_vector_0.15"),
    names_to = "Functional_role",
    values_to = "Score"
  ) %>%
  ggplot(aes(Age, Score, color = Functional_role)) +
  geom_hline(yintercept = 0, color = "red") +
  geom_jitter(height = 0.05, alpha = 0.1) +
  geom_smooth(linewidth = 2, method = "lm", alpha = .3) +
  scale_x_continuous(breaks = seq(20, 90, 5)) +
  scale_y_continuous(breaks = seq(0, 0.5, 0.1)) +
  coord_cartesian(ylim = c(0, 0.5)) +
  scale_color_brewer(palette = "PuOr") +
  theme_pubclean(base_size = 18) +
  theme(plot.title.position = "plot") +
  labs(y = "Raw proportions") +
  facet_wrap(~Consensus_vector_0.15, scales = "free") +
  ggtitle("Evolution of modular functional roles across the adult lifespan")

TFP_ComStruct_CLR %>%
  dplyr::select(Subj_ID, Age, Connector:Peripheral, Consensus_vector_0.15) %>%
  pivot_longer(
    cols = !c("Subj_ID", "Age", "Consensus_vector_0.15"),
    names_to = "Functional_role",
    values_to = "Score"
  ) %>%
  ggplot(aes(Age, Score, color = Functional_role)) +
  geom_hline(yintercept = 0, color = "red") +
  geom_jitter(height = 0.05, alpha = 0.1) +
  geom_smooth(linewidth = 2, method = "gam", alpha = .3) +
  scale_x_continuous(breaks = seq(20, 90, 5)) +
  scale_y_continuous(breaks = seq(-1, 1, 0.1)) +
  coord_cartesian(ylim = c(-0.2, 0.2)) +
  scale_color_brewer(palette = "PuOr") +
  theme_pubclean(base_size = 18) +
  theme(plot.title.position = "plot") +
  labs(y = "Centered log-ratios") +
  facet_wrap(~Consensus_vector_0.15) +
  ggtitle("Evolution of modular functional roles across the adult lifespan")

TFP_ComStruct_CLR %>%
  dplyr::select(Subj_ID, Age, Consensus_vector_0.15, Eloc, Eglob) %>%
  pivot_longer(
    cols = !c("Subj_ID", "Age", "Consensus_vector_0.15"),
    names_to = "Efficiencies",
    values_to = "value"
  ) %>%
  group_by(Consensus_vector_0.15, Efficiencies) %>% 
  mutate(value = as.numeric(scale(value))) %>% 
  ggplot(aes(Age, value, color = Efficiencies)) +
  geom_hline(yintercept = 0, color = "red") +
  geom_jitter(height = 0.05, alpha = 0.1) +
  geom_smooth(linewidth = 2, method = "lm", alpha = .3) +
  scale_x_continuous(breaks = seq(20, 90, 10)) +
  scale_y_continuous(breaks = seq(-0.5, 0.7, 0.2)) +
  coord_cartesian(ylim = c(-0.5, 0.7)) +
  scale_color_brewer(palette = "Dark2") +
  theme_pubclean(base_size = 18) +
  theme(plot.title.position = "plot") +
  labs(y = "Normalized Efficiencies") +
  facet_wrap(~Consensus_vector_0.15) +
  ggtitle("Evolution of the Balance I/S across the adult lifespan")


TFP_ComStruct_imputed_ILR <- TFP_ComStruct_imputed %>% 
  filter(Consensus_vector_0.15 != "RS-NET 5 (VMM)") %>% 
  # Integration
  mutate(Integration = ((3/4)^0.5)*log(Connector/(Provincial*Satellite*Peripheral)^(1/3))) %>%
  # Segregation
  mutate(Segregation = ((3/4)^0.5)*log(Provincial/(Connector*Satellite*Peripheral)^(1/3))) %>% 
  # Peripherisation
  mutate(Peripherisation = ((3/4)^0.5)*log(Peripheral/(Provincial*Satellite*Connector)^(1/3))) %>% 
  # Internodal role reconfiguration
  mutate(Polyvalent_interfaces = (((2/3)^0.5)*log((Super_Bridge)^(1)/((Global_Bridge*Local_Bridge)^(1/2)))))


TFP_ComStruct_imputed_ILR %>%
  filter(Consensus_vector_0.15 != 5) %>% 
  group_by(Consensus_vector_0.15) %>%
  mutate(Integration_scaled = as.numeric(scale(Integration))) %>%
  ungroup() %>%
  ggplot(aes(Age, Integration_scaled, color = forcats::fct_rev(Consensus_vector_0.15))) +
  geom_hline(yintercept = 0, color = "red") +
  geom_jitter(height = 0.05, alpha = 0.05, size = 4) +
  geom_smooth(linewidth = 2, method = "gam", alpha = 0) +
  scale_x_continuous(breaks = seq(20, 90, 10)) +
  scale_y_continuous(breaks = seq(-0.4, 0.4, 0.1)) +
  coord_cartesian(ylim = c(-0.4, 0.4)) +
  scale_color_brewer(palette = "Dark2") +
  theme_classic2(base_size = 18) +
  theme(plot.title.position = "plot") +
  labs(y = "Integraion\n (z-scored ILR)",
       color = "RSN") +
  ggtitle("Modulation of the Integration mechanism by RSN across the adult lifespan")

TFP_ComStruct_imputed_ILR %>%
  filter(Consensus_vector_0.15 != 5) %>% 
  group_by(Consensus_vector_0.15) %>%
  mutate(Segregation_scaled = as.numeric(scale(Segregation))) %>%
  ungroup() %>%
  ggplot(aes(Age, Segregation_scaled, color = forcats::fct_rev(Consensus_vector_0.15))) +
  geom_hline(yintercept = 0, color = "red") +
  geom_jitter(height = 0.05, alpha = 0.05, size = 4) +
  geom_smooth(linewidth = 2, method = "gam", alpha = 0) +
  scale_x_continuous(breaks = seq(20, 90, 10)) +
  scale_y_continuous(breaks = seq(-0.4, 0.4, 0.1)) +
  coord_cartesian(ylim = c(-0.4, 0.4)) +
  scale_color_brewer(palette = "Dark2") +
  theme_classic2(base_size = 18) +
  theme(plot.title.position = "plot") +
  labs(y = "Segregation\n (z-scored ILR)",
       color = "RSN") +
  ggtitle("Modulation of the Segregation mechanism by RSN across the adult lifespan")

TFP_ComStruct_imputed_ILR %>%
  filter(Consensus_vector_0.15 != 5) %>% 
  group_by(Consensus_vector_0.15) %>%
  mutate(Peripherisation_scaled = as.numeric(scale(Peripherisation))) %>%
  ungroup() %>%
  ggplot(aes(Age, Peripherisation_scaled, color = forcats::fct_rev(Consensus_vector_0.15))) +
  geom_hline(yintercept = 0, color = "red") +
  geom_jitter(height = 0.05, alpha = 0.05, size = 4) +
  geom_smooth(linewidth = 2, method = "gam", alpha = 0) +
  scale_x_continuous(breaks = seq(20, 90, 10)) +
  scale_y_continuous(breaks = seq(-0.4, 0.4, 0.1)) +
  coord_cartesian(ylim = c(-0.4, 0.4)) +
  scale_color_brewer(palette = "Dark2") +
  theme_classic2(base_size = 18) +
  theme(plot.title.position = "plot") +
  labs(y = "Peripherisation\n (z-scored ILR)",
       color = "RSN") +
  ggtitle("Modulation of the Peripherisation mechanism by RSN across the adult lifespan")

################################################################################
# DESCRIPTIVES
################################################################################

# data_TFP_analysis <- TFP_General %>% 
#   mutate(Age_group = ifelse(Age <= 39, "Young", ifelse(Age > 59, "Old", "Middle")))
# 
# data_TFP_analysis$Age_group <- factor(data_TFP_analysis$Age_group, levels = c("Young", "Middle", "Old"))
# 
# data_TFP_analysis %>%
#   group_by(Age_group) %>%
#   get_summary_stats("Age", type = "full")
# 
# data_TFP_analysis %>%
#   group_by(Age_group) %>%
#   count(Gender) %>%
#   mutate(n = prop.table(n))
# 
# 
# a <- gghistogram(data_TFP_analysis %>% subset(Age_group == "Young"),
#   x = "Age", y = "..density..", bins = 10,
#   fill = "purple", add_density = TRUE
# ) + theme_pubclean()
# b <- gghistogram(data_TFP_analysis %>% subset(Age_group == "Middle"),
#   x = "Age", y = "..density..", bins = 10,
#   fill = "purple", add_density = TRUE
# )
# c <- gghistogram(data_TFP_analysis %>% subset(Age_group == "Old"),
#   x = "Age", y = "..density..", bins = 15,
#   fill = "purple", add_density = TRUE
# ) + theme_pubclean()
# 
# 
# Rmisc::multiplot(a, b, c)
# 
# 
# # data_pca <- data_TFP_analysis %>% dplyr::select(-c("Subj_ID", "Gender", "Age", "Age_group", "Balance_eff"))
# # pc <- FactoMineR::PCA(data_pca)
# # pc$var$contrib
# # fviz_pca_ind(pc, geom.ind = "point", pointshape = 21,
# #              axes = c(1, 2),
# #              pointsize = 2,
# #              alpha.ind = .2,
# #              fill.ind = data_TFP_analysis$Age_group,
# #              col.ind = "black",
# #              palette = "jco",
# #              addEllipses = TRUE,
# #              label = "var",
# #              col.var = "black",
# #              repel = TRUE,
# #              legend.title = "Age group") +
# #   ggtitle("2D PCA-plot from the 6-graph-metric dataset") +
# #   theme(plot.title = element_text(hjust = 0.5))
# 
# 
# 
# 
# ################################################################################
# # Box plot ---------------------------------------------------------------------
# ################################################################################
# 
# data_box <- data_TFP_analysis %>%
#   pivot_longer(
#     cols = !c("Subj_ID", "Gender", "Age", "Age_group", "Eglob", "Eloc"),
#     names_to = "Metrics",
#     values_to = "Metric_value"
#   ) %>%
#   filter(Metrics != "Balance_eff")
# 
# 
# data_box$Metrics <- factor(data_box$Metrics, levels = c(
#   "Connector", "Provincial", "Satellite", "Peripheral",
#   "Global_Bridge", "Local_Bridge", "Super_Bridge", "Not_a_Bridge"
# ))
# 
# 
# # data_box_sig <- data_box %>%
# #   group_by(Metrics) %>%
# #   rstatix::t_test(Metric_value ~ Age_group) %>%
# #   adjust_pvalue(method = "fdr") %>%
# #   add_significance("p.adj") %>%
# #   add_xy_position(x = "Metrics")
# 
# ggplot(data_box, aes(x = Metrics, y = Metric_value)) +
#   geom_boxplot(aes(fill = Age_group)) +
#   scale_fill_brewer(palette = "YlOrRd") +
#   scale_y_continuous(limits = c(0, 50),
#                      breaks = seq(0, 50, 5)) +
#   theme_pubclean() 
# #   stat_pvalue_manual(data_box_sig,
# #   label = "p.adj.signif",
# #   tip.length = 0.01,
# #   hide.ns = TRUE,
# #   bracket.nudge.y = -5
# # ) 
# 
# # data_box_eff_size <- data_box %>%
# #   group_by(Metrics) %>%
# #   rstatix::cohens_d(Metric_value ~ Age_group, comparison = list(c("Young", "Old")), paired = FALSE, hedges.correction = TRUE) %>%
# #   mutate(effsize = effsize * (-1))
# # filter(magnitude != "negligible")
# 
# # ggdotchart(
# #   data_box_eff_size,
# #   x = "Metrics", y = "effsize",
# #   ylab = "Cohen's d (hedge corrected)",
# #   palette = "jco",
# #   add = "segment", position = position_dodge(0.3),
# #   sorting = "descending",
# #   rotate = TRUE, legend = "none"
# 
# 
# ################################################################################
# # Topologico-functional profile across clusters  -------------------------------
# ################################################################################
# 
# # Log ratio of the percentage of each metric of a cluster to the geometric mean of all individuals
# # equivalent to CLR-transform, preserves unit-sum constraint and removes value-range restriction
# 
# geometric_all <- data_TFP_analysis %>%
#   summarize_at(vars(Connector:Super_Bridge), funs(geomMeanExtension(., epsilon = 1e-1)))
# 
# Radar_functional_role_geometric_age <- data_TFP_analysis %>%
#   dplyr::select(-Age) %>%
#   filter(grepl("Young|Middle|Old", Age_group)) %>%
#   group_by(Age_group) %>%
#   summarize_at(vars(Connector:Super_Bridge), funs(geomMeanExtension(., epsilon = 1e-1))) %>%
#   mutate(Connector = log(Connector / geometric_all$Connector)) %>%
#   mutate(Provincial = log(Provincial / geometric_all$Provincial)) %>%
#   mutate(Satellite = log(Satellite / geometric_all$Satellite)) %>%
#   mutate(Peripheral = log(Peripheral / geometric_all$Peripheral)) %>%
#   mutate(Global_Bridge = log(Global_Bridge / geometric_all$Global_Bridge)) %>%
#   mutate(Local_Bridge = log(Local_Bridge / geometric_all$Local_Bridge)) %>%
#   mutate(Super_Bridge = log(Super_Bridge / geometric_all$Super_Bridge)) %>%
#   mutate(Not_a_Bridge = log(Not_a_Bridge / geometric_all$Not_a_Bridge)) %>%
#   remove_rownames() %>%
#   column_to_rownames(var = "Age_group")
# 
# 
# palette <- RColorBrewer::brewer.pal(3, "YlOrRd")
# 
# radarplotting_overlap(Radar_functional_role_geometric_age, 0.2, -0.2, 1, 1,
#   alpha = 0.1, label_size = 1, round = FALSE,
#   title_fill = "Lifespan Topologico-functional profile (log ratio of geometric means)",
#   palette = palette
# )
# 
# legend(
#   x = "topright",
#   legend = rownames(Radar_functional_role_geometric_age), horiz = TRUE,
#   bty = "n", pch = 20, col = palette,
#   text.col = "black", cex = 1, pt.cex = 2
# )
# 
# Radar_functional_role_age <- data_TFP_analysis %>%
#   dplyr::select(-Age) %>%
#   filter(grepl("Young|Middle|Old", Age_group)) %>%
#   group_by(Age_group) %>%
#   summarize_at(vars(Connector:Super_Bridge), funs(mean(.))) %>%
#   remove_rownames() %>%
#   column_to_rownames(var = "Age_group")
# 
# radarplotting_overlap(Radar_functional_role_age, 50, 0, 1, 1,
#   alpha = 0.1, label_size = 1,
#   title_fill = "Profile expressed in relative proportions",
#   palette = palette
# )
# 
# legend(
#   x = "topright",
#   legend = rownames(Radar_functional_role_age), horiz = TRUE,
#   bty = "n", pch = 20, col = palette,
#   text.col = "black", cex = 1, pt.cex = 2
# )
# 
# 
# ################################################################################
# # QUALITY CHECK -------------
# ################################################################################
# # Distribution of hubs across RSNs for each cluster for the individual hubs ----
# 
# data_cluster_selection("Young", "Old")
# 
# Radar_hub_RSN <-
#   tmp_cluster_final %>%
#   group_by(Subj_ID, Age_group, `1st_network`) %>%
#   summarise(n = n()) %>%
#   # mutate(freq = n / sum(n)) %>%
#   # dplyr::select(-n) %>%
#   group_by(Age_group, `1st_network`) %>%
#   summarize_at(vars(n), mean) %>%
#   spread(`1st_network`, n) %>%
#   remove_rownames() %>%
#   column_to_rownames("Age_group")
# # mutate_at(vars(everything()), funs(. * 100))
# 
# radarplotting_overlap(Radar_hub_RSN, 30, 0, 1, 1,
#   alpha = 0.3, label_size = 1,
#   title_fill = "Distribution of hubs regions across RSNs",
#   palette = RColorBrewer::brewer.pal(8, "Dark2")
# )
# 
# legend(
#   x = "topright",
#   legend = rownames(Radar_hub_RSN), horiz = TRUE,
#   bty = "n", pch = 20, col = RColorBrewer::brewer.pal(8, "Dark2"),
#   text.col = "black", cex = 1, pt.cex = 2
# )
# # 
# # # Distribution of hubs across communities for each Age_group for the individual hubs ----
# # 
# # Radar_hub_community <- tmp_cluster_final %>%
# #   group_by(Subj_ID, Age_group, Consensus_vector_0.15) %>%
# #   summarise(n = n()) %>%
# #   mutate(freq = n / sum(n)) %>%
# #   dplyr::select(-n) %>%
# #   group_by(Age_group, Consensus_vector_0.15) %>%
# #   summarize_at(vars(freq), mean) %>%
# #   spread(Consensus_vector_0.15, freq) %>%
# #   remove_rownames() %>%
# #   column_to_rownames("Age_group") %>%
# #   mutate_at(vars(everything()), funs(. * 100))
# # 
# # radarplotting_overlap(Radar_hub_community, 50, 0, 1, 1,
# #   alpha = 0.2, label_size = 1,
# #   title_fill = "Distribution of hubs regions across communities",
# #   palette = RColorBrewer::brewer.pal(8, "Dark2")
# # )
# # 
# # legend(
# #   x = "topright",
# #   legend = rownames(Radar_hub_community), horiz = TRUE,
# #   bty = "n", pch = 20, col = RColorBrewer::brewer.pal(8, "Dark2"),
# #   text.col = "black", cex = 1, pt.cex = 2
# # )
# 
# 
# ################################################################################
# # ILR TRANSFORMATION 
# ################################################################################
# 
# data_coda_modular <- data_TFP_analysis %>%
#   dplyr::select(Connector, Satellite, Provincial, Peripheral) %>% 
#   acomp(.) %>% 
#   # Preserves the ratios between non-zero components
#   cmultRepl(., output = "prop")
# 
# data_coda_interareal <- data_TFP_analysis %>%
#   dplyr::select(Global_Bridge, Local_Bridge, Super_Bridge, Not_a_Bridge) %>%
#   acomp(.) %>%
#   cmultRepl(., output = "prop")
# 
# data_imputed <- cbind(data_coda_modular, data_coda_interareal,
#                       data_TFP_analysis %>% dplyr::select(Subj_ID, Age, Age_group, Gender, Eglob, Eloc, Balance_eff))
# 
# data_imputed <- data_imputed %>%  
#   mutate(ilr_modular = (((1/2)^0.5)*log(Connector/((Provincial*Satellite)^0.5)))) %>% 
#   mutate(ilr_modular = as.numeric(scale(ilr_modular))) %>%
#   mutate(Eglob = as.numeric(scale(Eglob))) %>% 
#   mutate(Age = as.numeric(scale(Age))) %>% 
#   mutate(ilr_interareal = (((2/3)^0.5)*log(Global_Bridge/((Super_Bridge*Not_a_Bridge)^0.5))))
# 
# ################################################################################
# # MEDIATION ANALYSIS
# ################################################################################
# 
# library(robustbase)
# # STEP 1: Total effect DV ~ IV
# mod1 <- lmrob(Eglob~ Age, data = data_imputed, method = "MM")
# summary(mod1)
# # plot(mod1)
# data_flexplot <- data_imputed %>% filter(!(grepl("469|378|159", Subj_ID)))
# flexplot(Eglob~Age, data_flexplot, method = "lm")
# 
# # STEP 2: Indirect effect Mediator ~ IV
# mod2 <- lmrob(ilr_modular~ Age, data = data_imputed, method = "MM")
# summary(mod2)
# # plot(mod2)
# data_flexplot_mod2 <- data_imputed %>% filter(!(grepl("538|332|159", Subj_ID)))
# flexplot(ilr_modular~ Age, data_flexplot_mod2, method = "polynomial")
# 
# # Step 3: DV ~ IV + mediator
# mod3 <- lmrob(Eglob~ ilr_modular + Age, data = data_imputed, method = "MM")
# summary(mod3)
# # plot(mod3)
# data_flexplot_mod3 <- data_imputed %>% filter(!(grepl("469|242|378", Subj_ID)))
# flexplot(Eglob~ ilr_modular, data_flexplot_mod3, method = "polynomial")
# 
# psych::mediate(Eglob~Age + (ilr_modular), data = data_imputed) %>% summary()
# med <- robmed::test_mediation(Eglob~m(ilr_modular) + Age, data = data_imputed, robust = "MM")
# summary(med)
# robmed::ellipse_plot(med)
# 
# # The effect of Age on global efficiency was fully mediated
# # via the balance between Connector and Provincial/Satellite hubs.
# # We tested the significance of the indirect effect using bootstrapping procedures.
# # Unstandardized indirect effects were computed for each of 1000 boostrapped sample
# 
# ################################################################################
# ################################################################################
# # source("outliersFunction.R")
# # outliers(Eglob ~ ilr_modular, data_imputed)
# 
# data_imputed_plot <- data_imputed 
# #define datasets
# x1 <-  data_imputed_plot %>% filter(Age_group == "Young") %>% dplyr::select(ilr_modular) %>% as.matrix()
# y1 <-  data_imputed_plot %>% filter(Age_group == "Young") %>% dplyr::select(Eglob) %>% as.matrix()
# 
# x2 <-  data_imputed_plot %>% filter(Age_group == "Middle") %>% dplyr::select(ilr_modular) %>% as.matrix()
# y2 <-  data_imputed_plot %>% filter(Age_group == "Middle") %>% dplyr::select(Eglob) %>% as.matrix()
# 
# x3 <-  data_imputed_plot %>% filter(Age_group == "Old") %>% dplyr::select(ilr_modular) %>% as.matrix()
# y3 <-  data_imputed_plot %>% filter(Age_group == "Old") %>% dplyr::select(Eglob) %>% as.matrix()
# 
# #create scatterplot of x1 vs. y1
# plot(x1, y1, col='blue', pch=19, cex=1,
#      xlab='Connector vs Provincial & Satellite', ylab='Global efficiency', 
#      main='Modular balance and Global efficiency (proportional) ')
# 
# points(x2, y2, col='orange', pch=15, cex=1)
# points(x3, y3, col='red', pch=17, cex=1)
# #add legend
# legend("bottomleft", legend=c('Young', 'Middle', 'Old'), pch=c(19, 15, 17), col=c('blue', 'orange', 'red'))
# 
# #add polynomial curve to plot
# mod <- robustbase::lmrob(Eglob ~ ilr_modular + Age, data_imputed_plot, method = "MM")
# summary(mod)
# pred <- predict(mod)
# abline(robustbase::lmrob(Eglob ~ ilr_modular + Age, data_imputed_plot, method = "MM"), lwd = 3)
# # ix <- sort(data_imputed_plot$ilr_modular %>% as.matrix(), index.return=T)$ix
# # lines(data_imputed_plot$ilr_modular[ix], pred[ix], col='black', lwd=2)
# 
# text(0, 0.44, "Eglob = 0.51 - 0.02 * b1; p-value < 2e-16 (with MM estimator)")
# 
# 
