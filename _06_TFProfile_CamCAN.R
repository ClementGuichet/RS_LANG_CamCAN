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
# library(rio)
# TFP_General <- rio::import("TFP_General_PT_627subj.csv") %>% dplyr::select(-V1)

############################################################################
# Bayesian non-parametric multiplicative replacement
############################################################################
library(zCompositions)
library(compositions)

tmp_coda_modular <- TFP_General %>%
  dplyr::select(Connector, Satellite, Provincial, Peripheral) %>% 
  acomp(.) %>% 
  # Preserves the ratios between non-zero components
  cmultRepl(., output = "prop")

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
  scale_y_continuous(breaks = seq(-0.2, 0.2, 0.1)) +
  coord_cartesian(ylim = c(-0.2, 0.2)) +
  scale_color_brewer(palette = "PuOr") +
  geom_vline(xintercept = 52, color = "red", linewidth = 1.5, alpha = 1) +
  theme_pubclean(base_size = 18) +
  theme(plot.title.position = "plot") +
  labs(y = "Centered log-ratios") +
  
  ggtitle("Evolution of modular functional roles across adult lifespan")


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
  scale_y_continuous(breaks = seq(-0.2, 0.2, 0.1)) +
  coord_cartesian(ylim = c(-0.2, 0.2)) +
  scale_color_brewer(palette = "PuOr") +
  geom_vline(xintercept = 52, color = "red", linewidth = 1.5, alpha = 1) +
  theme_pubclean(base_size = 18) +
  theme(plot.title.position = "plot") +
  labs(y = "Centered log-ratios") +
  
  ggtitle("Evolution of internodal functional roles across adult lifespan")
  # facet_wrap(~Functional_role)

############################################################################
# DEFINING ILRs (ISOMETRIC LOG-RATIOS)
############################################################################

TFP_General_imputed_stats <- TFP_General_imputed %>% 
  # Inter-modular integration
  mutate(ilr_modular_inter = (((1/2)^0.5)*log(Connector/(Provincial)))) %>%
  # Intra-module integration
  mutate(ilr_modular_intra = (((2/3)^0.5)*log((Connector*Peripheral)^(1/2)/(Satellite)))) %>% 
  # Internodal role reconfiguration
  mutate(ilr_internodal = (((4/4)^0.5)*log((Global_Bridge*Local_Bridge)^(1/2)/((Super_Bridge*Not_a_Bridge)^(1/2)))))
      


mgcv::gam(ilr_modular_inter~Age + s(Age), data = TFP_General_imputed_stats) %>% 
  summary()

TFP_General_imputed_stats %>%  
  dplyr::select(Subj_ID, Age, ilr_modular_inter) %>% 
  ggplot(aes(Age, ilr_modular_inter)) +
  geom_hline(yintercept = 0, color = "red") +
  geom_jitter(height = 0.05, alpha = 0.1, size = 4) +
  geom_smooth(linewidth = 2, method = "gam", alpha = .3) +
  scale_x_continuous(breaks = seq(20, 90, 5)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_color_brewer(palette = "PuOr") +
  geom_vline(xintercept = 52, color = "red", linewidth = 1.5, alpha = 1) +
  theme_pubclean(base_size = 18) +
  theme(plot.title.position = "plot") +
  labs(y = "Balance Connector vs. Provincial hubs") +
  ggtitle("Inter-modular integration mechanism across adult lifespan")

mgcv::gam(ilr_modular_intra~Age + s(Age), data = TFP_General_imputed_stats) %>% 
  summary()

TFP_General_imputed_stats %>%  
  dplyr::select(Subj_ID, Age, ilr_modular_intra) %>% 
  ggplot(aes(Age, ilr_modular_intra)) +
  geom_hline(yintercept = 0, color = "red") +
  geom_jitter(height = 0.05, alpha = 0.1, size = 4) +
  geom_smooth(linewidth = 2, method = "gam", alpha = .3) +
  scale_x_continuous(breaks = seq(20, 90, 5)) +
  scale_y_continuous(breaks = seq(-1, 0, 0.2)) +
  coord_cartesian(ylim = c(-1, 0)) +
  scale_color_brewer(palette = "PuOr") +
  geom_vline(xintercept = 52, color = "red", linewidth = 1.5, alpha = 1) +
  theme_pubclean(base_size = 18) +
  theme(plot.title.position = "plot") +
  labs(y = "Balance Connector & Peripheral\n vs. Satellite hubs") +
  ggtitle("Intra-modular integration mechanism across adult lifespan")

mgcv::gam(ilr_internodal~Age + s(Age), data = TFP_General_imputed_stats) %>% 
  summary()

TFP_General_imputed_stats %>%  
  dplyr::select(Subj_ID, Age, ilr_internodal) %>% 
  ggplot(aes(Age, ilr_internodal)) +
  geom_hline(yintercept = 0, color = "red") +
  geom_jitter(height = 0.05, alpha = 0.1) +
  geom_smooth(linewidth = 2, method = "gam", alpha = .3) +
  scale_x_continuous(breaks = seq(20, 90, 5)) +
  scale_y_continuous(breaks = seq(-1, 0, 0.2)) +
  coord_cartesian(ylim = c(-1, 0)) +
  scale_color_brewer(palette = "PuOr") +
  geom_vline(xintercept = 52, color = "red", linewidth = 1.5, alpha = 1) +
  theme_pubclean(base_size = 18) +
  theme(plot.title.position = "plot") +
  labs(y = "ILR_modular_inter") +
  ggtitle("Evolution of the internodal balance across adult lifespan")

################################################################################
# Decomposition at the RSN-level
################################################################################
list_TFP_RSN <- TFP_RSN %>% 
  group_by(`1st_network`) %>% group_split()

list_tmp <- list()
list_raw_imputed <- list()
for (i in 1:length(list_TFP_RSN)) {
  library(zCompositions)
  library(compositions)
  
  tmp_raw <- rbindlist(list_TFP_RSN[i]) %>% arrange(Subj_ID) 
  
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

TFP_RSN_imputed <- rbindlist(list_raw_imputed)
TFP_RSN_CLR <- rbindlist(list_tmp)


RSN <- "DMN"

TFP_RSN_CLR %>% 
  filter(grepl(RSN, `1st_network`)) %>% 
  # dplyr::select(Subj_ID, Age, Connector:Peripheral) %>%
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
  coord_cartesian(ylim = c(-0.2, 0.2)) +
  scale_color_brewer(palette = "PuOr") +
  theme_pubclean(base_size = 18) +
  theme(plot.title.position = "plot") +
  labs(y = "Centered log-ratios") +
  ggtitle("Evolution of modular functional roles across adult lifespan")


################################################################################
# GRADIENT ANALYSIS
################################################################################

geometric_all <- TFP_RSN_imputed %>% 
  summarize_at(vars(Connector:Not_a_Bridge), funs(geometricmean(.)))

Gradient_stats <- TFP_RSN_imputed %>% 
  mutate(Age_group = ifelse(Age <= 39, "Young", ifelse(Age <= 59, "Middle", "Old"))) %>% 
  # Define ILR with raw components
   # Inter-modular integration
  mutate(ilr_modular_inter = (((1/2)^0.5)*log(Connector/(Provincial)))) %>%
  # Intra-module integration
  mutate(ilr_modular_intra = (((2/3)^0.5)*log((Connector*Peripheral)^(1/2)/(Satellite)))) %>% 
  # Internodal role reconfiguration
  mutate(ilr_internodal = (((4/4)^0.5)*log((Global_Bridge*Local_Bridge)^(1/2)/((Super_Bridge*Not_a_Bridge)^(1/2))))) %>% 
  # Standardized to CLR space
  mutate(Connector = log(Connector / geometric_all$Connector)) %>%
  mutate(Provincial = log(Provincial / geometric_all$Provincial)) %>%
  mutate(Satellite = log(Satellite / geometric_all$Satellite)) %>%
  mutate(Peripheral = log(Peripheral / geometric_all$Peripheral))

Gradient_stats$`1st_network` <- factor(Gradient_stats$`1st_network`) %>% 
  fct_reorder(Gradient_stats$G1, .desc = FALSE)

ggplot(Gradient_stats %>% group_by(`1st_network`) %>% 
         summarize_at(vars(G1), mean), aes(`1st_network`, G1)) +
  geom_col(aes(fill = G1)) +
  scale_fill_distiller(palette = "RdBu", direction = -1) +
  scale_y_continuous(breaks = seq(-8, 8, 2)) +
  coord_flip() +
  geom_hline(yintercept = 0, color = "red") +
  theme_classic2(base_size = 18)

Gradient_GAMM <- Gradient_stats %>%  
  dplyr::select(Subj_ID, Age, Age_group, ilr_modular_inter, ilr_modular_intra, ilr_internodal, G1) %>% 
  mutate(G1_binned = ifelse(G1 < -4, "-8 to -4", 
                            ifelse(G1 < 0, "-4 to 0",
                                   ifelse(G1 < 4, "0 to 4", ">4")))) 

mod1 <- gam(ilr_modular_inter~s(Age) + s(G1, k = 20) + 
              ti(Age, G1, bs = "tp"),
            data = Gradient_GAMM, method = "REML")

summary(mod1)
AIC(mod1)
gam.check(mod1)
concurvity(mod1, full = TRUE)

par(mfrow = c(1, 1))
vis.gam(mod1, view = c("Age", "G1"),
        plot.type = "persp", theta = 55, phi = 20,
        n.grid = 30, lwd = 0.4, 
        color = "topo", ticktype = "detailed")

Gradient_GAMM %>% 
  group_by(G1_binned) %>% 
  mutate(ilr_modular_inter_scaled = as.numeric(scale(ilr_modular_inter))) %>% 
  ungroup() %>% 
  ggplot(aes(Age, ilr_modular_inter_scaled, color = G1_binned)) +
  geom_hline(yintercept = 0, color = "red") +
  geom_jitter(height = 0.05, alpha = 0.05, size = 4) +
  geom_smooth(linewidth = 2, method = "gam", alpha = .1) +
  scale_x_continuous(breaks = seq(20, 90, 5)) +
  scale_y_continuous(breaks = seq(-0.2, 0.2, 0.1)) +
  coord_cartesian(ylim = c(-0.2, 0.2)) +
  scale_color_brewer(palette = "PuOr") +
  theme_pubclean(base_size = 18) +
  theme(plot.title.position = "plot") +
  labs(y = "Balance Connector vs. Provincial hubs") +
  ggtitle("Modulation of the inter-modular integration mechanism by G1\n across the adult lifespan")


mod1 <- gam(ilr_modular_intra~s(G1, by = Age, k = 20),
            data = Gradient_GAMM, method = "REML")

summary(mod1)
par(mfrow = c(1, 1))
vis.gam(mod1, view = c("Age", "G1"),
        theta = 35, n.grid = 30, lwd = 0.4, color = "topo", ticktype = "detailed")

Gradient_GAMM %>% 
  group_by(G1_binned) %>% 
  mutate(ilr_modular_intra_scaled = as.numeric(scale(ilr_modular_intra))) %>% 
  ungroup() %>% 
  ggplot(aes(Age, ilr_modular_intra_scaled, color = G1_binned)) +
  geom_hline(yintercept = 0, color = "red") +
  geom_jitter(height = 0.05, alpha = 0.05, size = 4) +
  geom_smooth(linewidth = 2, method = "lm", alpha = .1) +
  scale_x_continuous(breaks = seq(20, 90, 5)) +
  scale_y_continuous(breaks = seq(-0.2, 0.2, 0.1)) +
  coord_cartesian(ylim = c(-0.2, 0.2)) +
  scale_color_brewer(palette = "PuOr") +
  theme_pubclean(base_size = 18) +
  theme(plot.title.position = "plot") +
  labs(y = "Balance Connector vs. Provincial hubs") +
  ggtitle("Modulation of the intra-modular integration mechanism by G1\n across the adult lifespan")


mod1 <- gam(ilr_internodal~s(G1, by = Age, k = 20), data = Gradient_GAMM, method = "REML")
summary(mod1)

par(mfrow = c(1, 1))
vis.gam(mod1, view = c("Age", "G1"),
        theta = 35, n.grid = 30, lwd = 0.4, color = "topo", ticktype = "detailed")

Gradient_GAMM %>% 
  group_by(G1_binned) %>% 
  mutate(ilr_internodal_scaled = as.numeric(scale(ilr_internodal))) %>% 
  ungroup() %>% 
  ggplot(aes(Age, ilr_internodal_scaled, color = G1_binned)) +
  geom_hline(yintercept = 0, color = "red") +
  geom_jitter(height = 0.05, alpha = 0.05, size = 4) +
  geom_smooth(linewidth = 2, method = "lm", alpha = .1) +
  scale_x_continuous(breaks = seq(20, 90, 5)) +
  scale_y_continuous(breaks = seq(-0.2, 0.2, 0.1)) +
  coord_cartesian(ylim = c(-0.2, 0.2)) +
  scale_color_brewer(palette = "PuOr") +
  theme_pubclean(base_size = 18) +
  theme(plot.title.position = "plot") +
  labs(y = "Balance Connector vs. Provincial hubs") +
  ggtitle("Modulation of the internodal integration mechanism by G1\n across the adult lifespan")


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
