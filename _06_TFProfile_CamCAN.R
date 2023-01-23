##########################################################################################
# Script for Age-related analyses and visualization for the topologico-functional profile

# Written by CG
# 13-12-2022
##########################################################################################
library(rstatix)
library(ggpubr)
library(factoextra)

source("_05_Hub_detection_CamCAN.R")
source("_radarplotting_function.R")
source("_geometricmeanCruz.R")

################################################################################
# Investigating the evolution of graph-based metrics ----
################################################################################

TFP_General %>%
  dplyr::select(Subj_ID, Age, Connector, Provincial, Peripheral, Satellite) %>%
  pivot_longer(
    cols = !c("Subj_ID", "Age"),
    names_to = "Functional_role",
    values_to = "Score"
  ) %>%
  ggplot(aes(Age, Score, color = Functional_role)) +
  geom_jitter(height = 0.05, alpha = 0.1) +
  geom_smooth(linewidth = 2, method = "gam") +
  theme_pubclean() +
  ggtitle("Evolution of modular functional roles across adult lifespan")

TFP_General %>%
  dplyr::select(Subj_ID, Age, Global_Bridge, Local_Bridge, Super_Bridge, Not_a_Bridge) %>%
  pivot_longer(
    cols = !c("Subj_ID", "Age"),
    names_to = "Functional_role",
    values_to = "Score"
  ) %>%
  ggplot(aes(Age, Score, color = Functional_role)) +
  geom_jitter(height = 0.05, alpha = 0.1) +
  geom_smooth(linewidth = 2, method = "gam") +
  theme_pubclean() +
  ggtitle("Evolution of interareal functional roles across adult lifespan")

################################################################################

data_TFP_analysis <- TFP_General %>%
  filter(Age != "NaN") %>%
  plyr::rename(c("gender_text" = "Gender")) %>%
  mutate(Age_group = ifelse(Age <= 39, "Young", ifelse(Age > 59, "Old", "Middle")))

data_TFP_analysis$Age_group <- factor(data_TFP_analysis$Age_group, levels = c("Young", "Middle", "Old"))

data_TFP_analysis %>%
  group_by(Age_group) %>%
  get_summary_stats("Age", type = "full")

data_TFP_analysis %>%
  group_by(Age_group) %>%
  count(Gender) %>%
  mutate(n = prop.table(n))


a <- gghistogram(data_TFP_analysis %>% subset(Age_group == "Young"),
  x = "Age", y = "..density..", bins = 10,
  fill = "purple", add_density = TRUE
) + theme_pubclean()
b <- gghistogram(data_TFP_analysis %>% subset(Age_group == "Middle"),
  x = "Age", y = "..density..", bins = 10,
  fill = "purple", add_density = TRUE
)
c <- gghistogram(data_TFP_analysis %>% subset(Age_group == "Old"),
  x = "Age", y = "..density..", bins = 15,
  fill = "purple", add_density = TRUE
) + theme_pubclean()


Rmisc::multiplot(a, b, c)


# data_pca <- data_TFP_analysis %>% dplyr::select(-c("Subj_ID", "Gender", "Age", "Age_group", "Balance_eff"))
# pc <- FactoMineR::PCA(data_pca)
# pc$var$contrib
# fviz_pca_ind(pc, geom.ind = "point", pointshape = 21,
#              axes = c(1, 2),
#              pointsize = 2,
#              alpha.ind = .2,
#              fill.ind = data_TFP_analysis$Age_group,
#              col.ind = "black",
#              palette = "jco",
#              addEllipses = TRUE,
#              label = "var",
#              col.var = "black",
#              repel = TRUE,
#              legend.title = "Age group") +
#   ggtitle("2D PCA-plot from the 6-graph-metric dataset") +
#   theme(plot.title = element_text(hjust = 0.5))


# Final dataframe with only the subjects of chosen clusters, their hub regions and the RSNs -----
data_cluster_selection <- function(cluster1, cluster2) {
  # Retain only the rows specific of the two clusters
  tmp_cluster_0 <- data_TFP_analysis %>% subset(Age_group == cluster1 | Age_group == cluster2)
  # Get the associated Resting-state networks
  tmp_cluster_1 <- filter(data_functional_role, Subj_ID %in% tmp_cluster_0$Subj_ID)
  # Hub region specific to each subject yielded by hub detection procedure
  data_hub_selection_per_subject <- rbindlist(Hub_selection)
  # Select the subjects from the clusters
  data_hub_selection_cluster <- filter(
    data_hub_selection_per_subject,
    Subj_ID %in% tmp_cluster_1$Subj_ID
  )
  tmp_cluster_final <<- merge(data_hub_selection_cluster, tmp_cluster_0 %>%
    dplyr::select(Subj_ID, Age_group),
  by = "Subj_ID"
  )
}
data_cluster_selection("Young", "Old")





################################################################################
# Box plot ---------------------------------------------------------------------
################################################################################

data_box <- data_TFP_analysis %>%
  pivot_longer(
    cols = !c("Subj_ID", "Gender", "Age", "Age_group"),
    names_to = "Metrics",
    values_to = "Metric_value"
  ) %>%
  filter(Metrics != "Balance_eff")


data_box$Metrics <- factor(data_box$Metrics, levels = c(
  "Connector", "Provincial", "Satellite", "Peripheral",
  "Global_Bridge", "Local_Bridge", "Super_Bridge", "Not_a_Bridge"
))


# data_box_sig <- data_box %>%
#   group_by(Metrics) %>%
#   rstatix::t_test(Metric_value ~ Age_group) %>%
#   adjust_pvalue(method = "fdr") %>%
#   add_significance("p.adj") %>%
#   add_xy_position(x = "Metrics")

ggplot(data_box, aes(x = Metrics, y = Metric_value)) +
  geom_boxplot(aes(fill = Age_group)) +
  scale_fill_brewer(palette = "YlOrRd") +
  theme_pubr()
# stat_pvalue_manual(data_box_sig,
#   label = "p.adj.signif",
#   tip.length = 0.01,
#   hide.ns = TRUE,
#   bracket.nudge.y = -5
# )

# data_box_eff_size <- data_box %>%
#   group_by(Metrics) %>%
#   rstatix::cohens_d(Metric_value ~ Age_group, comparison = list(c("Young", "Old")), paired = FALSE, hedges.correction = TRUE) %>%
#   mutate(effsize = effsize * (-1))
# # filter(magnitude != "negligible")
#
# ggdotchart(
#   data_box_eff_size,
#   x = "Metrics", y = "effsize",
#   ylab = "Cohen's d (hedge corrected)",
#   palette = "jco",
#   add = "segment", position = position_dodge(0.3),
#   sorting = "descending",
#   rotate = TRUE, legend = "none"
# )



library(lme4)
source("PRE.R")
mod <- lm(Balance_eff ~ Age, data_TFP_analysis)
mod %>% summary(.)
effectsize::eta_squared(mod, partial = TRUE, alternative = "greater")

cor.test(TFP_General$Age, TFP_General$Balance_eff)
# Balance Integration/Segregation
plot(TFP_General$Age, TFP_General$Balance_eff, pch = 19, col = "darkblue")
# Regression line
abline(lm(TFP_General$Balance_eff ~ TFP_General$Age), col = "red", lwd = 3)
# Pearson correlation
text(paste(
  "Correlation between Age and I/S Balance (t(642) = -4.55, p < .001):",
  round(cor.test(TFP_General$Age, TFP_General$Balance_eff)$estimate, 2)
), x = 35, y = 6)


cor_efficiency <- data_functional_role %>%
  dplyr::select(Subj_ID, Age, Eglob, Eloc) %>%
  group_by(Age) %>%
  summarize_at(vars(Eglob, Eloc), mean) %>%
  mutate(Age_quadratic = Age^2)

cor.test(cor_efficiency$Age, cor_efficiency$Eloc)
# Balance Integration/Segregation
plot(cor_efficiency$Age, cor_efficiency$Eloc, pch = 19, col = "darkblue")
# Regression line
abline(lm(cor_efficiency$Eloc ~ cor_efficiency$Age), col = "red", lwd = 3)
# Pearson correlation
text(paste(
  "Correlation between Age and Local efficiency:",
  round(cor.test(cor_efficiency$Age, cor_efficiency$Eloc)$estimate, 2)
), x = 35, y = 0.675)

################################################################################
# Topologico-functional profile across clusters  -------------------------------
################################################################################

# Log ratio of the percentage of each metric of a cluster to the geometric mean of all individuals
# equivalent to CLR-transform, preserves unit-sum constraint and removes value-range restriction

geometric_all <- data_TFP_analysis %>%
  summarize_at(vars(Connector:Not_a_Bridge), funs(geomMeanExtension(., epsilon = 1e-1)))

Radar_functional_role_geometric_age <- data_TFP_analysis %>%
  dplyr::select(-Age) %>%
  filter(grepl("Young|Middle|Old", Age_group)) %>%
  group_by(Age_group) %>%
  summarize_at(vars(Connector:Not_a_Bridge), funs(geomMeanExtension(., epsilon = 1e-1))) %>%
  mutate(Connector = log(Connector / geometric_all$Connector)) %>%
  mutate(Provincial = log(Provincial / geometric_all$Provincial)) %>%
  mutate(Satellite = log(Satellite / geometric_all$Satellite)) %>%
  mutate(Peripheral = log(Peripheral / geometric_all$Peripheral)) %>%
  mutate(Global_Bridge = log(Global_Bridge / geometric_all$Global_Bridge)) %>%
  mutate(Local_Bridge = log(Local_Bridge / geometric_all$Local_Bridge)) %>%
  mutate(Super_Bridge = log(Super_Bridge / geometric_all$Super_Bridge)) %>%
  mutate(Not_a_Bridge = log(Not_a_Bridge / geometric_all$Not_a_Bridge)) %>%
  remove_rownames() %>%
  column_to_rownames(var = "Age_group")


palette <- RColorBrewer::brewer.pal(3, "YlOrRd")

radarplotting_overlap(Radar_functional_role_geometric_age, 0.4, -0.4, 1, 1,
  alpha = 0.1, label_size = 1, round = FALSE,
  title_fill = "Lifespan Topologico-functional profile (log ratio of geometric means)",
  palette = palette
)

legend(
  x = "topright",
  legend = rownames(Radar_functional_role_geometric_age), horiz = TRUE,
  bty = "n", pch = 20, col = palette,
  text.col = "black", cex = 1, pt.cex = 2
)

Radar_functional_role_age <- data_TFP_analysis %>%
  dplyr::select(-Age) %>%
  filter(grepl("Young|Middle|Old", Age_group)) %>%
  group_by(Age_group) %>%
  summarize_at(vars(Connector:Super_Bridge), funs(mean(.))) %>%
  remove_rownames() %>%
  column_to_rownames(var = "Age_group")

radarplotting_overlap(Radar_functional_role_age, 50, 10, 1, 1,
  alpha = 0.1, label_size = 1,
  title_fill = "Profile expressed in relative proportions",
  palette = palette
)

legend(
  x = "topright",
  legend = rownames(Radar_functional_role_age), horiz = TRUE,
  bty = "n", pch = 20, col = palette,
  text.col = "black", cex = 1, pt.cex = 2
)

# Globalisation & Stabilization ----


# glob <- data_TFP_analysis %>%
#   dplyr::select(-Age) %>%
#   filter(grepl("Young|Middle", Age_group)) %>%
#   group_by(Age_group) %>%
#   summarize_at(vars(Connector:Super_Bridge), funs(geomMeanExtension(., epsilon = 1e-1))) %>%
#   t(.) %>%
#   as.data.frame() %>%
#   janitor::row_to_names(1) %>%
#   mutate_at(vars(everything()), funs(as.numeric(.))) %>%
#   mutate(Globalisation = log(Middle / Young)) %>%
#   dplyr::select(Globalisation) %>%
#   t(.) %>%
#   as.data.frame() %>%
#   dplyr::select(-Satellite)
# 
# radarplotting(glob, 0.4, -0.4, 1, 1,
#   alpha = 0.1, label_size = 1, round = FALSE,
#   palette = "orange"
# )
# 
# periph <- data_TFP_analysis %>%
#   dplyr::select(-Age) %>%
#   filter(grepl("Old|Middle", Age_group)) %>%
#   group_by(Age_group) %>%
#   summarize_at(vars(Connector:Super_Bridge), funs(geomMeanExtension(., epsilon = 1e-1))) %>%
#   t(.) %>%
#   as.data.frame() %>%
#   janitor::row_to_names(1) %>%
#   mutate_at(vars(everything()), funs(as.numeric(.))) %>%
#   mutate(Stabilisation = log(Old / Middle)) %>%
#   dplyr::select(Stabilisation) %>%
#   t(.) %>%
#   as.data.frame() %>%
#   dplyr::select(-Satellite)
# 
# radarplotting(periph, 0.4, -0.4, 1, 1,
#   alpha = 0.1, label_size = 1, round = FALSE,
#   palette = "red"
# )


################################################################################
# QUALITY CHECK -------------
################################################################################
# Distribution of hubs across RSNs for each cluster for the individual hubs ----

# First averaging per subject then per clusters because grand mean is not equal to mean of means
# with unequal sample size i.e., subjects have a different total number of hubs

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
# 
# # Distribution of hubs across communities for each Age_group for the individual hubs ----
# 
# Radar_hub_community <- tmp_cluster_final %>%
#   group_by(Subj_ID, Age_group, Consensus_vector_0.15) %>%
#   summarise(n = n()) %>%
#   mutate(freq = n / sum(n)) %>%
#   dplyr::select(-n) %>%
#   group_by(Age_group, Consensus_vector_0.15) %>%
#   summarize_at(vars(freq), mean) %>%
#   spread(Consensus_vector_0.15, freq) %>%
#   remove_rownames() %>%
#   column_to_rownames("Age_group") %>%
#   mutate_at(vars(everything()), funs(. * 100))
# 
# radarplotting_overlap(Radar_hub_community, 50, 0, 1, 1,
#   alpha = 0.2, label_size = 1,
#   title_fill = "Distribution of hubs regions across communities",
#   palette = RColorBrewer::brewer.pal(8, "Dark2")
# )
# 
# legend(
#   x = "topright",
#   legend = rownames(Radar_hub_community), horiz = TRUE,
#   bty = "n", pch = 20, col = RColorBrewer::brewer.pal(8, "Dark2"),
#   text.col = "black", cex = 1, pt.cex = 2
# )




