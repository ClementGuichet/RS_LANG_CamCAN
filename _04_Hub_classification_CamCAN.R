##########################################################################################
# Script for functional hub classification with the data processed

# Written by CG
# 26-11-2022
##########################################################################################
library(jsonlite) # for working with json files
library(emmeans) # for post-hoc tests
library(data.table) # for working with lists
library(FactoMineR)

################################################################################
# Import processed data---------------------------------------------------------
source("_01_DataManip_CamCAN.R")

################################################################################
# ~~~~~~~~~~~ Hub classification ~~~~~~~~~~~ -----------------------------------

# High Betweenness centrality = global bridge
# High Flow centrality = local bridge
# High Participation coefficient (based on consensus group-level modular decomposition after 1000 iterations)
# High Within-z (based on consensus group-level modular decomposition after 1000 iterations)

# High_zPC/High z = connector
# High_zPC/low_z = satellite
# Low_zPC/High_z = provincial
# Low_zPC/Low_z = peripheral

# 72 Subjects, 131 Regions, one threshold = 0.15




################################################################################
# Consensual modular normalization ----

# !!!!!!!!!! PC & Wz have been generated using GraphVar consensus affiliation vector !!!!!!!!!!!!!!!

# PC_consensus <- as.data.frame(fromJSON("Participation_coefficient_consensus.json")) %>%
#   mutate(Subj_ID = rep(seq_len(645))) %>%
#   pivot_longer(
#     cols = !c("Subj_ID"),
#     names_to = "Region",
#     values_to = "PC_cons"
#   )
#
# Within_module_z_consensus <- as.data.frame(fromJSON("Within_module_z_score_consensus.json")) %>%
#   mutate(Subj_ID = rep(seq_len(645))) %>%
#   pivot_longer(
#     cols = !c("Subj_ID"),
#     names_to = "Region",
#     values_to = "Within_module_z_cons"
#   )
#
# nodal_metrics_cons <- cbind(PC_consensus, Within_module_z_cons = Within_module_z_consensus$Within_module_z_cons) %>%
#   dplyr::select(-Region)
#
# # Make sure dataframe is ordered identically to nodal_metrics
# data_full_thresholded <- data_full %>%
#   subset(threshold == "0.15") %>%
#   arrange(Subj_ID, Region)
#
# data_bind_PC_Wz <- cbind(data_full_thresholded,
#   PC_cons = nodal_metrics_cons$PC_cons, Within_module_z_cons = nodal_metrics_cons$Within_module_z_cons
# )
#
# 
#
# data_functional_role <- data_bind_PC_Wz %>%
#   # Normalizing at the connectomic level
#   # group_by(Subj_ID) %>%
#   mutate(zK = as.numeric(scale(degree))) %>%
#   mutate(zBT = as.numeric(scale(Betweenness))) %>%
#   mutate(zFlow = as.numeric(scale(Flow_coeff))) %>%
#   mutate(Bridgeness = ifelse(zBT > 0 & zFlow < 0, "Global_Bridge",
#     ifelse(zFlow > 0 & zBT < 0, "Local_Bridge",
#       ifelse(zBT > 0 & zFlow > 0, "Super_Bridge",
#         ifelse(zBT < 0 & zFlow < 0, "Not_a_Bridge", 0)
#       )
#     )
#   )) %>%
#   # Normalizing at the community level with the affiliation vector from consensus clustering
#   group_by(Consensus_vector_0.15) %>%
#   mutate(zPC_cons = as.numeric(scale(PC_cons))) %>%
#   # To resolve scaling issue
#   mutate(zPC_cons = ifelse(zPC_cons == "NaN", 0, zPC_cons)) %>%
#   mutate(Hub_consensus = ifelse(zPC_cons > 1e-5 & Within_module_z_cons > 1e-5, "Connector",
#     ifelse(zPC_cons > 1e-5 & Within_module_z_cons < 1e-5, "Satellite",
#       ifelse(zPC_cons < 1e-5 & Within_module_z_cons > 1e-5, "Provincial",
#         ifelse(zPC_cons < 1e-5 & Within_module_z_cons < 1e-5, "Peripheral", "Isolate")
#       )
#     )
#   )) %>%
#   relocate(Subj_ID, .after = "Hub_consensus") %>%
#   arrange(Subj_ID, Region) %>%
#   ungroup()
#


# Consensus modular normalization with Age-group-specific community structures ----
# YOUNG ----

PC_consensus <- as.data.frame(fromJSON("young_PC.json")) %>%
  mutate(Subj_ID = rep(seq_len(169))) %>%
  pivot_longer(
    cols = !c("Subj_ID"),
    names_to = "Region",
    values_to = "PC_cons"
  )

Within_module_z_consensus <- as.data.frame(fromJSON("young_Wz.json")) %>%
  mutate(Subj_ID = rep(seq_len(169))) %>%
  pivot_longer(
    cols = !c("Subj_ID"),
    names_to = "Region",
    values_to = "Within_module_z_cons"
  )

nodal_metrics_cons <- cbind(PC_consensus, Within_module_z_cons = Within_module_z_consensus$Within_module_z_cons) %>%
  dplyr::select(-Region)

# Make sure dataframe is ordered identically to nodal_metrics
data_full_thresholded <- data_full %>%
  subset(threshold == "0.15") %>%
  filter(Age_group == "Young") %>%
  arrange(Subj_ID, Region)

data_bind_PC_Wz <- cbind(data_full_thresholded,
  PC_cons = nodal_metrics_cons$PC_cons, Within_module_z_cons = nodal_metrics_cons$Within_module_z_cons
)

data_young <- data_bind_PC_Wz %>%
  # Normalizing at the community level with the affiliation vector from consensus clustering
  group_by(Consensus_young) %>%
  mutate(zPC_cons = as.numeric(scale(PC_cons))) %>%
  # To resolve scaling issue
  mutate(zPC_cons = ifelse(zPC_cons == "NaN", 0, zPC_cons)) %>%
  arrange(Subj_ID, Region) %>%
  ungroup()

# MIDDLE ----

PC_consensus <- as.data.frame(fromJSON("middle_PC.json")) %>%
  mutate(Subj_ID = rep(seq_len(202))) %>%
  pivot_longer(
    cols = !c("Subj_ID"),
    names_to = "Region",
    values_to = "PC_cons"
  )

Within_module_z_consensus <- as.data.frame(fromJSON("middle_Wz.json")) %>%
  mutate(Subj_ID = rep(seq_len(202))) %>%
  pivot_longer(
    cols = !c("Subj_ID"),
    names_to = "Region",
    values_to = "Within_module_z_cons"
  )

nodal_metrics_cons <- cbind(PC_consensus, Within_module_z_cons = Within_module_z_consensus$Within_module_z_cons) %>%
  dplyr::select(-Region)

# Make sure dataframe is ordered identically to nodal_metrics
data_full_thresholded <- data_full %>%
  subset(threshold == "0.15") %>%
  filter(Age_group == "Middle") %>%
  arrange(Subj_ID, Region)

data_bind_PC_Wz <- cbind(data_full_thresholded,
  PC_cons = nodal_metrics_cons$PC_cons, Within_module_z_cons = nodal_metrics_cons$Within_module_z_cons
)

data_middle <- data_bind_PC_Wz %>%
  # Normalizing at the community level with the affiliation vector from consensus clustering
  group_by(Consensus_middle) %>%
  mutate(zPC_cons = as.numeric(scale(PC_cons))) %>%
  # To resolve scaling issue
  mutate(zPC_cons = ifelse(zPC_cons == "NaN", 0, zPC_cons)) %>%
  arrange(Subj_ID, Region) %>%
  ungroup()


# OLD ----

PC_consensus <- as.data.frame(fromJSON("old_PC.json")) %>%
  mutate(Subj_ID = rep(seq_len(274))) %>%
  pivot_longer(
    cols = !c("Subj_ID"),
    names_to = "Region",
    values_to = "PC_cons"
  )

Within_module_z_consensus <- as.data.frame(fromJSON("old_Wz.json")) %>%
  mutate(Subj_ID = rep(seq_len(274))) %>%
  pivot_longer(
    cols = !c("Subj_ID"),
    names_to = "Region",
    values_to = "Within_module_z_cons"
  )

nodal_metrics_cons <- cbind(PC_consensus, Within_module_z_cons = Within_module_z_consensus$Within_module_z_cons) %>%
  dplyr::select(-Region)

# Make sure dataframe is ordered identically to nodal_metrics
data_full_thresholded <- data_full %>%
  subset(threshold == "0.15") %>%
  filter(Age_group == "Old") %>%
  arrange(Subj_ID, Region)

data_bind_PC_Wz <- cbind(data_full_thresholded,
  PC_cons = nodal_metrics_cons$PC_cons, Within_module_z_cons = nodal_metrics_cons$Within_module_z_cons
)


data_old <- data_bind_PC_Wz %>%
  # Normalizing at the community level with the affiliation vector from consensus clustering
  group_by(Consensus_old) %>%
  mutate(zPC_cons = as.numeric(scale(PC_cons))) %>%
  # To resolve scaling issue
  mutate(zPC_cons = ifelse(zPC_cons == "NaN", 0, zPC_cons)) %>%
  arrange(Subj_ID, Region) %>%
  ungroup()


# Putting it all together

data_functional_role <- rbind(data_young, data_middle, data_old) %>%
  # Normalizing at the connectomic level
  mutate(zK = as.numeric(scale(degree))) %>%
  mutate(zBT = as.numeric(scale(Betweenness))) %>%
  mutate(zFlow = as.numeric(scale(Flow_coeff))) %>%
  mutate(Bridgeness = ifelse(zBT > 0 & zFlow < 0, "Global_Bridge",
    ifelse(zFlow > 0 & zBT < 0, "Local_Bridge",
      ifelse(zBT > 0 & zFlow > 0, "Super_Bridge",
        ifelse(zBT < 0 & zFlow < 0, "Not_a_Bridge", 0)
      )
    )
  )) %>%
  # 1e-5 to avoid nodes with Wz = 0 to be classified as Connector or Provincial
  # 0 indicates that it forms its own module mathematically speaking
  mutate(Hub_consensus = ifelse(zPC_cons > 1e-5 & Within_module_z_cons > 1e-5, "Connector",
    ifelse(zPC_cons > 1e-5 & Within_module_z_cons < 1e-5, "Satellite",
      ifelse(zPC_cons < 1e-5 & Within_module_z_cons > 1e-5, "Provincial",
        ifelse(zPC_cons < 1e-5 & Within_module_z_cons < 1e-5, "Peripheral", "Isolate")
      )
    )
  )) %>%
  relocate(Subj_ID, .after = "Hub_consensus") %>%
  arrange(Subj_ID, Region) %>%
  ungroup()

data_functional_role$Hub_consensus <- factor(data_functional_role$Hub_consensus, levels = c(
  "Connector", "Provincial", "Satellite", "Peripheral"
))

data_functional_role$Bridgeness <- factor(data_functional_role$Bridgeness, levels = c(
  "Global_Bridge", "Local_Bridge", "Super_Bridge", "Not_a_Bridge"
))

################################################################################
# Normalizing by individual modular structure ignore the shared information between individuals
# We must use a consensus partitioning to reduce noise and achieve better representativity (Bian et al., 2023)
# The following code was not implemented in the final analysis

# Individual modular normalization ----
# Participation coefficient & Wz
# listfile <- list.files(getwd(), pattern = "*.txt")
# listfile_PC <- listfile[grep("participation", listfile)]
#
# PC_ind <- ldply(listfile_PC, read.table, header = T, sep = "\t") %>%
#   plyr::rename(c("X" = "Subj_ID")) %>%
#   replace("Subj_ID", seq_len(645)) %>%
#   pivot_longer(
#     cols = !c("Subj_ID"),
#     names_to = "Region",
#     values_to = "PC_ind"
#   )
#
# listfile_Wz <- listfile[grep("module", listfile)]
#
# Wz_ind <- ldply(listfile_Wz, read.table, header = T, sep = "\t") %>%
#   plyr::rename(c("X" = "Subj_ID")) %>%
#   replace("Subj_ID", seq_len(645)) %>%
#   pivot_longer(
#     cols = !c("Subj_ID"),
#     names_to = "Region",
#     values_to = "Wz_ind"
#   )
#
# listfile_affiliation <- listfile[grep("louvain", listfile)]
#
# affiliation_vector_ind <- ldply(listfile_affiliation, read.table, header = T, sep = "\t") %>%
#   plyr::rename(c("X" = "Subj_ID")) %>%
#   replace("Subj_ID", seq_len(645)) %>%
#   pivot_longer(
#     cols = !c("Subj_ID"),
#     names_to = "Region",
#     values_to = "affiliation_vector_ind"
#   )
#
#
# data_modular_thresholded <- list(PC_ind, Wz_ind, affiliation_vector_ind) %>%
#   reduce(full_join, by = c("Subj_ID", "Region")) %>%
#   mutate(threshold = rep(.15)) %>%
#   group_by(Subj_ID) %>%
#   group_split()
#
# data_modular_list <- list()
# for (i in 1:length(data_modular_thresholded)) {
#   tmp <- rbindlist(data_modular_thresholded[i]) %>% as.data.frame()
#   tmp_bis <- tmp %>%
#     group_by(affiliation_vector_ind) %>%
#     mutate(zPC_ind = as.numeric(scale(PC_ind))) %>%
#     # To resolve scaling issue
#     mutate(zPC_ind = ifelse(zPC_ind == "NaN", 0, zPC_ind))
#   data_modular_list[[i]] <- tmp_bis
# }
#
# data_full_thresholded <- data_full %>%
#   subset(threshold == "0.15") %>%
#   arrange(Subj_ID, Region)
#
# data_functional_role <- rbindlist(data_modular_list) %>% as.data.frame() %>%
#   dplyr::select(-threshold) %>%
#   merge(., data_full_thresholded, by = c("Subj_ID", "Region")) %>%
#   mutate(Hub_consensus = ifelse(zPC_ind > 1e-5 & Wz_ind > 1e-5, "Connector",
#                                 ifelse(zPC_ind > 1e-5 & Wz_ind < 1e-5, "Satellite",
#                                        ifelse(zPC_ind < 1e-5 & Wz_ind > 1e-5, "Provincial",
#                                               ifelse(zPC_ind < 1e-5 & Wz_ind < 1e-5, "Peripheral", "Isolate")
#                                        )
#                                 )
#   )) %>%
#   # Normalizing at the connectomic level
#   mutate(zK = as.numeric(scale(degree))) %>%
#   mutate(zBT = as.numeric(scale(Betweenness))) %>%
#   mutate(zFlow = as.numeric(scale(Flow_coeff))) %>%
#   mutate(Bridgeness = ifelse(zBT > 0 & zFlow < 0, "Global_Bridge",
#                              ifelse(zFlow > 0 & zBT < 0, "Local_Bridge",
#                                     ifelse(zBT > 0 & zFlow > 0, "Super_Bridge",
#                                            ifelse(zBT < 0 & zFlow < 0, "Not_a_Bridge", 0)
#                                     )
#                              )
#   ))


# # Multilevel PCA ----
# source("_multilevel_pca.R")
#
# validation_mPCA <- data_functional_role %>% dplyr::select(Subj_ID, zPC_cons, Within_module_z_cons, zBT, zFlow)
#
# mpca <- multilevel_pca(scale(validation_mPCA %>% as.matrix()), id = validation_mPCA$Subj_ID,
#                        twoway = TRUE, cov.method = "m1", pve = 0.9)
# mpca$evectors
# mpca$varofTot
#
#
# # validation_PCA <- validation_mPCA %>% group_by(Subj_ID) %>%
# #   summarize_at(vars(everything()), funs(mean))
# #
# # results <- FactoMineR::PCA(scale(validation_PCA %>% dplyr::select(-Subj_ID)))
# # results$eig
# # results$var
#
# library(mixOmics)
#
# pca.result_multi <- mixOmics::pca(validation_mPCA %>% dplyr::select(-Subj_ID),
#                             multilevel = validation_mPCA$Subj_ID, scale = TRUE)
#
# plotIndiv(pca.result_multi, ind.names = validation_mPCA$Subj_ID)
# plotVar(pca.result_multi)
