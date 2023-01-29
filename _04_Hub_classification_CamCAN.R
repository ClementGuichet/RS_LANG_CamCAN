##########################################################################################
# Script for functional hub classification with the data processed

# Written by CG - 2023
##########################################################################################
library(jsonlite)
library(data.table)
################################################################################
# Import processed data---------------------------------------------------------
source("_01_DataManip_CamCAN.R")

################################################################################
# ~~~~~~~~~~~ Hub classification ~~~~~~~~~~~ -----------------------------------
# 
# High Betweenness centrality = global bridge
# High Flow centrality = local bridge
# High Participation coefficient (based on consensus group-level modular decomposition after 1000 iterations)
# High Within-z (based on consensus group-level modular decomposition after 1000 iterations)

# High_zPC/High z = connector
# High_zPC/low_z = satellite
# Low_zPC/High_z = provincial
# Low_zPC/Low_z = peripheral

# 628 Subjects, 131 Regions, one threshold = 0.15 or OMST

################################################################################
################################################################################
################################################################################

# Consensus modular normalization with Age-group-specific community structures ----

Hub_classification_procedure <- function(filtering_scheme = NULL) {
  if (filtering_scheme == "OMST") {
    setwd(paste0(getwd(), "/OMST"))
    # YOUNG ----
    
    PC_consensus <- as.data.frame(fromJSON("Young_PC_norm.json")) %>%
      mutate(Subj_ID = rep(seq_len(167))) %>%
      pivot_longer(
        cols = !c("Subj_ID"),
        names_to = "Region",
        values_to = "PC_cons"
      )
    
    Within_module_z_consensus <- as.data.frame(fromJSON("Young_Wz.json")) %>%
      mutate(Subj_ID = rep(seq_len(167))) %>%
      pivot_longer(
        cols = !c("Subj_ID"),
        names_to = "Region",
        values_to = "Within_module_z_cons"
      )
    
    nodal_metrics_cons <- cbind(PC_consensus, Within_module_z_cons = Within_module_z_consensus$Within_module_z_cons) %>%
      dplyr::select(-Region)
    
    # Make sure dataframe is ordered identically to nodal_metrics
    data_full_thresholded <- data_full %>%
      # subset(threshold == "0.15") %>%
      filter(Age_group == "Young") %>%
      arrange(Subj_ID, Region)
    
    data_bind_PC_Wz <- cbind(data_full_thresholded,
                             PC_cons = nodal_metrics_cons$PC_cons, Within_module_z_cons = nodal_metrics_cons$Within_module_z_cons
    )
    
    data_young <- data_bind_PC_Wz 
    
    
    # MIDDLE ----
    
    PC_consensus <- as.data.frame(fromJSON("Middle_PC_norm.json")) %>% 
      mutate(Subj_ID = rep(seq_len(201))) %>%
      pivot_longer(
        cols = !c("Subj_ID"),
        names_to = "Region",
        values_to = "PC_cons"
      )
    
    Within_module_z_consensus <- as.data.frame(fromJSON("Middle_Wz.json")) %>%
      mutate(Subj_ID = rep(seq_len(201))) %>%
      pivot_longer(
        cols = !c("Subj_ID"),
        names_to = "Region",
        values_to = "Within_module_z_cons"
      )
    
    nodal_metrics_cons <- cbind(PC_consensus, Within_module_z_cons = Within_module_z_consensus$Within_module_z_cons) %>%
      dplyr::select(-Region)
    
    # Make sure dataframe is ordered identically to nodal_metrics
    data_full_thresholded <- data_full %>%
      # subset(threshold == "0.15") %>%
      filter(Age_group == "Middle") %>%
      arrange(Subj_ID, Region)
    
    data_bind_PC_Wz <- cbind(data_full_thresholded,
                             PC_cons = nodal_metrics_cons$PC_cons, Within_module_z_cons = nodal_metrics_cons$Within_module_z_cons
    )
    
    data_middle <- data_bind_PC_Wz
    
    
    # OLD ----
    
    PC_consensus <- as.data.frame(fromJSON("Old_PC_norm.json")) %>%
      mutate(Subj_ID = rep(seq_len(260))) %>%
      pivot_longer(
        cols = !c("Subj_ID"),
        names_to = "Region",
        values_to = "PC_cons"
      )
    
    Within_module_z_consensus <- as.data.frame(fromJSON("Old_Wz.json")) %>%
      mutate(Subj_ID = rep(seq_len(260))) %>%
      pivot_longer(
        cols = !c("Subj_ID"),
        names_to = "Region",
        values_to = "Within_module_z_cons"
      )
    
    nodal_metrics_cons <- cbind(PC_consensus, Within_module_z_cons = Within_module_z_consensus$Within_module_z_cons) %>%
      dplyr::select(-Region)
    
    # Make sure dataframe is ordered identically to nodal_metrics
    data_full_thresholded <- data_full %>%
      # subset(threshold == "0.15") %>%
      filter(Age_group == "Old") %>%
      arrange(Subj_ID, Region)
    
    data_bind_PC_Wz <- cbind(data_full_thresholded,
                             PC_cons = nodal_metrics_cons$PC_cons, Within_module_z_cons = nodal_metrics_cons$Within_module_z_cons
    )
    
    
    data_old <- data_bind_PC_Wz
    
    
    # Putting it all together ----
    
    data_functional_role_tmp <<- rbind(data_young, data_middle, data_old) %>% 
      group_by(Subj_ID) %>% 
      mutate(zBT = as.numeric(scale(Betweenness))) %>%
      mutate(zFlow = as.numeric(scale(Flow_coeff))) %>%
      mutate(Bridgeness = ifelse(zBT > 0 & zFlow < 0, "Global_Bridge",
                                 ifelse(zFlow > 0 & zBT < 0, "Local_Bridge",
                                        ifelse(zBT > 0 & zFlow > 0, "Super_Bridge",
                                               ifelse(zBT < 0 & zFlow < 0, "Not_a_Bridge", 0)
                                        )
                                 )
      )) %>% 
      mutate(zPC_cons = as.numeric(scale(PC_cons))) %>%
      # mutate(zPC_cons = ifelse(zPC_cons == "NaN", 0, zPC_cons)) %>%
      # 1e-5 to avoid nodes with Wz = 0 to be classified as Connector or Provincial
      # 0 indicates that it forms its own module mathematically speaking
      mutate(Hub_consensus = ifelse(zPC_cons >= 0 & Within_module_z_cons >= 1e-5, "Connector",
                                    ifelse(zPC_cons >= 0 & Within_module_z_cons < 1e-5, "Satellite",
                                           ifelse(zPC_cons < 0 & Within_module_z_cons >= 1e-5, "Provincial",
                                                  ifelse(zPC_cons < 0 & Within_module_z_cons < 1e-5, "Peripheral", "Isolate")
                                           )
                                    )
      )) %>%
      relocate(Subj_ID, .after = "Hub_consensus") %>%
      arrange(Subj_ID, Region) %>%
      ungroup()
    
    data_functional_role_tmp$Hub_consensus <<- factor(data_functional_role_tmp$Hub_consensus, levels = c(
      "Connector", "Provincial", "Satellite", "Peripheral"
    ))
    
    data_functional_role_tmp$Bridgeness <<- factor(data_functional_role_tmp$Bridgeness, levels = c(
      "Global_Bridge", "Local_Bridge", "Super_Bridge", "Not_a_Bridge"
    ))
    
    
    setwd(str_replace(getwd(), "\\/OMST", ""))
  } else if (filtering_scheme == "proportional") {
    setwd(paste0(getwd(), "/data_graphvar_T1"))
    # YOUNG ----
    
    PC_consensus <- as.data.frame(fromJSON("Young_PC_norm.json")) %>%
      mutate(Subj_ID = rep(seq_len(167))) %>%
      pivot_longer(
        cols = !c("Subj_ID"),
        names_to = "Region",
        values_to = "PC_cons"
      )
    
    Within_module_z_consensus <- as.data.frame(fromJSON("Young_Wz.json")) %>%
      mutate(Subj_ID = rep(seq_len(167))) %>%
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
    
    data_young <- data_bind_PC_Wz 
    
    
    # MIDDLE ----
    
    PC_consensus <- as.data.frame(fromJSON("Middle_PC_norm.json")) %>% 
      mutate(Subj_ID = rep(seq_len(201))) %>%
      pivot_longer(
        cols = !c("Subj_ID"),
        names_to = "Region",
        values_to = "PC_cons"
      )
    
    Within_module_z_consensus <- as.data.frame(fromJSON("Middle_Wz.json")) %>%
      mutate(Subj_ID = rep(seq_len(201))) %>%
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
    
    data_middle <- data_bind_PC_Wz
    
    
    # OLD ----
    
    PC_consensus <- as.data.frame(fromJSON("Old_PC_norm.json")) %>%
      mutate(Subj_ID = rep(seq_len(260))) %>%
      pivot_longer(
        cols = !c("Subj_ID"),
        names_to = "Region",
        values_to = "PC_cons"
      )
    
    Within_module_z_consensus <- as.data.frame(fromJSON("Old_Wz.json")) %>%
      mutate(Subj_ID = rep(seq_len(260))) %>%
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
    
    
    data_old <- data_bind_PC_Wz
    
    
    # Putting it all together ----
    
    data_functional_role <<- rbind(data_young, data_middle, data_old) %>% 
      group_by(Subj_ID) %>% 
      mutate(zBT = as.numeric(scale(Betweenness))) %>%
      mutate(zFlow = as.numeric(scale(Flow_coeff))) %>%
      mutate(Bridgeness = ifelse(zBT >= 0 & zFlow < 0, "Global_Bridge",
                                 ifelse(zFlow >= 0 & zBT < 0, "Local_Bridge",
                                        ifelse(zBT >= 0 & zFlow >= 0, "Super_Bridge",
                                               ifelse(zBT < 0 & zFlow < 0, "Not_a_Bridge", 0)
                                        )
                                 )
      )) %>% 
      mutate(zPC_cons = as.numeric(scale(PC_cons))) %>%
      # mutate(zPC_cons = ifelse(zPC_cons == "NaN", 0, zPC_cons)) %>%
      # 1e-5 to avoid nodes with Wz = 0 to be classified as Connector or Provincial
      # 0 indicates that it forms its own module mathematically speaking
      mutate(Hub_consensus = ifelse(zPC_cons >= 0 & Within_module_z_cons >= 1e-5, "Connector",
                                    ifelse(zPC_cons >= 0 & Within_module_z_cons < 1e-5, "Satellite",
                                           ifelse(zPC_cons < 0 & Within_module_z_cons >= 1e-5, "Provincial",
                                                  ifelse(zPC_cons < 0 & Within_module_z_cons < 1e-5, "Peripheral", "Isolate")
                                           )
                                    )
      )) %>%
      relocate(Subj_ID, .after = "Hub_consensus") %>%
      arrange(Subj_ID, Region) %>%
      ungroup()
    
    data_functional_role$Hub_consensus <<- factor(data_functional_role$Hub_consensus, levels = c(
      "Connector", "Provincial", "Satellite", "Peripheral"
    ))
    
    data_functional_role$Bridgeness <<- factor(data_functional_role$Bridgeness, levels = c(
      "Global_Bridge", "Local_Bridge", "Super_Bridge", "Not_a_Bridge"
    ))
    
    
    setwd(str_replace(getwd(), "\\/data_graphvar_T1", ""))
  }
}

Hub_classification_procedure("OMST")

################################################################################
################################################################################
################################################################################

# Consensual modular normalization ----

# !!!!!!!!!! PC & Wz have been generated using GraphVar consensus affiliation vector !!!!!!!!!!!!!!!

# PC_consensus <- as.data.frame(fromJSON("All_PC.json")) %>%
#   mutate(Subj_ID = rep(seq_len(628))) %>%
#   pivot_longer(
#     cols = !c("Subj_ID"),
#     names_to = "Region",
#     values_to = "PC_cons"
#   )
# 
# Within_module_z_consensus <- as.data.frame(fromJSON("All_Wz.json")) %>%
#   mutate(Subj_ID = rep(seq_len(628))) %>%
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
#   group_by(Subj_ID) %>%
#   mutate(zBT = as.numeric(scale(Betweenness))) %>%
#   mutate(zFlow = as.numeric(scale(Flow_coeff))) %>%
#   mutate(Bridgeness = ifelse(zBT >= 0 & zFlow < 0, "Global_Bridge",
#     ifelse(zFlow >= 0 & zBT < 0, "Local_Bridge",
#       ifelse(zBT >= 0 & zFlow >= 0, "Super_Bridge",
#         ifelse(zBT < 0 & zFlow < 0, "Not_a_Bridge", 0)
#       )
#     )
#   )) %>%
#   # Normalizing at the community level with the affiliation vector from consensus clustering
#   group_by(Consensus_vector_0.15) %>%
#   mutate(zPC_cons = as.numeric(scale(PC_cons))) %>%
#   # To resolve scaling issue
#   mutate(zPC_cons = ifelse(zPC_cons == "NaN", 0, zPC_cons)) %>%
#   mutate(Hub_consensus = ifelse(zPC_cons >= 0 & Within_module_z_cons >= 1e-5, "Connector",
#     ifelse(zPC_cons >= 0 & Within_module_z_cons < 1e-5, "Satellite",
#       ifelse(zPC_cons < 0 & Within_module_z_cons >= 1e-5, "Provincial",
#         ifelse(zPC_cons < 0 & Within_module_z_cons < 1e-5, "Peripheral", "Isolate")
#       )
#     )
#   )) %>%
#   relocate(Subj_ID, .after = "Hub_consensus") %>%
#   arrange(Subj_ID, Region) %>%
#   ungroup()
# 
# data_functional_role$Hub_consensus <- factor(data_functional_role$Hub_consensus, levels = c(
#   "Connector", "Provincial", "Satellite", "Peripheral"
# ))
# 
# data_functional_role$Bridgeness <- factor(data_functional_role$Bridgeness, levels = c(
#   "Global_Bridge", "Local_Bridge", "Super_Bridge", "Not_a_Bridge"
# ))
# 
# 
# setwd(str_replace(getwd(), "\\/data_graphvar_T1", ""))
