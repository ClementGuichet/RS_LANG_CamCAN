##########################################################################################
# Script for individual hub detection

# Written by CG
# 26-11-2022
##########################################################################################
library(data.table)
library(ggpubr)

source("_04_Hub_classification_CamCAN.R")

# What are the graph-based biomarkers of health aging ?
# Does the topological_functional profile, that is the proportion of the hubs' functional role, change throughout life ?

################################################################################
# ~~~~~~~~~~~ Hub Detection Procedure ~~~~~~~~~~~
################################################################################

# Method ~ Detect top % regions for each metric ------------------------------

Top_metric_Age_ind <- data_functional_role %>%
  group_by(Subj_ID, Region) %>%
  # summarize_at(vars(zK, Wz_ind, zPC_ind, zBT, zFlow), mean) %>%
  summarize_at(vars(zK, Within_module_z_cons, zPC_cons, zBT, zFlow), mean) %>%
  # mutate(across(degree:PC, ~ rank(-.x), .names = "{.col}_rank")) %>%
  pivot_longer(
    cols = !c("Subj_ID", "Region"),
    names_to = "Metric_name",
    values_to = "Metric_value"
  ) %>%
  ungroup() %>%
  group_by(Subj_ID, Metric_name) %>%
  group_split() %>%
  # Choosing 20% to achieve a middle ground between too much versus too little regions
  # Must have sufficient regions in each RSN for statistical analyses
  map_dfr(. %>% slice_max(Metric_value, n = 131 * 0.2) %>%
    mutate(rank = rep(seq(1:length(Region))))) %>%
  group_by(Subj_ID) %>%
  group_split()


Hub_selection <- list()
FR_list <- list()
for (i in 1:length(Top_metric_Age_ind)) {
  Hub_df <- rbindlist(Top_metric_Age_ind[i]) %>% distinct(Region, .keep_all = TRUE)
  # Here I subset the rows specific to each subject and their Hub regions
  tmp <- data_functional_role %>%
    filter(Region %in% Hub_df$Region) %>%
    filter(Subj_ID == i) %>%
    dplyr::select(Subj_ID, Region, `1st_network`, Consensus_vector_0.15, Hub_consensus, Bridgeness, zFlow, zBT)

  # Here I compute the proportion of functional roles regarding centrality and information flow
  # with the hubs specific to each subject
  FR_ind_hub <- tmp %>%
    group_by(Hub_consensus) %>%
    summarize(n = n()) %>%
    mutate(freq = n / sum(n)) %>%
    dplyr::select(-n) %>%
    spread(Hub_consensus, freq)
  FR_ind_bridge <- tmp %>%
    group_by(Bridgeness) %>%
    summarize(n = n()) %>%
    mutate(freq = n / sum(n)) %>%
    dplyr::select(-n) %>%
    spread(Bridgeness, freq)

  FR_ind <- cbind(FR_ind_hub, FR_ind_bridge)
  
  # Hub region of each subject
  Hub_selection[[i]] <- tmp
  # Topologico-functional profiles
  FR_list[[i]] <- FR_ind
}

# Adding a ratio score denoting propensity for segregation over integration
data_cluster_efficiency <- data_functional_role %>%
  dplyr::select(Subj_ID, Eglob, Eloc) %>%
  group_by(Subj_ID) %>%
  summarize_at(vars(Eglob, Eloc), mean) %>%
  mutate(Balance_eff = (Eloc - Eglob) / (Eloc + Eglob)) %>%
  dplyr::select(-c(Subj_ID, Eglob, Eloc)) %>%
  mutate_at(vars(Balance_eff), funs(. * 100))

################################################################################
# Putting everything together
# TPF at the Subject-level
TFP_General <- cbind(
  rbindlist(FR_list, fill = TRUE) %>%
    mutate_all(., ~ replace(., is.na(.), 0)) %>% mutate_at(vars(everything()), funs(. * 100)),
  data_functional_role %>% group_by(Subj_ID, gender_text) %>% summarise_at(vars(Age), mean) %>% arrange(Subj_ID),
  Balance_eff = data_cluster_efficiency$Balance_eff
)


# Putting everything together
# TPF at the Subject-level and the RSN-level
TFP_RSN <- cbind(
  rbindlist(FR_list, fill = TRUE) %>%
    mutate_all(., ~ replace(., is.na(.), 0)) %>% mutate_at(vars(everything()), funs(. * 100)),
  data_functional_role %>% group_by(Subj_ID, gender_text, `1st_network`) %>% summarise_at(vars(Age), mean) %>% arrange(Subj_ID, `1st_network`),
  Balance_eff = data_cluster_efficiency$Balance_eff
)


