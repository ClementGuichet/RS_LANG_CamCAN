##########################################################################################
# Script for individual hub detection

# Written by CG - 2023
##########################################################################################
source("_04_Hub_classification_CamCAN.R")

################################################################################
# ~~~~~~~~~~~ Hub Detection Procedure ~~~~~~~~~~~

Hub_detection_procedure <- function(filtering_scheme = NULL, percentage_hub_regions) {
  if (filtering_scheme == "OMST") {
    # Removing subjects with outlying global - cost efficiency
    setwd(paste0(getwd(), "/OMST"))
    cost_OMST <- as.data.frame(fromJSON("globalcosteff_max.json")) %>%
      mutate(Subj_ID = rep(seq_len(628))) %>% 
      rename(cost_OMST = `fromJSON("globalcosteff_max.json")`)
    setwd(str_replace(getwd(), "\\/OMST", ""))
    
    performance::check_outliers(cost_OMST[1])
    data_functional_role <- data_functional_role_tmp %>%
      filter(Subj_ID != "129") %>%
      mutate(helper_vector = rep(seq_len(627), each = 131)) 
    
    # Method ~ Detect top % regions for each metric ------------------------------
    Hub_detection <- data_functional_role %>% 
      dplyr::select(helper_vector, Region, 
                    degree, Within_module_z_cons, zPC_cons, zBT, zFlow,
                    `1st_network`, Consensus_OMST_all) %>%
      # mutate(across(degree:PC, ~ rank(-.x), .names = "{.col}_rank")) %>%
      pivot_longer(
        cols = !c("helper_vector","Region", "1st_network", "Consensus_OMST_all"),
        names_to = "Metric_name",
        values_to = "Metric_value"
      ) %>% 
      group_by(helper_vector, Metric_name) %>%
      group_split() %>%
      # Choosing 20% to achieve a middle ground between too much versus too little regions
      # Must have sufficient regions in each RSN for statistical analyses
      map_dfr(. %>% slice_max(Metric_value, n = 131 * percentage_hub_regions) %>%
                mutate(rank = rep(seq(1:length(Region))))) %>%
      group_by(helper_vector) %>%
      group_split()
    
    
    Hub_selection <<- list()
    FR_list <- list()
    for (i in 1:length(Hub_detection)) {
      Hub_df <- rbindlist(Hub_detection[i]) %>% distinct(Region, .keep_all = TRUE)
      # Here I subset the rows specific to each subject and their Hub regions
      tmp <- data_functional_role %>%
        filter(Region %in% Hub_df$Region) %>%
        filter(helper_vector == i) %>%
        dplyr::select(Subj_ID, helper_vector, Region, `1st_network`, Consensus_OMST_all, Hub_consensus, Bridgeness, zFlow, zBT)
      
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
      Hub_selection[[i]] <<- tmp
      # Topologico-functional profiles
      FR_list[[i]] <- FR_ind
    }
    
    # Adding a ratio score denoting the propensity for functional segregation over integration
    data_cluster_efficiency <- data_functional_role %>%
      dplyr::select(Subj_ID, Eglob, Eloc) %>%
      group_by(Subj_ID) %>%
      summarize_at(vars(Eglob, Eloc), mean) %>%
      mutate(Balance_eff = (Eloc - Eglob) / (Eloc + Eglob))
    
    ################################################################################
    # Putting everything together
    # TPF at the Subject-level
    TFP_General <<- cbind(
      rbindlist(FR_list, fill = TRUE) %>%
        mutate_all(., ~ replace(., is.na(.), 0)) %>% mutate_at(vars(everything()), funs(. * 100)),
      data_functional_role %>% group_by(Subj_ID, gender_text) %>% summarise_at(vars(Age), mean) %>% arrange(Subj_ID),
      Balance_eff = data_cluster_efficiency$Balance_eff,
      Eglob = data_cluster_efficiency$Eglob,
      Eloc = data_cluster_efficiency$Eloc
    )
    
    
    # Putting everything together
    # TPF at the Subject-level and the RSN-level
    TFP_RSN <<- cbind(
      rbindlist(FR_list, fill = TRUE) %>%
        mutate_all(., ~ replace(., is.na(.), 0)) %>% mutate_at(vars(everything()), funs(. * 100)),
      data_functional_role %>% group_by(Subj_ID, gender_text, `1st_network`) %>% summarise_at(vars(Age), mean) %>% arrange(Subj_ID, `1st_network`),
      Balance_eff = data_cluster_efficiency$Balance_eff
    )
  } else if (filtering_scheme == "proportional") {
    # Method ~ Detect top % regions for each metric ------------------------------
    Hub_detection <- data_functional_role %>% 
      dplyr::select(Subj_ID, Region, 
                    degree, Within_module_z_cons, zPC_cons, zBT, zFlow,
                    `1st_network`, Consensus_vector_0.15) %>%
      # mutate(across(degree:PC, ~ rank(-.x), .names = "{.col}_rank")) %>%
      pivot_longer(
        cols = !c("Subj_ID","Region", "1st_network", "Consensus_vector_0.15"),
        names_to = "Metric_name",
        values_to = "Metric_value"
      ) %>% 
      group_by(Subj_ID, Metric_name) %>%
      group_split() %>%
      # Choosing 20% to achieve a middle ground between too much versus too little regions
      # Must have sufficient regions in each RSN for statistical analyses
      map_dfr(. %>% slice_max(Metric_value, n = 131 * percentage_hub_regions) %>%
                mutate(rank = rep(seq(1:length(Region))))) %>%
      group_by(Subj_ID) %>%
      group_split()
    
    
    Hub_selection <<- list()
    FR_list <- list()
    for (i in 1:length(Hub_detection)) {
      Hub_df <- rbindlist(Hub_detection[i]) %>% distinct(Region, .keep_all = TRUE)
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
      Hub_selection[[i]] <<- tmp
      # Topologico-functional profiles
      FR_list[[i]] <- FR_ind
    }
    
    # Adding a ratio score denoting the propensity for functional segregation over integration
    data_cluster_efficiency <- data_functional_role %>%
      dplyr::select(Subj_ID, Eglob, Eloc) %>%
      group_by(Subj_ID) %>%
      summarize_at(vars(Eglob, Eloc), mean) %>%
      mutate(Balance_eff = (Eloc - Eglob) / (Eloc + Eglob))
    
    ################################################################################
    # Putting everything together
    # TPF at the Subject-level
    TFP_General <<- cbind(
      rbindlist(FR_list, fill = TRUE) %>%
        mutate_all(., ~ replace(., is.na(.), 0)) %>% mutate_at(vars(everything()), funs(. * 100)),
      data_functional_role %>% group_by(Subj_ID, gender_text) %>% summarise_at(vars(Age), mean) %>% arrange(Subj_ID),
      Balance_eff = data_cluster_efficiency$Balance_eff,
      Eglob = data_cluster_efficiency$Eglob,
      Eloc = data_cluster_efficiency$Eloc
    )
    
    
    # Putting everything together
    # TPF at the Subject-level and the RSN-level
    TFP_RSN <<- cbind(
      rbindlist(FR_list, fill = TRUE) %>%
        mutate_all(., ~ replace(., is.na(.), 0)) %>% mutate_at(vars(everything()), funs(. * 100)),
      data_functional_role %>% group_by(Subj_ID, gender_text, `1st_network`) %>% summarise_at(vars(Age), mean) %>% arrange(Subj_ID, `1st_network`),
      Balance_eff = data_cluster_efficiency$Balance_eff
    )
  }
}


Hub_detection_procedure("OMST", .2)
################################################################################

