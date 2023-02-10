##########################################################################################
# Script for individual hub detection

# Written by CG - 2023
##########################################################################################
source("_04_Hub_classification_CamCAN.R")
source("_geometricmeanCruz.R")

################################################################################
# ~~~~~~~~~~~ Hub Detection Procedure ~~~~~~~~~~~

Hub_detection_procedure <- function(filtering_scheme = NULL, percentage_hub_regions) {
  if (filtering_scheme == "OMST") {
    # Removing subjects with outlying global - cost efficiency
    setwd(paste0(getwd(), "/OMST"))
    cost_OMST <<- as.data.frame(fromJSON("globalcosteff_max.json")) %>%
      mutate(Subj_ID = rep(seq_len(628))) %>% 
      rename(cost_OMST = `fromJSON("globalcosteff_max.json")`)
    setwd(str_replace(getwd(), "\\/OMST", ""))
    
    performance::check_outliers(cost_OMST[1])
    data_functional_role <<- data_functional_role_tmp %>%
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
        dplyr::select(Subj_ID, Age, Region, `1st_network`, Consensus_vector_0.15, MODULAR, INTERNODAL)
      
      # Here I compute the proportion of functional roles regarding centrality and information flow
      # with the hubs specific to each subject
      FR_ind_modular <- tmp %>%
        group_by(MODULAR) %>%
        summarize(n = n()) %>%
        mutate(freq = n / sum(n)) %>%
        dplyr::select(-n) %>%
        spread(MODULAR, freq)
      FR_ind_internodal <- tmp %>%
        group_by(INTERNODAL) %>%
        summarize(n = n()) %>%
        mutate(freq = n / sum(n)) %>%
        dplyr::select(-n) %>%
        spread(INTERNODAL, freq)
      
      FR_ind <- cbind(FR_ind_modular, FR_ind_internodal)
      
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
      mutate(Balance_eff = (Eglob - Eloc) / (Eloc + Eglob))
    
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
    ) %>% 
      # Remove subject with fragmented graph
      filter(Subj_ID %in% LLC_filter$Subj_ID) %>% 
      plyr::rename(c("gender_text" = "Gender"))
    
    # Final dataframe with each subject's hub regions and the RSNs -----
    # Hub region specific to each subject yielded by hub detection procedure
    data_hub_selection_per_subject <- rbindlist(Hub_selection)
    # Select the subjects from the clusters
    tmp_cluster_final <<- filter(
      data_hub_selection_per_subject,
      Subj_ID %in% TFP_General$Subj_ID
    ) %>% 
      mutate(Age_decade = ifelse(Age <= 29, 25, 
                                 ifelse(Age <= 39, 35,
                                        ifelse(Age <= 49, 45,
                                               ifelse(Age <= 59, 55,
                                                      ifelse(Age <= 69, 65,
                                                             ifelse(Age <= 79, 75, 85)))))))
    
    ################################################################################
    # TPF at the Subject-level & RSN-level
    epsilon <- 1e-1
    
    trajectory_modular <- tmp_cluster_final %>%
      group_by(`1st_network`, Region, Subj_ID, Age, MODULAR) %>%
      summarise(n = n()) %>%
      mutate(freq = n / sum(n)) %>%
      spread(MODULAR, freq) %>%
      dplyr::select(-n) %>%
      mutate_all(., ~ replace(., is.na(.), 0)) %>%
      group_by(`1st_network`, Subj_ID, Age) %>%
      summarize_at(vars(Connector, Provincial, Satellite, Peripheral), mean) %>%
      ungroup() %>%
      pivot_longer(cols = !c("1st_network","Subj_ID", "Age"), names_to = "Functional_role", values_to = "freq") %>% 
      group_by(`1st_network`, Subj_ID, Age, Functional_role) %>%
      summarise_at(vars(freq), funs(geomMeanExtension(., epsilon = epsilon))) %>% 
      spread(Functional_role, freq)
    
    trajectory_internodal <- tmp_cluster_final %>%
      group_by(`1st_network`, Region, Subj_ID, Age, INTERNODAL) %>%
      summarise(n = n()) %>%
      mutate(freq = n / sum(n)) %>%
      spread(INTERNODAL, freq) %>%
      dplyr::select(-n) %>%
      mutate_all(., ~ replace(., is.na(.), 0)) %>%
      group_by(`1st_network`, Subj_ID, Age) %>%
      summarize_at(vars(Global_Bridge, Local_Bridge, Super_Bridge, Not_a_Bridge), mean) %>%
      ungroup() %>%
      pivot_longer(cols = !c("1st_network", "Subj_ID",  "Age"), names_to = "Functional_role", values_to = "freq") %>% 
      group_by(`1st_network`, Subj_ID, Age, Functional_role) %>%
      summarise_at(vars(freq), funs(geomMeanExtension(., epsilon = epsilon))) %>% 
      spread(Functional_role, freq)
    
    TFP_RSN <<- merge(trajectory_modular, trajectory_internodal, by = c("Subj_ID", "1st_network", "Age")) %>% 
      mutate_at(vars(Connector:Super_Bridge), funs(as.numeric(. * 100)))
    
    # TPF at the Age-level & RSN-level
    epsilon <- 1e-1
    
    trajectory_modular <- tmp_cluster_final %>%
      group_by(`1st_network`, Region, Subj_ID, Age_decade, MODULAR) %>%
      summarise(n = n()) %>%
      mutate(freq = n / sum(n)) %>%
      spread(MODULAR, freq) %>%
      dplyr::select(-n) %>%
      mutate_all(., ~ replace(., is.na(.), 0)) %>%
      group_by(`1st_network`, Subj_ID, Age_decade) %>%
      summarize_at(vars(Connector, Provincial, Satellite, Peripheral), mean) %>%
      ungroup() %>%
      pivot_longer(cols = !c("1st_network","Subj_ID", "Age_decade"), names_to = "Functional_role", values_to = "freq") %>% 
      group_by(`1st_network`, Age_decade, Functional_role) %>%
      summarise_at(vars(freq), funs(geomMeanExtension(., epsilon = epsilon))) %>% 
      spread(Functional_role, freq)
    
    trajectory_internodal <- tmp_cluster_final %>%
      group_by(`1st_network`, Region, Subj_ID, Age_decade, INTERNODAL) %>%
      summarise(n = n()) %>%
      mutate(freq = n / sum(n)) %>%
      spread(INTERNODAL, freq) %>%
      dplyr::select(-n) %>%
      mutate_all(., ~ replace(., is.na(.), 0)) %>%
      group_by(`1st_network`, Subj_ID, Age_decade) %>%
      summarize_at(vars(Global_Bridge, Local_Bridge, Super_Bridge, Not_a_Bridge), mean) %>%
      ungroup() %>%
      pivot_longer(cols = !c("1st_network", "Subj_ID",  "Age_decade"), names_to = "Functional_role", values_to = "freq") %>% 
      group_by(`1st_network`, Age_decade, Functional_role) %>%
      summarise_at(vars(freq), funs(geomMeanExtension(., epsilon = epsilon))) %>% 
      spread(Functional_role, freq)
    
    TFP_RSN_Age_decade <<- merge(trajectory_modular, trajectory_internodal, by = c("1st_network", "Age_decade")) %>% 
      mutate_at(vars(Connector:Super_Bridge), funs(as.numeric(. * 100)))
  }
}
Hub_detection_procedure("proportional", .2)
################################################################################

