################################################################################
# LOG-RATIO of the Topologico-functional profile at the RSN level --------------
################################################################################
source("_06_TFProfile_CamCAN.R")
source("_radarplotting_function.R")
source("_geometricmeanCruz.R")

palette <- RColorBrewer::brewer.pal(3, "YlOrRd")

# Ratio of the proportion of each functional role within each RSN between the two selected clusters
interaction_Age_FuncRole_RSN <- function(max, min, max2, min2, alpha, RSN_modular, RSN_interareal, epsilon, round) {
  # Get the associated Resting-state networks for all age groups
  # Hub region specific to each subject yielded by hub detection procedure
  data_hub_selection_per_subject <- rbindlist(Hub_selection)
  tmp_cluster_final <<- merge(data_hub_selection_per_subject, data_TFP_analysis %>%
                                dplyr::select(Subj_ID, Age_group),
                              by = "Subj_ID"
  )
  tmp_cluster_final$`1st_network` <- factor(tmp_cluster_final$`1st_network`,levels = 
                                              c("Auditory", "Language",
                                                "DMN", "FPN", 
                                                "SMN", "Visual_1", "VMM",
                                                "DAN", "PMM",
                                                "CON", "Visual_2"))

  ##############################################################################
  # MODULAR 
  ##############################################################################

  # Proportion of each functional role within each RSN for each subject
  trajectory_modular <<- tmp_cluster_final %>%
    group_by(`1st_network`, Region, Subj_ID, Age_group, Hub_consensus) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n)) %>%
    spread(Hub_consensus, freq) %>%
    dplyr::select(-n) %>%
    mutate_all(., ~ replace(., is.na(.), 0)) %>%
    group_by(`1st_network`, Subj_ID, Age_group) %>%
    summarize_at(vars(Connector, Provincial, Satellite, Peripheral), mean) %>% 
    group_by(`1st_network`, Age_group) %>%
    summarise_at(vars(Connector, Provincial, Satellite, Peripheral), funs(geomMeanExtension(., epsilon = epsilon)))

  geometric_trajectory_modular <<- tmp_cluster_final %>%
    group_by(`1st_network`, Region, Subj_ID, Age_group, Hub_consensus) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n)) %>%
    spread(Hub_consensus, freq) %>%
    dplyr::select(-n) %>%
    mutate_all(., ~ replace(., is.na(.), 0)) %>%
    group_by(`1st_network`, Subj_ID, Age_group) %>%
    summarize_at(vars(Connector, Provincial, Satellite, Peripheral), mean) %>%
    group_by(`1st_network`) %>%
    summarise_at(vars(Connector, Provincial, Satellite, Peripheral), funs(geomMeanExtension(., epsilon = epsilon)))
  
  delta_proportion_modular <<- trajectory_modular %>% group_by(Age_group) %>% 
    mutate(Connector = log(Connector / geometric_trajectory_modular$Connector)) %>%
    mutate(Provincial = log(Provincial / geometric_trajectory_modular$Provincial)) %>%
    mutate(Satellite = log(Satellite / geometric_trajectory_modular$Satellite)) %>%
    mutate(Peripheral = log(Peripheral / geometric_trajectory_modular$Peripheral)) %>%
    pivot_longer(cols = !c("1st_network", "Age_group"), names_to = "Hub_consensus", values_to = "delta_freq")
    

  if (RSN_modular == "All") {
    Radar_RSN_modular <-
      delta_proportion_modular %>% 
      spread(`1st_network`, delta_freq) %>% 
      group_by(Hub_consensus) %>% group_split()
  } else {
    Radar_RSN_modular <-
      delta_proportion_modular %>%
      filter(grepl(RSN_modular, `1st_network`)) %>%
      spread(`1st_network`, delta_freq) %>% 
      group_by(Hub_consensus) %>% group_split()
  }

  list_plot <- list()
  for (i in 1:length(Radar_RSN_modular)) {
    tmp_plot <- rbindlist(Radar_RSN_modular[i])
    title <- tmp_plot %>% dplyr::select(Hub_consensus) %>% distinct(.)
    
    radar_plot <- tmp_plot %>% dplyr::select(-Hub_consensus) %>% 
      remove_rownames() %>%
      column_to_rownames(var = "Age_group")
    
    radarplotting_overlap(radar_plot, max, min, 1, 1,
                          alpha = .1, label_size = 1, round = round,
                          title_fill = paste0(title),
                          palette = palette
    )
    
    # legend(
    #   x = "bottomleft", legend = rownames(radar_plot), horiz = TRUE,
    #   bty = "n", pch = 20, col = palette,
    #   text.col = "black", cex = 1, pt.cex = 2
    # )
  }
  
  ##############################################################################
  # INTERAREAL 
  ##############################################################################
  
  # Proportion of each functional role within each RSN for each subject
  trajectory_interareal <<- tmp_cluster_final %>%
    group_by(`1st_network`, Region, Subj_ID, Age_group, Bridgeness) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n)) %>%
    spread(Bridgeness, freq) %>%
    dplyr::select(-n) %>%
    mutate_all(., ~ replace(., is.na(.), 0)) %>%
    group_by(`1st_network`, Subj_ID, Age_group) %>%
    summarize_at(vars(Global_Bridge, Local_Bridge, Super_Bridge, Not_a_Bridge), mean) %>% 
    group_by(`1st_network`, Age_group) %>%
    summarise_at(vars(Global_Bridge, Local_Bridge, Super_Bridge, Not_a_Bridge), funs(geomMeanExtension(., epsilon = epsilon)))
  
  geometric_trajectory_interareal <<- tmp_cluster_final %>%
    group_by(`1st_network`, Region, Subj_ID, Age_group, Bridgeness) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n)) %>%
    spread(Bridgeness, freq) %>%
    dplyr::select(-n) %>%
    mutate_all(., ~ replace(., is.na(.), 0)) %>%
    group_by(`1st_network`, Subj_ID, Age_group) %>%
    summarize_at(vars(Global_Bridge, Local_Bridge, Super_Bridge, Not_a_Bridge), mean) %>%
    group_by(`1st_network`) %>%
    summarise_at(vars(Global_Bridge, Local_Bridge, Super_Bridge, Not_a_Bridge), funs(geomMeanExtension(., epsilon = epsilon)))
  
  delta_proportion_interareal <<- trajectory_interareal %>% group_by(Age_group) %>% 
    mutate(Global_Bridge = log(Global_Bridge / geometric_trajectory_interareal$Global_Bridge)) %>%
    mutate(Local_Bridge = log(Local_Bridge / geometric_trajectory_interareal$Local_Bridge)) %>%
    mutate(Super_Bridge = log(Super_Bridge / geometric_trajectory_interareal$Super_Bridge)) %>%
    mutate(Not_a_Bridge = log(Not_a_Bridge / geometric_trajectory_interareal$Not_a_Bridge)) %>%
    pivot_longer(cols = !c("1st_network", "Age_group"), names_to = "Bridgeness", values_to = "delta_freq")
  
  
  if (RSN_interareal == "All") {
    Radar_RSN_interareal <-
      delta_proportion_interareal %>% 
      spread(`1st_network`, delta_freq) %>% 
      group_by(Bridgeness) %>% group_split()
  } else {
    Radar_RSN_interareal <-
      delta_proportion_interareal %>%
      filter(grepl(RSN_interareal, `1st_network`)) %>%
      spread(`1st_network`, delta_freq) %>% 
      group_by(Bridgeness) %>% group_split()
  }
  
  list_plot <- list()
  for (i in 1:length(Radar_RSN_interareal)) {
    tmp_plot <- rbindlist(Radar_RSN_interareal[i])
    title <- tmp_plot %>% dplyr::select(Bridgeness) %>% distinct(.)
    
    radar_plot <- tmp_plot %>% dplyr::select(-Bridgeness) %>% 
      remove_rownames() %>%
      column_to_rownames(var = "Age_group")
    
    radarplotting_overlap(radar_plot, max2, min2, 1, 1,
                          alpha = alpha, label_size = 1, round = round,
                          title_fill = paste0(title),
                          palette = palette
    )
    
    # legend(
    #   x = "bottomleft", legend = rownames(radar_plot), horiz = TRUE,
    #   bty = "n", pch = 20, col = palette,
    #   text.col = "black", cex = 1, pt.cex = 2
    # )
  }
}

  
# Choose the clusters
# Max and min value on radar plot grid
# mfrow settings
# alpha level
# RSN to be displayed -- either "All" or specify the RSN with | for separation
# Choose round = FALSE if scores have decimal points

interaction_Age_FuncRole_RSN(0.5, -0.5, 0.5, -0.5, 0.1,
                             RSN_modular = "Auditory|Language|CON|DMN|FPN|SMN|Visual_1|DAN|PMM|Visual_2",
                             RSN_interareal = "Auditory|Language|CON|DMN|FPN|SMN|Visual_1|DAN|PMM|Visual_2",
                             epsilon = 1e-1,
                             round = FALSE
)


################################################################################
# For ratio between two age group
################################################################################

# trajectory_b <<- tmp_cluster_final %>%
#   group_by(`1st_network`, Region, Subj_ID, Age_group, Bridgeness) %>%
#   summarise(n = n()) %>%
#   mutate(freq = n / sum(n)) %>%
#   spread(Bridgeness, freq) %>%
#   dplyr::select(-n) %>%
#   mutate_all(., ~ replace(., is.na(.), 0)) %>%
#   group_by(`1st_network`, Subj_ID, Age_group) %>%
#   summarize_at(vars(Global_Bridge, Local_Bridge, Super_Bridge, Not_a_Bridge), mean) %>%
#   ungroup() %>%
#   pivot_longer(cols = !c("1st_network", "Subj_ID", "Age_group"), names_to = "Bridgeness", values_to = "freq") %>%
#   group_by(`1st_network`, Age_group, Bridgeness) %>%
#   summarise_at(vars(freq), funs(geomMeanExtension(., epsilon = epsilon)))
# 
# delta_proportion_b <<- trajectory_b %>%
#   group_by(`1st_network`, Bridgeness) %>%
#   # Compute the log ratios
#   mutate(delta_freq = log(freq / lag(freq))) %>%
#   filter(!(grepl("NaN|-Inf", delta_freq))) %>%
#   dplyr::select(-c(freq, Age_group)) %>%
#   na.omit()
# 
# if (RSN_interareal == "All") {
#   Radar_functional_role_RSN_delta_b <<-
#     delta_proportion_b %>%
#     dplyr::select(`1st_network`, Bridgeness, delta_freq) %>%
#     spread(`1st_network`, delta_freq) %>%
#     remove_rownames() %>%
#     column_to_rownames(var = "Bridgeness")
# } else {
#   Radar_functional_role_RSN_delta_b <<-
#     delta_proportion_b %>%
#     dplyr::select(`1st_network`, Bridgeness, delta_freq) %>%
#     filter(grepl(RSN_interareal, `1st_network`)) %>%
#     spread(`1st_network`, delta_freq) %>%
#     remove_rownames() %>%
#     column_to_rownames(var = "Bridgeness")
# }
# 
# 
# radarplotting_overlap(Radar_functional_role_RSN_delta_b, max2, min2, 1, 1,
#                       alpha = alpha, label_size = 1, round = FALSE,
#                       title_fill = paste("Log ratios of geometric mean proportions of interareal roles"),
#                       palette = RColorBrewer::brewer.pal(8, "Dark2")
# )
# 
# legend(
#   x = "bottomleft", legend = rownames(Radar_functional_role_RSN_delta_b), horiz = TRUE,
#   bty = "n", pch = 20, col = RColorBrewer::brewer.pal(8, "Dark2"),
#   text.col = "black", cex = 1, pt.cex = 2
# )

################################################################################
################################################################################
################################################################################
# BOOTSTRAPPING for RSN modulation
################################################################################
################################################################################
################################################################################


# set.seed(123)
# bootstrap_ci <- function(n_boot, cluster1, cluster2, RSN_modular, RSN_interareal, epsilon) {
#   data_cluster_selection(cluster1, cluster2)
#   
#   data_boot <- tmp_cluster_final %>%
#     group_by(`1st_network`, Region, Subj_ID, Age_group, Hub_consensus) %>%
#     summarise(n = n()) %>%
#     mutate(freq = n / sum(n)) %>%
#     spread(Hub_consensus, freq) %>%
#     dplyr::select(-n) %>%
#     mutate_all(., ~ replace(., is.na(.), 0)) %>%
#     group_by(`1st_network`, Subj_ID, Age_group) %>%
#     summarize_at(vars(Connector:Satellite), mean) %>%
#     ungroup() %>%
#     filter(grepl(RSN_modular, `1st_network`)) %>%
#     # Draw samples with replacement for every RSN, and Age_group
#     group_by(Age_group, `1st_network`) %>%
#     group_split()
#   
#   list_boot <- list()
#   for (i in 1:length(data_boot)) {
#     tmp <- rbindlist(data_boot[i]) %>% mutate(n_data_split = rep(i, times = length(nrow(.))))
#     tmp_bis <- lapply(1:n_boot, function(i) tmp[sample(nrow(tmp), replace = TRUE), ])
#     list_boot_bis <- list()
#     for (j in 1:length(tmp_bis)) {
#       list_boot_bis[[j]] <- rbindlist(tmp_bis[j]) %>% mutate(n_data_boot = rep(j, times = length(nrow(.))))
#     }
#     delta_boot_a <- rbindlist(list_boot_bis) %>%
#       as.data.frame() %>%
#       pivot_longer(cols = !c("n_data_boot", "n_data_split", "1st_network", "Subj_ID", "Age_group"), names_to = "Hub_consensus", values_to = "freq") %>%
#       group_by(n_data_boot, n_data_split, `1st_network`, Age_group, Hub_consensus) %>%
#       summarise_at(vars(freq), funs(geomMeanExtension(., epsilon = epsilon)))
#     
#     list_boot[[i]] <- delta_boot_a
#   }
#   
#   # 1000 samples for each Age_group*1st_network combination
#   # So that for each RSN, and within each Age_group, 1000 resamples of the same size are taken
#   
#   resamples <- rbindlist(list_boot) %>% as.data.frame()
#   
#   log_ratio_boot_a <- resamples %>%
#     group_by(n_data_boot, `1st_network`, Hub_consensus) %>%
#     # Compute the log ratios
#     mutate(delta_freq = log(freq / lag(freq))) %>%
#     filter(!(grepl("NaN|-Inf", delta_freq))) %>%
#     dplyr::select(-c(freq, Age_group)) %>%
#     na.omit()
#   
#   summary_bootstrap_a <<- log_ratio_boot_a %>%
#     group_by(`1st_network`, Hub_consensus) %>%
#     get_summary_stats(delta_freq, type = "full") %>%
#     mutate(mean_025 = mean - ci) %>%
#     mutate(mean_0975 = mean + ci) %>%
#     dplyr::select(`1st_network`, Hub_consensus, mean, mean_025, mean_0975)
#   
#   write.csv(summary_bootstrap_a, "summary_bootstrap_modular.csv")
#   
#   ##############################################################################
#   ##############################################################################
#   
#   data_boot_bis <- tmp_cluster_final %>%
#     group_by(`1st_network`, Region, Subj_ID, Age_group, Bridgeness) %>%
#     summarise(n = n()) %>%
#     mutate(freq = n / sum(n)) %>%
#     spread(Bridgeness, freq) %>%
#     dplyr::select(-n) %>%
#     mutate_all(., ~ replace(., is.na(.), 0)) %>%
#     group_by(`1st_network`, Subj_ID, Age_group) %>%
#     summarize_at(vars(Global_Bridge:Super_Bridge), mean) %>%
#     ungroup() %>%
#     filter(grepl(RSN_interareal, `1st_network`)) %>%
#     # Draw samples with replacement for every RSN, and Age_group
#     group_by(Age_group, `1st_network`, .add = TRUE) %>%
#     group_split()
#   
#   list_boot <- list()
#   for (i in 1:length(data_boot_bis)) {
#     tmp <- rbindlist(data_boot_bis[i]) %>% mutate(n_data_split = rep(i, times = length(nrow(.))))
#     tmp_bis <- lapply(1:n_boot, function(i) tmp[sample(nrow(tmp), replace = TRUE), ])
#     list_boot_bis <- list()
#     for (j in 1:length(tmp_bis)) {
#       list_boot_bis[[j]] <- rbindlist(tmp_bis[j]) %>% mutate(n_data_boot = rep(j, times = length(nrow(.))))
#     }
#     delta_boot_b <<- rbindlist(list_boot_bis) %>%
#       as.data.frame() %>%
#       pivot_longer(
#         cols = !c("n_data_boot", "n_data_split", "1st_network", "Subj_ID", "Age_group"),
#         names_to = "Bridgeness",
#         values_to = "freq"
#       ) %>%
#       group_by(n_data_boot, n_data_split, `1st_network`, Age_group, Bridgeness) %>%
#       summarise_at(vars(freq), funs(geomMeanExtension(., epsilon = epsilon)))
#     
#     list_boot[[i]] <- delta_boot_b
#   }
#   
#   # 1000 samples for each Age_group*1st_network*Region combination
#   # So that for each <- <- <- <-  region for each RSN, and within each Age_group, 1000 resamples of the same size are taken
#   
#   resamples <- rbindlist(list_boot) %>% as.data.frame()
#   
#   log_ratio_boot_b <- resamples %>%
#     group_by(n_data_boot, `1st_network`, Bridgeness) %>%
#     mutate(delta_freq = log(freq / lag(freq))) %>%
#     filter(!(grepl("NaN|-Inf", delta_freq))) %>%
#     dplyr::select(-c(freq, Age_group)) %>%
#     na.omit()
#   
#   summary_bootstrap_b <<- log_ratio_boot_b %>%
#     group_by(`1st_network`, Bridgeness) %>%
#     get_summary_stats(delta_freq, type = "full") %>%
#     mutate(mean_025 = mean - ci) %>%
#     mutate(mean_0975 = mean + ci) %>%
#     dplyr::select(`1st_network`, Bridgeness, mean, mean_025, mean_0975)
#   
#   write.csv(summary_bootstrap_b, "summary_bootstrap_interareal.csv")
#   ##############################################################################
#   ##############################################################################
# }
# 
# bootstrap_ci(
#   n_boot = 1000, "Young", "Old",
#   RSN_modular = "Auditory|CON|DMN|FPN|Language|SMN|PMM",
#   RSN_interareal = "Auditory|CON|DMN|FPN|Language|SMN|PMM",
#   epsilon = 1e-1
# )

