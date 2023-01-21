################################################################################
# LOG-RATIO of the Topologico-functional profile at the RSN level --------------
################################################################################

source("_06_TFProfile_CamCAN.R")
source("_radarplotting_function.R")
source("_geometricmeanCruz.R")

# Ratio of the proportion of each functional role within each RSN between the two selected clusters
interaction_Age_FuncRole_RSN <- function(cluster1, cluster2, max, min, max2, min2, alpha, RSN_modular, RSN_interareal, epsilon, round) {
  data_cluster_selection(cluster1, cluster2)

  ##############################################################################
  # For modular roles
  ##############################################################################

  # Proportion of each functional role within each RSN for each subject
  trajectory_a <<- tmp_cluster_final %>%
    group_by(`1st_network`, Region, Subj_ID, Age_group, Hub_consensus) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n)) %>%
    spread(Hub_consensus, freq) %>%
    dplyr::select(-n) %>%
    mutate_all(., ~ replace(., is.na(.), 0)) %>%
    group_by(`1st_network`, Subj_ID, Age_group) %>%
    summarize_at(vars(Connector, Provincial, Satellite, Peripheral), mean) %>%
    ungroup() %>%
    pivot_longer(
      cols = !c("1st_network", "Subj_ID", "Age_group"),
      names_to = "Hub_consensus", values_to = "freq"
    ) %>%
    group_by(`1st_network`, Age_group, Hub_consensus) %>%
    summarise_at(vars(freq), funs(geomMeanExtension(., epsilon = epsilon)))

  delta_proportion_a <<- trajectory_a %>%
    group_by(`1st_network`, Hub_consensus) %>%
    # Compute the log ratios
    mutate(delta_freq = log(freq / lag(freq))) %>%
    filter(!(grepl("NaN|-Inf", delta_freq))) %>%
    dplyr::select(-c(freq, Age_group)) %>%
    na.omit()

  if (RSN_modular == "All") {
    Radar_functional_role_RSN_delta_a <<-
      delta_proportion_a %>%
      dplyr::select(`1st_network`, Hub_consensus, delta_freq) %>%
      spread(`1st_network`, delta_freq) %>%
      remove_rownames() %>%
      column_to_rownames(var = "Hub_consensus")
  } else {
    Radar_functional_role_RSN_delta_a <<-
      delta_proportion_a %>%
      dplyr::select(`1st_network`, Hub_consensus, delta_freq) %>%
      filter(grepl(RSN_modular, `1st_network`)) %>%
      spread(`1st_network`, delta_freq) %>%
      remove_rownames() %>%
      column_to_rownames(var = "Hub_consensus")
  }

  radarplotting_overlap(Radar_functional_role_RSN_delta_a, max, min, 1, 1,
    alpha = alpha, label_size = 1, round = round,
    title_fill = paste("Log ratio of geometric mean proportions of modular roles"),
    palette = RColorBrewer::brewer.pal(8, "Dark2")
  )

  legend(
    x = "bottomleft", legend = rownames(Radar_functional_role_RSN_delta_a), horiz = TRUE,
    bty = "n", pch = 20, col = RColorBrewer::brewer.pal(8, "Dark2"),
    text.col = "black", cex = 1, pt.cex = 2
  )



  ##############################################################################
  # For interareal roles
  ##############################################################################

  trajectory_b <<- tmp_cluster_final %>%
    group_by(`1st_network`, Region, Subj_ID, Age_group, Bridgeness) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n)) %>%
    spread(Bridgeness, freq) %>%
    dplyr::select(-n) %>%
    mutate_all(., ~ replace(., is.na(.), 0)) %>%
    group_by(`1st_network`, Subj_ID, Age_group) %>%
    summarize_at(vars(Global_Bridge, Local_Bridge, Super_Bridge, Not_a_Bridge), mean) %>%
    ungroup() %>%
    pivot_longer(cols = !c("1st_network", "Subj_ID", "Age_group"), names_to = "Bridgeness", values_to = "freq") %>%
    group_by(`1st_network`, Age_group, Bridgeness) %>%
    summarise_at(vars(freq), funs(geomMeanExtension(., epsilon = epsilon)))

  delta_proportion_b <<- trajectory_b %>%
    group_by(`1st_network`, Bridgeness) %>%
    # Compute the log ratios
    mutate(delta_freq = log(freq / lag(freq))) %>%
    filter(!(grepl("NaN|-Inf", delta_freq))) %>%
    dplyr::select(-c(freq, Age_group)) %>%
    na.omit()

  if (RSN_interareal == "All") {
    Radar_functional_role_RSN_delta_b <<-
      delta_proportion_b %>%
      dplyr::select(`1st_network`, Bridgeness, delta_freq) %>%
      spread(`1st_network`, delta_freq) %>%
      remove_rownames() %>%
      column_to_rownames(var = "Bridgeness")
  } else {
    Radar_functional_role_RSN_delta_b <<-
      delta_proportion_b %>%
      dplyr::select(`1st_network`, Bridgeness, delta_freq) %>%
      filter(grepl(RSN_interareal, `1st_network`)) %>%
      spread(`1st_network`, delta_freq) %>%
      remove_rownames() %>%
      column_to_rownames(var = "Bridgeness")
  }


  radarplotting_overlap(Radar_functional_role_RSN_delta_b, max2, min2, 1, 1,
    alpha = alpha, label_size = 1, round = FALSE,
    title_fill = paste("Log ratios of geometric mean proportions of interareal roles"),
    palette = RColorBrewer::brewer.pal(8, "Dark2")
  )

  legend(
    x = "bottomleft", legend = rownames(Radar_functional_role_RSN_delta_b), horiz = TRUE,
    bty = "n", pch = 20, col = RColorBrewer::brewer.pal(8, "Dark2"),
    text.col = "black", cex = 1, pt.cex = 2
  )
}


# Choose the clusters
# Max and min value on radar plot grid
# mfrow settings
# alpha level
# RSN to be displayed -- either "All" or specify the RSN with | for separation
# Choose round = FALSE if scores have decimal points

interaction_Age_FuncRole_RSN("Young", "Old", 1, -2, 0.5, -1, 0.1,
  RSN_modular = "DMN|FPN|Visual_1",
  RSN_interareal = "All",
  epsilon = 1e-1,
  round = FALSE
)

################################################################################
################################################################################
################################################################################
# Group-based multivariate trajectory analysis ----
################################################################################
################################################################################
################################################################################

# Get the TFP for the first temporal segment
interaction_Age_FuncRole_RSN("Young", "Middle", 1, -2, 0.5, -1, 0.1,
  RSN_modular = "All",
  RSN_interareal = "All",
  epsilon = 1e-1,
  round = FALSE
)

traj_YM <- trajectory_a %>%
  spread(Hub_consensus, freq) %>%
  merge(., trajectory_b %>% spread(Bridgeness, freq),
    by = c("Age_group", "1st_network")
  )

# Get the TFP for the second temporal segment
interaction_Age_FuncRole_RSN("Middle", "Old", 1, -2, 0.5, -1, 0.1,
  RSN_modular = "All",
  RSN_interareal = "All",
  epsilon = 1e-1,
  round = FALSE
)
traj_MO <- trajectory_a %>%
  spread(Hub_consensus, freq) %>%
  merge(., trajectory_b %>% spread(Bridgeness, freq) %>% filter(Age_group == "Old"),
    by = c("Age_group", "1st_network")
  )


traj_global <- rbind(traj_YM, traj_MO) %>%
  arrange(Age_group) %>%
  as.data.frame()

# Normalization procedure
geometric_all <- traj_global %>%
  group_by(`1st_network`) %>%
  summarize_at(vars(Connector:Super_Bridge), funs(geomMeanExtension(., epsilon = 1e-5)))

gbmt_rsn <- pca_clust %>%
  # Log ratio to the geometric mean of all subjects for each compositions for each RSN
  group_by(Age_group, `1st_network`) %>%
  summarize_at(vars(Connector:Super_Bridge), funs(geomMeanExtension(., epsilon = 1e-5))) %>%
  mutate(Connector = log(Connector / geometric_all$Connector)) %>%
  mutate(Provincial = log(Provincial / geometric_all$Provincial)) %>%
  mutate(Satellite = log(Satellite / geometric_all$Satellite)) %>%
  mutate(Peripheral = log(Peripheral / geometric_all$Peripheral)) %>%
  mutate(Global_Bridge = log(Global_Bridge / geometric_all$Global_Bridge)) %>%
  mutate(Local_Bridge = log(Local_Bridge / geometric_all$Local_Bridge)) %>%
  mutate(Super_Bridge = log(Super_Bridge / geometric_all$Super_Bridge)) %>%
  mutate(Not_a_Bridge = log(Not_a_Bridge / geometric_all$Not_a_Bridge)) %>%
  mutate(Age_group = ifelse(Age_group == "Young", 1,
    ifelse(Age_group == "Middle", 2,
      3
    )
  ))


# https://journals.sagepub.com/doi/10.1177/0962280216673085
library(gbmt)
set.seed(4)
mod <- gbmt::gbmt(
  x.names = c(
    "Connector", "Provincial", "Satellite", "Peripheral",
    "Global_Bridge", "Local_Bridge", "Super_Bridge", "Not_a_Bridge"
  ),
  unit = "1st_network",
  time = "Age_group",
  scaling = 0,
  data = gbmt_rsn %>% as.data.frame(),
  nstart = 10,
  ng = 5,
  d = 2
)

# Restart 10 - EM iteration 1. Log likelihood: 1697.9793 
# > mod$assign.list
# $`1`
# [1] "DAN"      "SMN"      "Visual_1"
# 
# $`2`
# [1] "Auditory" "FPN"     
# 
# $`3`
# [1] "DMN"      "Language"
# 
# $`4`
# [1] "PMM" "VMM"
# 
# $`5`
# [1] "CON"      "Visual_2"
# 
# > mod$ic
# aic       bic      caic     ssbic      hqic 
# -1958.443 -1582.820 -1331.820 -2365.743 -1832.058 

mod$assign.list
mod$ic


plot(mod,
  x.names = c("Connector", "Provincial", "Satellite", "Peripheral"),
  bands = F,
  cex.legend = 1,
  equal = T
)
plot(mod,
  x.names = c("Global_Bridge", "Local_Bridge", "Super_Bridge", "Not_a_Bridge"),
  bands = F,
  cex.legend = 1,
  equal = T
)


mod$reg$`1`$Not_a_Bridge %>% summary(.)
mod$beta

################################################################################
################################################################################
################################################################################
# BOOTSTRAPPING for RSN modulation
################################################################################
################################################################################
################################################################################


set.seed(123)
bootstrap_ci <- function(n_boot, cluster1, cluster2, RSN_modular, RSN_interareal, epsilon) {
  data_cluster_selection(cluster1, cluster2)
  
  data_boot <- tmp_cluster_final %>%
    group_by(`1st_network`, Region, Subj_ID, Age_group, Hub_consensus) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n)) %>%
    spread(Hub_consensus, freq) %>%
    dplyr::select(-n) %>%
    mutate_all(., ~ replace(., is.na(.), 0)) %>%
    group_by(`1st_network`, Subj_ID, Age_group) %>%
    summarize_at(vars(Connector:Satellite), mean) %>%
    ungroup() %>%
    filter(grepl(RSN_modular, `1st_network`)) %>%
    # Draw samples with replacement for every RSN, and Age_group
    group_by(Age_group, `1st_network`) %>%
    group_split()
  
  list_boot <- list()
  for (i in 1:length(data_boot)) {
    tmp <- rbindlist(data_boot[i]) %>% mutate(n_data_split = rep(i, times = length(nrow(.))))
    tmp_bis <- lapply(1:n_boot, function(i) tmp[sample(nrow(tmp), replace = TRUE), ])
    list_boot_bis <- list()
    for (j in 1:length(tmp_bis)) {
      list_boot_bis[[j]] <- rbindlist(tmp_bis[j]) %>% mutate(n_data_boot = rep(j, times = length(nrow(.))))
    }
    delta_boot_a <- rbindlist(list_boot_bis) %>%
      as.data.frame() %>%
      pivot_longer(cols = !c("n_data_boot", "n_data_split", "1st_network", "Subj_ID", "Age_group"), names_to = "Hub_consensus", values_to = "freq") %>%
      group_by(n_data_boot, n_data_split, `1st_network`, Age_group, Hub_consensus) %>%
      summarise_at(vars(freq), funs(geomMeanExtension(., epsilon = epsilon)))
    
    list_boot[[i]] <- delta_boot_a
  }
  
  # 1000 samples for each Age_group*1st_network combination
  # So that for each RSN, and within each Age_group, 1000 resamples of the same size are taken
  
  resamples <- rbindlist(list_boot) %>% as.data.frame()
  
  log_ratio_boot_a <- resamples %>%
    group_by(n_data_boot, `1st_network`, Hub_consensus) %>%
    # Compute the log ratios
    mutate(delta_freq = log(freq / lag(freq))) %>%
    filter(!(grepl("NaN|-Inf", delta_freq))) %>%
    dplyr::select(-c(freq, Age_group)) %>%
    na.omit()
  
  summary_bootstrap_a <<- log_ratio_boot_a %>%
    group_by(`1st_network`, Hub_consensus) %>%
    get_summary_stats(delta_freq, type = "full") %>%
    mutate(mean_025 = mean - ci) %>%
    mutate(mean_0975 = mean + ci) %>%
    dplyr::select(`1st_network`, Hub_consensus, mean, mean_025, mean_0975)
  
  write.csv(summary_bootstrap_a, "summary_bootstrap_modular.csv")
  
  ##############################################################################
  ##############################################################################
  
  data_boot_bis <- tmp_cluster_final %>%
    group_by(`1st_network`, Region, Subj_ID, Age_group, Bridgeness) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n)) %>%
    spread(Bridgeness, freq) %>%
    dplyr::select(-n) %>%
    mutate_all(., ~ replace(., is.na(.), 0)) %>%
    group_by(`1st_network`, Subj_ID, Age_group) %>%
    summarize_at(vars(Global_Bridge:Super_Bridge), mean) %>%
    ungroup() %>%
    filter(grepl(RSN_interareal, `1st_network`)) %>%
    # Draw samples with replacement for every RSN, and Age_group
    group_by(Age_group, `1st_network`, .add = TRUE) %>%
    group_split()
  
  list_boot <- list()
  for (i in 1:length(data_boot_bis)) {
    tmp <- rbindlist(data_boot_bis[i]) %>% mutate(n_data_split = rep(i, times = length(nrow(.))))
    tmp_bis <- lapply(1:n_boot, function(i) tmp[sample(nrow(tmp), replace = TRUE), ])
    list_boot_bis <- list()
    for (j in 1:length(tmp_bis)) {
      list_boot_bis[[j]] <- rbindlist(tmp_bis[j]) %>% mutate(n_data_boot = rep(j, times = length(nrow(.))))
    }
    delta_boot_b <<- rbindlist(list_boot_bis) %>%
      as.data.frame() %>%
      pivot_longer(
        cols = !c("n_data_boot", "n_data_split", "1st_network", "Subj_ID", "Age_group"),
        names_to = "Bridgeness",
        values_to = "freq"
      ) %>%
      group_by(n_data_boot, n_data_split, `1st_network`, Age_group, Bridgeness) %>%
      summarise_at(vars(freq), funs(geomMeanExtension(., epsilon = epsilon)))
    
    list_boot[[i]] <- delta_boot_b
  }
  
  # 1000 samples for each Age_group*1st_network*Region combination
  # So that for each <- <- <- <-  region for each RSN, and within each Age_group, 1000 resamples of the same size are taken
  
  resamples <- rbindlist(list_boot) %>% as.data.frame()
  
  log_ratio_boot_b <- resamples %>%
    group_by(n_data_boot, `1st_network`, Bridgeness) %>%
    mutate(delta_freq = log(freq / lag(freq))) %>%
    filter(!(grepl("NaN|-Inf", delta_freq))) %>%
    dplyr::select(-c(freq, Age_group)) %>%
    na.omit()
  
  summary_bootstrap_b <<- log_ratio_boot_b %>%
    group_by(`1st_network`, Bridgeness) %>%
    get_summary_stats(delta_freq, type = "full") %>%
    mutate(mean_025 = mean - ci) %>%
    mutate(mean_0975 = mean + ci) %>%
    dplyr::select(`1st_network`, Bridgeness, mean, mean_025, mean_0975)
  
  write.csv(summary_bootstrap_b, "summary_bootstrap_interareal.csv")
  ##############################################################################
  ##############################################################################
}

bootstrap_ci(
  n_boot = 1000, "Young", "Old",
  RSN_modular = "Auditory|CON|DMN|FPN|Language|SMN|PMM",
  RSN_interareal = "Auditory|CON|DMN|FPN|Language|SMN|PMM",
  epsilon = 1e-1
)

################################################################################
################################################################################
################################################################################

# Trying with GBMT at the subject level but this is unfeasible considering cross-sectional design
# 
# data_coda_modular <- TFP_RSN %>%
#   dplyr::select(Connector, Satellite, Provincial, Peripheral) %>%
#   acomp(.) %>%
#   cmultRepl(., output = "prop")
#
# data_coda_interareal <- TFP_RSN %>%
#   dplyr::select(Global_Bridge, Local_Bridge, Super_Bridge, Not_a_Bridge) %>%
#   acomp(.) %>%
#   # Bayesian non-parametric multiplicative replacement that preserves the ratios between non-zero components
#   cmultRepl(., output = "prop")
#
# gbmt_subj <- cbind(data_coda_modular, data_coda_interareal, TFP_RSN %>%
#                      dplyr::select(c("Subj_ID", "Age", "1st_network")))
#
# geometric_all <- gbmt_subj %>%
#   group_by(`1st_network`) %>%
#   summarize_at(vars(Connector:Not_a_Bridge), funs(geomMeanExtension(., epsilon = 1e-1)))
#
# gbmt_subj <- gbmt_subj %>%
#   # Log ratio to the geometric mean of all subjects for each compositions for each RSN
#   group_by(Age, `1st_network`) %>%
#   summarize_at(vars(Connector:Not_a_Bridge), funs(geomMeanExtension(., epsilon = 1e-1))) %>%
#   mutate(Connector = log(Connector / geometric_all$Connector)) %>%
#   mutate(Provincial = log(Provincial / geometric_all$Provincial)) %>%
#   mutate(Satellite = log(Satellite / geometric_all$Satellite)) %>%
#   mutate(Peripheral = log(Peripheral / geometric_all$Peripheral)) %>%
#   mutate(Global_Bridge = log(Global_Bridge / geometric_all$Global_Bridge)) %>%
#   mutate(Local_Bridge = log(Local_Bridge / geometric_all$Local_Bridge)) %>%
#   mutate(Super_Bridge = log(Super_Bridge / geometric_all$Super_Bridge)) %>%
#   mutate(Not_a_Bridge = log(Not_a_Bridge / geometric_all$Not_a_Bridge))
#
# library(gbmt)
# set.seed(4)
# mod <- gbmt::gbmt(
#   x.names = c("Connector", "Provincial", "Peripheral", "Satellite",
#               "Global_Bridge", "Local_Bridge", "Super_Bridge", "Not_a_Bridge"),
#   unit = "1st_network",
#   time = "Age",
#   scaling = 0,
#   data = gbmt_subj %>% as.data.frame(),
#   nstart = 10,
#   ng = 6,
#   d = 4
# )
#
# mod$assign.list
# mod$ic
#
#
# plot(mod, x.names = c("Connector", "Provincial", "Peripheral", "Satellite"),
#      bands = F,
#      cex.legend = 1,
#      equal = F)
# plot(mod, x.names = c("Global_Bridge", "Local_Bridge", "Super_Bridge", "Not_a_Bridge"),
#      bands = F,
#      cex.legend = 1,
#      equal = F)

