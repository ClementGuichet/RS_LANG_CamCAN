source("_06_TFProfile_CamCAN.R")

# ******************************************************************************
# Most likely trajectory
# ******************************************************************************

# Pick the appropriate dataframe - modular or interareal
# Pick the desired RSN - to be specified in a grepl format e.g., "DMN|FPN|Language"
# Select the z-scored probability threshold for trajectories. 1.5 seems to work quite well

trajectory <- function(composition, list_RSN, threshold) {
  # Get the associated Resting-state networks for all age groups
  # Hub region specific to each subject yielded by hub detection procedure
  data_hub_selection_per_subject <- rbindlist(Hub_selection)
  tmp_cluster_final <<- merge(data_hub_selection_per_subject, data_TFP_analysis %>%
    dplyr::select(Subj_ID, Age_group),
  by = "Subj_ID"
  )
  ##############################################################################
  
  if (composition == "modular") {
    modular <<- tmp_cluster_final %>%
      group_by(`1st_network`, Region, Age_group, Subj_ID, Hub_consensus) %>%
      summarise(n = n()) %>%
      mutate(freq = n / sum(n)) %>%
      dplyr::select(-n) %>%
      spread(Hub_consensus, freq) %>%
      mutate_all(., ~ replace(., is.na(.), 0)) %>%
      group_by(`1st_network`, Region, Age_group) %>%
      summarize_at(vars(Connector, Provincial, Satellite, Peripheral), mean) %>%
      ungroup() %>%
      arrange(Region) %>%
      pivot_longer(
        cols = !c("1st_network", "Region", "Age_group"),
        names_to = "Hub_consensus", values_to = "freq"
      )

    modular$Age_group <- factor(modular$Age_group, levels = c("Young", "Middle", "Old"))

    if (list_RSN != "All") {
      Cond_PMF <- modular %>%
        filter(grepl(list_RSN, `1st_network`)) %>%
        spread(Age_group, freq) %>%
        arrange(Hub_consensus) %>%
        group_by(Region) %>%
        group_split()
    } else {
      Cond_PMF <- modular %>%
        spread(Age_group, freq) %>%
        arrange(Hub_consensus) %>%
        group_by(Region) %>%
        group_split()
    }
  } else {
    interareal <<- tmp_cluster_final %>%
      group_by(`1st_network`, Region, Age_group, Subj_ID, Bridgeness) %>%
      summarise(n = n()) %>%
      mutate(freq = n / sum(n)) %>%
      dplyr::select(-n) %>%
      spread(Bridgeness, freq) %>%
      mutate_all(., ~ replace(., is.na(.), 0)) %>%
      group_by(`1st_network`, Region, Age_group) %>%
      summarize_at(vars(Global_Bridge, Local_Bridge, Super_Bridge, Not_a_Bridge), mean) %>%
      ungroup() %>%
      arrange(Region) %>%
      pivot_longer(
        cols = !c("1st_network", "Region", "Age_group"),
        names_to = "Bridgeness", values_to = "freq"
      )

    interareal$Age_group <- factor(interareal$Age_group, levels = c("Young", "Middle", "Old"))

    if (list_RSN != "All") {
      Cond_PMF <- interareal %>%
        filter(grepl(list_RSN, `1st_network`)) %>%
        spread(Age_group, freq) %>%
        arrange(Bridgeness) %>%
        group_by(Region) %>%
        group_split()
    } else {
      Cond_PMF <- interareal %>%
        spread(Age_group, freq) %>%
        arrange(Bridgeness) %>%
        group_by(Region) %>%
        group_split()
    }
  }

  ##############################################################################
  ##############################################################################
  # Convert an adjacency dataframe to a 2-column dataframe
  adjacency_to_2col <- function(data) {
    crossdata <- lapply(rownames(data), function(x) sapply(colnames(data), function(y) list(x, y, data[x, y])))
    crossdatatmp <- matrix(unlist(crossdata), nrow = 3)
    crossdatamat <- t(crossdatatmp)
    if (colnames(data)[5] == "YM") {
      colnames(crossdatamat) <- c("Young", "Middle", "Value")
    } else {
      colnames(crossdatamat) <- c("Middle", "Old", "Value")
    }
    crossdatadf <- as.data.frame(crossdatamat, stringsAsFactors = F)
    crossdatadf[, 3] <- as.numeric(crossdatadf[, 3] %>% na.omit())
    return(crossdatadf %>% na.omit())
  }
  ##############################################################################
  ##############################################################################

  # For each region, find the most probable trajectory
  Cond_PMF_list <- list()
  for (i in 1:length(Cond_PMF)) {
    tmp <<- rbindlist(Cond_PMF[i])
    
    ############################################################################
    # First segment - Young to Middle
    ############################################################################

    # Beware of the order here, first young then old, make sure of the factor levels
    # Outer compute all pairwise probabilities between each of the 4 functional roles
    # Returning a 4 by 4 matrix
    
    outer_young_to_middle <- outer(tmp[, 4] %>% as.matrix(), tmp[, 5] %>% as.matrix()) %>% as.data.frame()
    if (composition == "modular") {
      rownames(outer_young_to_middle) <- unlist(tmp$Hub_consensus)
      colnames(outer_young_to_middle) <- unlist(tmp$Hub_consensus)
    } else {
      rownames(outer_young_to_middle) <- unlist(tmp$Bridgeness)
      colnames(outer_young_to_middle) <- unlist(tmp$Bridgeness)
    }

    outer_young_to_middle <- outer_young_to_middle %>% mutate(YM = "YM")
    crossdatadf <- adjacency_to_2col(outer_young_to_middle)
    outer_young_to_middle_2 <- crossdatadf

    tmp_young_to_middle <- cbind(tmp %>% slice(1:nrow(outer_young_to_middle_2)) %>%
      dplyr::select(`1st_network`, Region), outer_young_to_middle_2) %>%
      plyr::rename(c("Value" = "Value_YM"))

    ############################################################################
    # Second segment - Middle to Old
    ############################################################################
    
    # Returns a matrix of size Most likely YM by 4 functional role for Old
    outer_middle_to_old <- outer(tmp_young_to_middle[, 5] %>% as.matrix(), tmp[, 6] %>% as.matrix()) %>% as.data.frame()
    if (composition == "modular") {
      # Create a vector to bypass duplicate rows not allowed
      # This essentially keeps trajectories if they are tied
      # e.g., Provincial-Connector & Provincial-Provincial
      
      bypass_duplicate_row <- tmp_young_to_middle %>%
        mutate(helper_vector = rep(seq(1:nrow(.)))) %>%
        unite(., Young, "Young", "helper_vector", remove = FALSE) %>%
        unite(., Middle, "Middle", "helper_vector", remove = FALSE) %>%
        dplyr::select(-helper_vector)

      rownames(outer_middle_to_old) <- unlist(bypass_duplicate_row$Middle)
      colnames(outer_middle_to_old) <- unlist(tmp$Hub_consensus)
    } else {
      bypass_duplicate_row <- tmp_young_to_middle %>%
        mutate(helper_vector = rep(seq(1:nrow(.)))) %>%
        unite(., Young, "Young", "helper_vector", remove = FALSE) %>%
        unite(., Middle, "Middle", "helper_vector", remove = FALSE) %>%
        dplyr::select(-helper_vector)

      rownames(outer_middle_to_old) <- unlist(bypass_duplicate_row$Middle)
      colnames(outer_middle_to_old) <- unlist(tmp$Bridgeness)
    }

    outer_middle_to_old <- outer_middle_to_old %>% mutate(MO = "MO")
    # Normalize the trajectory probabilities for each region
    crossdatadf <- adjacency_to_2col(outer_middle_to_old) %>% mutate(Value_norm = as.numeric(scale(Value)))

    ############################################################################
    # Probabilistic Threshold selection
    ############################################################################
    max <- crossdatadf %>% dplyr::slice_max(Value_norm, n = 1)
    if (max$Value_norm[1] >= threshold) {
      outer_middle_to_old_2 <- crossdatadf %>% filter(Value_norm >= threshold)
    } else {
      outer_middle_to_old_2 <- crossdatadf %>% dplyr::slice_max(Value_norm, n = 1)
    }

    tmp_middle_to_old <- cbind(tmp %>% slice(1:nrow(outer_middle_to_old_2)) %>%
      dplyr::select(`1st_network`, Region), outer_middle_to_old_2) %>%
      plyr::rename(c("Value" = "Value_MO"))

    full_trajectory <- tmp_middle_to_old %>%
      merge(., bypass_duplicate_row %>% dplyr::select(c("Young", "Middle", "Value_YM")), by = c("Middle"))

    Cond_PMF_list[[i]] <- full_trajectory
  }

  if (composition == "modular") {
    Cond_PMF_final <<- rbindlist(Cond_PMF_list) %>%
      # Giving the original labels back
      # They had an index number from the helper vector
      mutate(Young = ifelse(base::startsWith(Young, "Connector"), "Connector",
        ifelse(base::startsWith(Young, "Provincial"), "Provincial",
          ifelse(base::startsWith(Young, "Peripheral"), "Peripheral",
            ifelse(base::startsWith(Young, "Satellite"), "Satellite", Young)
          )
        )
      )) %>%
      mutate(Middle = ifelse(base::startsWith(Middle, "Connector"), "Connector",
        ifelse(base::startsWith(Middle, "Provincial"), "Provincial",
          ifelse(base::startsWith(Middle, "Peripheral"), "Peripheral",
            ifelse(base::startsWith(Middle, "Satellite"), "Satellite", Middle)
          )
        )
      )) %>%
      # Identify each region with a unique label
      mutate(helper_vector = rep(seq(nrow(.)))) %>%
      unite(., Region, c("Region", "helper_vector"), remove = FALSE) %>%
      # Transform  back to an alluvial format
      pivot_longer(
        cols = c("Young", "Middle", "Old"),
        names_to = "Age_group",
        values_to = "Hub_consensus"
      )

    Cond_PMF_final$Age_group <- factor(Cond_PMF_final$Age_group, levels = c("Young", "Middle", "Old"))
  } else {
    Cond_PMF_final <<- rbindlist(Cond_PMF_list) %>%
      mutate(Young = ifelse(base::startsWith(Young, "Global_Bridge"), "Global_Bridge",
        ifelse(base::startsWith(Young, "Local_Bridge"), "Local_Bridge",
          ifelse(base::startsWith(Young, "Super_Bridge"), "Super_Bridge",
            ifelse(base::startsWith(Young, "Not_a_Bridge"), "Not_a_Bridge", Young)
          )
        )
      )) %>%
      mutate(Middle = ifelse(base::startsWith(Middle, "Global_Bridge"), "Global_Bridge",
        ifelse(base::startsWith(Middle, "Local_Bridge"), "Local_Bridge",
          ifelse(base::startsWith(Middle, "Super_Bridge"), "Super_Bridge",
            ifelse(base::startsWith(Middle, "Not_a_Bridge"), "Not_a_Bridge", Middle)
          )
        )
      )) %>%
      # Identify each region with a unique label
      # This is to deal with tied trajectories
      # Essentially, we keep the same label for the functional role but in the dataframe,
      # tied trajectories are denoted by different name for the region so that it can be used
      # for alluvial visualization 
      
      mutate(helper_vector = rep(seq(nrow(.)))) %>%
      unite(., Region, c("Region", "helper_vector"), remove = FALSE) %>%
      
      # Transform into a long alluvial format
      pivot_longer(
        cols = c("Young", "Middle", "Old"),
        names_to = "Age_group",
        values_to = "Bridgeness"
      )

    Cond_PMF_final$Age_group <- factor(Cond_PMF_final$Age_group, levels = c("Young", "Middle", "Old"))
  }

  ############################################################################
  # ALLUVIAL PLOT
  ############################################################################
  library(ggalluvial)
  if (composition == "modular") {
    display_cluster <- Cond_PMF_final %>%
      group_by(Age_group, Hub_consensus) %>%
      summarize(s = n()) %>%
      arrange(Age_group, desc(Hub_consensus)) %>%
      .$Hub_consensus

    alluvial_cluster <- ggplot(
      Cond_PMF_final,
      aes(x = Age_group, stratum = Hub_consensus, alluvium = Region, fill = Hub_consensus)
    ) +
      geom_flow(alpha = .7, curve_type = "arctangent", width = .2, na.rm = TRUE) +
      geom_stratum(alpha = .8) +
      scale_x_discrete(expand = c(.1, .1)) +
      # geom_text(stat = "stratum", label = display_percentage, nudge_x = -0.05) +
      geom_text(stat = "stratum", label = display_cluster) +
      scale_fill_brewer(palette = "Oranges", direction = 1) +
      labs(title = paste0("Most probable topological reconfiguration trajectory (z-scored probability threshold: ", threshold, ")")) +
      labs(x = "Age group") +
      theme_pubclean()

    alluvial_cluster
  } else {
    display_cluster <- Cond_PMF_final %>%
      group_by(Age_group, Bridgeness) %>%
      summarize(s = n()) %>%
      arrange(Age_group, desc(Bridgeness)) %>%
      .$Bridgeness

    alluvial_cluster <- ggplot(
      Cond_PMF_final,
      aes(x = Age_group, stratum = Bridgeness, alluvium = Region, fill = Bridgeness)
    ) +
      geom_flow(alpha = .7, curve_type = "arctangent", width = .2, na.rm = TRUE) +
      geom_stratum(alpha = .8) +
      scale_x_discrete(expand = c(.1, .1)) +
      # geom_text(stat = "stratum", label = display_percentage, nudge_x = -0.05) +
      geom_text(stat = "stratum", label = display_cluster) +
      scale_fill_brewer(palette = "Oranges", direction = 1) +
      labs(title = paste0("Most probable topological reconfiguration trajectory (z-scored probability threshold: ", threshold, ")")) +
      labs(x = "Age group") +
      theme_pubclean()

    alluvial_cluster
  }
}

trajectory("modular", "All", 2.5)

############################################################################
# Easter egg :)
#############################################################################

# seq(from=-10, to=10, by = 0.05) %>%
#   expand.grid(x=., y=.) %>%
#   ggplot(aes(x=(x+pi*sin(y)), y=(y+pi*sin(x)))) +
#   geom_point(alpha=.1, shape=20, size=1, color="black")+
#   theme_void()
