################################################################################
# Written by CG
# 23-11-2022

# What is the proportion of each CAB-NP RSN in each of the 4 LANG Nets ?
# What is the proportion of each CAB-NP RSN in the community structure yielded by rs-data ?
# How does the community structure reconfigures between states ?

# I compute the proportion using the results yielded by our network overlap calculator from the ImCalc method
# The certainty factor corresponds to the percentage of overlap between a region and the highest corresponding network
################################################################################
library(ggalluvial)
##########################################################################################
# Import processed data------------------------------------------------------------------
source("_01_DataManip_CamCAN.R")
source("_NMI&AMI_functions.R")
source("_radarplotting_function.R")

# Define palette for visualization
custom_palette <- c(
  "Auditory" = "#FF99FF",
  "Language" = "#FF6600",
  "CON" = "#9900CC",
  "DMN" = "#FF0033",
  "FPN" = "#FFCC33",
  "SMN" = "#0099CC",
  "DAN" = "#33FF66",
  "Visual_1" = "#CCCCCC",
  "Visual_2" = "#CCCCCC",
  "PMM" = "#CC0033",
  "VMM" = "#CC0033",
  "NaN" = "white",
  "Multi" = "white",
  "Multi/SM" = "white"
)

################################################################################
# State reconfiguration between task & rs-fMRI across the 131 LANG ROIs
################################################################################

# Keeping a 131 row-dataframe
data_131 <- data_full_per_region 
# %>% subset(threshold == "0.15")

# NMI & AMI
NMI_func(factor(data_131$LANG_Net_assign), factor(data_131$Consensus_OMST_all))
AMI_func(factor(data_131$LANG_Net_assign), factor(data_131$Consensus_OMST_all))

# Contingency table
addmargins(table(data_131$Consensus_OMST_all, data_131$LANG_Net_assign))

# Create alluvial diagram between the two community structures

data_alluvial_community <- data_131 %>%
  dplyr::select(LANG_Net_assign, Consensus_OMST_all, Region) %>%
  mutate(LANG_Net_assign = ifelse(LANG_Net_assign == "1", "Encoding-Decoding",
    ifelse(LANG_Net_assign == "2", "Control-Executive",
      ifelse(LANG_Net_assign == "3", "Abstract-Conceptual",
        "Sensorimotor"
      )
    )
  )) %>%
  mutate(Consensus_OMST_all = ifelse(Consensus_OMST_all == "1", "RS-NET 1",
    ifelse(Consensus_OMST_all == "2", "RS-NET 2",
      ifelse(Consensus_OMST_all == "3", "RS-NET 3",
        ifelse(Consensus_OMST_all == "4", "RS-NET 4",
          "RS-NET 5"
        )
      )
    )
  )) %>%
  plyr::rename(c("Consensus_OMST_all" = "RS-Nets")) %>%
  plyr::rename(c("LANG_Net_assign" = "LANG Nets")) %>%
  pivot_longer(
    cols = c("LANG Nets", "RS-Nets"),
    names_to = "Community_structure",
    values_to = "Communities"
  )

display_percentage <- data_alluvial_community %>%
  group_by(Community_structure, Communities) %>%
  summarize(s = n()) %>%
  group_by(Community_structure) %>%
  mutate(s = scales::percent(s / sum(s), accuracy = 0.1)) %>%
  arrange(Community_structure, desc(Communities)) %>%
  .$s

display_communities <- data_alluvial_community %>%
  group_by(Community_structure, Communities) %>%
  summarize(s = n()) %>%
  arrange(Community_structure, desc(Communities)) %>%
  .$Communities


alluvial_community <- ggplot(
  data_alluvial_community,
  aes(x = Community_structure, stratum = Communities, alluvium = Region, fill = Communities)
) +
  geom_flow(alpha = .7, curve_type = "arctangent", width = .2, na.rm = TRUE) +
  geom_stratum(alpha = .8) +
  scale_x_discrete(expand = c(.1, .1)) +
  # geom_text(stat = "stratum", label = display_percentage, nudge_x = -0.05) +
  geom_text(stat = "stratum", label = display_communities) +
  scale_fill_brewer(palette = "Oranges", direction = 1) +
  labs(title = "Reconfiguration of community structure across the 131 LANG ROIs between task-fMRI & rs-fMRI") +
  theme_pubclean()

alluvial_community

