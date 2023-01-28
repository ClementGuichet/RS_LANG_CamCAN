source("_06_TFProfile_CamCAN.R")
library(introdataviz)

##########################################################################################
# Bootstrapping procedure to compute confidence intervals of log-ratios
##########################################################################################

# The geometric mean of each behavior (%) in all clusters is computed
# The log-ratio of geometric means relative to that of the group is computed
# First, 1000 virtual data sets are drawn with replacement from the source population and of the same size.
# That, for each of Age group, we take 1000 resample with replacement
# For each resample, the log-ratio of the geometric mean is computed
# The resulting distribution of 1000 log-ratios are averaged to calculate bootstrapped mean,
# and the 2.5th and 97.5th percentiles are selected as upper and lower limits of 95% confidence intervals of the bootstrapped mean.

# If the confidence interval contains the value ‘0’, no difference between the two groups for this particular cluster/RSN combination is concluded.
# Only combinations for which the intervals are outside 0 are considered responsible for the cluster differences.

set.seed(123)

bootstrap_ci <- function(n_boot, epsilon) {
  data_boot <- data_TFP_analysis %>% 
    # Draw samples with replacement for every Age_group
    group_by(Age_group) %>%
    group_split()
  
  list_boot <- list()
  for (i in 1:length(data_boot)) {
    # Number the group being sampled
    tmp <- rbindlist(data_boot[i]) %>% mutate(n_data_split = rep(i, times = length(nrow(.))))
    tmp_bis <- lapply(1:n_boot, function(i) tmp[sample(nrow(tmp), replace = TRUE), ])
    list_boot_bis <- list()
    for (j in 1:length(tmp_bis)) {
      # Number the boot sample
      list_boot_bis[[j]] <- rbindlist(tmp_bis[j]) %>% mutate(n_data_boot = rep(j, times = length(nrow(.))))
    }
    delta_boot_a <<- rbindlist(list_boot_bis) %>%
      as.data.frame()
    list_boot[[i]] <- delta_boot_a
  }
  resamples <<- rbindlist(list_boot) %>% as.data.frame()
  
  # Overall geom mean
  geometric_all_boot <<- resamples %>%
    group_by(n_data_boot) %>% 
    summarize_at(vars(Connector:Not_a_Bridge), funs(geomMeanExtension(., epsilon = epsilon)))
  # Geom Mean per Age group
  geometric_group_boot <<- resamples %>% group_by(n_data_boot, Age_group) %>% 
    summarize_at(vars(Connector:Not_a_Bridge), funs(geomMeanExtension(., epsilon = epsilon)))
  
  # Log ratios
  final_results <<- geometric_group_boot %>% group_by(Age_group) %>% 
    mutate(Connector = log(Connector / geometric_all_boot$Connector)) %>% 
    mutate(Provincial = log(Provincial / geometric_all_boot$Provincial)) %>%
    mutate(Satellite = log(Satellite / geometric_all_boot$Satellite)) %>%
    mutate(Peripheral = log(Peripheral / geometric_all_boot$Peripheral)) %>%
    mutate(Global_Bridge = log(Global_Bridge / geometric_all_boot$Global_Bridge)) %>%
    mutate(Local_Bridge = log(Local_Bridge / geometric_all_boot$Local_Bridge)) %>%
    mutate(Super_Bridge = log(Super_Bridge / geometric_all_boot$Super_Bridge)) %>%
    mutate(Not_a_Bridge = log(Not_a_Bridge /  geometric_all_boot$Not_a_Bridge))
  
  summary_bootstrap <<- final_results %>%
    pivot_longer(
      cols = !c("Age_group", "n_data_boot"),
      names_to = "Functional_role",
      values_to = "Log_ratio_geom_mean"
    )
  
  write.csv(summary_bootstrap, "summary_bootstrap.csv")
}

bootstrap_ci(1e4, 1e-1)

##########################################################################################
# Bootstrap results visualization
##########################################################################################

plot_ci <- function(functional_role, rain_height) {
  if (functional_role == "modular") {
    selection <- c("Connector", "Provincial", "Satellite", "Peripheral")
  } else if (functional_role == "interareal") {
    selection <- c("Global_Bridge", "Local_Bridge", "Super_Bridge", "Not_a_Bridge")
  }
  rain_height <- rain_height
  summary_bootstrap <- read.csv("summary_bootstrap.csv")
  summary_bootstrap$Age_group <- factor(summary_bootstrap$Age_group, levels = c("Young", "Middle", "Old"))
  
  ggplot(summary_bootstrap %>% filter(Functional_role == selection), aes(x = "", y = Log_ratio_geom_mean, fill = Age_group)) +
    # clouds
    introdataviz::geom_flat_violin(trim=FALSE, alpha = 0.4,
                                   position = position_nudge(x = rain_height+.05)) +
    # rain
    geom_point(aes(colour = Age_group), size = 2, alpha = .5, show.legend = FALSE, 
               position = position_jitter(width = rain_height, height = 0)) +
    # boxplots
    geom_boxplot(width = rain_height, alpha = 0.4, show.legend = FALSE, 
                 outlier.shape = NA,
                 position = position_nudge(x = -rain_height*2)) +
    # mean and SE point in the cloud
    stat_summary(fun.data = mean_cl_normal, mapping = aes(color = Age_group), show.legend = FALSE,
                 position = position_nudge(x = rain_height * 3)) +
    # adjust layout
    scale_x_discrete(name = "", expand = c(rain_height*3, 0, 0, 0.7)) +
    scale_y_continuous(name = "Log ratios to the geometric mean",
                       breaks = seq(-0.2, 0.2, 0.1), 
                       limits = c(-0.2, 0.2)) +
    geom_hline(yintercept = 0, color = "red", linewidth = 1) +
    coord_flip() +
    facet_wrap(~Functional_role,
               nrow = 3) +
    # custom colours and theme
    scale_fill_brewer(palette = "YlOrRd", name = "Age group") +
    scale_colour_brewer(palette = "YlOrRd") +
    theme_minimal() +
    theme(panel.grid.major.y = element_blank(),
          legend.position = c(0.5, 0.5),
          legend.background = element_rect(fill = "white", color = "white"))
}

plot_ci("modular", 0.075)
plot_ci("interareal", 0.075)
