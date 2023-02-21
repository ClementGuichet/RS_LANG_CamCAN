source("_06_TFProfile_CamCAN.R")
source("_radarplotting_function.R")
source("_geometricmeanCruz.R")

library(gbmt)
################################################################################
################################################################################
################################################################################
# Group-based multivariate trajectory analysis ----
################################################################################
################################################################################
################################################################################

list_TFP_RSN_Age_group <- TFP_RSN_Age_group %>% 
  group_by(`1st_network`) %>% group_split()

list_tmp_Age_group <- list()
for (i in 1:length(list_TFP_RSN_Age_group)) {
  library(zCompositions)
  
  tmp_raw <- rbindlist(list_TFP_RSN_Age_group[i]) %>% arrange(Age_group) 
  
  tmp_coda_modular <- tmp_raw %>%
    dplyr::select(Connector, Satellite, Provincial, Peripheral)
  
  if (min(tmp_coda_modular) == 0) {
    tmp_coda_modular_bis <- tmp_coda_modular %>% 
      acomp(.) %>% 
      cmultRepl(., output = "prop")
  } else {
    tmp_coda_modular_bis <- tmp_coda_modular
  }
    
  
  tmp_coda_internodal <- tmp_raw %>%
    dplyr::select(Global_Bridge, Local_Bridge, Super_Bridge, Not_a_Bridge)
  
  if (min(tmp_coda_internodal) == 0) {
    tmp_coda_internodal_bis <- tmp_coda_internodal %>% 
      acomp(.) %>%
      cmultRepl(., output = "prop")
  } else {
    tmp_coda_internodal_bis <- tmp_coda_internodal
  }
    
  tmp_raw_imputed <- cbind(tmp_raw %>% dplyr::select(-c(Connector, Satellite, Provincial, Peripheral,
                                                        Global_Bridge, Local_Bridge, Super_Bridge, Not_a_Bridge)),
                           tmp_coda_modular_bis, 
                           tmp_coda_internodal_bis)
  
  tmp_geometric_all <- tmp_raw_imputed %>% 
    summarize_at(vars(Connector:Not_a_Bridge), funs(geometricmean(.)))
  
  tmp_final <- tmp_raw_imputed %>% 
    mutate(Connector = log(Connector / tmp_geometric_all$Connector)) %>%
    mutate(Provincial = log(Provincial / tmp_geometric_all$Provincial)) %>%
    mutate(Satellite = log(Satellite / tmp_geometric_all$Satellite)) %>%
    mutate(Peripheral = log(Peripheral / tmp_geometric_all$Peripheral)) %>% 
    mutate(Global_Bridge = log(Global_Bridge / tmp_geometric_all$Global_Bridge)) %>%
    mutate(Local_Bridge = log(Local_Bridge / tmp_geometric_all$Local_Bridge)) %>%
    mutate(Super_Bridge = log(Super_Bridge / tmp_geometric_all$Super_Bridge)) %>%
    mutate(Not_a_Bridge = log(Not_a_Bridge / tmp_geometric_all$Not_a_Bridge)) 
  
  list_tmp_Age_group[[i]] <- tmp_final
}

GBMT_TFP_RSN_Age_group <- rbindlist(list_tmp_Age_group)

# https://journals.sagepub.com/doi/10.1177/0962280216673085
set.seed(1)
mod <- gbmt::gbmt(
  x.names = c(
    "Connector", "Satellite", "Provincial", "Peripheral",
    "Global_Bridge", "Local_Bridge", "Super_Bridge", "Not_a_Bridge"
  ),
  unit = "1st_network",
  time = "Age_group",
  scaling = 0,
  data = GBMT_TFP_RSN_Age_group %>% as.data.frame(),
  nstart = 10,
  ng = 5,
  d =2
)

mod$assign.list
mod$ic

plot(mod,
     x.names = c("Connector", "Satellite", "Provincial", "Peripheral"),
     bands = F,
     cex.legend = 1,
     equal = T
     
)

# Stability of the solution ----
tmp_loglik <- list()
stable_groups <- numeric(0)
stability <- function(start, i) {
  for (i in start:i) {
    base::print(i)
    set.seed(i)
    mod <- gbmt::gbmt(
      x.names = c(
        "Connector", "Satellite", "Provincial", "Peripheral",
        "Global_Bridge", "Local_Bridge", "Super_Bridge", "Not_a_Bridge"
      ),
      unit = "1st_network",
      time = "Age_group",
      scaling = 0,
      data = GBMT_TFP_RSN_Age_group %>% as.data.frame(),
      nstart = 1e2,
      ng = 5,
      d = 2
    )
    tmp_grouping <- mod$assign
    tmp_loglik[[i]] <- mod$logLik
    
    stable_groups <<- cbind(stable_groups, tmp_grouping) %>% as.data.frame()
  }
}
stability(1, 100)


consensus_grouping <- stable_groups %>%
  janitor::clean_names() %>% 
  mutate_at(vars(everything()), funs(as.numeric(.))) %>%
  as.data.frame()

write.csv(consensus_grouping, "consensus_GBMT_grouping_20022023_PT_100iter.csv")

consensus_grouping <- read.csv("consensus_GBMT_grouping_20022023_PT_100iter.csv") %>%
  as.data.frame() %>% remove_rownames() %>% tibble::column_to_rownames("X")


consClustering <- function(X, K) {
  df <- consensus_grouping
  X <- consensus_grouping
  n <- nrow(X)
  
  # Empty matrices
  norm_pairs <- matrix(0, nrow = n, ncol = n, dimnames = list(rownames(X), rownames(X)))
  pairs <- matrix(0, nrow = n, ncol = n, dimnames = list(rownames(X), rownames(X)))
  
  for (i in 1:n) {
    base::print(i)
    for (j in 1:n) {
      base::print(j)
      # Looking at two observations and their cluster assignments
      pair <- df[c(i, j), df[i, ] != 0 & df[j, ] != 0, drop = FALSE]
      
      # Total times the two obs were clustered similarly in the same iter
      aggr <- sum(pair[1, ] == pair[2, ]) 
      
      # Total times both obs existed in iter
      denom <- ncol(pair)
      
      # If neither of them appeared in iter together
      if (denom == 0) denom <- 1
      
      # Regular pairs just counts number of times clustered together
      pairs[i, j] <- aggr
      # Normalized pairs divides all entries by number of iterations
      norm_pairs[i, j] <- aggr / denom
    }
  }
  
  # Clustering consensus matrix
  dist <- as.dist(1 - norm_pairs)
  hier.out <- hclust(dist, method = "ward.D2")
  hier.clust <- cutree(hier.out, k = 5)
  # normalized matrix reordered
  reorder.npairs <- norm_pairs[order(hier.clust), order(hier.clust)]
  # regular, counting matrix reordered
  reorder.pairs <- pairs[order(hier.clust), order(hier.clust)]
  
  return(list(
    norm.mat = reorder.npairs,
    clusterings = hier.clust,
    count.mat = reorder.pairs
  ))
}

output <- consClustering(consensus_grouping, 5)
final_partition <- output$clusterings %>% as.data.frame()

computeDistribution <- function(cons.mat) {
  upper_t <- cons.mat[upper.tri(cons.mat, diag = FALSE)] %>% unlist()
  ggplot() +
    geom_density(mapping = aes(upper_t), fill = "steelblue") +
    theme_light() +
    labs(
      x = "Consensus Index Value", y = "Density",
      title = "Consensus Distribution"
    )
}
computeDistribution(output$norm.mat)

heatmap(output$count.mat, verbose = TRUE)

plot_cluster <- cluster::agnes(1-norm_pairs, method = "ward", metric = "euclidian")
factoextra::fviz_dend(plot_cluster,
                      cex = 0.8,
                      k = 5,
                      palette = "jco",
                      rect = T,
                      color_labels_by_k = TRUE,
                      main = "Ward Hierarchical clustering  (100 iter PC_normed, Hub selection = .2,
                      100 iter)"
)

# LDA ----
# set.seed(1)
# data_LDA <- gbmt_rsn
# sample <- sample(c(TRUE, FALSE), nrow(data_LDA), replace = TRUE, prob = c(.8, .2))
# train <- data_LDA[sample, ]
# test <- data_LDA[!sample, ]
# 
# 
# model <- lda(Age_group~., data = train)
# model
# 
# predicted <- predict(model, test)
# predicted$class
# 
# mean(predicted$class==test$Age_group)
# 
# lda_coef <- model$scaling %>% as.data.frame() %>% mutate(names = rownames(.)) %>% arrange(LD1, LD2)
# #define data to plot
# lda_plot <- cbind(train, predict(model)$x)
# 
# #create plot
# ggplot(lda_plot, aes(LD1, LD2)) +
#   geom_point(aes(color = Age_group), size = 4) +
#   theme_pubclean()

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


