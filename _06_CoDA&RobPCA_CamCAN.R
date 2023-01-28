# Justification for using the graph-metrics
# CoDA & Robust PCA with DetMCD estimator

# Written by CG
# 26-11-2022
##########################################################################################
pacman::p_load("zCompositions", "compositions", "cluster", "robCompositions", "vegan", "mvoutlier",
               "FactoMineR", "factoextra")

source("_05_Hub_detection_CamCAN.R")

###############################################################################
# Non-parametric Bayesian multiplicative replacement
###############################################################################

data_coda_modular <- TFP_General %>%
  dplyr::select(Connector, Satellite, Provincial, Peripheral) %>% 
  acomp(.) %>%
  # Preserves the ratios between non-zero components
  cmultRepl(., output = "prop")

data_coda_interareal <- TFP_General %>%
  dplyr::select(Global_Bridge, Local_Bridge, Super_Bridge, Not_a_Bridge) %>%
  acomp(.) %>%
  cmultRepl(., output = "prop")

# The compositional data set is expressed in isometric logratio coordinates.
# Then, robust principal component analysis is performed with a robust covariance MCD-estimator.
# ILR - moving the D-part compositional data from the simplex to a D-1 dimensional real space
# Orthonormal basis recommended by Egozcue et al. (2003)

ilrV <- function(x) {
  x.ilr <- matrix(NA, nrow = nrow(x), ncol = ncol(x) - 1)
  for (i in 1:ncol(x.ilr)) {
    x.ilr[, i] <- sqrt((i) / (i + 1)) * log(((apply(as.matrix(x[, 1:i]), 1, prod))^(1 / i)) / (x[, i + 1]))
  }
  return(x.ilr)
}
data_ilr <- cbind(ilrV(data_coda_modular), ilrV(data_coda_interareal))

# Correlation in a symmetric Procrustes rotation: 0.8901
vegan::protest(
  TFP_General %>%
    dplyr::select(
      Connector, Satellite, Provincial, Peripheral,
      Global_Bridge, Local_Bridge, Super_Bridge, Not_a_Bridge
    ) %>%
    as.matrix(),
  data_ilr
)

################################################################################
################################################################################

# PCA of ILR-transformed data because a non-singular covariance matrix is needed/ robust covariance estimation need a full-rank matrix
# CLR removes the value-range restriction but not the unit-sum constraint which makes PCA sensitive to outliers ----
set.seed(4)
# library(PPcovMcd)
cv <- robustbase::covMcd(data_ilr, nsamp = "deterministic")
pcaIlr <- princomp(data_ilr, covmat = cv)
pcaIlr$scores
biplot(pcaIlr, scale = 0)
plot(cv$mah)


# CLR-backtransform for interpretability of compositional variability and compositional biplot
data_CODA <- cbind(data_coda_modular, data_coda_interareal) %>% as.data.frame()

# Back-transform to clr fo rinterpretability of compositional variability
robCODA <- robCompositions::pcaCoDa(data_CODA,
                                    method = "robust", solve = "eigen",
                                    mult_comp = list(
                                      c(1, 2, 3, 4), c(5, 6, 7, 8)
                                    )
)


summary(robCODA)
biplot(robCODA,
       scale = 0,
       xlabs = rep("o", 627),
       col = c("red", "black"),
       xlim = c(-2, 2),
       cex = c(0.5, 1),
       choices = c(1, 2)
)



plot(robCODA, type = "l")
robCODA$scores
robCODA$loadings

# Outlier detection
par(mfrow = c(2, 1))
res <- mvoutlier.CoDa(data_CODA)
plot(res, which = "parallel", onlyout = TRUE)
plot(outCoDa(data_CODA))
