library(infotheo)

# EMI, NMI & AMI functions-----------------------------------------
EMI_func <- function(s1, s2, l1, l2, n) { # Expected mutual information
  s_emi <- 0
  for (i in 1:l1) {
    for (j in 1:l2) {
      min_nij <- max(1, s1[i] + s2[j] - n)
      max_nij <- min(s1[i], s2[j])
      n.ij <- seq(min_nij, max_nij)
      t1 <- (n.ij / n) * log((n.ij * n) / (s1[i] * s2[j]))
      t2 <- exp(lfactorial(s1[i]) + lfactorial(s2[j]) + lfactorial(n - s1[i]) + lfactorial(n - s2[j]) - lfactorial(n) - lfactorial(n.ij) - lfactorial(s1[i] - n.ij) - lfactorial(s2[j] - n.ij) - lfactorial(n - s1[i] - s2[j] + n.ij))
      emi <- sum(t1 * t2)
      s_emi <- s_emi + emi
    }
  }
  return(s_emi)
}

NMI_func <- function(v1, v2) {
  library(infotheo)
  s1 <- tabulate(v1)
  s2 <- tabulate(v2)
  l1 <- length(s1)
  l2 <- length(s2)
  N <- length(v1)
  tij <- table(v1, v2, useNA = "always")
  mi <- mutinformation(v1, v2) # function for Mutual Information from package infotheo
  h1 <- -sum(s1 * log(s1 / N)) / N
  h2 <- -sum(s2 * log(s2 / N)) / N
  nmi <- mi / max(h1, h2) # NMI Normalized MI
  emi <- EMI_func(s1, s2, l1, l2, N) # EMI Expected MI
  return(nmi)
}

AMI_func <- function(v1, v2) {
  library(infotheo)
  s1 <- tabulate(v1)
  s2 <- tabulate(v2)
  l1 <- length(s1)
  l2 <- length(s2)
  N <- length(v1)
  tij <- table(v1, v2, useNA = "always")
  mi <- mutinformation(v1, v2) # function for Mutual Information from package infotheo
  h1 <- -sum(s1 * log(s1 / N)) / N
  h2 <- -sum(s2 * log(s2 / N)) / N
  emi <- EMI_func(s1, s2, l1, l2, N) # EMI Expected MI
  ami <- (mi - emi) / (max(h1, h2) - emi) # AMI Adjusted MI
  return(ami)
}
