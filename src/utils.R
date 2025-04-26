library(dplyr)
library(tidyr)
library(mvtnorm) 
library(MASS)    
library(RiskPortfolios) 
library(Matrix)  
library(here)

make_pd <- function(S)            
  if (min(eigen(S, TRUE)$values) < 1e-8) as.matrix(nearPD(S)$mat) else S

cov_fun <- function(R) make_pd(cov(R, use = "pairwise.complete.obs"))

sharpe  <- function(w, R) { w <- w/sum(w); p <- drop(R %*% w); mean(p)/sd(p) }

draw_iid   <- function(R) R[sample.int(nrow(R), replace = TRUE), ]

draw_block <- function(R, l)
  R[tseries::tsbootstrap(1:nrow(R), nb = 1, statistic = NULL,
                         b = l, type = "stationary"), , drop = FALSE]
draw_gauss <- function(R) mvtnorm::rmvnorm(nrow(R), colMeans(R), cov(R))

prep_data <- function(Rfull, n_assets = 25) {
  as.matrix(Rfull[, 1:n_assets])
}
