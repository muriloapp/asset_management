library("here")
library(RiskPortfolios)


compute_in_sample_metrics <- function(weights, returns_matrix) {
  port_ret <- rowSums(returns_matrix * matrix(weights, nrow = nrow(returns_matrix), 
                                              ncol = ncol(returns_matrix), byrow = TRUE))
  list(mean_return = mean(port_ret), vol = sd(port_ret)
  )
}

analyze_dim_impact <- function(rets, dims = c(5,10,25), num_subsets = 5) {
  # rets: T x N matrix (T = time, N = total assets)
  # dims: integer vector, e.g. c(5, 10, 25)
  # num_subsets: how many random sub-universes to sample for each d
  Sigma_full <- cov(rets)
  n_assets   <- ncol(rets)
  
  # 2 portfolios (minvol, maxdiv) per subset, per dimension:
  n_rows <- length(dims) * num_subsets * 2
  
  results <- data.frame(
    dimension   = integer(n_rows),
    subset_id   = integer(n_rows),
    portfolio   = character(n_rows),
    mean_return = numeric(n_rows),
    vol         = numeric(n_rows),
    stringsAsFactors = FALSE
  )
  
  row_i <- 1
  portfolio_types <- c("minvol", "maxdiv")
  
  for (d in dims) {
    for (s in seq_len(num_subsets)) {
      idx <- sample(seq_len(n_assets), d)
      Sigma_sub <- Sigma_full[idx, idx, drop = FALSE]
      rets_sub  <- rets[, idx, drop = FALSE]
      
      for (ptype in portfolio_types) {
        w <- optimalPortfolio(
          Sigma_sub,
          control = list(type = ptype, constraint = "lo")
        )
  
        perf <- compute_in_sample_metrics(w, rets_sub)
        # store results
        results$dimension[row_i]   = d
        results$subset_id[row_i]   = s
        results$portfolio[row_i]   = ptype
        results$mean_return[row_i] = perf$mean_return
        results$vol[row_i]         = perf$vol
        
        row_i <- row_i + 1
      }
    }
  }
  return(results)
}


build_scaled_cov <- function(corr_mat, sd_vec, factor) {
  # Build scaled covariance from a correlation matrix and a vector of stdev * factor
  # If sd_vec is length d, scale it by 'factor'
  sd_scaled <- sd_vec * factor
  D <- diag(sd_scaled)
  D %*% corr_mat %*% D
}

portfolio_performance_matrix <- function(w, Sigma_scaled, mean_vec) {
  # Compute portfolio mean and volatility using matrix approach
  # w: vector of weights (d)
  # Sigma_scaled: d x d
  # mean_vec: vector of average returns (d)
  var_p <- as.numeric(t(w) %*% Sigma_scaled %*% w)
  vol_p <- sqrt(var_p)
  mu_p  <- sum(w * mean_vec)
  list(vol = vol_p, mean_return = mu_p)
}

analyze_variance_impact <- function(rets, dims = c(5,10,25), num_subsets = 3,
                                    scale_factors = c(0.5,1,1.5,2)) {
  # rets: T x N matrix
  # dims: vector of chosen portfolio sizes
  # num_subsets: how many random sub-universes for each dimension
  # scale_factors: e.g. c(0.5, 1, 1.5, 2)
  
  n_assets <- ncol(rets)
  # 2 portfolios * length(scale_factors) * num_subsets * length(dims) total rows
  n_rows <- 2 * length(scale_factors) * num_subsets * length(dims)
  
  results <- data.frame(
    dimension    = integer(n_rows),
    subset_id    = integer(n_rows),
    scale_factor = numeric(n_rows),
    portfolio    = character(n_rows),
    vol          = numeric(n_rows),
    mean_return  = numeric(n_rows),
    stringsAsFactors = FALSE
  )
  
  row_i <- 1
  portfolio_types <- c("minvol", "maxdiv")
  
  for (d in dims) {
    for (s in seq_len(num_subsets)) {
      idx <- sample(seq_len(n_assets), d)
      rets_sub <- rets[, idx, drop = FALSE]
      
      # Estimate correlation & stdev from the sub-universe
      Sigma_sub <- cov(rets_sub)
      sd_sub    <- sqrt(diag(Sigma_sub))
      corr_sub  <- diag(1/sd_sub) %*% Sigma_sub %*% diag(1/sd_sub)
      
      mean_sub  <- colMeans(rets_sub)
      
      for (f in scale_factors) {
        Sigma_scaled <- build_scaled_cov(corr_sub, sd_sub, f)
        
        for (ptype in portfolio_types) {
          w <- optimalPortfolio(
            Sigma_scaled,
            control = list(type = ptype, constraint = "lo")
          )
          
          perf <- portfolio_performance_matrix(w, Sigma_scaled, mean_sub)
          
          # Store results
          results$dimension[row_i]    = d
          results$subset_id[row_i]    = s
          results$scale_factor[row_i] = f
          results$portfolio[row_i]    = ptype
          results$vol[row_i]          = perf$vol
          results$mean_return[row_i]  = perf$mean_return
          
          row_i <- row_i + 1
        }
      }
    }
  }
  
  results
}

equicorr_matrix <- function(d, rho) {
  # Build an equicorrelation matrix (d x d) with correlation 'rho'
  R <- matrix(rho, nrow = d, ncol = d)
  diag(R) <- 1
  R
}

build_equi_cov <- function(sd_vec, rho) {
  # Build covariance from eqicorr matrix & stdev vector
  d <- length(sd_vec)
  R <- equicorr_matrix(d, rho)
  D <- diag(sd_vec)
  Sigma_rho <- D %*% R %*% D
  Sigma_rho
}

analyze_correlation_impact <- function(rets, dims = c(5,10,25), num_subsets = 100,
                                       corr_values = seq(0.1, 0.9, by=0.1)) {
  # rets: T x N
  # dims: vector of portfolio sizes, e.g. c(5,10)
  # num_subsets: number of random draws for each dimension
  # corr_values: vector of correlation values, e.g. seq(0.1,1,0.2)
  
  n_assets <- ncol(rets)
  # 2 portfolios (minvol, maxdiv) * length(corr_values) * num_subsets * length(dims)
  n_rows <- 2 * length(corr_values) * num_subsets * length(dims)
  
  results <- data.frame(
    dimension   = integer(n_rows),
    subset_id   = integer(n_rows),
    correlation = numeric(n_rows),
    portfolio   = character(n_rows),
    vol         = numeric(n_rows),
    mean_return = numeric(n_rows),
    stringsAsFactors = FALSE
  )
  
  row_i <- 1
  portfolio_types <- c("minvol", "maxdiv")
  
  for (d in dims) {
    for (s in seq_len(num_subsets)) {
      # pick random subset of size d
      idx <- sample(seq_len(n_assets), d)
      rets_sub <- rets[, idx, drop = FALSE]
      
      # from rets_sub, get stdev & mean
      Sigma_sub <- cov(rets_sub)
      sd_sub    <- sqrt(diag(Sigma_sub))
      mean_sub  <- colMeans(rets_sub)
      
      for (rho in corr_values) {
        # Build equicorr covariance
        Sigma_rho <- build_equi_cov(sd_sub, rho)
        
        for (ptype in portfolio_types) {
          # Solve for portfolio weights
          w <- optimalPortfolio(
            Sigma_rho,
            control = list(type = ptype, constraint = "lo")  # long-only
          )

          perf <- portfolio_performance_matrix(w, Sigma_rho, mean_sub)
          
          results$dimension[row_i]   = d
          results$subset_id[row_i]   = s
          results$correlation[row_i] = rho
          results$portfolio[row_i]   = ptype
          results$vol[row_i]         = perf$vol
          results$mean_return[row_i] = perf$mean_return
          
          row_i <- row_i + 1
        }
      }
    }
  }
  
  results
}
