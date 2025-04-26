library(dplyr)
library(tidyr)
library(ggplot2)
library(mvtnorm) 
library(MASS)    
library(RiskPortfolios) 
library(Matrix)  
library(here)
library(fitHeavyTail) 

source(here('src','utils.R'))

# 1.

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

# 2.1; 2.2

compute_in_sample_metrics_sharpe <- function(weights, returns_matrix, risk_free_rate = 0) {
  if(ncol(returns_matrix) != length(weights)) { stop("Dim mismatch: weights vs returns.")}
  weights <- weights / sum(weights)
  port_ret <- returns_matrix %*% weights
  mean_ret <- mean(port_ret)
  vol <- sd(port_ret)
  sharpe <- ifelse(vol == 0, NA, (mean_ret - risk_free_rate) / vol) 
  list(mean = mean_ret, vol = vol, sharpe = sharpe)
}

run_bootstrap_simulation <- function(rets_orig, model_type = "gaussian",
                                     B = 1000, seed = 123) {
  set.seed(seed)
  T_obs <- nrow(rets_orig)
  n_assets <- ncol(rets_orig)
  portfolio_types <- c("minvol", "maxdiv")
  params <- list()
  # Calibration Gaussian or Student-t
  if (model_type == "gaussian") {
    message("Calibrating Gaussian model...")
    params$mu <- colMeans(rets_orig)
    params$Sigma <- cov(rets_orig)
    eigen_check <- eigen(params$Sigma, symmetric = TRUE, only.values = TRUE)$values
    if(any(eigen_check <= 1e-8)) {
      message("Using nearPD.")
      params$Sigma <- as.matrix(Matrix::nearPD(params$Sigma)$mat)
    }
    message("Gaussian calibration done.")
  } else if (model_type == "student_t") {
    message("Calibrating Multivariate Student-t model")
    fit_t <- fitHeavyTail::fit_mvt(rets_orig) 
    params$mu <- fit_t$mu; params$Sigma <- fit_t$scatter; params$nu <- fit_t$nu
    message(paste("Student-t calibration done. Est. nu =", round(params$nu, 2)))
    eigen_check <- eigen(params$Sigma, symmetric = TRUE, only.values = TRUE)$values
    if(any(eigen_check <= 1e-8)) {
      message("Using nearPD.")
      params$Sigma <- as.matrix(Matrix::nearPD(params$Sigma)$mat)
    }
  }
  
  # Storage 
  results_list <- vector("list", B * length(portfolio_types))
  res_counter <- 1
  
  # Simulation
  for (b in 1:B) {
    if (b %% (B/10) == 0) print(paste("Model:", model_type, "- Sim:", b, "/", B))
    rets_sim <- if (model_type == "gaussian") {
      MASS::mvrnorm(n = T_obs, mu = params$mu, Sigma = params$Sigma)
    } else {
      mvtnorm::rmvt(n = T_obs, delta = params$mu, sigma = params$Sigma, df = params$nu)
    }
    Sigma_sim <- cov(rets_sim)
    eigen_sim <- eigen(Sigma_sim, symmetric = TRUE, only.values = TRUE)$values
    if(any(eigen_sim <= 1e-8)) {
      Sigma_sim <- as.matrix(Matrix::nearPD(Sigma_sim)$mat)
    }
    
    for (ptype in portfolio_types) {
      w <- tryCatch({
        optimalPortfolio(Sigma_sim, control = list(type = ptype, constraint = "lo"))
      }, error = function(e) {
        warning(paste("optimalPortfolio failed: Model=", model_type, ", Sim=", b, ", Type=", ptype, ". Error:", e$message), call. = FALSE, immediate. = TRUE)
        NULL 
      })
      
      if (!is.null(w)) {
        perf <- compute_in_sample_metrics_sharpe(w, rets_sim)
        results_list[[res_counter]] <- data.frame(
          simulation_id = b, portfolio = ptype,
          mean = perf$mean, vol = perf$vol, sharpe = perf$sharpe
        )
        res_counter <- res_counter + 1
      }
    } 
  } 
  
  results_df <- dplyr::bind_rows(results_list) 
  
  return(results_df)
}

# 2.3

metrics <- function(w, R) {
  w <- w / sum(w)
  pr <- drop(R %*% w)
  c(mean = mean(pr),
    vol  = sd(pr),
    SR   = mean(pr) / sd(pr))
}

one_draw <- function(R, portfolios = c("minvol","maxdiv")) {
  Rb <- R[sample.int(nrow(R), replace = TRUE), ]      
  
  # Covariance estimators
  S_sample <- make_pd(cov(Rb))
  S_LW <- covEstimation(rets, control = list(type = 'lw'))       
  S_F <- covEstimation(rets, control = list(type = 'factor', K = 3))
  cov_set <- list(sample = S_sample, LW = S_LW, FM = S_F)
  
  # Build portfolios & evaluate 
  res <- vector("list", length(portfolios) * length(cov_set))
  k <- 1
  for (ptype in portfolios)
    for (cname in names(cov_set)) {
      w <- tryCatch(
        optimalPortfolio(cov_set[[cname]],
                         control = list(type = ptype, constraint = "lo")),
        error = function(e) NULL)
      if (!is.null(w)) {
        perf <- metrics(w, R)         
        res[[k]] <- data.frame(portfolio = ptype,
                               cov_est  = cname,
                               t(perf))
        k <- k + 1
      }
    }
  dplyr::bind_rows(res)
}

# Main bootstrap loop 
run_boot <- function(R, B = 500, seed = 1) {
  set.seed(seed)
  out <- vector("list", B)
  for (b in seq_len(B)) {
    if (b %% (B/10) == 0) message("bootstrap ", b, "/", B)
    out[[b]] <- one_draw(R)
  }
  bind_rows(out) |>
    mutate(cov_est = factor(cov_est, levels = c("sample","LW","FM")))
}

# 2.4

boot_portfolio_stats <- function(rets,
                                 schemes   = c("iid","block","gauss"),
                                 ports     = c("minvol","maxdiv"),
                                 B         = 100,
                                 block_len = 10) {
  
  weight_store <- list()
  metric_store <- list()
  
  for (sc in schemes) {
    draw_fun <- switch(sc,
                       iid   = draw_iid,
                       block = function(R) draw_block(R, l = block_len),
                       gauss = draw_gauss)
    
    for (pt in ports) {
      tag     <- paste(sc, pt, sep = "_")
      Wmat    <- matrix(NA, B, ncol(rets))
      metrics <- data.frame(mean = numeric(0),
                            vol  = numeric(0),
                            sharpe = numeric(0))
      
      for (b in seq_len(B)) {
        Rb <- na.omit(draw_fun(rets));  if (nrow(Rb) < 2) next
        w  <- tryCatch(
          optimalPortfolio(make_pd(cov(Rb)),
                           control = list(type = pt, constraint = "lo")),
          error = function(e) NA)
        if (anyNA(w)) next
        
        pr <- drop(Rb %*% w)
        Wmat[b, ] <- w
        metrics   <- rbind(metrics,
                           data.frame(mean   = mean(pr),
                                      vol    = sd(pr),
                                      sharpe = mean(pr)/sd(pr)))
      }
      
      weight_store[[tag]] <- Wmat
      metric_store[[tag]] <- transform(metrics, tag = tag)
    }
  }
  
  list(weights_raw = weight_store,
       metrics_raw = metric_store)
}

# 3.1; 3.2

get_cov <- function(R, method = c("sample", "lw", "factor3")) {
  method <- match.arg(method)
  switch(method,
         sample  = make_pd(cov(R)),
         lw      = RiskPortfolios::covEstimation(R, control = list(type = "lw")),
         factor3 = RiskPortfolios::covEstimation(R, control = list(type = "factor", K = 3))
  )
}

get_weights <- function(S, port_type = c("minvol","maxdiv")) {
  optimalPortfolio(S, control = list(type = match.arg(port_type), constraint = "lo"))
}

oos_one_step <- function(Rtrain, Rnext, port_type, cov_method) {
  S <- get_cov(Rtrain, cov_method)
  w <- get_weights(S, port_type)
  sum(w * Rnext) # one-week OOS return
}

rolling_oos <- function(R25,
                        window_len = 104,
                        port_types = c("minvol","maxdiv"),
                        cov_methods = c("sample","lw","factor3"))
{
  stopifnot(nrow(R25) > window_len + 1)
  
  out <- data.frame()       # store results per week × portfolio × cov
  for (t in window_len:(nrow(R25) - 1)) {
    Rtrain <- R25[(t - window_len + 1):t, ]
    Rnext  <- R25[t + 1, ]
    
    for (p in port_types)
      for (cm in cov_methods) {
        ret <- oos_one_step(Rtrain, Rnext, p, cm)
        out <- rbind(out,
                     data.frame(week = t + 1,
                                portfolio = p,
                                cov_est   = cm,
                                ret       = ret))
      }
  }
  out
}

# 3.3

avg_weights <- function(R, draw_fun, B, port_type) {
  W <- matrix(NA, B, ncol(R))
  for (b in 1:B) {
    Rb <- draw_fun(R)
    S  <- make_pd(cov(Rb))
    W[b,] <- tryCatch(
      optimalPortfolio(S, control = list(type = port_type,
                                         constraint = "lo")),
      error = function(e) NA)
  }
  colMeans(W, na.rm = TRUE)
}


rolling_oos_resamp <- function(R25, window_len = 104,
                               B = 100, block_len = 10,
                               schemes = c("iid","block","gauss"),
                               ports   = c("minvol","maxdiv"))
{
  out <- data.frame()
  for (t in window_len:(nrow(R25) - 1)) {
    
    Rtrain <- R25[(t - window_len + 1):t, ]
    Rnext  <- R25[t + 1, ]
    
    for (sc in schemes) {
      draw_fun <- switch(sc,
                         iid   = draw_iid,
                         block = function(R) draw_block(R, block_len),
                         gauss = draw_gauss)
      
      for (pt in ports) {
        w_bar <- avg_weights(Rtrain, draw_fun, B, pt)
        ret   <- sum(w_bar * Rnext)
        out   <- rbind(out,
                       data.frame(week = t + 1,
                                  resample = sc,
                                  portfolio = pt,
                                  ret = ret))
      }
    }
  }
  out
}

