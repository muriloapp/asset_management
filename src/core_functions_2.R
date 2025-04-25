library(dplyr)
library(tidyr)
library(ggplot2)
library(mvtnorm) 
library(MASS)    
library(RiskPortfolios) 
library(Matrix)  
library(here)
library(fitHeavyTail) 


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

make_pd <- function(S)            
  if (min(eigen(S, TRUE)$values) < 1e-8) as.matrix(nearPD(S)$mat) else S


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
        perf <- metrics(w, R)          # <-- evaluated on *original* rets
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

make_pd <- function(S) {
  if (any(!is.finite(S))) return(diag(mean(diag(S), na.rm = TRUE), ncol(S)))
  ev <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) <= 1e-8) as.matrix(nearPD(S)$mat) else S
}
cov_fun <- function(R) make_pd(cov(R, use = "pairwise.complete.obs"))
sharpe  <- function(w, R) { w <- w/sum(w); p <- drop(R %*% w); mean(p)/sd(p) }

draw_iid   <- function(R) R[sample.int(nrow(R), replace = TRUE), ]
draw_block_rows <- function(R, l)
  R[tseries::tsbootstrap(1:nrow(R), nb = 1, statistic = NULL,
                         b = l, type = "stationary"), , drop = FALSE]
draw_gauss <- function(R) mvtnorm::rmvnorm(nrow(R), colMeans(R), cov(R))


boot_portfolio_stats <- function(rets,
                                 schemes   = c("iid","block","gauss"),
                                 ports     = c("minvol","maxdiv"),
                                 B         = 100,
                                 block_len = 20) {
  
  weight_store <- list()
  metric_store <- list()
  
  for (sc in schemes) {
    draw_fun <- switch(sc,
                       iid   = draw_iid,
                       block = function(R) draw_block_rows(R, l = block_len),
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























