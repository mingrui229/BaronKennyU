#' @keywords internal
trans <- function(x){
  n <- length(x)
  y <- rep(1, n+1)
  for (i in c(1: n)){
    y[i] <- y[i] * cos(x[i])
    y[(i+1):(n+1)] <- y[(i+1):(n+1)] * sin(x[i])
  }
  return(y)
}

#' @keywords internal
direct_compute_T <- function(y, m, a, c){
  if (is.null(c)){
    y <- as.matrix(lm(y~1)$residuals)
    m <- as.matrix(lm(m~1)$residuals)
    a <- as.matrix(lm(a~1)$residuals)
  }else{
    y <- as.matrix(lm(y~1+c)$residuals)
    m <- as.matrix(lm(m~1+c)$residuals)
    a <- as.matrix(lm(a~1+c)$residuals)
  }
  T0 <- lm(y~a+m)$coefficients[2]
  c1 <- var(lm(y~a+m)$residuals) / var(lm(a~m)$residuals)
  c2 <- var(lm(a~m)$residuals) / var(a)
  T1 <- -sqrt(c1) * sqrt(c2)
  T2 <- (sqrt(c1) / sqrt(c2) / sqrt(var(a))) %*% cov(a, m) %*% solve(cov(m)) %*%
    expm::sqrtm(cov(as.matrix(lm(m~a)$residuals)))
  return(c(T0, T1, T2))
}

#' @keywords internal
direct_point_solve <- function(T0, T1, T2, rho_y, rho_m, rho_a, algorithm){
  dim_m <- length(T2)

  if (dim_m > 1){
    fun <- function(x){
      x1 <- x[1]
      x2 <- trans(x[2: dim_m])
      x3 <- x[dim_m+1]
      f <- T0 + T1 * sqrt(rho_y * rho_a / (1-rho_a)) * x1 / cos(x3) + sum(T2*x2) * sqrt(rho_y) * tan(x3)
      return((2*(T0>=0)-1) * f)
    }
    x0 <- c(runif(1, -1, 1), runif(dim_m-1, 0, 2*pi), runif(1, 0, atan(sqrt(rho_m/(1-rho_m)))))
    reT0 <- nloptr::nloptr(x0=x0,
                   eval_f = fun,
                   lb = c(-1, rep(0, dim_m-1), 0),
                   ub = c(1, rep(2*pi, dim_m-1), atan(sqrt(rho_m/(1-rho_m)))),
                   opts = list(algorithm = algorithm,
                               xtol_rel= 1.0e-8,
                               maxeval = 100000))
  }else{
    fun <- function(x){
      x1 <- x[1]
      x3 <- x[2]
      f1 <- T0 + T1 * sqrt(rho_y * rho_a / (1-rho_a)) * x1 / cos(x3) + T2 * sqrt(rho_y) * tan(x3)
      f2 <- T0 + T1 * sqrt(rho_y * rho_a / (1-rho_a)) * x1 / cos(x3) - T2 * sqrt(rho_y) * tan(x3)
      return(min((2*(T0>=0)-1) * f1, (2*(T0>=0)-1) * f2))
    }
    x0 <- c(runif(1, -1, 1), runif(1, 0, atan(sqrt(rho_m/(1-rho_m)))))
    reT0 <- nloptr::nloptr(x0=x0,
                   eval_f = fun,
                   lb = c(-1, 0),
                   ub = c(1, atan(sqrt(rho_m/(1-rho_m)))),
                   opts = list(algorithm = algorithm,
                               xtol_rel= 1.0e-8,
                               maxeval = 100000))
  }

  return(reT0$objective)
}

#' @keywords internal
direct_t_solve <- function(T0, T1, T2, T0_bootstrap, T1_bootstrap, T2_bootstrap,
                           rho_y, rho_m, rho_a, algorithm){
  dim_m <- length(T2)

  if (dim_m > 1){
    fun <- function(x){
      x1 <- x[1]
      x2 <- trans(x[2: dim_m])
      x3 <- x[dim_m+1]
      f <- T0 + T1 * sqrt(rho_y * rho_a / (1-rho_a)) * x1 / cos(x3) + sum(T2*x2) * sqrt(rho_y) * tan(x3)
      f_bootstrap <- T0_bootstrap + T1_bootstrap * sqrt(rho_y * rho_a / (1-rho_a)) * x1 / cos(x3) + (T2_bootstrap %*% x2) * sqrt(rho_y) * tan(x3)
      return((2*(T0>=0)-1) * f / sd(f_bootstrap))
    }
    x0 <- c(runif(1, -1, 1), runif(dim_m-1, 0, 2*pi), runif(1, 0, atan(sqrt(rho_m/(1-rho_m)))))
    reT0 <- nloptr::nloptr(x0=x0,
                   eval_f = fun,
                   lb = c(-1, rep(0, dim_m-1), 0),
                   ub = c(1, rep(2*pi, dim_m-1), atan(sqrt(rho_m/(1-rho_m)))),
                   opts = list(algorithm = algorithm,
                               xtol_rel= 1.0e-8,
                               maxeval = 100000))
  }else{
    fun <- function(x){
      x1 <- x[1]
      x3 <- x[2]
      f1 <- T0 + T1 * sqrt(rho_y * rho_a / (1-rho_a)) * x1 / cos(x3) + T2 * sqrt(rho_y) * tan(x3)
      f2 <- T0 + T1 * sqrt(rho_y * rho_a / (1-rho_a)) * x1 / cos(x3) - T2 * sqrt(rho_y) * tan(x3)
      f1_bootstrap <- T0_bootstrap + T1_bootstrap * sqrt(rho_y * rho_a / (1-rho_a)) * x1 / cos(x3) + T2_bootstrap * sqrt(rho_y) * tan(x3)
      f2_bootstrap <- T0_bootstrap + T1_bootstrap * sqrt(rho_y * rho_a / (1-rho_a)) * x1 / cos(x3) - T2_bootstrap * sqrt(rho_y) * tan(x3)
      return(min((2*(T0>=0)-1) * f1 / sd(f1_bootstrap), (2*(T0>=0)-1) * f2 / sd(f2_bootstrap)))
    }
    x0 <- c(runif(1, -1, 1), runif(1, 0, atan(sqrt(rho_m/(1-rho_m)))))
    reT0 <- nloptr::nloptr(x0=x0,
                   eval_f = fun,
                   lb = c(-1, 0),
                   ub = c(1, atan(sqrt(rho_m/(1-rho_m)))),
                   opts = list(algorithm = algorithm,
                               xtol_rel= 1.0e-8,
                               maxeval = 100000))
  }

  return(reT0$objective)
}

#' Compute the robustness values for the direct effect
#' @param y a vector of the outcome variable
#' @param m a matrix of the mediator variables
#' @param a a vector of the exposure variable
#' @param c a matrix of the covariates
#' @param randomized if the exposure is randomized
#' @param rho_values the values of rho's to be considered for the grid search
#' @param algorithm the optimization algorithm in the R package "nloptr"
#' @return a summary table, the worst point estimates and the worst t statistics given each value of rho.
#' @examples
#' set.seed(1234)
#' library(MASS)
#' n = 200
#' a = rnorm(n)
#' c = as.matrix(rnorm(n))
#' m = mvrnorm(n, mu = c(0, 0), Sigma = matrix(c(1,0.2,0.2,0.8), 2, 2))
#' m = m + cbind(2.5 * a, 1.5 * a)
#' y = a + 0.2 * m[, 1] + 0.25 * m[, 2] + rnorm(n)
#' result = bku_rv_direct(y, m, a, c)
#' @importFrom stats cor cov lm runif sd var
#' @importFrom MASS mvrnorm
#' @importFrom expm sqrtm
#' @importFrom nloptr nloptr
#' @export
bku_rv_direct <- function(y, m, a, c = NULL, randomized = F, rho_values = c(1:99)/100, algorithm = "NLOPT_GN_ORIG_DIRECT_L"){
  n <- length(y)
  m <- as.matrix(m)
  Ts <- direct_compute_T(y, m, a, c)
  T0 <- Ts[1]
  T1 <- Ts[2]
  T2 <- Ts[-c(1,2)]
  sample_index <- replicate(1000, sample(n, n, replace = T))
  Ts_bootstrap <- apply(sample_index, 2, function(x) direct_compute_T(y[x], m[x, ], a[x], c[x, ]))
  T0_bootstrap <- Ts_bootstrap[1, ]
  T1_bootstrap <- Ts_bootstrap[2, ]
  T2_bootstrap <- t(Ts_bootstrap[-c(1,2), ])
  point_rho <- matrix(nrow = length(rho_values), ncol = 4)
  t_rho <- matrix(nrow = length(rho_values), ncol = 4)
  for (i in c(1: length(rho_values))){
    point_rho[i, 1] <- rho_values[i]
    point_rho[i, 2] <- rho_values[i]
    point_rho[i, 3] <- rho_values[i] * (1-randomized)
    t_rho[i, 1] <- rho_values[i]
    t_rho[i, 2] <- rho_values[i]
    t_rho[i, 3] <- rho_values[i] * (1-randomized)
    point_rho[i, 4] <- direct_point_solve(T0, T1, T2, point_rho[i, 1], point_rho[i, 2], point_rho[i, 3], algorithm)
    t_rho[i, 4] <- direct_t_solve(T0, T1, T2, T0_bootstrap, T1_bootstrap, T2_bootstrap, t_rho[i, 1], t_rho[i, 2], t_rho[i, 3], algorithm)
  }
  rho_values <- c(0, rho_values)
  point_rv <- rho_values[min(which(point_rho[, 4] < 0))]
  t_rv <- rho_values[min(which(t_rho[, 4] < 1.96))]
  colnames(point_rho) <- c("rho_y", "rho_m", "rho_a", "worst point estimate")
  colnames(t_rho) <- c("rho_y", "rho_m", "rho_a", "worst t statistic")
  point_rho[, 4] <- point_rho[, 4] * (2*(T0>=0)-1)
  t_rho[, 4] <- t_rho[, 4] * (2*(T0>=0)-1)
  point_obs <- T0
  se_obs <- sd(T0_bootstrap)
  t_obs <- point_obs / se_obs
  summary_tab <- c(point_obs, se_obs, t_obs, point_rv, t_rv)
  names(summary_tab) <- c("Est.", "Std. Error", "t value", "R.V. for Est.", "R.V. for 95% C.I.")
  return(list("summary_table" = summary_tab, "worst_point" = point_rho, "worst_t" = t_rho))
}

#' Compute the sensitivity bound for the direct effect with prespecified sensitivity parameters
#' @param y a vector of the outcome variable
#' @param m a matrix of the mediator variables
#' @param a a vector of the exposure variable
#' @param c a matrix of the covariates
#' @param Ry a number of the R parameter for the outcome-confounder correlation
#' @param Rm a vector of the R parameter for the mediator-confounder correlation
#' @param Ra a number of the R parameter for the exposure-confounder correlation
#' @return the point estimate, the standard error, and the t statistic
#' @examples
#' set.seed(1234)
#' library(MASS)
#' n = 200
#' a = rnorm(n)
#' c = as.matrix(rnorm(n))
#' m = mvrnorm(n, mu = c(0, 0), Sigma = matrix(c(1,0.2,0.2,0.8), 2, 2))
#' m = m + cbind(2.5 * a, 1.5 * a)
#' y = a + 0.2 * m[, 1] + 0.25 * m[, 2] + rnorm(n)
#' result = bku_direct(y, m, a, c, 0.1, c(0.2, 0.1), 0.05)
#' @importFrom stats cor cov lm runif sd var
#' @importFrom MASS mvrnorm
#' @importFrom expm sqrtm
#' @export
bku_direct <- function(y, m, a, c = NULL, Ry, Rm, Ra){
  n <- length(y)
  m <- as.matrix(m)
  Ts <- direct_compute_T(y, m, a, c)
  T0 <- Ts[1]
  T1 <- Ts[2]
  T2 <- Ts[-c(1,2)]
  sample_index <- replicate(1000, sample(n, n, replace = T))
  Ts_bootstrap <- apply(sample_index, 2, function(x) direct_compute_T(y[x], m[x, ], a[x], c[x, ]))
  T0_bootstrap <- Ts_bootstrap[1, ]
  T1_bootstrap <- Ts_bootstrap[2, ]
  T2_bootstrap <- t(Ts_bootstrap[-c(1,2), ])
  phi1 <- (Ry * Ra) / sqrt(1 - Ra^2) / sqrt(1 - sqrt(sum(Rm^2)))
  phi2 <- (Ry * Rm) / sqrt(1 - sqrt(sum(Rm^2)))
  point_est <- T0 + T1 * phi1 + sum(T2 * phi2)
  if (length(T2) > 1){
    std_err <- sd(T0_bootstrap + T1_bootstrap * phi1 + T2_bootstrap %*% phi2)
  }else{
    std_err <- sd(T0_bootstrap + T1_bootstrap * phi1 + T2_bootstrap * phi2)
  }
  t_stat <- point_est / std_err
  return(c("point_est" = point_est, "std_err" = std_err, "t_stat" = t_stat))

}

#' @keywords internal
indirect_compute_T <- function(y, m, a, c){
  if (is.null(c)){
    y <- as.matrix(lm(y~1)$residuals)
    m <- as.matrix(lm(m~1)$residuals)
    a <- as.matrix(lm(a~1)$residuals)
  }else{
    y <- as.matrix(lm(y~1+c)$residuals)
    m <- as.matrix(lm(m~1+c)$residuals)
    a <- as.matrix(lm(a~1+c)$residuals)
  }
  beta1obs <- matrix(lm(m~a+0)$coefficients)
  theta3obs <- matrix(lm(y~a+m+0)$coefficients[-1])
  var_y_res <- c(var(lm(y~a+m+0)$residuals))
  cov_m_res <- cov(as.matrix(lm(m~a+0)$residuals))
  var_a_res <- c(var(a))

  T3 <- sum(beta1obs * theta3obs)
  T4 <- sqrt(var_y_res) / sqrt(var_a_res)
  T5 <- -expm::sqrtm(cov_m_res) %*% theta3obs / sqrt(var_a_res)
  T6 <- -solve(expm::sqrtm(cov_m_res)) %*% beta1obs * sqrt(var_y_res)
  return(c(T3, T4, T5, T6))
}

#' @keywords internal
indirect_point_solve <- function(T3, T4, T5, T6, rho_y, rho_m, rho_a, algorithm){
  dim_m <- length(T5)

  if (dim_m > 1){
    fun <- function(x){
      x1 <- trans(x[1: (dim_m-1)])
      x2 <- x[dim_m]
      x3 <- x[dim_m+1]
      f <- T3 + T4 * sqrt(rho_m * rho_a / (1-rho_a)) * sqrt(rho_y * rho_m / (1-rho_m)) * x2 * x3 +
        sqrt(rho_m * rho_a / (1-rho_a)) * sum(T5 * x1) * x2 +
        sqrt(rho_y * rho_m / (1-rho_m)) * sum(T6 * x1) * x3
      return((2*(T3>=0)-1) * f)
    }
    x0 <- c(runif(dim_m-1, 0, 2*pi), runif(2, -1, 1))
    reT3 <- nloptr::nloptr(x0=x0,
                   eval_f = fun,
                   lb = c(rep(0, dim_m-1), -1, -1),
                   ub = c(rep(2*pi, dim_m-1), 1, 1),
                   opts = list(algorithm = algorithm,
                               xtol_rel= 1.0e-8,
                               maxeval = 100000))
  }else{
    fun <- function(x){
      x2 <- x[1]
      x3 <- x[2]
      f1 <- T3 + T4 * sqrt(rho_m * rho_a / (1-rho_a)) * sqrt(rho_y * rho_m / (1-rho_m)) * x2 * x3 +
        sqrt(rho_m * rho_a / (1-rho_a)) * T5 * x2 +
        sqrt(rho_y * rho_m / (1-rho_m)) * T6 * x3
      f2 <- T3 + T4 * sqrt(rho_m * rho_a / (1-rho_a)) * sqrt(rho_y * rho_m / (1-rho_m)) * x2 * x3 -
        sqrt(rho_m * rho_a / (1-rho_a)) * T5 * x2 -
        sqrt(rho_y * rho_m / (1-rho_m)) * T6 * x3
      return(min((2*(T3>=0)-1) * f1, (2*(T3>=0)-1) * f2))
    }
    x0 <- c(runif(2, -1, 1))
    reT3 <- nloptr::nloptr(x0 = x0,
                   eval_f = fun,
                   lb = c(-1, -1),
                   ub = c(1, 1),
                   opts = list(algorithm = algorithm,
                               xtol_rel= 1.0e-8,
                               maxeval = 100000))
  }

  return(reT3$objective)
}

#' @keywords internal
indirect_point_solve_vu <- function(T3, T4, T5, T6, rho_y, rho_m, rho_a, algorithm){
  dim_m <- length(T5)

  if (dim_m > 1){
    fun <- function(x){
      x1 <- trans(x[1: (dim_m-1)]) * x[2*dim_m-1]
      x2 <- trans(x[dim_m: (2*dim_m-2)]) * x[2*dim_m]
      f <- T3 + T4 * sqrt(rho_m * rho_a / (1-rho_a)) * sqrt(rho_y * rho_m / (1-rho_m)) * sum(x1*x2) +
        sqrt(rho_m * rho_a / (1-rho_a)) * sum(T5 * x1) +
        sqrt(rho_y * rho_m / (1-rho_m)) * sum(T6 * x2)
      return((2*(T3>=0)-1) * f)
    }
    x0 <- c(runif(2*dim_m-2, 0, 2*pi), runif(2, 0, 1))
    reT3 <- nloptr::nloptr(x0 = x0,
                   eval_f = fun,
                   lb = c(rep(0, 2*dim_m-2), 0, 0),
                   ub = c(rep(2*pi, 2*dim_m-2), 1, 1),
                   opts = list(algorithm = algorithm,
                               xtol_rel= 1.0e-8,
                               maxeval = 100000))
  }else{
    fun <- function(x){
      x1 <- x[1]
      x2 <- x[2]
      f <- T3 + T4 * sqrt(rho_m * rho_a / (1-rho_a)) * sqrt(rho_y * rho_m / (1-rho_m)) * sum(x1*x2) +
        sqrt(rho_m * rho_a / (1-rho_a)) * sum(T5 * x1) +
        sqrt(rho_y * rho_m / (1-rho_m)) * sum(T6 * x2)
      return((2*(T3>=0)-1) * f)
    }
    x0 <- c(runif(2, -1, 1))
    reT3 <- nloptr::nloptr(x0=x0,
                   eval_f = fun,
                   lb = c(-1, -1),
                   ub = c(1, 1),
                   opts = list(algorithm = algorithm,
                               xtol_rel= 1.0e-8,
                               maxeval = 100000))
  }

  return(reT3$objective)
}

#' @keywords internal
indirect_t_solve <- function(T3, T4, T5, T6, T3_bootstrap, T4_bootstrap, T5_bootstrap, T6_bootstrap,
                             rho_y, rho_m, rho_a, algorithm){
  dim_m <- length(T5)

  if (dim_m > 1){
    fun <- function(x){
      x1 <- trans(x[1: (dim_m-1)])
      x2 <- x[dim_m]
      x3 <- x[dim_m+1]
      f <- T3 + T4 * sqrt(rho_m * rho_a / (1-rho_a)) * sqrt(rho_y * rho_m / (1-rho_m)) * x2 * x3 +
        sqrt(rho_m * rho_a / (1-rho_a)) * sum(T5 * x1) * x2 +
        sqrt(rho_y * rho_m / (1-rho_m)) * sum(T6 * x1) * x3
      f_bootstrap <- T3_bootstrap + T4_bootstrap * sqrt(rho_m * rho_a / (1-rho_a)) * sqrt(rho_y * rho_m / (1-rho_m)) * x2 * x3 +
        sqrt(rho_m * rho_a / (1-rho_a)) * (T5_bootstrap %*% x1) * x2 +
        sqrt(rho_y * rho_m / (1-rho_m)) * (T6_bootstrap %*% x1) * x3
      return((2*(T3>=0)-1) * f / sd(f_bootstrap))
    }
    x0 <- c(runif(dim_m-1, 0, 2*pi), runif(2, -1, 1))
    reT3 <- nloptr::nloptr(x0=x0,
                   eval_f = fun,
                   lb = c(rep(0, dim_m-1), -1, -1),
                   ub = c(rep(2*pi, dim_m-1), 1, 1),
                   opts = list(algorithm = algorithm,
                               xtol_rel= 1.0e-8,
                               maxeval = 100000))
  }else{
    fun <- function(x){
      x2 <- x[1]
      x3 <- x[2]
      f1 <- T3 + T4 * sqrt(rho_m * rho_a / (1-rho_a)) * sqrt(rho_y * rho_m / (1-rho_m)) * x2 * x3 +
        sqrt(rho_m * rho_a / (1-rho_a)) * T5 * x2 +
        sqrt(rho_y * rho_m / (1-rho_m)) * T6 * x3
      f2 <- T3 + T4 * sqrt(rho_m * rho_a / (1-rho_a)) * sqrt(rho_y * rho_m / (1-rho_m)) * x2 * x3 -
        sqrt(rho_m * rho_a / (1-rho_a)) * T5 * x2 -
        sqrt(rho_y * rho_m / (1-rho_m)) * T6 * x3
      f1_bootstrap <- T3_bootstrap + T4_bootstrap * sqrt(rho_m * rho_a / (1-rho_a)) * sqrt(rho_y * rho_m / (1-rho_m)) * x2 * x3 +
        sqrt(rho_m * rho_a / (1-rho_a)) * T5_bootstrap * x2 +
        sqrt(rho_y * rho_m / (1-rho_m)) * T6_bootstrap * x3
      f2_bootstrap <- T3_bootstrap + T4_bootstrap * sqrt(rho_m * rho_a / (1-rho_a)) * sqrt(rho_y * rho_m / (1-rho_m)) * x2 * x3 -
        sqrt(rho_m * rho_a / (1-rho_a)) * T5_bootstrap * x2 -
        sqrt(rho_y * rho_m / (1-rho_m)) * T6_bootstrap * x3
      return(min((2*(T3>=0)-1) * f1 / sd(f1_bootstrap), (2*(T3>=0)-1) * f2 / sd(f2_bootstrap)))
    }
    x0 <- c(runif(2, -1, 1))
    reT3 <- nloptr::nloptr(x0=x0,
                   eval_f = fun,
                   lb = c(-1, -1),
                   ub = c(1, 1),
                   opts = list(algorithm = algorithm,
                               xtol_rel= 1.0e-8,
                               maxeval = 100000))
  }

  return(reT3$objective)
}

#' @keywords internal
indirect_t_solve_vu <- function(T3, T4, T5, T6, T3_bootstrap, T4_bootstrap, T5_bootstrap, T6_bootstrap,
                                rho_y, rho_m, rho_a, algorithm){
  dim_m <- length(T5)

  if (dim_m > 1){
    fun <- function(x){
      x1 <- trans(x[1: (dim_m-1)]) * x[2*dim_m-1]
      x2 <- trans(x[dim_m: (2*dim_m-2)]) * x[2*dim_m]
      f <- T3 + T4 * sqrt(rho_m * rho_a / (1-rho_a)) * sqrt(rho_y * rho_m / (1-rho_m)) * sum(x1*x2) +
        sqrt(rho_m * rho_a / (1-rho_a)) * sum(T5 * x1) +
        sqrt(rho_y * rho_m / (1-rho_m)) * sum(T6 * x2)
      f_bootstrap <- T3_bootstrap + T4_bootstrap * sqrt(rho_m * rho_a / (1-rho_a)) * sqrt(rho_y * rho_m / (1-rho_m)) * sum(x1*x2) +
        sqrt(rho_m * rho_a / (1-rho_a)) * (T5_bootstrap %*% x1) +
        sqrt(rho_y * rho_m / (1-rho_m)) * (T6_bootstrap %*% x2)
      return((2*(T3>=0)-1) * f / sd(f_bootstrap))
    }
    x0 <- c(runif(2*dim_m-2, 0, 2*pi), runif(2, 0, 1))
    reT3 <- nloptr::nloptr(x0=x0,
                   eval_f = fun,
                   lb = c(rep(0, 2*dim_m-2), 0, 0),
                   ub = c(rep(2*pi, 2*dim_m-2), 1, 1),
                   opts = list(algorithm = algorithm,
                               xtol_rel= 1.0e-8,
                               maxeval = 100000))
  }else{
    fun <- function(x){
      x1 <- x[1]
      x2 <- x[2]
      f <- T3 + T4 * sqrt(rho_m * rho_a / (1-rho_a)) * sqrt(rho_y * rho_m / (1-rho_m)) * sum(x1*x2) +
        sqrt(rho_m * rho_a / (1-rho_a)) * sum(T5 * x1) +
        sqrt(rho_y * rho_m / (1-rho_m)) * sum(T6 * x2)
      f_bootstrap <- T3_bootstrap + T4_bootstrap * sqrt(rho_m * rho_a / (1-rho_a)) * sqrt(rho_y * rho_m / (1-rho_m)) * sum(x1*x2) +
        sqrt(rho_m * rho_a / (1-rho_a)) * T5_bootstrap * x1 +
        sqrt(rho_y * rho_m / (1-rho_m)) * T6_bootstrap * x2
      return((2*(T3>=0)-1) * f / sd(f_bootstrap))
    }
    x0 <- c(runif(2, -1, 1))
    reT3 <- nloptr::nloptr(x0 = x0,
                   eval_f = fun,
                   lb = c(-1, -1),
                   ub = c(1, 1),
                   opts = list(algorithm = algorithm,
                               xtol_rel= 1.0e-8,
                               maxeval = 100000))
  }

  return(reT3$objective)
}

#' Compute the robustness values for the indirect effect
#' @param y a vector of the outcome variable
#' @param m a matrix of the mediator variables
#' @param a a vector of the exposure variable
#' @param c a matrix of the covariates
#' @param randomized if the exposure is randomized
#' @param rho_values the values of rho's to be considered for the grid search
#' @param vector_u if the unmeasured confounder is a vector
#' @param algorithm the optimization algorithm in the R package "nloptr"
#' @return a summary table, the worst point estimates and the worst t statistics given each value of rho.
#' @examples
#' set.seed(1234)
#' library(MASS)
#' n = 200
#' a = rnorm(n)
#' c = as.matrix(rnorm(n))
#' m = mvrnorm(n, mu = c(0, 0), Sigma = matrix(c(1,0.2,0.2,0.8), 2, 2))
#' m = m + cbind(2.5 * a, 1.5 * a)
#' y = a + 0.2 * m[, 1] + 0.25 * m[, 2] + rnorm(n)
#' result = bku_rv_indirect(y, m, a, c)
#' @importFrom stats cor cov lm runif sd var
#' @importFrom MASS mvrnorm
#' @importFrom expm sqrtm
#' @importFrom nloptr nloptr
#' @export
bku_rv_indirect <- function(y, m, a, c = NULL, randomized = F, rho_values = c(1:99)/100, vector_u = F, algorithm = "NLOPT_GN_ORIG_DIRECT_L"){
  n <- length(y)
  m <- as.matrix(m)
  dim_m <- dim(m)[2]
  Ts <- indirect_compute_T(y, m, a, c)
  T3 <- Ts[1]
  T4 <- Ts[2]
  T5 <- Ts[3: (2+dim_m)]
  T6 <- Ts[(3+dim_m): (2+2*dim_m)]
  sample_index <- replicate(1000, sample(n, n, replace = T))
  Ts_bootstrap <- apply(sample_index, 2, function(x) indirect_compute_T(y[x], m[x, ], a[x], c[x, ]))
  T3_bootstrap <- Ts_bootstrap[1, ]
  T4_bootstrap <- Ts_bootstrap[2, ]
  T5_bootstrap <- t(Ts_bootstrap[3: (2+dim_m), ])
  T6_bootstrap <- t(Ts_bootstrap[(3+dim_m): (2+2*dim_m), ])
  point_rho <- matrix(nrow = length(rho_values), ncol = 4)
  t_rho <- matrix(nrow = length(rho_values), ncol = 4)
  for (i in c(1: length(rho_values))){
    point_rho[i, 1] <- rho_values[i]
    point_rho[i, 2] <- rho_values[i]
    point_rho[i, 3] <- rho_values[i] * (1-randomized)
    t_rho[i, 1] <- rho_values[i]
    t_rho[i, 2] <- rho_values[i]
    t_rho[i, 3] <- rho_values[i] * (1-randomized)
    if (vector_u){
      point_rho[i, 4] <- indirect_point_solve_vu(T3, T4, T5, T6, point_rho[i, 1], point_rho[i, 2], point_rho[i, 3], algorithm)
      t_rho[i, 4] <- indirect_t_solve_vu(T3, T4, T5, T6, T3_bootstrap, T4_bootstrap, T5_bootstrap, T6_bootstrap, t_rho[i, 1], t_rho[i, 2], t_rho[i, 3], algorithm)
    }else{
      point_rho[i, 4] <- indirect_point_solve(T3, T4, T5, T6, point_rho[i, 1], point_rho[i, 2], point_rho[i, 3], algorithm)
      t_rho[i, 4] <- indirect_t_solve(T3, T4, T5, T6, T3_bootstrap, T4_bootstrap, T5_bootstrap, T6_bootstrap, t_rho[i, 1], t_rho[i, 2], t_rho[i, 3], algorithm)
    }
  }
  rho_values <- c(0, rho_values)
  point_rv <- rho_values[min(which(point_rho[, 4] < 0))]
  t_rv <- rho_values[min(which(t_rho[, 4] < 1.96))]
  colnames(point_rho) <- c("rho_y", "rho_m", "rho_a", "worst point estimate")
  colnames(t_rho) <- c("rho_y", "rho_m", "rho_a", "worst t statistic")
  point_rho[, 4] <- point_rho[, 4] * (2*(T3>=0)-1)
  t_rho[, 4] <- t_rho[, 4] * (2*(T3>=0)-1)
  point_obs <- T3
  se_obs <- sd(T3_bootstrap)
  t_obs <- point_obs / se_obs
  summary_tab <- c(point_obs, se_obs, t_obs, point_rv, t_rv)
  names(summary_tab) <- c("Est.", "Std. Error", "t value", "R.V. for Est.", "R.V. for 95% C.I.")
  return(list("summary_table" = summary_tab, "worst_point" = point_rho, "worst_t" = t_rho))
}

#' Compute the sensitivity bound for the indirect effect with prespecified sensitivity parameters
#' @param y a vector of the outcome variable
#' @param m a matrix of the mediator variables
#' @param a a vector of the exposure variable
#' @param c a matrix of the covariates
#' @param Ry a number of the R parameter for the outcome-confounder correlation
#' @param Rm a vector of the R parameter for the mediator-confounder correlation
#' @param Ra a number of the R parameter for the exposure-confounder correlation
#' @return the point estimate, the standard error, and the t statistic
#' @examples
#' set.seed(1234)
#' library(MASS)
#' n = 200
#' a = rnorm(n)
#' c = as.matrix(rnorm(n))
#' m = mvrnorm(n, mu = c(0, 0), Sigma = matrix(c(1,0.2,0.2,0.8), 2, 2))
#' m = m + cbind(2.5 * a, 1.5 * a)
#' y = a + 0.2 * m[, 1] + 0.25 * m[, 2] + rnorm(n)
#' result = bku_indirect(y, m, a, c, 0.1, c(0.2, 0.1), 0.05)
#' @importFrom stats cor cov lm runif sd var
#' @importFrom MASS mvrnorm
#' @importFrom expm sqrtm
#' @export
bku_indirect <- function(y, m, a, c = NULL, Ry, Rm, Ra){
  n <- length(y)
  m <- as.matrix(m)
  dim_m <- dim(m)[2]
  Ts <- indirect_compute_T(y, m, a, c)
  T3 <- Ts[1]
  T4 <- Ts[2]
  T5 <- Ts[3: (2+dim_m)]
  T6 <- Ts[(3+dim_m): (2+2*dim_m)]
  sample_index <- replicate(1000, sample(n, n, replace = T))
  Ts_bootstrap <- apply(sample_index, 2, function(x) indirect_compute_T(y[x], m[x, ], a[x], c[x, ]))
  T3_bootstrap <- Ts_bootstrap[1, ]
  T4_bootstrap <- Ts_bootstrap[2, ]
  T5_bootstrap <- t(Ts_bootstrap[3: (2+dim_m), ])
  T6_bootstrap <- t(Ts_bootstrap[(3+dim_m): (2+2*dim_m), ])
  phi3 <- Ra / sqrt(1 - Ra^2) * Rm
  phi4 <- Ry / sqrt(1 - sum(Rm^2)) * Rm
  point_est <- T3 + T4 * sum(phi3*phi4) + sum(T5 * phi3) + sum(T6 * phi4)
  if (dim_m > 1){
    std_err <- sd(T3_bootstrap + T4_bootstrap * sum(phi3*phi4) + T5_bootstrap %*% phi3 + T6_bootstrap %*% phi4)
  }else{
    std_err <- sd(T3_bootstrap + T4_bootstrap * sum(phi3*phi4) + T5_bootstrap * phi3 + T6_bootstrap * phi4)
  }
  t_stat <- point_est / std_err
  return(c("point_est" = point_est, "std_err" = std_err, "t_stat" = t_stat))

}

#' Compute the robustness values for both direct and indirect effects
#' @param y a vector of the outcome variable
#' @param m a matrix of the mediator variables
#' @param a a vector of the exposure variable
#' @param c a matrix of the covariates
#' @param randomized if the exposure is randomized
#' @param rho_values the values of rho's to be considered for the grid search
#' @param vector_u if the unmeasured confounder is a vector
#' @param algorithm the optimization algorithm in the R package "nloptr"
#' @return a summary table, the worst point estimates and the worst t statistics given each value of rho.
#' @examples
#' set.seed(1234)
#' library(MASS)
#' n = 200
#' a = rnorm(n)
#' c = as.matrix(rnorm(n))
#' m = mvrnorm(n, mu = c(0, 0), Sigma = matrix(c(1,0.2,0.2,0.8), 2, 2))
#' m = m + cbind(2.5 * a, 1.5 * a)
#' y = a + 0.2 * m[, 1] + 0.25 * m[, 2] + rnorm(n)
#' result = bku_rv_indirect(y, m, a, c)
#' @importFrom stats cor cov lm runif sd var
#' @importFrom MASS mvrnorm
#' @importFrom expm sqrtm
#' @importFrom nloptr nloptr
#' @export
bku_rv <- function(y, m, a, c = NULL, randomized = F, rho_values = c(1:99)/100, vector_u = F, algorithm = "NLOPT_GN_ORIG_DIRECT_L"){
  n <- length(y)
  m <- as.matrix(m)
  dim_m <- dim(m)[2]
  Ts1 <- direct_compute_T(y, m, a, c)
  T0 <- Ts1[1]
  T1 <- Ts1[2]
  T2 <- Ts1[-c(1,2)]

  Ts2 <- indirect_compute_T(y, m, a, c)
  T3 <- Ts2[1]
  T4 <- Ts2[2]
  T5 <- Ts2[3: (2+dim_m)]
  T6 <- Ts2[(3+dim_m): (2+2*dim_m)]

  sample_index <- replicate(1000, sample(n, n, replace = T))
  Ts1_bootstrap <- apply(sample_index, 2, function(x) direct_compute_T(y[x], m[x, ], a[x], c[x, ]))
  T0_bootstrap <- Ts1_bootstrap[1, ]
  T1_bootstrap <- Ts1_bootstrap[2, ]
  T2_bootstrap <- t(Ts1_bootstrap[-c(1,2), ])

  Ts2_bootstrap <- apply(sample_index, 2, function(x) indirect_compute_T(y[x], m[x, ], a[x], c[x, ]))
  T3_bootstrap <- Ts2_bootstrap[1, ]
  T4_bootstrap <- Ts2_bootstrap[2, ]
  T5_bootstrap <- t(Ts2_bootstrap[3: (2+dim_m), ])
  T6_bootstrap <- t(Ts2_bootstrap[(3+dim_m): (2+2*dim_m), ])

  direct_point_rho <- matrix(nrow = length(rho_values), ncol = 4)
  direct_t_rho <- matrix(nrow = length(rho_values), ncol = 4)
  indirect_point_rho <- matrix(nrow = length(rho_values), ncol = 4)
  indirect_t_rho <- matrix(nrow = length(rho_values), ncol = 4)

  for (i in c(1: length(rho_values))){
    direct_point_rho[i, 1] <- rho_values[i]
    direct_point_rho[i, 2] <- rho_values[i]
    direct_point_rho[i, 3] <- rho_values[i] * (1-randomized)
    direct_t_rho[i, 1] <- rho_values[i]
    direct_t_rho[i, 2] <- rho_values[i]
    direct_t_rho[i, 3] <- rho_values[i] * (1-randomized)
    direct_point_rho[i, 4] <- direct_point_solve(T0, T1, T2, direct_point_rho[i, 1], direct_point_rho[i, 2], direct_point_rho[i, 3], algorithm)
    direct_t_rho[i, 4] <- direct_t_solve(T0, T1, T2, T0_bootstrap, T1_bootstrap, T2_bootstrap, direct_t_rho[i, 1], direct_t_rho[i, 2], direct_t_rho[i, 3], algorithm)
  }

  for (i in c(1: length(rho_values))){
    indirect_point_rho[i, 1] <- rho_values[i]
    indirect_point_rho[i, 2] <- rho_values[i]
    indirect_point_rho[i, 3] <- rho_values[i] * (1-randomized)
    indirect_t_rho[i, 1] <- rho_values[i]
    indirect_t_rho[i, 2] <- rho_values[i]
    indirect_t_rho[i, 3] <- rho_values[i] * (1-randomized)
    if (vector_u){
      indirect_point_rho[i, 4] <- indirect_point_solve_vu(T3, T4, T5, T6, indirect_point_rho[i, 1], indirect_point_rho[i, 2], indirect_point_rho[i, 3], algorithm)
      indirect_t_rho[i, 4] <- indirect_t_solve_vu(T3, T4, T5, T6, T3_bootstrap, T4_bootstrap, T5_bootstrap, T6_bootstrap, indirect_t_rho[i, 1], indirect_t_rho[i, 2], indirect_t_rho[i, 3], algorithm)
    }else{
      indirect_point_rho[i, 4] <- indirect_point_solve(T3, T4, T5, T6, indirect_point_rho[i, 1], indirect_point_rho[i, 2], indirect_point_rho[i, 3], algorithm)
      indirect_t_rho[i, 4] <- indirect_t_solve(T3, T4, T5, T6, T3_bootstrap, T4_bootstrap, T5_bootstrap, T6_bootstrap, indirect_t_rho[i, 1], indirect_t_rho[i, 2], indirect_t_rho[i, 3], algorithm)
    }
  }
  rho_values <- c(0, rho_values)

  direct_point_rv <- rho_values[min(which(direct_point_rho[, 4] < 0))]
  direct_t_rv <- rho_values[min(which(direct_t_rho[, 4] < 1.96))]
  indirect_point_rv <- rho_values[min(which(indirect_point_rho[, 4] < 0))]
  indirect_t_rv <- rho_values[min(which(indirect_t_rho[, 4] < 1.96))]

  colnames(direct_point_rho) <- c("rho_y", "rho_m", "rho_a", "worst point estimate")
  colnames(direct_t_rho) <- c("rho_y", "rho_m", "rho_a", "worst t statistic")

  direct_point_rho[, 4] <- direct_point_rho[, 4] * (2*(T0>=0)-1)
  direct_t_rho[, 4] <- direct_t_rho[, 4] * (2*(T0>=0)-1)
  indirect_point_rho[, 4] <- indirect_point_rho[, 4] * (2*(T3>=0)-1)
  indirect_t_rho[, 4] <- indirect_t_rho[, 4] * (2*(T3>=0)-1)

  direct_point_obs <- T0
  direct_se_obs <- sd(T0_bootstrap)
  direct_t_obs <- direct_point_obs / direct_se_obs
  indirect_point_obs <- T3
  indirect_se_obs <- sd(T3_bootstrap)
  indirect_t_obs <- indirect_point_obs / indirect_se_obs
  summary_tab <- rbind(c(direct_point_obs, direct_se_obs, direct_t_obs, direct_point_rv, direct_t_rv),
                       c(indirect_point_obs, indirect_se_obs, indirect_t_obs, indirect_point_rv, indirect_t_rv))
  colnames(summary_tab) <- c("Est.", "Std. Error", "t value", "R.V. for Est.", "R.V. for 95% C.I.")
  rownames(summary_tab) <- c("Direct effect", "Indirect effect")
  return(list("summary_table" = summary_tab, "direct_worst_point" = direct_point_rho, "direct_worst_t" = direct_t_rho,
              "indirect_worst_point" = indirect_point_rho, "indirect_worst_t" = indirect_t_rho))
}

#' @keywords internal
fb_compute_T <- function(y, m, a, cj, cmj){
  dim_m <- dim(as.matrix(m))[2]

  T1 <- cor(lm(a~cmj)$residuals, lm(cj~cmj)$residuals)
  T2 <- solve(expm::sqrtm(cov(as.matrix(lm(m~a+cj+cmj)$residuals)))) %*%
    expm::sqrtm(cov(as.matrix(lm(m~a+cmj)$residuals)))
  T3 <- solve(expm::sqrtm(cov(as.matrix(lm(m~a+cmj)$residuals)))) %*%
    cov(as.matrix(lm(m~a+cmj)$residuals), lm(cj~a+cmj)$residuals) %*%
    (1/sqrt(var(lm(cj~a+cmj)$residuals)))
  T4 <- cor(lm(y~a+m+cmj)$residuals, lm(cj~a+m+cmj)$residuals)

  T5 <- lm(y~a+m+cj+cmj)$coefficients[2]
  T6 <- sqrt(var(lm(y~a+m+cj+cmj)$residuals)) / sqrt(var(lm(a~m+cj+cmj)$residuals))
  T7 <- solve(expm::sqrtm(cov(as.matrix(lm(m~cj+cmj)$residuals)))) %*%
    cov(as.matrix(lm(m~cj+cmj)$residuals), lm(a~cj+cmj)$residuals) %*%
    (1/sqrt(var(lm(a~cj+cmj)$residuals)))
  T8 <- solve(expm::sqrtm(cov(as.matrix(lm(m~cj+cmj)$residuals)))) %*%
    expm::sqrtm(cov(as.matrix(lm(m~a+cj+cmj)$residuals)))
  T9 <- sqrt(var(lm(a~m+cj+cmj)$residuals) / var(lm(a~cj+cmj)$residuals))

  T10 <- as.matrix(lm(m~a+cj+cmj)$coefficients)[2, ]
  T11 <- lm(y~m+a+cj+cmj)$coefficients[2: (dim_m+1)]
  T12 <- expm::sqrtm(cov(as.matrix(lm(m~a+cj+cmj)$residuals))) / sqrt(var(lm(a~cj+cmj)$residuals))
  T13 <- solve(expm::sqrtm(cov(as.matrix(lm(m~a+cj+cmj)$residuals)))) * sqrt(var(lm(y~a+m+cj+cmj)$residuals))

  return(list(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13))
}

#' @keywords internal
fb_compute_T_resample <- function(y, m, a, cj, cmj){
  n <- length(y)
  sample_index <- sample(n, n, replace = T)
  y <- y[sample_index]
  m <- as.matrix(m)[sample_index, ]
  a <- a[sample_index]
  cj <- cj[sample_index]
  cmj <- as.matrix(cmj)[sample_index, ]
  return(fb_compute_T(y, m, a, cj, cmj))
}

#' @keywords internal
fb_point <- function(T_point, ry, rm, ra, rc){
  Ra <- (ra - T_point[[1]] * rc) / sqrt(1 - T_point[[1]]^2) / sqrt(1 - rc^2)
  R1 <- (rc - T_point[[1]] * ra) / sqrt(1 - ra^2) / sqrt(1 - T_point[[1]]^2)
  Rm <- T_point[[2]] %*% (rm - T_point[[3]] * R1) / sqrt(1 - R1^2)
  R2 <- (rc - sum(T_point[[3]] * rm)) / sqrt(1 - sum(T_point[[3]]^2)) / sqrt(1 - sum(rm^2))
  Ry <- (ry - T_point[[4]] * R2) / sqrt(1 - T_point[[4]]^2) / sqrt(1 - R2^2)

  point_direct <- T_point[[5]] - Ry * T_point[[6]] *
    (Ra - t(T_point[[7]]) %*% (sqrt(1-Ra^2) * T_point[[8]] %*% Rm + Ra * T_point[[7]])) /
    sqrt(1-Ra^2) / T_point[[9]] / sqrt(1 - sum(Rm^2))
  point_indirect <- t(T_point[[10]] - Ra / sqrt(1-Ra^2) * T_point[[12]] %*% Rm) %*%
    (T_point[[11]] - Ry / sqrt(1-sum(Rm^2)) * T_point[[13]] %*% Rm)
  return(list(point_direct, point_indirect))
}

#' @keywords internal
fb_bootstrap <- function(T_bootstrap, ry, rm, ra, rc){
  Ra <- sapply(1: 1000, function(x)
    (ra - T_bootstrap[1, ][[x]] * rc) / sqrt(1 - T_bootstrap[1, ][[x]]^2) / sqrt(1 - rc^2),
    simplify = F)

  R1 <- sapply(1: 1000, function(x)
    (rc - T_bootstrap[1, ][[x]] * ra) / sqrt(1 - ra^2) / sqrt(1 - T_bootstrap[1, ][[x]]^2),
    simplify = F)

  Rm <- sapply(1: 1000, function(x)
    T_bootstrap[2, ][[x]] %*% (rm - T_bootstrap[3, ][[x]] * R1[[x]]) / sqrt(1 - R1[[x]]^2),
    simplify = F)

  R2 <- sapply(1: 1000, function(x)
    (rc - sum(T_bootstrap[3, ][[x]] * rm)) / sqrt(1 - sum(T_bootstrap[3, ][[x]]^2)) / sqrt(1 - sum(rm^2)),
    simplify = F)

  Ry <- sapply(1: 1000, function(x)
    (ry - T_bootstrap[4, ][[x]] * R2[[x]]) / sqrt(1 - T_bootstrap[4, ][[x]]^2) / sqrt(1 - R2[[x]]^2),
    simplify = F)

  direct_bootstrap <- sapply(1: 1000, function(x)
    T_bootstrap[5, ][[x]] - Ry[[x]] * T_bootstrap[6, ][[x]] *
      (Ra[[x]] - t(T_bootstrap[7, ][[x]]) %*% (sqrt(1-Ra[[x]]^2) * T_bootstrap[8, ][[x]] %*% Rm[[x]] + Ra[[x]] * T_bootstrap[7, ][[x]])) /
      sqrt(1-Ra[[x]]^2) / T_bootstrap[9, ][[x]] / sqrt(1 - sum(Rm[[x]]^2)))

  indirect_bootstrap <- sapply(1: 1000, function(x)
    t(T_bootstrap[10, ][[x]] - Ra[[x]] / sqrt(1-Ra[[x]]^2) * T_bootstrap[12, ][[x]] %*% Rm[[x]]) %*%
      (T_bootstrap[11, ][[x]] - Ry[[x]] / sqrt(1-sum(Rm[[x]]^2)) * T_bootstrap[13, ][[x]] %*% Rm[[x]]))

  return(list(direct_bootstrap, indirect_bootstrap))
}

#' @keywords internal
fb_direct_point_solve <- function(T_point, rho_y, rho_m, rho_a, algorithm = "NLOPT_GN_ORIG_DIRECT_L"){
  dim_m <- dim(T_point[[2]])[1]

  if (dim_m > 1){
    fun <- function(x){
      return((2*(T_point[[5]]>=0)-1)*fb_point(T_point, x[1], x[2] * trans(x[4: (dim_m+2)]), x[3], 0)[[1]])
    }
    x0 <- c(runif(1, -sqrt(rho_y), sqrt(rho_y)),
            runif(1, -sqrt(rho_m), sqrt(rho_m)),
            runif(1, -sqrt(rho_a), sqrt(rho_a)),
            runif(dim_m-1, 0, 2*pi))
    res0 <- nloptr::nloptr(x0=x0,
                   eval_f = fun,
                   lb = c(-sqrt(rho_y), -sqrt(rho_m), -sqrt(rho_a), rep(0, dim_m-1)),
                   ub = c(sqrt(rho_y), sqrt(rho_m), sqrt(rho_a), rep(2*pi, dim_m-1)),
                   opts = list(algorithm = algorithm,
                               xtol_rel= 1.0e-8,
                               maxeval = 100000))
  }else{
    fun <- function(x){
      return((2*(T_point[[5]]>=0)-1)*fb_point(T_point, x[1], x[2], x[3], 0)[[1]])
    }
    x0 <- c(runif(1, -sqrt(rho_y), sqrt(rho_y)),
            runif(1, -sqrt(rho_m), sqrt(rho_m)),
            runif(1, -sqrt(rho_a), sqrt(rho_a)))
    res0 <- nloptr::nloptr(x0=x0,
                   eval_f = fun,
                   lb = c(-sqrt(rho_y), -sqrt(rho_m), -sqrt(rho_a)),
                   ub = c(sqrt(rho_y), sqrt(rho_m), sqrt(rho_a)),
                   opts = list(algorithm = algorithm,
                               xtol_rel= 1.0e-8,
                               maxeval = 100000))
  }
  return(res0$objective)
}

#' @keywords internal
fb_indirect_point_solve <- function(T_point, rho_y, rho_m, rho_a, algorithm = "NLOPT_GN_ORIG_DIRECT_L"){
  dim_m <- dim(T_point[[2]])[1]

  if (dim_m > 1){
    fun <- function(x){
      return((2*(sum(T_point[[10]]*T_point[[11]])>=0)-1)*fb_point(T_point, x[1], x[2] * trans(x[4: (dim_m+2)]), x[3], 0)[[2]])
    }
    x0 <- c(runif(1, -sqrt(rho_y), sqrt(rho_y)),
            runif(1, -sqrt(rho_m), sqrt(rho_m)),
            runif(1, -sqrt(rho_a), sqrt(rho_a)),
            runif(dim_m-1, 0, 2*pi))
    res0 <- nloptr::nloptr(x0=x0,
                   eval_f = fun,
                   lb = c(-sqrt(rho_y), -sqrt(rho_m), -sqrt(rho_a), rep(0, dim_m-1)),
                   ub = c(sqrt(rho_y), sqrt(rho_m), sqrt(rho_a), rep(2*pi, dim_m-1)),
                   opts = list(algorithm = algorithm,
                               xtol_rel= 1.0e-8,
                               maxeval = 100000))
  }else{
    fun <- function(x){
      return((2*(sum(T_point[[10]]*T_point[[11]])>=0)-1)*fb_point(T_point, x[1], x[2], x[3], 0)[[2]])
    }
    x0 <- c(runif(1, -sqrt(rho_y), sqrt(rho_y)),
            runif(1, -sqrt(rho_m), sqrt(rho_m)),
            runif(1, -sqrt(rho_a), sqrt(rho_a)))
    res0 <- nloptr::nloptr(x0=x0,
                   eval_f = fun,
                   lb = c(-sqrt(rho_y), -sqrt(rho_m), -sqrt(rho_a)),
                   ub = c(sqrt(rho_y), sqrt(rho_m), sqrt(rho_a)),
                   opts = list(algorithm = algorithm,
                               xtol_rel= 1.0e-8,
                               maxeval = 100000))
  }
  return(res0$objective)
}

#' @keywords internal
fb_direct_t_solve <- function(T_point, T_bootstrap, rho_y, rho_m, rho_a, algorithm = "NLOPT_GN_ORIG_DIRECT_L"){
  dim_m <- dim(T_point[[2]])[1]

  if (dim_m > 1){
    fun <- function(x){
      return((2*(T_point[[5]]>=0)-1)*fb_point(T_point, x[1], x[2] * trans(x[4: (dim_m+2)]), x[3], 0)[[1]] /
               sd(fb_bootstrap(T_bootstrap, x[1], x[2] * trans(x[4: (dim_m+2)]), x[3], 0)[[1]]))
    }
    x0 <- c(runif(1, -sqrt(rho_y), sqrt(rho_y)),
            runif(1, -sqrt(rho_m), sqrt(rho_m)),
            runif(1, -sqrt(rho_a), sqrt(rho_a)),
            runif(dim_m-1, 0, 2*pi))
    res0 <- nloptr::nloptr(x0=x0,
                   eval_f = fun,
                   lb = c(-sqrt(rho_y), -sqrt(rho_m), -sqrt(rho_a), rep(0, dim_m-1)),
                   ub = c(sqrt(rho_y), sqrt(rho_m), sqrt(rho_a), rep(2*pi, dim_m-1)),
                   opts = list(algorithm = algorithm,
                               xtol_rel= 1.0e-8,
                               maxeval = 100000))
  }else{
    fun <- function(x){
      return((2*(T_point[[5]]>=0)-1)*fb_point(T_point, x[1], x[2], x[3], 0)[[1]] /
               sd(fb_bootstrap(T_bootstrap, x[1], x[2], x[3], 0)[[1]]))
    }
    x0 <- c(runif(1, -sqrt(rho_y), sqrt(rho_y)),
            runif(1, -sqrt(rho_m), sqrt(rho_m)),
            runif(1, -sqrt(rho_a), sqrt(rho_a)))
    res0 <- nloptr::nloptr(x0=x0,
                   eval_f = fun,
                   lb = c(-sqrt(rho_y), -sqrt(rho_m), -sqrt(rho_a)),
                   ub = c(sqrt(rho_y), sqrt(rho_m), sqrt(rho_a)),
                   opts = list(algorithm = algorithm,
                               xtol_rel= 1.0e-8,
                               maxeval = 100000))
  }
  return(res0$objective)
}

#' @keywords internal
fb_indirect_t_solve <- function(T_point, T_bootstrap, rho_y, rho_m, rho_a, algorithm = "NLOPT_GN_ORIG_DIRECT_L"){
  dim_m <- dim(T_point[[2]])[1]

  if (dim_m > 1){
    fun <- function(x){
      return((2*(sum(T_point[[10]]*T_point[[11]])>=0)-1)*fb_point(T_point, x[1], x[2] * trans(x[4: (dim_m+2)]), x[3], 0)[[2]] /
               sd(fb_bootstrap(T_bootstrap, x[1], x[2] * trans(x[4: (dim_m+2)]), x[3], 0)[[2]]))
    }
    x0 <- c(runif(1, -sqrt(rho_y), sqrt(rho_y)),
            runif(1, -sqrt(rho_m), sqrt(rho_m)),
            runif(1, -sqrt(rho_a), sqrt(rho_a)),
            runif(dim_m-1, 0, 2*pi))
    res0 <- nloptr::nloptr(x0=x0,
                   eval_f = fun,
                   lb = c(-sqrt(rho_y), -sqrt(rho_m), -sqrt(rho_a), rep(0, dim_m-1)),
                   ub = c(sqrt(rho_y), sqrt(rho_m), sqrt(rho_a), rep(2*pi, dim_m-1)),
                   opts = list(algorithm = algorithm,
                               xtol_rel= 1.0e-8,
                               maxeval = 100000))
  }else{
    fun <- function(x){
      return((2*(sum(T_point[[10]]*T_point[[11]])>=0)-1)*fb_point(T_point, x[1], x[2], x[3], 0)[[2]] /
               sd(fb_bootstrap(T_bootstrap, x[1], x[2], x[3], 0)[[2]]))
    }
    x0 <- c(runif(1, -sqrt(rho_y), sqrt(rho_y)),
            runif(1, -sqrt(rho_m), sqrt(rho_m)),
            runif(1, -sqrt(rho_a), sqrt(rho_a)))
    res0 <- nloptr::nloptr(x0=x0,
                   eval_f = fun,
                   lb = c(-sqrt(rho_y), -sqrt(rho_m), -sqrt(rho_a)),
                   ub = c(sqrt(rho_y), sqrt(rho_m), sqrt(rho_a)),
                   opts = list(algorithm = algorithm,
                               xtol_rel= 1.0e-8,
                               maxeval = 100000))
  }
  return(res0$objective)
}


#' Compute the sensitivity bound with formal benchmarking
#' @param y a vector of the outcome variable
#' @param m a matrix of the mediator variables
#' @param a a vector of the exposure variable
#' @param c a matrix of the covariates
#' @param j a number of index to specify the refernece covaraite cj
#' @param ky_bound a number of the upper bound for the relative strength of the outcome-confounder correlation
#' @param km_bound a number of the upper bound for the relative strength of the mediator-confounder correlation
#' @param ka_bound a number of the upper bound for the relative strength of the exposure-confounder correlation
#' @return the worst point estimate for the direct effect, the worst point estimate for the indirect effect, the worst t statistic for the direct effect, and the worst t statistic for the indirect effect
#' @examples
#' set.seed(1234)
#' library(MASS)
#' n = 200
#' a = rnorm(n)
#' c = cbind(rnorm(n), rnorm(n))
#' m = mvrnorm(n, mu = c(0, 0), Sigma = matrix(c(1,0.2,0.2,0.8), 2, 2))
#' m = m + cbind(2.5 * a, 1.5 * a)
#' y = a + 0.2 * m[, 1] + 0.25 * m[, 2] + rnorm(n)
#' result = bku_fb(y, m, a, c, 1, 0.2, 0.2, 0.2)
#' @importFrom stats cor cov lm runif sd var
#' @importFrom MASS mvrnorm
#' @importFrom expm sqrtm
#' @importFrom nloptr nloptr
#' @export
bku_fb <- function(y, m, a, c, j, ky_bound, km_bound, ka_bound){
  cj <- c[, j]
  cmj <- c[, -j]
  rho_a <- (1 - var(lm(a~cj+cmj)$residuals) / var(lm(a~cmj)$residuals)) * ka_bound
  rho_m <- (1 - var(lm(cj~a+m+cmj)$residuals) / var(lm(cj~a+cmj)$residuals)) * km_bound
  rho_y <- (1 - var(lm(y~cj+a+m+cmj)$residuals) / var(lm(y~a+m+cmj)$residuals)) * ky_bound
  T_point <- fb_compute_T(y, m, a, cj, cmj)
  direct_point <- fb_direct_point_solve(T_point, rho_y, rho_m, rho_a)
  indirect_point <- fb_indirect_point_solve(T_point, rho_y, rho_m, rho_a)
  T_bootstrap <- replicate(1000, fb_compute_T_resample(y, m, a, cj, cmj))
  direct_t <- fb_direct_t_solve(T_point, T_bootstrap, rho_y, rho_m, rho_a)
  indirect_t <- fb_indirect_t_solve(T_point, T_bootstrap, rho_y, rho_m, rho_a)
  return(c(direct_point, indirect_point, direct_t, indirect_t))
}
