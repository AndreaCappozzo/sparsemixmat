

source("../example_inverse_sparse_MBC_MICHAEL.R")
source("R/utils.R")
source("R/initialize.R")
source("R/mstep_inverse_with_M_penalized.R")

library(Rcpp)
library(RcppArmadillo)
sourceCpp("src/w_array.cpp")
sourceCpp("src/mstep_obj.cpp")
sourceCpp("src/estep_calc.cpp")
sourceCpp("src/penalization_M_mat.cpp")

# FIXME atm the code seems to work: it performs penalized MBC of three-way data, where both the precision matrices omega and gamma and the mean matrices are penalized.
# Next things to do and to improve:
# create a wrapper function for model fitting with multiple penalties and K

data <- scale_matrix_data(X)
data_dim <- dim(data)
p <- data_dim[1]
q <- data_dim[2]
N <- data_dim[3]
penalty_sigma = 0.1 # FIXME we do not need this
penalty_theta = 0.1 # FIXME we do not need this
penalty_omega = 0.1 # 0.1
penalty_gamma = 0.1 # 0.2
# penalty_mu = 100
penalty_mu = 0
penalize_diag = rep(FALSE, 2)
control = EM_controls()
gamma_init <- array(1, dim = c(q, q, K))
gamma_tmp <- gamma_init


em_mix_mat <- function(data,
                       K = 2,
                       control = EM_controls()){
  call <- match.call()
  data_dim <- dim(data)
  p <- data_dim[1]
  q <- data_dim[2]
  N <- data_dim[3]
  dims <- list(p = p, q = q, N = N, K = K)

  # set penalization matrices
  if (!is.matrix(penalty_sigma)) penalty_sigma <- matrix(penalty_sigma, nrow = p, ncol = p)
  if (penalize_diag[1] == FALSE) diag(penalty_sigma) <- 0
  #
  if (!is.matrix(penalty_theta)) penalty_theta <- matrix(penalty_theta, nrow = q, ncol = q)
  if (penalize_diag[2] == FALSE) diag(penalty_theta) <- 0
  #
  if (!is.matrix(penalty_omega)) penalty_omega <- matrix(penalty_omega, nrow = p, ncol = p)
  if (penalize_diag[1] == FALSE) diag(penalty_omega) <- 0
  #
  if (!is.matrix(penalty_gamma)) penalty_gamma <- matrix(penalty_gamma, nrow = q, ncol = q)
  if (penalize_diag[2] == FALSE) diag(penalty_gamma) <- 0

  # mu
  if (!is.matrix(penalty_mu)) penalty_mu <- matrix(penalty_mu, nrow = p, ncol = q)

  # store EM parameters
  tol <- control$tol
  max_iter <- control$max_iter
  type_start = control$type_start
  n_subset_start = control$n_subset_start
  n_random_start = control$n_random_start

  # initialization of z -----------------------------------------------------------------
  # TODO: remove it and include in main wrapper for model selection
  hc_init <- if ( type_start == "hc" ) {
    mclust::hcVVV( data = matrix(data, N, p*q, byrow = TRUE) )
  } else NULL

  sigma <- omega <- array(diag(p), dim = c(p, p, K))   # start with identity matrix
  theta <- gamma <- array(diag(q), dim = c(q, q, K))
  # init <- initialize(data, type_start, hc_init, gamma, dims)
  init <- initialize(data, type_start = "random", hc_init, sigma, omega, theta, gamma, dims)
  z <- init$z
  tau <- init$parameters$tau
  mu <- init$parameters$mu
  sigma <- init$parameters$sigma
  theta <- init$parameters$theta
  omega <- init$parameters$omega
  gamma <- init$parameters$gamma
  data_cent <- init$data_cent
  det_sigma <- init$det_sigma
  det_theta <- init$det_theta

  # EM ----------------------------------------------------------------------------------
  crit <- TRUE
  iter <- 0
  loglik <- loglik_prev <- loglik_pen <- loglik_pen_prev <- -.Machine$double.xmax
  # sanity check
  LLK <- c()

  while ( crit ) {

    iter <- iter + 1

    # M step -----------
    # out_mstep <-
    #   mstep_inverse_sparse_M(
    #     data = data,
    #     z = z,
    #     penalty_omega = penalty_omega,
    #     penalty_gamma = penalty_gamma,
    #     penalty_mu = penalty_mu,
    #     gamma = gamma,
    #     omega = omega,
    #     mu = mu,
    #     control = control,
    #     dims = dims
    #   )

    out_mstep <-
      mstep_inverse(
        data = data,
        z = z,
        penalty_omega = penalty_omega,
        penalty_gamma = penalty_gamma,
        sigma = sigma,
        omega = omega,
        theta = theta,
        gamma = gamma,
        control = control,
        dims = dims
      )

    # TODO: remove
    out_mstep <- list(parameters = list(tau = tau,
                                        mu = mu,
                                        sigma = sigma,
                                        theta = theta,
                                        omega = omega,
                                        gamma = gamma),
                      data_cent = data_cent,
                      det_sigma = det_sigma,
                      det_theta = det_theta)


    tau <- out_mstep$parameters$tau
    mu <- out_mstep$parameters$mu
    sigma <- out_mstep$parameters$sigma
    theta <- out_mstep$parameters$theta
    omega <- out_mstep$parameters$omega
    gamma <- out_mstep$parameters$gamma
    data_cent <- out_mstep$data_cent
    det_sigma <- out_mstep$det_sigma
    det_theta <- out_mstep$det_theta

    # E step -----------
    out_estep <- estep_calc(data_cent, z, mu, tau,
                            sigma, theta, omega, gamma,
                            det_sigma, det_theta)
    z <- out_estep$z


   # Loglik -----------

    loglik <- out_estep$loglik

    pen_value_omega <- sum( sweep(abs(omega), c(1,2), penalty_omega, "*") )
    pen_value_gamma <- sum( sweep(abs(gamma), c(1,2), penalty_gamma, "*") )
    pen_value_mu <- sum( sweep(abs(mu), c(1,2), penalty_mu, "*") )

    loglik_pen <- loglik - (pen_value_omega + pen_value_gamma + pen_value_mu)
    err <- abs(loglik_pen - loglik_pen_prev) / (1 + abs(loglik_pen))
    loglik_pen_prev <- loglik_pen
    LLK <- c(LLK, loglik_pen)

    crit <- ( err > tol[1] & iter < max_iter[1] )

  }

  out_mstep <- mstep_inverse_sparse_M(
    data = data,
    z = z,
    penalty_omega = penalty_omega,
    penalty_gamma = penalty_gamma,
    penalty_mu = penalty_mu,
    gamma = gamma,
    omega = omega,
    mu = mu,
    control = control,
    dims = dims
  )

  tau <- out_mstep$parameters$tau
  mu <- out_mstep$parameters$mu
  sigma <- out_mstep$parameters$sigma
  theta <- out_mstep$parameters$theta
  omega <- out_mstep$parameters$omega
  gamma <- out_mstep$parameters$gamma
  data_cent <- out_mstep$data_cent
  det_sigma <- out_mstep$det_sigma
  det_theta <- out_mstep$det_theta

}
