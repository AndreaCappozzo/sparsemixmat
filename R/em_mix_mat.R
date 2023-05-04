#' @export
em_mix_mat <- function(data,
                       K = 2,
                       penalty_omega,
                       penalty_gamma,
                       penalty_mu,
                       penalize_diag,
                       hc_init,
                       data_dim,
                       penalty_factor,
                       type_penalty_mu,
                       penalization_M_mat_coord_ascent,
                       control = EM_controls()) {
  
  call <- match.call()
  p <- data_dim[1]
  q <- data_dim[2]
  N <- data_dim[3]
  
  dims <- list(p = p, q = q, N = N, K = K)

  if (!is.matrix(penalty_omega)) penalty_omega <- matrix(penalty_omega, nrow = p, ncol = p)
  if (penalize_diag[1] == FALSE) diag(penalty_omega) <- 0
  #
  if (!is.matrix(penalty_gamma)) penalty_gamma <- matrix(penalty_gamma, nrow = q, ncol = q)
  if (penalize_diag[2] == FALSE) diag(penalty_gamma) <- 0

  # mu (varies depending on lasso or group-lasso)
  if (!is.matrix(penalty_mu) & type_penalty_mu == "lasso"){
    penalty_mu <- matrix(penalty_mu, nrow = p, ncol = q)
  } else if (!is.matrix(penalty_mu) & type_penalty_mu == "group-lasso") {
    penalty_mu <- rep(penalty_mu,p) *penalty_factor
  }

  # store EM parameters
  tol <- control$tol
  max_iter <- control$max_iter
  type_start <-  control$type_start
  n_subset_start <-  control$n_subset_start
  n_random_start <-  control$n_random_start

  # initialization of z -----------------------------------------------------------------
  
  omega <- array(diag(p), dim = c(p, p, K))   # start with identity matrix
  gamma <- array(diag(q), dim = c(q, q, K))
  
  init <-
    initialize(
      data = data,
      type_start = type_start,
      hc_init = hc_init,
      omega = omega,
      gamma = gamma,
      dims = dims,
      control=control,
      penalty_omega=penalty_omega,
      penalty_gamma=penalty_gamma,
      penalty_mu=penalty_mu,
      penalization_M_mat_coord_ascent=penalization_M_mat_coord_ascent
    )
  
  z <- init$z
  tau <- init$parameters$tau
  mu <- init$parameters$mu
  sigma <- init$parameters$sigma
  psi <- init$parameters$psi
  omega <- init$parameters$omega
  gamma <- init$parameters$gamma
  data_cent <- init$data_cent
  det_sigma <- init$det_sigma
  det_psi <- init$det_psi

  # ME ----------------------------------------------------------------------------------
  crit <- TRUE
  iter <- 0
  loglik <- loglik_prev <- loglik_pen <- loglik_pen_prev <- -.Machine$double.xmax
  # sanity check
  LLK <- c()

  while ( crit ) {

    iter <- iter + 1

    # M step -----------
    
    out_mstep <-
      mstep_inverse_sparse_M(
        data = data,
        z = z,
        penalty_omega = penalty_omega,
        penalty_gamma = penalty_gamma,
        penalty_mu = penalty_mu,
        omega = omega,
        gamma = gamma,
        control = control,
        dims = dims,
        penalization_M_mat_coord_ascent = penalization_M_mat_coord_ascent
      )

    tau <- out_mstep$parameters$tau
    mu <- out_mstep$parameters$mu
    sigma <- out_mstep$parameters$sigma
    psi <- out_mstep$parameters$psi
    omega <- out_mstep$parameters$omega
    gamma <- out_mstep$parameters$gamma
    data_cent <- out_mstep$data_cent
    det_sigma <- out_mstep$det_sigma
    det_psi <- out_mstep$det_psi

    # E step -----------
    out_estep <- estep_calc(data_cent, z, mu, tau,
                            sigma, psi, omega, gamma,
                            det_sigma, det_psi)
    z <- out_estep$z


   # Loglik -----------

    loglik <- out_estep$loglik

    pen_value_omega <- sum( sweep(abs(omega), c(1,2), penalty_omega, "*") )
    pen_value_gamma <- sum( sweep(abs(gamma), c(1,2), penalty_gamma, "*") )
    pen_value_mu <-
      pen_mu_f(
        type_penalty_mu = type_penalty_mu,
        mu = mu,
        penalty_mu = penalty_mu,
        K=K
      )

    loglik_pen <- loglik - (pen_value_omega + pen_value_gamma + pen_value_mu)
    err <- abs(loglik_pen - loglik_pen_prev) / (1 + abs(loglik_pen))
    loglik_pen_prev <- loglik_pen
    LLK <- c(LLK, loglik_pen)

    crit <- ( err > tol[1] & iter < max_iter[1] )

  }

  out_mstep <-
    mstep_inverse_sparse_M(
      data = data,
      z = z,
      penalty_omega = penalty_omega,
      penalty_gamma = penalty_gamma,
      penalty_mu = penalty_mu,
      omega = omega,
      gamma = gamma,
      control = control,
      dims = dims,
      penalization_M_mat_coord_ascent = penalization_M_mat_coord_ascent
    )
  
  # Compute bic
  
  n_par_pro <- K - 1
  n_par_mu <- sum(!out_mstep$parameters$mu==0) # non-zero values for mu matrices
  n_par_omega <-
    p * K + sum(apply(out_mstep$parameters$omega, 3, function(A)
      A[upper.tri(A)] != 0)) # non-zero values for omega matrices
  n_par_gamma <-
    q * K + sum(apply(out_mstep$parameters$gamma, 3, function(A)
      A[upper.tri(A)] != 0)) # non-zero values for gamma matrices
  
  # FIXME do we need to recompute the loglik here? 
  bic_final <- 2*loglik-(n_par_pro+n_par_mu+n_par_omega+n_par_gamma)*log(N)
  
  OUT <-
    list(
      loglik = loglik,
      loglik_pen = loglik_pen,
      K=K,
      parameters = out_mstep$parameters,
      z = z,
      classification = mclust::map(z),
      bic=bic_final,
      n_params = c(
        pro = n_par_pro,
        mu = n_par_mu,
        omega = n_par_omega,
        gamma = n_par_gamma
      ),
      penalty = list(penalty_omega = penalty_omega,
                     penalty_gamma = penalty_gamma,
                     penalty_mu=penalty_mu),
      type_penalty_mu=type_penalty_mu,
      LLK_trace = LLK, # FIXME to be deleted in the final version
      iter = iter
    )
  
  OUT
}
