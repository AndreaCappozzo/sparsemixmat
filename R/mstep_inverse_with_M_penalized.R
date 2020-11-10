
mstep_inverse <- function(data, z,
                          penalty_omega,
                          penalty_gamma,
                          sigma,
                          omega,
                          theta,
                          gamma,   # needed to initialize gamma
                          control,
                          dims,
                          initialization = FALSE)
  # M step for inverse covariance parameterization
{
  p <- dims$p
  q <- dims$q
  N <- dim(data)[3] # compute N in this way to account for a possible subset in the initialization step
  # N <- dims$N
  K <- dims$K

  # inner M step control parameters
  tol <- control$tol[2]
  max_iter <- control$max_iter[2]

  # update tau and mu
  out <- mean_w_array(data, z)
  Nk <- c(out$Nk)
  tau <- c(out$tau)
  mu <- out$mean
  data_cent <- out$data_cent

  # store parameters
  # if ( initialization ) {
  #   sigma[] <- omega[] <- 0
  #   theta[] <- 0
  # }
  # sigma <- omega <- array(0, dim = c(p, p, K))
  # theta <- array(0, dim = c(q, q, K))

  # start iterative algorithm
  crit <- TRUE
  iter <- 0
  obj <- obj_prev <- -.Machine$double.xmax
  obj_vec <- det_sigma <- det_theta <- rep(NA, K)

  # sanity check
  OBJ <- c()

  while ( crit ) {
    iter <- iter + 1
    start <- if ( iter < 2 & initialization ) "cold" else "warm"

    for ( k in 1:K ) {

      # estimate sparse omega
      # scattering matrix needed for graphical lasso problem for sigma - corresponds to S in eq (2.1) of \cite{Friedman2008}
      S_sigma <- cov_w_array(data_cent[k], z[,k, drop = FALSE], mu[,,k, drop = FALSE],
                             gamma[,,k, drop = FALSE], Nk[k], pbyp = TRUE)$covmat[,,1]

      # gl <- glassoFast::glassoFast(S = S_sigma, rho = penalty_omega)
      # gl <- glasso::glasso(s = S_sigma, rho = penalty_omega[1,2], nobs = tau[k]*N)
      gl <- glassoFast::glassoFast(S = S_sigma, rho = penalty_omega, start = start, w.init = sigma[,,k], wi.init = omega[,,k])

      omega[,,k] <- gl$wi

      # normalize omega so that omega[1,1,k] == 1
      omega[,,k] <- omega[,,k] / omega[1,1,k]
      #
      # normalize to det(omega) = 1
      # e <- eigen(omega[,,k], only.values = TRUE)$val
      # omega[,,k] <- omega[,,k] / exp( 1/p*sum(log(e)) )

      # estimate sparse gamma
      # scattering matrix needed for graphical lasso problem for theta - corresponds to S in eq (2.1) of \cite{Friedman2008}
      S_theta <- cov_w_array(data_cent[k], z[,k, drop = FALSE], mu[,,k, drop = FALSE],
                             omega[,,k, drop = FALSE], Nk[k], pbyp = FALSE)$covmat[,,1]

      # gl <- glassoFast::glassoFast(S = S_theta, rho = penalty_theta)
      # gl <- glasso::glasso(s = S_theta, rho = penalty_theta[1,2], nobs = tau[k]*N)
      gl <- glassoFast::glassoFast(S = S_theta, rho = penalty_theta, start = start, w.init = theta[,,k], wi.init = gamma[,,k])

      gamma[,,k] <- gl$wi

      # estimate covariance matrices
      sigma[,,k] <- solve(omega[,,k])
      theta[,,k] <- solve(gamma[,,k])

      # compute objective function for convergence
      m_obj <- mstep_obj(data_cent[k], z[,k, drop = FALSE], mu[,,k, drop = FALSE],
                         sigma[,,k, drop = FALSE], theta[,,k, drop = FALSE],
                         omega[,,k, drop = FALSE], gamma[,,k, drop = FALSE],tau = tau[k])
      obj_vec[k] <- m_obj$obj
      det_sigma[k] <- m_obj$det_sigma
      det_theta[k] <- m_obj$det_theta
    }

    # check convergence
    pen_value_omega <- sum( sweep(abs(omega), c(1,2), penalty_omega, "*") )
    pen_value_gamma <- sum( sweep(abs(gamma), c(1,2), penalty_gamma, "*") )

    obj <- sum(obj_vec) - (pen_value_omega + pen_value_gamma)
    OBJ <- c(OBJ, obj)
    err <- abs(obj - obj_prev)/(1 + abs(obj))
    obj_prev <- obj
    crit <- (err > tol) & (iter < max_iter)

  }

  return(
    list(parameters = list(tau = tau,
                           mu = mu,
                           sigma = sigma,
                           theta = theta,
                           omega = omega,
                           gamma = gamma),
         data_cent = data_cent,
         det_sigma = det_sigma,
         det_theta = det_theta,
         OBJ = OBJ)
  )
}







mstep_inverse_sparse_M <- function(data,
                                   z,
                                   penalty_omega,
                                   penalty_gamma,
                                   penalty_mu,
                                   gamma, # needed for the penalization of the mu matrices
                                   omega,
                                   mu, # the estimates of the previous m-step
                                   control,
                                   dims
)
  # M step for inverse covariance parameterization
{
  p <- dims$p
  q <- dims$q
  N <- dim(data)[3] # compute N in this way to account for a possible subset in the initialization step
  # N <- dims$N
  K <- dims$K

  # inner M step control parameters
  tol <- control$tol[2]
  max_iter <- control$max_iter[2]

  # update tau and sparse mu
  out <- w_array_given_mean(data, z,mu)
  Nk <- c(out$Nk)
  tau <- c(out$tau)
  mu <- out$mean
  data_cent <- out$data_cent

  mu_penalized <- array(dim = c(p,q,K))
  data_cent_penalized <- vector(mode = "list", length = K)

  for (k in 1:K) {
    mu_penalized[, , k] <-
      penalization_M_mat_cpp(
        data = data,
        data_cent = data_cent[[k]],
        z = z[, k],
        Nk = Nk[k],
        mu = mu[, , k],
        # from the previous m-step
        omega = omega[, , k],
        gamma = gamma[, , k],
        penalty_mu = penalty_mu,
        p = p,
        q = q,
        N = N
      )

    data_cent_penalized[[k]] <-
      array(apply(data, 3, function(A)
        A - mu_penalized[, , k]), dim = c(p, q, N))
  }

  # store parameters
  sigma <- array(0, dim = c(p, p, K))
  theta <- array(0, dim = c(q, q, K))
  # gamma <- array(1, dim = c(q, q, K))  # start to identity matrix

  # start iterative algo
  crit <- TRUE
  iter <- 0
  obj <- obj_prev <- -.Machine$double.xmax
  obj_vec <- det_sigma <- det_theta <- rep(NA, K)

  # sanity check
  OBJ <- c()

  while ( crit ) {
    iter <- iter + 1

    for ( k in 1:K ) {

      # estimate sparse mu

      # mu_penalized[, , k] <- # FIXME as it is this is wrong: since the m-step is not in closed form, do we need to recompute the penalized mu within the loop? Need to think about it
      #   penalization_M_mat_cpp(
      #     data = data,
      #     data_cent = data_cent[[k]],
      #     z = z[, k],
      #     Nk = Nk[k],
      #     mu = mu[, , k], # from the previous m-step
      #     omega = omega[, , k],
      #     gamma = gamma[, , k],
      #     penalty_mu = penalty_mu,
      #     p = p,
      #     q = q,
      #     N = N
      #   )
      #
      # data_cent_penalized[[k]] <- array(apply(data, 3, function(A) A-mu_penalized[,,k]), dim = c(p,q,N))

      # estimate sparse omega
      # scattering matrix needed for graphical lasso problem for sigma - corresponds to S in eq (2.1) of \cite{Friedman2008}
      S_sigma <- cov_w_array(data_cent_penalized[k], z[,k, drop = FALSE], mu_penalized[,,k, drop = FALSE],
                             gamma[,,k, drop = FALSE], Nk[k], pbyp = TRUE)$covmat[,,1]
      gl <- glassoFast::glassoFast(S = S_sigma, rho = penalty_omega)
      omega[,,k] <- gl$wi

      # normalize theta so that omega[1,1,k] == 1
      omega[,,k] <- omega[,,k] / omega[1,1,k]

      # estimate sparse gamma
      # scattering matrix needed for graphical lasso problem for theta - corresponds to S in eq (2.1) of \cite{Friedman2008}
      S_theta <- cov_w_array(data_cent_penalized[k], z[,k, drop = FALSE], mu_penalized[,,k, drop = FALSE],
                             omega[,,k, drop = FALSE], Nk[k], pbyp = FALSE)$covmat[,,1]
      gl <- glassoFast::glassoFast(S = S_theta, rho = penalty_gamma)
      gamma[,,k] <- gl$wi

      # estimate covariance matrices
      sigma[,,k] <- solve(omega[,,k])
      theta[,,k] <- solve(gamma[,,k])

      # compute objective function for convergence
      m_obj <- mstep_obj(data_cent_penalized[k], z[,k, drop = FALSE], mu_penalized[,,k, drop = FALSE],
                         sigma[,,k, drop = FALSE], theta[,,k, drop = FALSE],
                         omega[,,k, drop = FALSE], gamma[,,k, drop = FALSE],tau = tau[k])
      obj_vec[k] <- m_obj$obj
      det_sigma[k] <- m_obj$det_sigma
      det_theta[k] <- m_obj$det_theta
    }

    # #FIXME check theory
    # mu <- mu_penalized
    # data_cent <- data_cent_penalized

    # check convergence
    pen_value_omega <- sum( sweep(abs(omega), c(1,2), penalty_omega, "*") )
    pen_value_gamma <- sum( sweep(abs(gamma), c(1,2), penalty_gamma, "*") )

    obj <- sum(obj_vec) - (pen_value_omega +pen_value_gamma)
    OBJ <- c(OBJ, obj)
    err <- abs(obj - obj_prev)/(1 + abs(obj))
    obj_prev <- obj
    crit <- (err > tol) & (iter < max_iter)

  }

  # for (k in 1:K) {
  #   mu_penalized[, , k] <-
  #     penalization_M_mat_cpp(
  #       data = data,
  #       data_cent = data_cent[[k]],
  #       z = z[, k],
  #       Nk = Nk[k],
  #       mu = mu[, , k],
  #       # from the previous m-step
  #       omega = omega[, , k],
  #       gamma = gamma[, , k],
  #       penalty_mu = penalty_mu,
  #       p = p,
  #       q = q,
  #       N = N
  #     )
  #
  #   data_cent_penalized[[k]] <-
  #     array(apply(data, 3, function(A)
  #       A - mu_penalized[, , k]), dim = c(p, q, N))
  # }

  return(
    list(parameters = list(tau = tau,
                           mu = mu_penalized,
                           sigma = sigma,
                           theta = theta,
                           omega = omega,
                           gamma = gamma),
         data_cent = data_cent_penalized,
         det_sigma = det_sigma,
         det_theta = det_theta)
  )
}
