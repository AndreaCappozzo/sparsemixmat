# M step for inverse covariance matrices parametrization
#' @export
mstep_inverse_sparse_M <- function(data,
                                   z,
                                   penalty_omega,
                                   penalty_gamma,
                                   penalty_mu,
                                   omega,
                                   gamma, 
                                   control,
                                   dims){
  p <- dims$p
  q <- dims$q
  N <- dim(data)[3] # compute N in this way to account for a possible subset in the initialization step
  K <- dims$K
  
  # inner M step control parameters
  tol <- control$tol[2]
  max_iter <- control$max_iter[2]
  
  # Compute sample mean, Nk and tau
  out <- mean_w_array(data, z)
  Nk <- c(out$Nk)
  tau <- c(out$tau)
  mu <- out$mean
  data_cent <- out$data_cent
  
  if(any(penalty_mu!=0)) {
    # Compute sparse mean matrices via coordinate ascent algorithm
    for (k in 1:K) {
      out_penalized <-
        penalization_M_mat_coord_ascent(
          data = data,
          data_cent = data_cent[[k]],
          z = z[, k],
          Nk = Nk[k],
          mu = mu[, , k], # I initialize it with the sample mean at the current iteration
          omega = omega[, , k],
          gamma = gamma[, , k],
          penalty_mu = penalty_mu,
          p = p,
          q = q,
          N = N
        )
      mu[, , k] <- out_penalized$mu_penalized
      data_cent[[k]] <- out_penalized$data_cent_penalized
    }
  }
  
  # store parameters
  sigma <- array(0, dim = c(p, p, K))
  psi <- array(0, dim = c(q, q, K))
  
  # start iterative algorithm
  crit <- TRUE
  iter <- 0
  obj <- obj_prev <- -.Machine$double.xmax
  obj_vec <- det_sigma <- det_psi <- rep(NA, K)
  
  # sanity check
  OBJ <- c()
  
  while ( crit ) {
    iter <- iter + 1
    # start <- if ( iter < 2 & initialization ) "cold" else "warm" FIXME this seems to create issues: need to investigate it
    start <- "cold"
    
    for ( k in 1:K ) {
      
      # estimate sparse omega
      # scattering matrix needed for graphical lasso problem for sigma - corresponds to S in eq (2.1) of \cite{Friedman2008}
      S_sigma <- cov_w_array(data_cent[k], z[,k, drop = FALSE], mu[,,k, drop = FALSE],
                             gamma[,,k, drop = FALSE], Nk[k], pbyp = TRUE)$covmat[,,1]
      
      if(any(penalty_omega!=0)) {
      gl <- glassoFast::glassoFast(S = S_sigma, rho = (2*penalty_omega)/(Nk[k]*q), start = start, w.init = sigma[,,k], wi.init = omega[,,k])
      
      omega[,,k] <- gl$wi
      sigma[,,k] <- solve(omega[,,k])
      } else {
        sigma[,,k] <- S_sigma
        omega[,,k] <- solve(S_sigma)
      }
      
      # estimate sparse gamma
      # scattering matrix needed for graphical lasso problem for psi - corresponds to S in eq (2.1) of \cite{Friedman2008}
      S_psi <- cov_w_array(data_cent[k], z[,k, drop = FALSE], mu[,,k, drop = FALSE],
                             omega[,,k, drop = FALSE], Nk[k], pbyp = FALSE)$covmat[,,1]
      if(any(penalty_gamma!=0)) {
      gl <- glassoFast::glassoFast(S = S_psi, rho = (2*penalty_gamma)/(Nk[k]*p), start = start, w.init = psi[,,k], wi.init = gamma[,,k])
      
      gamma[,,k] <- gl$wi
      } else {
        gamma[,,k] <- solve(S_psi)
      }
      
      # normalize to det(gamma) = 1
      e <- eigen(gamma[,,k], only.values = TRUE)$val
      gamma[,,k] <- gamma[,,k] / exp( 1/q*sum(log(e)) )
      psi[,,k] <- solve(gamma[,,k])
      # compute objective function for convergence
      m_obj <- mstep_obj(data_cent[k], z[,k, drop = FALSE], mu[,,k, drop = FALSE],
                         sigma[,,k, drop = FALSE], psi[,,k, drop = FALSE],
                         omega[,,k, drop = FALSE], gamma[,,k, drop = FALSE],tau = tau[k])
      obj_vec[k] <- m_obj$obj
      det_sigma[k] <- m_obj$det_sigma
      det_psi[k] <- m_obj$det_psi
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
                           psi = psi,
                           omega = omega,
                           gamma = gamma),
         data_cent = data_cent,
         det_sigma = det_sigma,
         det_psi = det_psi)
  )
}
