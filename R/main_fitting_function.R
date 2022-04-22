# WRAPPER FUNCTION
# Function for fitting Sparse Mixture of matrix normals for Model-based Clustering  ---------------------
#' @export
fit_sparsemixmat <- function(data,
                         K = 2,
                         penalty_omega,
                         penalty_gamma,
                         penalty_mu,
                         type_penalty_mu = c("lasso", "group-lasso"),
                         penalize_diag,
                         control = EM_controls(),
                         verbose = interactive()) {
  
  # The best model is the one that maximizes the BIC
  
  # data <- scale_matrix_data(data) # work with standardized data
  # Depending on the penalty type for M_k, I define two distinct functions for updating M and for computing the objective function
  penalization_M_mat_coord_ascent <- penalization_M_mat_coord_ascent_f(type_penalty_mu = type_penalty_mu)
  
  call <- match.call()
  data_dim <- dim(data)
  
  p <- data_dim[1]
  q <- data_dim[2]
  N <- data_dim[3]
  
  all_hyperparameters <-
    expand.grid(K = K,
                penalty_omega = penalty_omega,
                penalty_gamma=penalty_gamma,
                penalty_mu=penalty_mu)
  
  n_different_models <- nrow(all_hyperparameters)
  
  models_container <-
    vector(mode = "list", length = n_different_models)
  
  if (verbose) {
    cat("Fitting:\n")
    utils::flush.console()
    pbar <- utils::txtProgressBar(min = 0,
                                  max = n_different_models,
                                  style = 3)
    on.exit(close(pbar))
    ipbar <- 0
  }
  
  # Initialization common for all models if performed via hierarchical clustering
  # FIXME atm only type_start == "hc" should be used
  hc_init <- if ( control$type_start == "hc" ) {
    # mclust::hcVVV( data = matrix(data, N, p*q, byrow = TRUE) )
    mclust::hcEII( data = matrix(data, N, p*q, byrow = TRUE) )
  } else NULL
  
  for (model in 1:n_different_models) {
    models_container[[model]] <-
      tryCatch(em_mix_mat(
        data = data,
        K = all_hyperparameters[model, "K"],
        penalty_omega = all_hyperparameters[model, "penalty_omega"],
        penalty_gamma = all_hyperparameters[model, "penalty_gamma"],
        penalty_mu=all_hyperparameters[model, "penalty_mu"],
        penalize_diag = penalize_diag,
        hc_init = hc_init,
        control = control,
        data_dim=data_dim,
        type_penalty_mu=type_penalty_mu,
        penalization_M_mat_coord_ascent=penalization_M_mat_coord_ascent
      ),error = function(e) {
        list(bic = NA)
      })
    
    if (verbose) {
      ipbar <- ipbar + 1
      utils::setTxtProgressBar(pbar, ipbar)
    }
  }
  
  models_BICS <- sapply(models_container, "[[", "bic")
  max_bic_model <- which.max(models_BICS)
  
  selected_model <- models_container[[max_bic_model]]
  selected_model$BIC <- cbind(all_hyperparameters, bic=models_BICS)
  
  if (verbose) {
    cat(
      "\nModel with K=",
      selected_model$K,
      ", penalty_omega=",
      all_hyperparameters[max_bic_model, "penalty_omega"],
      ", penalty_gamma=",
      all_hyperparameters[max_bic_model, "penalty_gamma"],
      ", penalty_mu=",
      all_hyperparameters[max_bic_model, "penalty_mu"],
      " returns the highest BIC.",
      sep = ""
    )
  }
  selected_model
  
}