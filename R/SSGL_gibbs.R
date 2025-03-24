########################################
########################################
## FUNCTION FOR IMPLEMENTING THE SSGL ##
## USING GIBBS SAMPLING               ##
########################################
########################################

# This function implements group-regularized regression models in the exponential
# dispersion family with the spike-and-slab group lasso (SSGL) penalty.

# INPUTS:
# Y = n x 1 vector of responses (y_1, ...., y_n) for training data
# X = n x p design matrix for training data, where ith row is (x_{i1},..., x_{ip})
# groups = G x 1 vector of group indices or factor level names for each of the p individual covariates.
# family = the exponential family. Currently allows "gaussian", "binomial", and "poisson".
# X_test = n_test x p design matrix for test data. If missing, then program computes in-sample
#          predictions on training data X
# group_weights = group-specific weights. Default is to use the square roots of the group sizes.
# lambda0 = spike hyperparameter. Default is 5
# lambda1 = Slab hyperparameter. Default is 1
# a = shape hyperparameter in B(a,b) prior on mixing proportion. Default is 1
# b = shape hyperparameter in B(a,b) prior on mixing proportion. Default is G for number of groups
# burn = number of MCMC samples to discard during the burn-in period. Default is 1000.
# n_mcmc = number of MCMC iterations to save after burn-in. Default is 2000.
# save_samples = Boolean variable whether to save all MCMC samples. Default is TRUE

# OUTPUT:
#	beta_hat = posterior mean estimate for beta coefficients
# Y_pred_hat = posterior mean predictions E(Y) based on test data X_test
#              If X_test was left blank or X_test=X, then in-sample predictions on X are returned.
# beta_lower = lower endpoint of the 95% credible interval for beta coefficients
# beta_upper = upper endpoint of the 95% credible interval for beta coefficients
# Y_pred_lower = lower endpoint of the 95% prediction interval for mean response values based on test data X_test
# Y_pred_upper = upper endpoint of the 95% prediction interval for mean response values based on test data X_test
# beta_samples = posterior samples of beta after burnin
# Y_pred_samples = posterior samples of Y_pred after burnin


SSGL_gibbs = function(Y, X, groups, family=c("gaussian","binomial","poisson"),
                      X_test, group_weights, lambda0=5, lambda1=1, 
                      a=1, b=length(unique(groups)),
                      burn=1000, n_mcmc=2000, save_samples=TRUE) {
  
  
  ##################
  ##################
  ### PRE-CHECKS ###
  ##################
  ##################

  ## Enumerate groups if not already done
  group_numbers = as.numeric(groups)
  
  ## Number of groups and covariates overall
  X = as.matrix(X)
  G = length(unique(group_numbers))
  n = dim(X)[1]
  p = dim(X)[2]

  ## If test data is missing, make it the same as the training data
  if(missing(X_test)) X_test = X
  n_test = dim(X_test)[1]
  X_test = as.matrix(X_test)
  
  ## Check that dimensions are conformable
  if(length(Y) != dim(X)[1])
    stop("Non-conformable dimensions of Y and X.")
  ## Check that X and X_test have the same number of columns
  if(dim(X_test)[2] != p)
    stop("X and X_test should have the same number of columns.")

  ## Set group-specific weights
  if(missing(group_weights))
    group_weights = sqrt(as.vector(table(group_numbers)))
  
  ## Coercion
  family <- match.arg(family)
  
  ## Check that the data can be used for the respective family
  if(family=="poisson"){
    if(any(Y<0))
      stop("All counts must be greater than or equal to zero.")
    if(any(Y-floor(Y)!=0))
      stop("All counts must be whole numbers.")
  }
  if(family=="binomial"){
    if(any(Y<0))
      stop("All binary responses must be either '0' or '1.'")
    if(any(Y>1))
      stop("All binary responses must be either '0' or '1.'")
    if(any(Y-floor(Y)!=0))
      stop("All binary responses must be either '0' or '1.'")
  }
  
  ## Check that weights are all greater than or equal to 0
  if(!missing(group_weights)){
    if(!all(group_weights>=0))
      stop("All group-specific weights should be nonnegative.")
  }

  ## Check hyperparameters to be safe
  if ((lambda1 <= 0) || (lambda0 <= 0) || (a <= 0) || (b <= 0))
    stop("Please make sure that all hyperparameters are strictly positive.")
  
  ## Check that burn and n_mcmc are reasonable
  if((burn<0) || (n_mcmc<=0))
    stop("Please specify valid values for burn and/or n_mcmc.")
  
  ## Rescale the lambda0s
  scaled_lambda0 = lambda0*group_weights

  ## Initialize values for beta, theta, gamma, xi, and omega
  beta = rep(0,p)               # initial beta
  theta = 0.5                   # initial theta
  latent_gamma = rep(0, G)      # initial gamma
  xi = rep(1,G)                 # to hold xi (vector)
  xi_diag = rep(1, p)           # to hold entries of block diagonal matrix D_xi
  omega = rep(1,n)              # to hold omega (vector)

  ## total iterations
  n_iter = burn+n_mcmc
  
  ## Other initializations
  if(family=="gaussian"){
    Z_tilde = Y
    X_tilde = X
    XtildeXtilde = crossprod(X,X)
    XtildeZtilde = crossprod(X,Z_tilde)
  } else if(family=="binomial"){
    kappa = Y-0.5
    # Initial Z
    Z = kappa/omega
  } else if(family=="poisson"){
    M = max(Y)+1
    kappa = 0.5*(Y-M)
    # Initial Z
    Z = kappa/omega+log(M)
  }
    
  ## Matrices and vectors to hold solutions
  beta_samples = matrix(0, nrow=p, ncol=n_iter)
  Y_pred_samples = matrix(0, nrow=n_test, ncol=n_iter)
  
  ################################
  ################################
  ### Gibbs sampling algorithm ###
  ################################
  ################################
  j = 0
  while (j < n_iter){
    
    ## Print iteration number
    j = j + 1
    if (j %% 100 == 0) {
      cat("Gibbs sampling iteration:", j, "\n")
    }
      
    ## Update beta
    if((family=="binomial") || (family=="poisson")){
      X_tilde = sqrt(omega)*X
      Z_tilde = sqrt(omega)*Z
      XtildeXtilde = crossprod(X_tilde,X_tilde)
      XtildeZtilde = crossprod(X_tilde,Z_tilde)
    }
    
    ## Sample beta
    if(p<=n){
      ridge_beta = XtildeXtilde + diag(1/xi_diag)
      inv_ridge_beta = Matrix::chol2inv(Matrix::chol(ridge_beta))
      mu_beta = crossprod(inv_ridge_beta, XtildeZtilde)  
      beta = MASS::mvrnorm(1, mu_beta, inv_ridge_beta)
      
    } else if(p>n) {
      
      ## Use the exact sampling algorithm of Bhattacharya et al. (2016)
      m = stats::rnorm(p, mean=rep(0,p), sd=sqrt(xi_diag))
      delta = stats::rnorm(n, mean=0, sd=1)
      v = crossprod(t(X_tilde),m)+delta
      v_star = Z_tilde-v
      K = crossprod(t(X_tilde), (xi_diag*t(X_tilde))) + diag(n)
      w = crossprod(t(Matrix::chol2inv(Matrix::chol(K))), v_star)
      beta = m + crossprod(t(xi_diag*t(X_tilde)), w)
    } 
    beta_samples[,j] = beta
    
    ## Update gamma and xi
    for(g in 1:G){
      m_g = group_weights[g]^2
      lambda0_g = scaled_lambda0[g]
      
      # Sample latent_gamma[g]
      pi1 = theta*lambda1^(m_g+1)*exp(-lambda1^2*xi[g]/2)
      pi2 = (1-theta)*lambda0_g^(m_g+1)*exp(-lambda0_g^2*xi[g]/2)
      if(pi1 == 0 & pi2 == 0){
        prop = 1
      } else {
        prop = pi1/(pi1+pi2)
      }
      latent_gamma[g] = stats::rbinom(1, 1, prop)
      # Set lambda_star 
      lambda_star_g = latent_gamma[g]*lambda1 + (1-latent_gamma[g])*lambda0_g
      # Sample xi[g]
      norm_term = sum(beta[which(groups==g)]^2)
      v <- max(norm_term, .Machine$double.eps)
      xi[g] = GIGrvg::rgig(1, lambda=0.5, chi=v, psi=lambda_star_g^2)
      
      # Store entries in xi_diag
      xi_diag[which(groups==g)] = xi[g]
    }

    ## Update theta
    sum_gamma = sum(latent_gamma)
    theta = stats::rbeta(1, a+sum_gamma, b+G-sum_gamma)
    
    ## Sample omegas and update Z if binomial or Poisson
    if(family=="binomial"){
      eta = crossprod(t(X),beta)
      omega = pmax(BayesLogit::rpg(n,1, eta), 0.03) # For numerical stability
      Z = kappa/omega
    }
    if(family=="poisson"){
      eta = crossprod(t(X),beta)
      omega = BayesLogit::rpg(n, M, eta-log(M))
      Z = kappa/omega+log(M)
    }
    
    ## Save MCMC samples for predictions
    if(family=="gaussian"){  
        Y_pred_samples[,j] = crossprod(t(X_test),beta)
    } else if(family=="binomial") {
        Y_pred_samples[,j] = 1/(1+exp(-(crossprod(t(X_test),beta))))
    } else if(family=="poisson"){
        Y_pred_samples[,j] = exp(crossprod(t(X_test),beta))
    }
    
  }

  ## Discard burnin and extract summary statistics
  beta_samples = beta_samples[,(burn+1):n_iter]
  Y_pred_samples = Y_pred_samples[,(burn+1):n_iter]
  rownames(beta_samples) = groups
  
  ## Posterior means
  beta_hat = rowMeans(beta_samples)
  Y_pred_hat = rowMeans(Y_pred_samples)
  
  # 95 percent CIs
  beta_intervals = apply(beta_samples, 1, function(x) stats::quantile(x, prob=c(.025,.975)))
  beta_lower = beta_intervals[1,]
  beta_upper = beta_intervals[2,]
  
  Y_pred_intervals = apply(Y_pred_samples, 1, function(x) stats::quantile(x, prob=c(.025,.975)))
  Y_pred_lower = Y_pred_intervals[1,]
  Y_pred_upper = Y_pred_intervals[2,]
  
  
  #####################
  #####################
  ### Return a list ###
  #####################
  #####################
  if(save_samples==TRUE){

      SSGL_MCMC_output <- list(beta_hat = beta_hat,
                               Y_pred_hat = Y_pred_hat,
                               beta_lower = beta_lower,
                               beta_upper = beta_upper,
                               Y_pred_lower = Y_pred_lower,
                               Y_pred_upper = Y_pred_upper,
                               beta_samples = beta_samples,
                               Y_pred_samples = Y_pred_samples)
      
  } else {
    SSGL_MCMC_output <- list(beta_hat = beta_hat,
                             Y_pred_hat = Y_pred_hat,
                             beta_lower = beta_lower,
                             beta_upper = beta_upper,
                             Y_pred_lower = Y_pred_lower,
                             Y_pred_upper = Y_pred_upper)
  }
  
  # Return list
  return(SSGL_MCMC_output)
}