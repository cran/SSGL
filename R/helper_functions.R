####################
# Helper functions #
####################

#' Prior density Psi. No need for normalizing constant C_d as it cancels out
#' @keywords internal
Psi = function(beta, lambda) {
  m = length(beta)
  dens = lambda^m * exp(-lambda*sqrt(sum(beta^2)))
  
  return(dens)
}

#' pStar function
#' #' @keywords internal
pStar = function(beta, lambda1, lambda0, theta) {
  Psi1 = Psi(beta=beta, lambda=lambda1)
  Psi0 = Psi(beta=beta, lambda=lambda0)
  
  ## if a coefficient is really large then both these will 
  ## numerically be zero because R can't handle such small numbers
  if ((theta*Psi1) == 0 & (1 - theta)*Psi0 == 0) {
    p = 1
  } else {
    p = (theta*Psi1) / (theta*Psi1 + (1 - theta)*Psi0)
  }
  return(p)
}

#' Lambda star function
#' @keywords internal
lambdaStar = function(beta, lambda1, lambda0, theta) {
  p = pStar(beta = beta, lambda1 = lambda1,
            lambda0 = lambda0, theta = theta)
  
  l = lambda1*p + lambda0*(1 - p)
  return(l)
}

#' EM algorithm for SB-GAM.
#' Here, lambda0 is a single tuning parameter
#' @keywords internal
SSGL_EM = function(Y, X, groups, 
                   family=c("gaussian","binomial","poisson"), 
                   n, G, a, b, group_weights, lambda0, lambda1, beta0_init, 
                   beta_init, theta_init, max_iter, tol){
  
  ## Coercion
  family <- match.arg(family)
  
  ## Initialize the following values
  difference = 100*tol
  counter = 0
  pstar_k = rep(0,G)
  lambdastar_k=rep(0,G) # To hold lambdastar for each group of coefficients
  
  ## Initialize parameters
  beta0 = beta0_init
  beta = beta_init
  theta = theta_init
  
  ## Update the parameters
  while( (difference > tol) & (counter < max_iter) ){
  
    ## Iterate counter  
    counter = counter+1
    ## Keep track of old beta
    beta_old = beta
      
    ##############
    ##############
    ### E-step ###
    ##############
    ##############
    for(k in 1:G){
      ## Which groups are active
      active = which(groups == k)
      ## Update pStar
      pstar_k[k] = pStar(beta_old[active], lambda1, lambda0, theta)
      # Update lambdastar_k for groups 1,...,G 
      lambdastar_k[k] = lambda1*pstar_k[k] + lambda0*(1-pstar_k[k])
    }
    
    ############## 
    ##############
    ### M-step ###
    ##############
    ##############
    
    ## Update theta
    theta = (a-1 + sum(pstar_k))/(a+b+G-2)
    
    ## Update beta0 and beta
    ## Note that grpreg solves is -(1/n)*loglik(beta0,beta) + pen(beta)
    ## so we have to multiply by 1/n in the penalty

    solve_obj = grpreg::grpreg(X, Y, group=groups, penalty="grLasso", family=family, 
                               lambda=1, group.multiplier=group_weights*(lambdastar_k/n))
    beta0 = solve_obj$beta[1]
    beta = solve_obj$beta[-1]

    ## Update diff
    diff = sum((beta-beta_old)^2)/(sum(beta_old^2)+1e-8)
  }
  
  ## Store beta0, beta, theta in a list
  SSGL_EM_output <- list(beta0 = beta0,
                         beta = beta,
                         theta = theta)
  
  # Return list
  return(SSGL_EM_output)
}