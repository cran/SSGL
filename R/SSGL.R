########################################
########################################
## FUNCTION FOR IMPLEMENTING THE SSGL ##
########################################
########################################

# This function implements group-regularized regression models for linear, logistic
# and Poisson regression with the spike-and-slab group lasso (SSGL) penalty.

# INPUTS:
# Y = n x 1 vector of responses (y_1, ...., y_n) for training data
# X = n x p design matrix for training data, where ith row is (x_{i1},..., x_{ip})
# groups = p x 1 vector of group indices or factor level names for each of the p columns of X.
# family = the exponential family. Currently allows "gaussian", "binomial", and "poisson". 
# X_test = n_test x p design matrix for test data. If missing, then program computes in-sample
#          predictions on training data X
# group_weights = group-specific weights. Default is to use the square roots of the group sizes.
# n_lambda0 = number of lambda0's in our grid. Default is 20.
# lambda0 = a grid of values for the spike hyperparameter. Note that the program automatically orders 
#           these in descending order, so the solution path is from least sparse to most sparse. 
#           If lambda0 is not specified, the program chooses a default equispaced grid from 100 to 5.
# lambda1 = a fixed small value for the slab hyperparameter. Default is 1
# a = shape hyperparameter in Beta(a,b) prior on mixing proportion. Default is 1
# b = shape hyperparameter in Beta(a,b) prior on mixing proportion. Default is G for number of groups
# max_iter = maximum number of iterations. Default is 100
# tol = convergence criteria. Default is 1e-6
# return_GIC = TRUE
# print_lambda0 = Boolean variable whether to print the current lambda0 in our grid. Default is TRUE

# OUTPUT:
# lambda0 = grid of lambda0's in descending order. 
#	beta = p x L matrix of estimated basis coefficients. The kth column in beta corresponds to the kth 
#        spike parameter in our descending lambda0 grid.
#	beta0 = L x 1 vector of estimated intercepts. The kth entry in beta0 corresponds to the kth spike parameter 
#         in our descending lambda0 grid.
# classifications = p x L matrix of classifications, where "1" indicates that the variable was selected 
#                   and "0" indicates that the variable was not selected. The kth column in classifications 
#                   corresponds to the kth spike parameter in our descending lambda0 grid.
# Y_pred = n_test x L matrix of predicted mean response values based on test data in X_test. If
#          X_test was left blank or X_test=X, then in-sample predictions on X are returned.
# GIC = Lx1 vector of GIC values. The lth entry of GIC corresponds to the lth entry in our lambda0 grid.
#       Not returned if return_GIC=FALSE. 
# lambda0_GIC_min = value of lambda0 that minimizes the GIC. Not returned if return_GIC=FALSE.
# min_GIC_index = the index of lambda0_min in lambda0. Not returned if return_GIC=FALSE.

SSGL = function(Y, X, groups, family=c("gaussian","binomial","poisson"),
                X_test, group_weights, n_lambda0=25, lambda0, lambda1=1, a=1, b=length(unique(groups)),
                max_iter=100, tol = 1e-6, return_GIC=TRUE, print_lambda0=TRUE) {
  
  
  ##################
  ##################
  ### PRE-CHECKS ###
  ##################
  ##################

  ## Enumerate groups if not already done
  groups = as.factor(groups)
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
  ## Number of lambdas
  if(n_lambda0 < 1)
    stop("The number of lambdas must be at least one.")
  ## If user specified lambda, check that all lambdas are greater than 0
  if(!missing(lambda0)) {
    n_lambda0 = length(lambda0) # Override n_lambda0 with the length of lambda0
    if (!all(lambda0>0))
      stop("All lambda0s should be strictly positive.")
  }  

  ## Default parameters for missing arguments
  L = n_lambda0
  if(missing(lambda1)) lambda1 = 1
  if(missing(a)) a = 1
  if(missing(b)) b = G
  ## Check hyperparameters to be safe
  if ((lambda1 <= 0) || (a <= 0) || (b <= 0))
    stop("Please make sure that all hyperparameters are strictly positive.")
  if(missing(lambda0)){ 
    if(family=="binomial"){
      # Compute max_lambda
      norms_vec = rep(0,G)
      for(g in 1:G){
        ind_g = which(group_numbers==g)
        X_g = X[,ind_g]
        norms_vec[g] = sqrt(sum((t(X_g)%*%(0.25*(Y-0.5)))^2))/sqrt(length(ind_g))  
      }
      max_lambda = min(10, max(abs(norms_vec)))
      max_lambda = max(max_lambda, lambda1+7)
      
      ## Create grid of lambdas
      if(n_lambda0==1){ 
        lambda0 = max_lambda/2
      } else if(n_lambda0 > 1){
        lambda0 = seq(from=max_lambda, to=lambda1+1, length=n_lambda0)
      }
    } 
    if(family!="binomial") {
      # Compute max_lambda
      norms_vec = rep(0, G)
      if(family=="gaussian"){
        mu=0
      } else {
        mu=1
      }
      for(g in 1:G){
        ind_g = which(group_numbers==g)
        X_g = X[,ind_g]
        norms_vec[g] = sqrt(sum((t(X_g)%*%(Y-mu))^2))/sqrt(length(ind_g))  
      }
      
      ## Create grid of lambdas
      if(n_lambda0==1){ 
        lambda0 = max_lambda/2
      } else if(n_lambda0 > 1) {
        max_lambda = min(100, max(abs(norms_vec)))
        max_lambda = max(max_lambda, lambda1+50)
        if(family=="poisson"){
          lambda0 = seq(from=max_lambda, to=lambda1+7, length=n_lambda0)
        } else if(family=="gaussian"){
          lambda0 = seq(from=max_lambda, to=lambda1+1, length=n_lambda0)
        }
      }
    }
  }
  
  # Sort lambda0 just in case
  lambda0 = sort(lambda0, decreasing=TRUE)
    
  ## Initialize values for beta0, beta, and theta
  beta0_init = 0                # initial beta0
  beta_init = rep(0,p)          # initial beta
  theta_init = 0.5              # initial theta
  
  ## Matrices and vectors to hold solutions
  beta = matrix(0, nrow=p, ncol=L)
  beta0 = rep(0,L)


  ####################
  ####################
  ### EM algorithm ###
  ####################
  ####################
  
  for(l in 1:L){

    ## To output iteration
    if(print_lambda0==TRUE)
        cat("lambda0 = ", lambda0[l], "\n")
    
    ## Run EM algorithm on training data
    EM_output = SSGL_EM(Y=Y, X=X, groups=group_numbers, family=family, 
                        n=n, G=G, a=a, b=b, group_weights=group_weights, 
                        lambda0=lambda0[l], lambda1=lambda1, 
                        beta0_init=beta0_init, beta_init=beta_init, 
                        theta_init=theta_init, max_iter=max_iter, tol=tol)
    
    ## Save values for the next pass and also store them as a 'warm start'
    ## to the next lambda0 in our grid
    beta[,l] = EM_output$beta
    beta0[l] = EM_output$beta0
    
    if(L>1) {  
      beta0_init = beta[l]
      beta_init = beta[,l]
      theta_init = EM_output$theta
    }
  }
  
  # Row names for beta
  beta = as.matrix(beta)
  rownames(beta) = groups
  
  ## Compute predictions on test data
  Y_pred = matrix(0, n_test, L) # To compute predictions on test set
  for(l in 1:L){
    if(family=="gaussian"){  
      Y_pred[,l] = beta0[l]+X_test%*%beta[,l]
    } else if(family=="binomial") {
      Y_pred[,l] = 1/(1+exp(-(beta0[l]+X_test%*%beta[,l])))
    } else if(family=="poisson"){
      Y_pred[,l] = exp(beta0[l]+X_test%*%beta[,l])
    }
  }
 
  ## Compute classifications
  classifications = matrix(0, nrow=G, ncol=L)
  
  for(l in 1:L){
    for (g in 1:G) {
      active = which(group_numbers == g)
      
      ## Update classifications
      if(!identical(as.numeric(beta[active,l]), rep(0,length(active))))
          classifications[g,l] = 1   
    }
  } 
  
  ## Row names for classifications
  row.names(classifications) = unique(groups)
  
  if(return_GIC==TRUE){
    ## Compute GIC
    SSGL_dev = rep(0, L)
    GIC = rep(0, L)
  
    for(l in 1:L){
      # Compute deviance first
      if(family=="gaussian"){
        mu = beta0[l]+X%*%beta[,l]
        SSGL_dev[l] = sum((Y-mu)^2)
      } 
      if(family=="binomial"){
        mu = 1/(1+exp(-(beta0[l]+X%*%beta[,l])))
        term1 = sum(log(mu[Y==1]))
        term2 = sum(log(1-mu[Y==0]))
        SSGL_dev[l] = -2*(term1+term2)
      }
      if(family=="poisson"){
        mu = exp(beta0[l]+X%*%beta[,l])
        term1 = rep(0, length(Y))
        term1[Y!=0] = Y[Y!=0]*log(Y[Y!=0]/mu[Y!=0])
        SSGL_dev[l] = 2*sum(term1-(Y-mu))
      }
    
      # Compute GIC
      GIC[l] = SSGL_dev[l]/n + (log(log(n))*log(p)*length(which(beta[,l]!=0)))/n
    }
    min_GIC = GIC[which.min(GIC)]
    min_GIC_index = max(which(GIC==min_GIC))
    lambda0_GIC_min = lambda0[min_GIC_index]
  }
  
  #####################
  #####################
  ### Return a list ###
  #####################
  #####################
  if(return_GIC==TRUE){
    SSGL_output <- list(lambda0=lambda0,
                        beta=beta,
                        beta0=beta0,
                        classifications=classifications,
                        Y_pred=Y_pred,
                        GIC=GIC,
                        lambda0_GIC_min=lambda0_GIC_min,
                        min_GIC_index=min_GIC_index) 
  } else {
    SSGL_output <- list(lambda0=lambda0,
                        beta=beta,
                        beta0=beta0,
                        classifications=classifications,
                        Y_pred=Y_pred) 
  }
  # Return list
  return(SSGL_output)
}