#################################################
#################################################
## FUNCTION FOR IMPLEMENTING GROUP-REGULARIZED ##
## NEGATIVE BINOMIAL REGRESSION MODELS         ##
#################################################
#################################################

# This function implements group-regularized negative binomial regression models. This
# program employs the least squares approximation approach of Wang and Leng (2007) and
# hence does not allow for the number of columns of X to be greater than sample size.
# The size parameter alpha is automatically computed using the MLE of the non-penalized
# problem.

# INPUTS:
# Y = n x 1 vector of responses (y_1, ...., y_n) for training data
# X = n x p design matrix for training data
# groups = p x 1 vector of group indices or factor level names for each of the p columns of X.
# X_test = n_test x J design matrix for test data. If missing, then the program sets X_test=X
#          and computes in-sample predictions on training data X. X_test must have the same 
#          number of columns as X, but not necessarily the same number of rows.
# nb_size = known size parameter alpha in NB(alpha, mu_i) distribution for the responses.
#           Default is nb_size=1
# penalty = group regularization method to apply. Options are "gLASSO" for group lasso,
#           "gSCAD" for group SCAD, and "gMCP" for group MCP. Negative binomial regression for
#           the SSGL penalty is available in the stand-alone SSGL function.
# group_weights = group-specific weights for the penalty. Default is to use the square roots of the 
#           group sizes
# taper = tapering term in group SCAD and group MCP controlling how rapidly the penalty 
#         tapers off. Default is taper=4 for group SCAD and taper=3 for group MCP. This is ignored 
#         if "gLASSO" is specified as the penalty.
# n_lambda = number of tuning parameters to use. Default is 100
# lambda = a grid of tuning parameters. If the user does not specify this, then the program
#          chooses a grid automatically
# max_iter = maximum number of iterations. Default is 10,000
# tol = convergence criteria. Default is 1e-4
# return_GIC = Boolean variable for whether or not to return the GIC. Default is return_GIC=TRUE

# OUTPUT:
# lambda = grid of L lambda's in descending order.
# beta = p x L matrix of regression coefficient estimates. The lth column of beta corresponds to the
#        lth entry in our lambda grid.
# beta0 = Lx1 vector of intercept estimates. The lth entry of beta0 corresponds to the lth entry in
#         our lambda grid.
# classifications = G x L matrix of group classifications, where G is the number of groups. "1" indicates
#                   that the group was selected and "0" indicates that the group was not selected. 
#                   The lth column in this matrix corresponds to the lth entry in our lambda grid.
# Y_pred = n_test x L matrix of predicted mean response values based on test data in X_test. If
#          X_test was left blank or X_test=X, then in-sample predictions on X are returned.
# GIC = Lx1 vector of GIC values. The lth entry of GIC corresponds to the lth entry in our lambda grid.
#       Not returned if return_GIC=FALSE. 
# lambda_min = value of lambda that minimizes the GIC. Not returned if return_GIC=FALSE.
# min_index = the index of lambda_min in lambda. Not returned if return_GIC=FALSE.

nb_grpreg = function(Y, X, groups, X_test, nb_size=1, penalty=c("gLASSO","gSCAD","gMCP"),
                     group_weights, taper, n_lambda=100, lambda, max_iter=10000, tol=1e-4,
                     return_GIC=TRUE) {

  ##################
  ##################
  ### PRE-CHECKS ###
  ##################
  ##################
  
  ## Coercion
  penalty = match.arg(penalty)
  if(penalty=="gLASSO"){
    penalty="grLasso"
  } else if(penalty=="gSCAD"){ 
    penalty="grSCAD"
  } else if(penalty=="gMCP"){
    penalty="grMCP"
  }
  
  ## Enumerate groups if not already done
  groups = as.factor(groups)
  group_numbers = as.numeric(groups)
  
  ## Number of groups and covariates overall
  G = length(unique(group_numbers))
  n = dim(X)[1]
  p = dim(X)[2]
  
  ## Check that p is less than or equal to n
  if(p > n) {
    stop("For group-regularized negative binomial regression, we require the
          total number of covariates to be less than or equal to sample size. 
          Consider reducing the number of covariates.")
  }
  
  ## Set test data as training data if test data not provided.
  X = as.matrix(X)
  if(missing(X_test)) X_test = X
  n_test = dim(X_test)[1]
  
  ## Set taper parameter if not specified.
  if(missing(taper)){
    if(penalty=="grSCAD") taper=4
    if(penalty=="grMCP") taper=3
  }
  
  ## Set group-specific weights
  if(missing(group_weights))
    group_weights = sqrt(as.vector(table(group_numbers)))

  ## Check that dimensions are conformable
  if(length(Y) != dim(X)[1])
    stop("Non-conformable dimensions of Y and X.")
  ## Check that X and X_test have the same number of columns
  if(dim(X_test)[2] != p)
    stop("X and X_test should have the same number of columns.")
  
  ## Check that taper parameter is greater than 2 for gSCAD and greater than 1 for gMCP
  if(penalty=="grSCAD"){
    if(taper<=2)
      stop("The taper parameter must be greater than 2 for the group SCAD penalty.")
  }
  if(penalty=="grMCP"){
    if(taper<=1)
      stop("The taper parameter must be greater than 1 for the group MCP penalty.")
  }
  
  ## Check that nb_size is strictly positive
  if (nb_size<=0)
    stop("Size parameter for negative binomial density must be strictly positive.")
  ## Check that the response variables are strictly nonnegative whole numbers
  if(any(Y<0))
    stop("All counts must be greater than or equal to zero.")
  if(!all(Y==floor(Y)))
    stop("All counts must be whole numbers.")
  ## Check that weights are all greater than or equal to 0
  if(!missing(group_weights)){
    if(!all(group_weights>=0))
      stop("All group-specific weights should be nonnegative.")
  }
  
  ## Number of lambdas
  if(n_lambda < 1)
      stop("The number of lambdas must be at least one.")
  ## If user specified lambda, check that all lambdas are greater than 0
  if(!missing(lambda)) {
    n_lambda = length(lambda) # Override n_lambda with the length of lambda

    if (!all(lambda>0))
      stop("All lambdas should be strictly positive.")
  }  
  
  ## If lambda is not specified
  if(missing(lambda)) {
    max_lambda = 1
    eps = .05
    
    ## Create grid of lambdas
    if(n_lambda==1){ 
      lambda = max_lambda*eps 
    } else if(n_lambda > 1) { 
      lambda = rep(0, n_lambda)
      lambda[1] = max_lambda
      lambda[n_lambda] = max_lambda*eps
      
      if(n_lambda >= 3){
        for(l in 2:(n_lambda-1)){
          ## equispaced lambdas on log scale
          loglambda = log(lambda[1])-(l-1)*((log(lambda[1])-log(lambda[n_lambda]))/(n_lambda-1))
          lambda[l] = exp(loglambda)
        }
      }
    }
  }
  
  # order the lambdas just in case
  lambda = sort(lambda, decreasing=TRUE)
  
  # Find MLE of regression coefficients
  glm_mod = stats::glm(Y~X, family=MASS::negative.binomial(nb_size))
  nb_MLE = glm_mod$coefficients # MLE of regression coefficients 
  
  # Find asymptotic variance-covariances
  MLE_hessian = stats::vcov(glm_mod, type="hessian")
  
  # Square root of Sigma_inv
  sqrt_Sigma_inv = pracma::sqrtm(solve(MLE_hessian))$B
  
  # New parameter Y_check
  Y_check = sqrt_Sigma_inv %*% nb_MLE
  
  # Fit group regression. Do not penalize the intercept
  
  groups_MLE = group_numbers + 1
  groups_MLE = c(1, groups_MLE)
  weights_MLE = c(0, group_weights)
  
  ## Solve for beta
  ## LSA approach
  if(n_lambda > 1){
      nb_group = grpreg::grpreg(X=sqrt_Sigma_inv, y=Y_check, group=groups_MLE, penalty=penalty, 
                                  nlambda=n_lambda, lambda=lambda, eps=tol, max.iter=max_iter,
                                  group.multiplier = weights_MLE)
  } else if (n_lambda == 1){
      nb_group = grpreg::grpreg(X=sqrt_Sigma_inv, y=Y_check, group=groups_MLE, penalty=penalty, 
                                  lambda=lambda, eps=tol, max.iter=max_iter,
                                  group.multiplier = weights_MLE)
  }
    
  
  ## Estimated intercepts
  beta0 = nb_group$beta[2,] 
  
  ## Estimated regression coefficients
  beta = as.matrix(nb_group$beta[-c(1,2),])
  colnames(beta) = NULL
  rownames(beta) = groups
  
  ## Compute predictions matrix on test data
  L = n_lambda
  Y_pred = matrix(0, n_test, L) # To compute predictions on test set
  for(l in 1:n_lambda){
      Y_pred[,l] = exp(beta0[l]+X_test%*%beta[,l])
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
  
  ## if GIC is desired
  if(return_GIC==TRUE){
    ## Compute GIC
    nbreg_dev = rep(0, L)
    GIC = rep(0, L)
    
    for(l in 1:L){
      # Compute deviance first
      mu = exp(beta0[l]+X%*%beta[,l])
      term1 = rep(0, length(Y))
      term1[Y!=0] = Y[Y!=0]*log(Y[Y!=0]/mu[Y!=0])
      nbreg_dev[l] = 2*sum(term1-(Y+nb_size)*log((1+Y/nb_size)/(1+mu/nb_size)))
      # Compute GIC
      GIC[l] = nbreg_dev[l]/n + (log(log(n))*log(p)*length(which(beta[,l]!=0)))/n
    }
    min_GIC = GIC[which.min(GIC)]
    min_index = max(which(GIC==min_GIC))
    lambda_min = lambda[min_index]
  }
  
  ## Return a list
  if(return_GIC==TRUE){
    nb_grpreg_output <- list(lambda=lambda,
                             beta=beta,
                             beta0=beta0,
                             classifications=classifications,
                             Y_pred=Y_pred,
                             GIC=GIC,
                             lambda_min=lambda_min,
                             min_index=min_index) 
  } else {
    nb_grpreg_output <- list(lambda=lambda,
                             beta=beta,
                             beta0=beta0,
                             classifications = classifications,
                             Y_pred=Y_pred)
  }
  # Return list
  return(nb_grpreg_output)
}