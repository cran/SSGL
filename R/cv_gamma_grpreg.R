############################################
############################################
## FUNCTION FOR IMPLEMENTING K-FOLD CROSS ##
## VALIDATION IN GROUP-REGULARIZED        ##
## GAMMA REGRESSION MODELS                ##
############################################
############################################

# This function implements K-fold cross-validation for group-regularized
# gamma regression models. We employ the least squares approximation 
# approach of Wang and Leng (2007), and hence, the program does not allow 
# for the number of columns of X to be greater than sample size.


# INPUTS:
# Y = n x 1 vector of responses (y_1, ...., y_n) for training data
# X = n x p design matrix for training data
# groups = p x 1 vector of group indices or factor level names for each of the p columns of X.
# gamma_shape = known shape parameter nu in Gamma(mu_i, nu) distribution for the responses.
#               Default is gamma_shape=1.
# penalty = group regularization method to apply. Options are "gLASSO" for group lasso,
#           "gSCAD" for group SCAD, and "gMCP" for group MCP. Gamma regression with
#           the SSGL penalty is available in the stand-alone SSGL function.
# n_folds = number of folds in K-fold cross-validation. Default is n_folds=10
# group_weights = group-specific weights for the penalty. Default is to use the square roots of the 
#                 group sizes
# taper = tapering term in group SCAD and group MCP controlling how rapidly the penalty 
#         tapers off. Default is taper=4 for group SCAD and taper=3 for group MCP. This is ignored 
#         if "gLASSO" is specified as the penalty.
# n_lambda = number of tuning parameters to use. Default is 100
# lambda = a grid of tuning parameters. If the user does not specify this, then the program
#          chooses a grid automatically
# max_iter = maximum number of iterations. Default is 10,000
# tol = convergence criteria. Default is 1e-4

# OUTPUT:
# lambda = grid of L lambda's in descending order.
# cve = Lx1 vector of mean cross-validation error across all K folds. The kth entry in cve corresponds
#       to the kth regularization parameter in our lambda grid. The CVE on each of the K validation sets 
#       is the mean gamma regression deviance evaluated on that set.
# cvse = Lx1 vector of standard errors for cross-validation error across all K folds. 
#        The kth entry in cvse corresponds to the kth regularization parameter in our lambda grid.
# lambda_min = value of lambda that minimizes mean cross-validation error.
# min_index = index of lambda_min in lambda

cv_gamma_grpreg = function(Y, X, groups, gamma_shape=1, penalty=c("gLASSO","gSCAD","gMCP"),
                           n_folds=10, group_weights, taper, n_lambda=100, lambda, max_iter=10000, tol=1e-4) {

  ##################
  ##################
  ### PRE-CHECKS ###
  ##################
  ##################
  
  ## Coercion
  penalty = match.arg(penalty)
  
  ## Enumerate groups if not already done
  groups = as.factor(groups)
  group_numbers = as.numeric(groups)
  
  ## Number of groups and covariates overall
  X = as.matrix(X)
  G = length(unique(group_numbers))
  n = dim(X)[1]
  p = dim(X)[2]
  
  ## Check that p is less than or equal to (n_folds-1)/n_folds * n
  if(p > (n_folds-1)/n_folds*n) {
    stop("For cross-validation in group-regularized gamma regression, we require the total
          number of covariates to be less than or equal to sample size*(n_folds-1)/n_folds. 
          Consider reducing the number of covariates.")
  }
  
  ## Set taper parameter if not specified.
  if(missing(taper)){
    if(penalty=="gSCAD") taper=4
    if(penalty=="gMCP") taper=3
  }
  
  ## Set group-specific weights
  if(missing(group_weights))
    group_weights = sqrt(as.vector(table(group_numbers)))
  
  ## Check that dimensions are conformable
  if(length(Y) != dim(X)[1])
    stop("Non-conformable dimensions of Y and X.")
  
  ## Check that taper parameter is greater than 2 for gSCAD and greater than 1 for gMCP
  if(penalty=="gSCAD"){
    if(taper<=2)
      stop("The taper parameter must be greater than 2 for the group SCAD penalty.")
  }
  if(penalty=="gMCP"){
    if(taper<=1)
      stop("The taper parameter must be greater than 1 for the group MCP penalty.")
  }
  
  ## Check that gamma.shape is strictly positive
  if (gamma_shape<=0)
    stop("Shape parameter for gamma density must be strictly positive.")
  ## Check that the response variables are strictly positive
  if(any(Y<=0))
    stop("All responses must be greater than zero.")
  ## Check that weights are all greater than or equal to 0
  if(!missing(group_weights)){
    if(!all(group_weights>=0))
      stop("All group-specific weights should be nonnegative.")
  }
  
  ## Number of lambdas
  if(n_lambda < 2)
      stop("For cross-validation, n_lambda should be at least 2.")
  ## If user specified lambda, check that all lambdas are greater than 0
  if(!missing(lambda)) {
    n_lambda = length(lambda) # Override n_lambda with the length of lambda

    if (!all(lambda>0))
      stop("All lambdas should be strictly positive.")
  }  
  
  ## If lambda is not specified
  if(missing(lambda)) {
    # Compute max_lambda
    norms_vec = rep(0, G)
    for(g in 1:G){
      ind_g = which(group_numbers==g)
      X_g = X[,ind_g]
      norms_vec[g] = sqrt(sum((t(X_g)%*%(Y-1))^2))/sqrt(length(ind_g))  
    }
    max_lambda = max(abs(norms_vec))/n
    eps = 0.001
    
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
  
  # Sort the lambdas just in case
  lambda = sort(lambda, decreasing=TRUE)
  
  ## Randomly shuffle the data
  new_order = sample(n)
  X_new = X[new_order,]
  Y_new = Y[new_order]
  
  ## Create K equally-sized folds
  folds = cut(seq(1,n), breaks=n_folds, labels=FALSE)
  
  ## To store the cross-validation error
  folds_cve = matrix(0, n_folds, n_lambda)
  
  for(k in 1:n_folds){
    
    ## Indices for validation set
    val_ind = which(folds==k,arr.ind=TRUE)
    
    Y_val = Y_new[val_ind] 
    X_val = X_new[val_ind, ]
    
    ## Indices for training set
    Y_train = Y_new[-val_ind]
    X_train = X_new[-val_ind, ]
    
    ## Train model on training set
    gamma_mod_train = gamma_grpreg(Y=Y_train, X=X_train, groups=group_numbers, penalty=penalty, 
                                   X_test=X_val, gamma_shape=gamma_shape,
                                   group_weights=group_weights, taper=taper, lambda=lambda,
                                   n_lambda=n_lambda, max_iter=max_iter, tol=tol, return_GIC=FALSE)
    
    ## Compute cross-validation error, which is the deviance
    for(l in 1:n_lambda){
      ## Store CVE for all lambdas in the kth row of folds_cve
      mu = gamma_mod_train$Y_pred[,l]
      folds_cve[k,l] = 2*sum(-log(Y_val/mu)+(Y_val-mu)/mu)
    }
  }
  
  ## Mean cross-validation error
  cve = colMeans(folds_cve)
  ## CVE standard error
  cvse = apply(folds_cve,2,stats::sd)
  ## Lambda which minimizes cve
  min_index = which.min(cve)
  lambda_min = lambda[min_index]
  
  ## Return a list
 cv_gamma_grpreg_output <- list(lambda=lambda,
                                cve=cve,
                                cvse=cvse,
                                lambda_min=lambda_min,
                                min_index=min_index)
  # Return list
  return(cv_gamma_grpreg_output)
}