#################################################
#################################################
## FUNCTION FOR IMPLEMENTING GROUP-REGULARIZED ##
## GAMMA REGRESSION MODELS                     ##
#################################################
#################################################

# This function implements K-fold cross-validation for group-regularized regression
# models in the exponential dispersion family with the spike-and-slab group lasso
# (SSGL) penalty.

# INPUTS:
# Y = n x 1 vector of responses (y_1, ...., y_n) for training data
# X = n x p design matrix for training data
# groups = p x 1 vector of group indices or factor level names for each of the p columns of X.
# family = the exponential disperison family. Currently allows "gaussian", "binomial", "poisson," 
#          "negativebinomial", and "gamma". Default is "gaussian" for non-binary responses and
#          "binomial" for binary responses, but "poisson" or "negativebinomial" may be more
#          appropriate for count data and "gamma" may be more appropriate for right-skewed,
#          strictly positive, continuous data.
# nb_size = known shape parameter for negative binomial regression. Default is 1
# gamma_shape = known shape parameter for gamma regression. Default is 1
# group_weights = group-specific weights. Default is to use the square roots of the group 
# n_folds = number of folds to use in K-fold cross-validation. Default is n_folds=5.
# n_lambda0 = number of lambda0's in our grid. Default is 20.
# lambda0 = a grid of values for the spike hyperparameter. Note that the program automatically orders 
#           these in descending order, so the solution path is from least sparse to most sparse. 
#           If lambda0 is not specified, the program chooses a default equispaced grid from 100 to 5.
# lambda1 = a fixed small value for the slab hyperparameter. Default is 1
# a = shape hyperparameter in Beta(a,b) prior on mixing proportion. Default is 1
# b = shape hyperparameter in Beta(a,b) prior on mixing proportion. Default is G for number of groups
# max_iter = maximum number of iterations. Default is 10,000
# tol = convergence criteria. Default is 1e-6
# print_fold = Boolean variable whether to print the current fold. Default is TRUE

# OUTPUT:
# lambda0 = grid of L spike hyperparameters lambda0's in descending order.
# cve = Lx1 vector of mean cross-validation error across all K folds. The kth entry in cve corresponds
#       to the kth regularization parameter in our lambda grid. The CVE on each of the K validation sets 
#       is the mean GLM deviance evaluated on that set.
# cvse = Lx1 vector of standard errors for cross-validation error across all K folds. 
#        The kth entry in cvse corresponds to the kth spike hyperparameter in our lambda0 grid.
# lambda0_min = value of lambda0 that minimizes mean cross-validation error.
# min_index = index of lambda0_min in lambda0

cv_SSGL = function(Y, X, groups, 
                   family=c("gaussian","binomial","poisson","negativebinomial","gamma"), 
                   nb_size=1, gamma_shape=1, group_weights, n_folds=5, n_lambda0=25,
                   lambda0, lambda1=1, a=1, b=dim(X)[2], 
                   max_iter=100, tol=1e-6, print_fold=TRUE) {

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

  ## Check that dimensions are conformable
  if(length(Y) != dim(X)[1])
    stop("Non-conformable dimensions of Y and X.")  
  
  ## Coercion
  if(missing(family)){
    if(all(Y==1 || Y==0)){ 
      family=="binomial"
    } else {
      family=="gaussian"
    }
  }
  family = match.arg(family)
  
  ## Check that the data can be used for the respective family
  if(family=="poisson" || family=="negativebinomial"){
    if(any(Y<0))
      stop("All counts y must be greater than or equal to zero.")
    if(any(Y-floor(Y)!=0))
      stop("All counts must be whole numbers.")
  }
  if(family=="negativebinomial"){
    ## Check that nb_size is strictly positive
    if (nb_size<=0)
      stop("Size parameter for negative binomial density must be strictly positive.")
    ## Check that p is less than or equal to n
    if(p > n) {
      stop("For group-regularized negative binomial regression, we require the
          total number of covariates to be less than or equal to sample size. 
          Consider reducing the number of covariates.")
    }
  }
  if(family=="binomial"){
    if(any(Y<0))
      stop("All binary responses must be either '0' or '1.'")
    if(any(Y>1))
      stop("All binary responses must be either '0' or '1.'")
    if(any(Y-floor(Y)!=0))
      stop("All binary responses must be either '0' or '1.'")
  }
  if(family=="gamma"){
    if(any(Y<=0))
      stop("All responses must be strictly positive.")
    if(gamma_shape<=0)
      stop("Shape parameter for gamma density must be strictly positive.")
    ## Check that p is less than or equal to n
    if(p > n) {
      stop("For group-regularized gamma regression, we require the total
          number of covariates to be less than or equal to sample size. 
          Consider reducing the number of covariates.")
    }
  }
  
  ## Set group-specific weights
  if(missing(group_weights))
    group_weights = sqrt(as.vector(table(group_numbers)))

  ## Check that weights are all greater than or equal to 0
  if(!missing(group_weights)){
    if(!all(group_weights>=0))
      stop("All group-specific weights should be nonnegative.")
  }
  
  ## Number of lambdas
  if(n_lambda0 < 2)
    stop("For cross-validation, n_lambda0 should be at least 2.")  
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
        if(family=="poisson" || family=="negativebinomial"){
          lambda0 = seq(from=max_lambda, to=lambda1+5, length=n_lambda0)
        } else if(family=="gaussian" || family=="gamma"){
          lambda0 = seq(from=max_lambda, to=lambda1+1, length=n_lambda0)
        }
      }
    }
  }
  
  # Sort lambda0 just in case
  lambda0 = sort(lambda0, decreasing=TRUE)
  
  ## Randomly shuffle the data
  new_order = sample(n)
  X_new = X[new_order,]
  Y_new = Y[new_order]
  
  ## Create K equally-sized folds
  folds = cut(seq(1,n), breaks=n_folds, labels=FALSE)
  
  ## To store the cross-validation error
  folds_cve = matrix(0, n_folds, n_lambda0)
  
  for(k in 1:n_folds){
    
    ## To output iteration
    if(print_fold==TRUE)
      cat("Fold number", k, "\n")
    
    ## Indices for validation set
    val_ind = which(folds==k,arr.ind=TRUE)
    
    Y_val = Y_new[val_ind] 
    X_val = X_new[val_ind, ]
    
    ## Indices for training set
    Y_train = Y_new[-val_ind]
    X_train = X_new[-val_ind, ]
    
    ## Train model on training set
    SSGL_mod_train = SSGL(Y=Y_train, X=X_train, groups=group_numbers, X_test=X_val, family=family, 
                          nb_size=nb_size, gamma_shape=gamma_shape, group_weights=group_weights,
                          n_lambda0=n_lambda0, lambda0=lambda0, lambda1=lambda1,
                          a=a, b=b, max_iter=max_iter, tol=tol, 
                          return_GIC=FALSE, print_lambda0=FALSE)
    
    # Store CVE, i.e. the deviance, for all lambdas in the kth row of folds_cve
    for(l in 1:n_lambda0){
        if(family=="gaussian"){
          folds_cve[k,l] = sum((Y_val-SSGL_mod_train$Y_pred[,l])^2)
        } 
        if(family=="binomial"){
          mu = SSGL_mod_train$Y_pred[,l]
          term1 = sum(log(mu[Y_val==1]))
          term2 = sum(log(1-mu[Y_val==0]))
          folds_cve[k,l] = -2*(term1+term2)
        }
        if(family=="poisson"){
          mu = SSGL_mod_train$Y_pred[,l]
          term1 = rep(0, length(Y_val))
          term1[Y_val!=0] = Y_val[Y_val!=0]*log(Y_val[Y_val!=0]/mu[Y_val!=0])
          folds_cve[k,l] = 2*sum(term1-(Y_val-mu))
        }
        if(family=="negativebinomial"){
          mu = SSGL_mod_train$Y_pred[,l]
          term1 = rep(0, length(Y_val))
          term1[Y_val!=0] = Y_val[Y_val!=0]*log(Y_val[Y_val!=0]/mu[Y_val!=0])
          folds_cve[k,l] = 2*sum(term1-(Y_val+nb_size)*log((1+Y_val/nb_size)/(1+mu/nb_size)))
        } 
        if(family=="gamma"){
          mu = SSGL_mod_train$Y_pred[,l]
          folds_cve[k,l] = 2*sum(-log(Y_val/mu)+(Y_val-mu)/mu)
        }
    } 
  }
  
  ## Mean cross-validation error
  cve = colMeans(folds_cve)
  ## CVE standard error
  cvse = apply(folds_cve,2,stats::sd)
  ## Lambda which minimizes cve
  min_index = which.min(cve)
  lambda0_min = lambda0[min_index]
  
  
  ## Return a list
 cv_SSGL_output <- list(lambda0=lambda0,
                        cve=cve,
                        cvse=cvse,
                        lambda0_min=lambda0_min,
                        min_index=min_index)
  # Return list
  return(cv_SSGL_output)
}