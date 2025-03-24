\name{SSGL}
\alias{SSGL}

\title{Spike-and-Slab Group Lasso for Group-Regularized Generalized Linear Models (GLMs)}

\description{The \code{SSGL} function implements maximum a posteriori (MAP) estimation for group-regularized GLMs with the spike-and-slab group lasso (SSGL) penalty of Bai et al. (2022) and Bai (2023). The identity link function is used for Gaussian regression, the logit link is used for binomial regression, and the log link is used for Poisson regression. If the covariates in each \eqn{x_i} are grouped according to known groups \eqn{g=1, ..., G}, then this function can estimate some of the \eqn{G} groups of coefficients as all zero, depending on the amount of regularization. 

This function only returns point estimates. Please refer to the \code{SSGL_gibbs} function if uncertainty quantification of the model parameters is desired. In general, we recommend using \code{SSGL} for estimation and variable selection and \code{SSGL_gibbs} for uncertainty quantification. 

 The \code{SSGL} function also has the option of returning the generalized information criterion (GIC) of Fan and Tang (2013) for each regularization parameter in the grid \code{lambda0}. The GIC can be used for model selection and serves as a useful alternative to cross-validation. The formula for the GIC and a given \eqn{\lambda_0} is
    
      \deqn{DIC(\lambda_0) = \frac{1}{n} Deviance_{\lambda_0} + a_n \times \nu),}

where \eqn{Deviance_{\lambda_0}} is the deviance computed with the estimate of \code{beta} based on spike hyperparameter \eqn{\lambda_0}, \eqn{\nu_0} is the number of nonzero elements in the estimated \code{beta}, and \eqn{a_n} is a sequence that diverges at a suitable rate relative to \eqn{n}. As recommended by Fan and Tang (2013), we set \eqn{a_n = \{\log(\log(n))\}\log(p)}. 

If cross-validation is preferred for tuning \eqn{\lambda_0}, please refer to the \code{SSGL_cv} function.}

\usage{
SSGL(Y, X, groups, family=c("gaussian","binomial","poisson"), 
     X_test, group_weights, n_lambda0=25, 
     lambda0, lambda1=1, a=1, b=length(unique(groups)), 
     max_iter=100, tol = 1e-6, return_GIC=TRUE, print_lambda0=TRUE) 
}

\arguments{
  \item{Y}{\eqn{n \times 1} vector of responses for training data.}
  \item{X}{\eqn{n \times p} design matrix for training data, where the \eqn{j}th column of \code{X} corresponds to the \eqn{j}th overall covariate.}
  \item{groups}{\eqn{p}-dimensional vector of group labels. The \eqn{j}th entry in \code{groups} should contain either the group number \emph{or} the factor level name that the feature in the \eqn{j}th column of \code{X} belongs to. \code{groups} must be either a vector of integers or factors.}
  \item{family}{exponential dispersion family of the response variables. Allows for \code{"gaussian"}, \code{"binomial"}, and \code{"poisson"}.}
  \item{X_test}{\eqn{n_{test} \times p} design matrix for test data to calculate predictions. \code{X_test} must have the \emph{same} number of columns as X, but not necessarily the same number of rows. If \emph{no} test data is provided or if in-sample predictions are desired, then the function automatically sets \code{X_test=X} in order to calculate \emph{in-sample} predictions.}
  \item{group_weights}{group-specific, nonnegative weights for the penalty. Default is to use the square roots of the group sizes.}
  \item{n_lambda0}{number of spike hyperparameters \eqn{L}. Default is \code{n_lambda0=25}.}
  \item{lambda0}{grid of \eqn{L} spike hyperparameters \eqn{\lambda_0}. The user may specify either a scalar or a vector. If the user does not provide this, the program chooses the grid automatically.}
  \item{lambda1}{slab hyperparameter \eqn{\lambda_1} in the SSGL prior. Default is \code{lambda1=1}.}
  \item{a}{shape hyperparameter for the \eqn{Beta(a,b)} prior on the mixing proportion in the SSGL prior. Default is \code{a=1}.}
  \item{b}{shape hyperparameter for the \eqn{Beta(a,b)} prior on the mixing proportion in the SSGL prior. Default is \code{b=length(unique(groups))}, i.e. the number of groups.}
  \item{max_iter}{maximum number of iterations in the algorithm. Default is \code{max_iter=100}.}
  \item{tol}{convergence threshold for algorithm. Default is \code{tol=1e-6}.}
  \item{return_GIC}{Boolean variable for whether or not to return the GIC. Default is \code{return_GIC=TRUE}.}
  \item{print_lambda0}{Boolean variable for whether or not to print the current value in \code{lambda0}. Default is \code{print_lambda0=TRUE}.}
}

\value{The function returns a list containing the following components:
    \item{lambda0}{\eqn{L \times 1} vector of spike hyperpameters \code{lambda0} used to fit the model. \code{lambda0} is displayed in descending order.}
  \item{beta}{\eqn{p \times L} matrix of estimated regression coefficients. The \eqn{k}th column in \code{beta} corresponds to the \eqn{k}th spike hyperparameter in \code{lambda0}.}
  \item{beta0}{\eqn{L \times 1} vector of estimated intercepts. The \eqn{k}th entry in \code{beta0} corresponds to the \eqn{k}th spike hyperparameter in \code{lambda0}.}
    \item{classifications}{\eqn{G \times L} matrix of classifications, where \eqn{G} is the number of groups. An entry of "1" indicates that the group was classified as nonzero, and an entry of "0" indicates that the group was classified as zero. The \eqn{k}th column of \code{classifications} corresponds to the \eqn{k}th spike hyperparameter in \code{lambda0}.}
  \item{Y_pred}{\eqn{n_{test} \times L} matrix of predicted mean response values \eqn{\mu_{test} = E(Y_{test})} based on the test data in \code{X_test} (or training data \code{X} if no argument was specified for \code{X_test}). The \eqn{k}th column in \code{Y_pred} corresponds to the predictions for the \eqn{k}th spike hyperparameter in \code{lambda0}.}
    \item{GIC}{\eqn{L \times 1} vector of GIC values. The \eqn{k}th entry of \code{GIC} corresponds to the \eqn{k}th entry in our \code{lambda0} grid. This is not returned if \code{return_GIC=FALSE}.}
  \item{lambda0_GIC_min}{The value in \code{lambda0} that minimizes \code{GIC}. This is not returned if \code{return_GIC=FALSE}.}
    \item{min_GIC_index}{The index of \code{lambda0_GIC_min} in \code{lambda0}. This is not returned if \code{return_GIC=FALSE}.}
}

\references{
Bai, R. (2023). "Bayesian group regularization in generalized linear models with a continuous spike-and-slab prior." \emph{arXiv pre-print arXiv:2007.07021}.

Bai, R., Moran, G. E., Antonelli, J. L., Chen, Y., and Boland, M.R. (2022). "Spike-and-slab group lassos for grouped regression and sparse generalized additive models." \emph{Journal of the American Statistical Association}, \bold{117}:184-197.

Fan, Y. and Tang, C. Y. (2013). "Tuning parameter selection in high dimensional penalized likelihood." \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, \bold{75}:531-552.
}

\examples{
## Generate data
set.seed(12345)
X = matrix(runif(100*10), nrow=100)
n = dim(X)[1]
groups = c("A","A","A","B","B","B","C","C","D","D")
groups = as.factor(groups)
beta_true = c(-2.5,1.5,1.5,0,0,0,2,-2,0,0)

## Generate responses from Gaussian distribution
Y = crossprod(t(X), beta_true) + rnorm(n)

## Generate test data
n_test = 50
X_test = matrix(runif(n_test*10), nrow=n_test)

## Fit SSGL model with 10 spike hyperparameters
## NOTE: If you do not specify lambda0, the program will automatically choose a suitable grid.
SSGL_mod = SSGL(Y, X, groups, family="gaussian", X_test, lambda0=seq(from=50,to=5,by=-5))

## Regression coefficient estimates
SSGL_mod$beta

## Predicted n_test-dimensional vectors mu=E(Y.test) based on test data, X_test. 
## The kth column of 'Y_pred' corresponds to the kth entry in 'lambda.'
SSGL_mod$Y_pred 

## Classifications of the 8 groups. The kth column of 'classifications'
## corresponds to the kth entry in 'lambda.'
SSGL_mod$classifications

## Plot lambda vs. GIC
plot(SSGL_mod$lambda0, SSGL_mod$GIC, type='l')

## Model selection with the lambda that minimizes GIC
SSGL_mod$lambda0_GIC_min
SSGL_mod$min_GIC_index 
SSGL_mod$classifications[, SSGL_mod$min_GIC_index]
SSGL_mod$beta[, SSGL_mod$min_GIC_index]

\donttest{
## Example with Poisson regression

## Generate data
set.seed(1234)
X = matrix(runif(100*10), nrow=100)
n = dim(X)[1]
groups = c("A","A","A","B","B","B","C","C","D","D")
groups = as.factor(groups)
beta_true = c(-2.5,1.5,1.5,0,0,0,2,-2,0,0)

## Generate count responses
eta = crossprod(t(X), beta_true)
Y = rpois(n, exp(eta))

## Generate test data
n_test = 50
X_test = matrix(runif(n_test*10), nrow=n_test)

## Fit SSGL model
SSGL_poisson_mod = SSGL(Y, X, groups, family="poisson")

## Regression coefficient estimates
SSGL_poisson_mod$beta

## Predicted n_test-dimensional vectors mu=E(Y.test) based on test data, X_test. 
## The kth column of 'Y_pred' corresponds to the kth entry in 'lambda.'
SSGL_poisson_mod$Y_pred 

## Classifications of the 8 groups. The kth column of 'classifications'
## corresponds to the kth entry in 'lambda.'
SSGL_poisson_mod$classifications

## Plot lambda vs. GIC
plot(SSGL_poisson_mod$lambda0, SSGL_poisson_mod$GIC, type='l')

## Model selection with the lambda that minimizes GIC
SSGL_poisson_mod$lambda0_GIC_min
SSGL_poisson_mod$min_GIC_index 
SSGL_poisson_mod$classifications[, SSGL_mod$min_GIC_index]
SSGL_poisson_mod$beta[, SSGL_mod$min_GIC_index]
}
}