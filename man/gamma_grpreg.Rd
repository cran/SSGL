\name{gamma_grpreg}
\alias{gamma_grpreg}

\title{Group-regularized Gamma Regression}

\description{This function implements group-regularized gamma regression with a known shape parameter \eqn{\nu} and the log link. In gamma regression, we assume that \eqn{y_i \sim Gamma(\mu_i, \nu)}, where
  
  \deqn{f(y_i | \mu_i, \nu ) = \frac{1}{\Gamma(\nu)} (\frac{\nu}{\mu_i})^{\nu} \exp(-\frac{\nu}{\mu_i}y_i) y_i^{\nu-1}, y > 0.}
  
  Then \eqn{E(y_i) = \mu_i}, and we relate \eqn{\mu_i} to a set of \eqn{p} covariates \eqn{x_i} through the log link,
  
  \deqn{\log(\mu_i) = \beta_0 + x_i^T \beta, i=1,..., n}
  
  If the covariates in each \eqn{x_i} are grouped according to known groups \eqn{g=1, ..., G}, then this function can estimate some of the \eqn{G} groups of coefficients as all zero, depending on the amount of regularization. Our implementation for regularized gamma regression is based on the least squares approximation approach of Wang and Leng (2007), and hence, the function does not allow the total number of covariates \eqn{p} to be greater than sample size.
  
  In addition, this function has the option of returning the generalized information criterion (GIC) of Fan and Tang (2013) for each regularization parameter in the grid \code{lambda}. The GIC can be used for model selection and serves as a useful alternative to cross-validation. 
}

\usage{
gamma_grpreg(Y, X, groups, X_test, gamma_shape=1, 
             penalty=c("gLASSO","gSCAD","gMCP"), 
             group_weights, taper, n_lambda=100, lambda, 
             max_iter=10000, tol=1e-4, return_GIC=TRUE)
}

\arguments{
  \item{Y}{\eqn{n \times 1} vector of strictly positive, continuous responses for training data.}
  \item{X}{\eqn{n \times p} design matrix for training data, where the \eqn{j}th column corresponds to the \eqn{j}th overall feature.}
  \item{groups}{\eqn{p}-dimensional vector of group labels. The \eqn{j}th entry in \code{groups} should contain either the group number \emph{or} the factor level name that the feature in the \eqn{j}th column of \code{X} belongs to. \code{groups} must be either a vector of integers or factors.}
  \item{X_test}{\eqn{n_{test} \times p} design matrix for test data to calculate predictions. \code{X_test} must have the \emph{same} number of columns as \code{X}, but not necessarily the same number of rows. If \emph{no} test data is provided or if in-sample predictions are desired, then the function automatically sets \code{X_test=X} in order to calculate \emph{in-sample} predictions.}
  \item{gamma_shape}{known shape parameter \eqn{\nu} in \eqn{Gamma(\mu_i,\nu)} distribution for the responses. Default is \code{gamma_shape=1}.}
  \item{penalty}{group regularization method to use on the groups of regression coefficients. The options are \code{"gLASSO"}, \code{"gSCAD"}, \code{"gMCP"}. To implement gamma regression with the SSGL penalty, use the \code{SSGL} function.}
  \item{group_weights}{group-specific, nonnegative weights for the penalty. Default is to use the square roots of the group sizes.}
  \item{taper}{tapering term \eqn{\gamma} in group SCAD and group MCP controlling how rapidly the penalty tapers off. Default is \code{taper=4} for group SCAD and \code{taper=3} for group MCP. Ignored if \code{"gLASSO"} is specified as the penalty.}
  \item{n_lambda}{number of regularization parameters \eqn{L}. Default is \code{n_lambda=100}.}
  \item{lambda}{grid of \eqn{L} regularization parameters. The user may specify either a scalar or a vector. If the user does not provide this, the program chooses the grid automatically.}
  \item{max_iter}{maximum number of iterations in the algorithm. Default is \code{max_iter=10000}.}
  \item{tol}{convergence threshold for algorithm. Default is \code{tol=1e-4}.}
  \item{return_GIC}{Boolean variable for whether or not to return the GIC. Default is \code{return_GIC=TRUE}.}
}

\value{The function returns a list containing the following components:
    \item{lambda}{\eqn{L \times 1} vector of regularization parameters \code{lambda} used to fit the model. \code{lambda} is displayed in descending order.}
    \item{beta}{\eqn{p \times L} matrix of estimated regression coefficients. The \eqn{k}th column in \code{beta} corresponds to the \eqn{k}th regularization parameter in \code{lambda}.}
  \item{beta0}{\eqn{L \times 1} vector of estimated intercepts. The \eqn{k}th entry in \code{beta0} corresponds to the \eqn{k}th regularization parameter in \code{lambda}.}
  \item{classifications}{\eqn{G \times L} matrix of classifications, where \eqn{G} is the number of groups. An entry of "1" indicates that the group was classified as nonzero, and an entry of "0" indicates that the group was classified as zero. The \eqn{k}th column of \code{classifications} corresponds to the \eqn{k}th regularization parameter in \code{lambda}.}
  \item{Y_pred}{\eqn{n_{test} \times L} matrix of predicted mean response values \eqn{\mu_{test} = E(Y_{test})} based on the \emph{test} data in \code{X_test} (or training data \code{X} if no argument was specified for \code{X_test}). The \eqn{k}th column in \code{Y_pred} corresponds to the predictions for the \eqn{k}th regularization parameter in \code{lambda}.}
  \item{GIC}{\eqn{L \times 1} vector of GIC values. The \eqn{k}th entry of \code{GIC} corresponds to the \eqn{k}th entry in our \code{lambda} grid. This is not returned if \code{return_GIC=FALSE}.}
  \item{lambda_min}{The value in \code{lambda} that minimizes \code{GIC}. This is not returned if \code{return_GIC=FALSE}.}
    \item{min_index}{The index of \code{lambda_min} in \code{lambda}. This is not returned if \code{return_GIC=FALSE}.}
}

\references{
  Breheny, P. and Huang, J. (2015). "Group descent algorithms for nonconvex penalized linear and logistic regression models with grouped predictors." \emph{Statistics and Computing}, \bold{25}:173-187.
  
  Fan, Y. and Tang, C. Y. (2013). "Tuning parameter selection in high dimensional penalized likelihood." \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, \bold{75}:531-552.
  
  Wang, H. and Leng, C. (2007). "Unified LASSO estimation by least squares approximation." \emph{Journal of the American Statistical Association}, \bold{102}:1039-1048.
  
  Yuan, M. and Lin, Y. (2006). "Model selection and estimation in regression with grouped variables." \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, \bold{68}:49-67.
}

\examples{
## Generate data
set.seed(1234)
X = matrix(runif(100*11), nrow=100)
n = dim(X)[1]
groups = c(1,1,1,2,2,2,3,3,4,5,5)
beta_true = c(-1,1,1,0,0,0,0,0,0,1.5,-1.5)

## Generate responses from gamma regression with known shape parameter 1
eta = crossprod(t(X), beta_true)
shape = 1
Y = rgamma(n, rate=shape/exp(eta), shape=shape)

## Generate test data
n_test = 50
X_test = matrix(runif(n_test*11), nrow=n_test)

## Fit gamma regression models with the group SCAD penalty
gamma_mod = gamma_grpreg(Y, X, groups, X_test, penalty="gSCAD")

## Tuning parameters used to fit models 
gamma_mod$lambda

## Predicted n_test-dimensional vectors mu=E(Y_test) based on test data, X_test. 
## The kth column of 'Y_pred' corresponds to the kth entry in 'lambda.'
gamma_mod$Y_pred

## Classifications of the 5 groups. The kth column of 'classifications'
# corresponds to the kth entry in 'lambda.'
gamma_mod$classifications

## Plot lambda vs. GIC
plot(gamma_mod$lambda, gamma_mod$GIC, type='l')

## Model selection with the lambda that minimizes GIC
gamma_mod$lambda_min
gamma_mod$min_index 
gamma_mod$classifications[, gamma_mod$min_index]
gamma_mod$beta[, gamma_mod$min_index]
}