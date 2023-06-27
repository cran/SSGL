\name{cv_gamma_grpreg}
\alias{cv_gamma_grpreg}

\title{Cross-validation for Group-Regularized Gamma Regression}

\description{This function implements \eqn{K}-fold cross-validation for group-regularized gamma regression with a known shape parameter \eqn{\nu} and the log link. The cross-validation error (CVE) and cross-validation standard error (CVSE) are computed using the deviance for gamma regression.

For a description of group-regularized gamma regression, see the description for the \code{gamma_grpreg} function. Our implementation is based on the least squares approximation approach of Wang and Leng (2007), and hence, the function does not allow the total number of covariates \eqn{p} to be greater than \eqn{\frac{K-1}{K} \times} sample size, where \eqn{K} is the number of folds.

Note that the \code{gamma_grpreg} function also returns the generalized information criterion (GIC) of Fan and Tang (2013) for each regularization parameter in \code{lambda}, and the GIC can also be used for model selection instead of cross-validation.
}

\usage{
cv_gamma_grpreg(Y, X, groups, gamma_shape=1, penalty=c("gLASSO","gSCAD","gMCP"),
                n_folds=10, group_weights, taper, n_lambda=100, lambda, 
                max_iter=10000, tol=1e-4)
}

\arguments{
  \item{Y}{\eqn{n \times 1} vector of strictly positive, continuous responses for training data.}
  \item{X}{\eqn{n \times p} design matrix for training data, where the \eqn{j}th column corresponds to the \eqn{j}th overall feature.}
  \item{groups}{\eqn{p}-dimensional vector of group labels. The \eqn{j}th entry in \code{groups} should contain either the group number \emph{or} the factor level name that the feature in the \eqn{j}th column of \code{X} belongs to. \code{groups} must be either a vector of integers or factors.}
  \item{gamma_shape}{known shape parameter \eqn{\nu} in \eqn{Gamma(\mu_i,\nu)} distribution for the responses. Default is \code{gamma_shape=1}.}
  \item{penalty}{group regularization method to use on the groups of regression coefficients. The options are \code{"gLASSO"}, \code{"gSCAD"}, \code{"gMCP"}. To implement cross-validation for gamma regression with the SSGL penalty, use the \code{cv_SSGL} function.}
  \item{n_folds}{number of folds \eqn{K} to use in \eqn{K}-fold cross-validation. Default is \code{n_folds=10}.}
  \item{group_weights}{group-specific, nonnegative weights for the penalty. Default is to use the square roots of the group sizes.}
  \item{taper}{tapering term \eqn{\gamma} in group SCAD and group MCP controlling how rapidly the penalty tapers off. Default is \code{taper=4} for group SCAD and \code{taper=3} for group MCP. Ignored if \code{"gLASSO"} is specified as the penalty.}
  \item{n_lambda}{number of regularization parameters \eqn{L}. Default is \code{n_lambda=100}.}
  \item{lambda}{grid of \eqn{L} regularization parameters. The user may specify either a scalar or a vector. If the user does not provide this, the program chooses the grid automatically.}
  \item{max_iter}{maximum number of iterations in the algorithm. Default is \code{max_iter=10000}.}
  \item{tol}{convergence threshold for algorithm. Default is \code{tol=1e-4}.}
}

\value{The function returns a list containing the following components:
    \item{lambda}{\eqn{L \times 1} vector of regularization parameters \code{lambda} used to fit the model. \code{lambda} is displayed in descending order.}
  \item{cve}{\eqn{L \times 1} vector of mean cross-validation error across all \eqn{K} folds. The \eqn{k}th entry in \code{cve} corresponds to the \eqn{k}th regularization parameter in \code{lambda}.}
  \item{cvse}{\eqn{L \times 1} vector of standard errors for cross-validation error across all \eqn{K} folds. The \eqn{k}th entry in \code{cvse} corresponds to the \eqn{k}th regularization parameter in \code{lambda}.}
  \item{lambda_min}{The value in \code{lambda} that minimizes mean cross-validation error \code{cve}.}
  \item{min_index}{The index of \code{lambda_min} in \code{lambda}.}
}

\references{
  Breheny, P. and Huang, J. (2015). "Group descent algorithms for nonconvex penalized linear and logistic regression models with grouped predictors." \emph{Statistics and Computing}, \bold{25}:173-187.
  
  Fan, Y. and Tang, C. Y. (2013). "Tuning parameter selection in high-dimensional penalized likelihood." \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, \bold{75}:531-552.
  
  Wang, H. and Leng, C. (2007). "Unified LASSO estimation by least squares approximation." \emph{Journal of the American Statistical Association}, \bold{102}:1039-1048.
  
  Yuan, M. and Lin, Y. (2006). "Model selection and estimation in regression with grouped variables." \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, \bold{68}:49-67.
}

\examples{
## Generate data
set.seed(12345)
X = matrix(runif(100*11), nrow=100)
n = dim(X)[1]
groups = c(1,1,1,2,2,2,3,3,4,5,5)
beta_true = c(-1,1,1,0,0,0,0,0,0,1.5,-1.5)

## Generate responses from gamma regression with known shape parameter 1
eta = crossprod(t(X), beta_true)
shape = 1
Y = rgamma(n, rate=shape/exp(eta), shape=shape)

## 10-fold cross-validation for group-regularized gamma regression
## with the group LASSO penalty
gamma_cv = cv_gamma_grpreg(Y, X, groups, penalty="gLASSO")

## Plot cross-validation curve
plot(gamma_cv$lambda, gamma_cv$cve, type="l", xlab="lambda", ylab="CVE")
## lambda which minimizes mean CVE
gamma_cv$lambda_min
## index of lambda_min in lambda
gamma_cv$min_index
}