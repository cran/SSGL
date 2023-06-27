\name{cv_nb_grpreg}
\alias{cv_nb_grpreg}

\title{Cross-validation for Group-Regularized Negative Binomial Regression}

\description{This function implements \eqn{K}-fold cross-validation for group-regularized negative binomial regression with a known size parameter \eqn{\alpha} and the log link. The cross-validation error (CVE) and cross-validation standard error (CVSE) are computed using the deviance for negative binomial regression.

For a description of group-regularized negative binomial regression, see the description for the \code{nb_grpreg} function. Our implementation is based on the least squares approximation approach of Wang and Leng (2007), and hence, the function does not allow the total number of covariates \eqn{p} to be greater than \eqn{\frac{K-1}{K} \times} sample size, where \eqn{K} is the number of folds.

Note that the \code{nb_grpreg} function also returns the generalized information criterion (GIC) of Fan and Tang (2013) for each regularization parameter in \code{lambda}, and the GIC can also be used for model selection instead of cross-validation.
}

\usage{
cv_nb_grpreg(Y, X, groups, nb_size=1, penalty=c("gLASSO","gSCAD","gMCP"),
            n_folds=10, group_weights, taper, n_lambda=100, lambda, 
            max_iter=10000, tol=1e-4) 
}

\arguments{
  \item{Y}{\eqn{n \times 1} vector of strictly nonnegative integer responses for training data.}
  \item{X}{\eqn{n \times p} design matrix for training data, where the \eqn{j}th column corresponds to the \eqn{j}th overall feature.}
  \item{groups}{\eqn{p}-dimensional vector of group labels. The \eqn{j}th entry in \code{groups} should contain either the group number \emph{or} the factor level name that the feature in the \eqn{j}th column of \code{X} belongs to. \code{groups} must be either a vector of integers or factors.}
  \item{nb_size}{known size parameter \eqn{\alpha} in \eqn{NB(\alpha,\mu_i)} distribution for the responses. Default is \code{nb_size=1}.}
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
  
  Fan, Y. and Tang, C. Y. (2013). "Tuning parameter selection in high dimensional penalized likelihood." \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, \bold{75}:531-552. 
  
  Wang, H. and Leng, C. (2007). "Unified LASSO estimation by least squares approximation." \emph{Journal of the American Statistical Association}, \bold{102}:1039-1048.
  
  Yuan, M. and Lin, Y. (2006). "Model selection and estimation in regression with grouped variables." \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, \bold{68}:49-67.
}

\examples{
## Generate data
set.seed(1234)
X = matrix(runif(100*14), nrow=100)
n = dim(X)[1]
groups = c(1,1,1,2,2,2,2,3,3,4,5,5,6,6)
beta_true = c(-1,1,1,0,0,0,0,-1,1,0,0,0,-1.5,1.5)

## Generate count responses from negative binomial regression
eta = crossprod(t(X), beta_true)
Y = rnbinom(n, size=1, mu=exp(eta))

## 10-fold cross-validation for group-regularized negative binomial
## regression with the group MCP penalty
nb_cv = cv_nb_grpreg(Y, X, groups, penalty="gMCP")

## Plot cross-validation curve
plot(nb_cv$lambda, nb_cv$cve, type="l", xlab="lambda", ylab="CVE")
## lambda which minimizes mean CVE
nb_cv$lambda_min 
## index of lambda_min in lambda
nb_cv$min_index
}