\name{SSGL_cv}
\alias{SSGL_cv}

\title{Cross-Validation for Spike-and-Slab Group Lasso in Group-Regularized Generalized Linear Models (GLMs)}

\description{The \code{SSGL_cv} function implements \eqn{K}-fold cross-validation for choosing the regularization parameter \eqn{\lambda_0} in group-regularized GLMs with the spike-and-slab group lasso (SSGL) penalty of Bai et al. (2022) and Bai (2023). The default is \eqn{K=10}. The identity link function is used for Gaussian regression, the logit link is used for binomial regression, and the log link is used for Poisson regression. 

Although you can choose \code{lambda0} from cross-validation with this function, it can be time-consuming to do so if the number of groups \eqn{G} and/or the number of total covariantes \eqn{p} is moderate to large. In this case, you may choose to set the argument \code{parallelize=TRUE}, which will perform \eqn{K}-fold cross-validation in parallel across the \eqn{K} folds.
If \eqn{K} cores are used, then this may offer a speed-up of roughly the order of \eqn{K}.

As an alternative to cross-validation, you can also simply use the \code{SSGL} function on your data and select the final model according to the \code{lambda0} which minimizes the generalized information criterion (GIC). See description of the \code{SSGL} function for more details.
}

\usage{
SSGL_cv(Y, X, groups, 
        family=c("gaussian","binomial","poisson"), 
        group_weights, n_folds=10, n_lambda0=25,
        lambda0, lambda1=1, a=1, b=length(unique(groups)), 
        max_iter=100, tol=1e-6, parallelize=FALSE, n_cores) 
}

\arguments{
  \item{Y}{\eqn{n \times 1} vector of responses for training data.}
  \item{X}{\eqn{n \times p} design matrix for training data, where the \eqn{j}th column corresponds to the \eqn{j}th overall feature.}
  \item{groups}{\eqn{p}-dimensional vector of group labels. The \eqn{j}th entry in \code{groups} should contain either the group number \emph{or} the factor level name that the feature in the \eqn{j}th column of \code{X} belongs to. \code{groups} must be either a vector of integers or factors.}
  \item{family}{exponential dispersion family of the response variables. Allows for \code{"gaussian"}, \code{"binomial"}, and \code{"poisson"}.}
  \item{group_weights}{group-specific, nonnegative weights for the penalty. Default is to use the square roots of the group sizes.}
  \item{n_folds}{number of folds \eqn{K} to use in \eqn{K}-fold cross-validation. Default is \code{n_folds=10}.}
  \item{n_lambda0}{number of spike hyperparameters \eqn{L}. Default is \code{n_lambda0=25}.}
  \item{lambda0}{grid of \eqn{L} spike hyperparameters \eqn{\lambda_0}. The user may specify either a scalar or a vector. If the user does not provide this, the program chooses the grid automatically.}
  \item{lambda1}{slab hyperparameter \eqn{\lambda_1} in the SSGL prior. Default is \code{lambda1=1}.}
  \item{a}{shape hyperparameter for the \eqn{Beta(a,b)} prior on the mixing proportion in the SSGL prior. Default is \code{a=1}.}
  \item{b}{shape hyperparameter for the \eqn{Beta(a,b)} prior on the mixing proportion in the SSGL prior. Default is \code{b=length(unique(groups))}, i.e. the number of groups.}
  \item{max_iter}{maximum number of iterations in the algorithm. Default is \code{max_iter=100}.}
  \item{tol}{convergence threshold for algorithm. Default is \code{tol=1e-6}.}
  \item{parallelize}{Boolean variable for whether or not to parallelize \eqn{K}-fold cross-validation across the \eqn{K} folds. If the number of group \eqn{G} and/or the number of predictors \eqn{p} is moderate or large, then it may be preferable to perform cross-validation in parallel. In this case, the user can set \code{parallelize=TRUE}.}
  \item{n_cores}{Number of cores to use for parallelization. If the user does not specify this, the function will use the minimum of either \eqn{K} or the number of available cores minus one.}
}

\value{The function returns a list containing the following components:
    \item{lambda0}{\eqn{L \times 1} vector of spike hyperparameters \code{lambda0} used to fit the model. \code{lambda0} is displayed in descending order.}
  \item{cve}{\eqn{L \times 1} vector of mean cross-validation error across all \eqn{K} folds. The \eqn{k}th entry in \code{cve} corresponds to the \eqn{k}th spike hyperparameter parameter in \code{lambda0}.}
  \item{cvse}{\eqn{L \times 1} vector of standard errors for cross-validation error across all \eqn{K} folds. The \eqn{k}th entry in \code{cvse} corresponds to the \eqn{k}th spike hyperparameter parameter in \code{lambda0}.}
  \item{lambda0_cve_min}{The value in \code{lambda0} that minimizes mean cross-validation error \code{cve}.}
  \item{min_cve_index}{The index of \code{lambda0_cve_min} in \code{lambda0}.}
}

\references{
Bai, R. (2023). "Bayesian group regularization in generalized linear models with a continuous spike-and-slab prior." \emph{arXiv pre-print arXiv:2007.07021}.

Bai, R., Moran, G. E., Antonelli, J. L., Chen, Y., and Boland, M.R. (2022). "Spike-and-slab group lassos for grouped regression and sparse generalized additive models." \emph{Journal of the American Statistical Association}, \bold{117}:184-197.
}

\examples{
## Generate data
set.seed(12345)
X = matrix(runif(50*6), nrow=50)
n = dim(X)[1]
groups = c(1,1,1,2,2,2)
beta_true = c(-2,1,1.5,0,0,0)

## Generate responses from Gaussian distribution
Y = crossprod(t(X), beta_true) + rnorm(n)

## K-fold cross-validation 
## NOTE: If you do not specify lambda0, the function will automatically choose a suitable grid.

ssgl_mods = SSGL_cv(Y, X, groups, family="gaussian", n_folds=5, lambda0=seq(from=16,to=4,by=-4))

## Plot cross-validation curve
plot(ssgl_mods$lambda0, ssgl_mods$cve, type="l", xlab="lambda0", ylab="CVE")

## lambda which minimizes mean CVE
ssgl_mods$lambda0_cve_min
ssgl_mods$min_cve_index

\donttest{
## Example with binary logistic regression

## Generate binary responses
set.seed(123)
X = matrix(runif(50*6), nrow=50)
n = dim(X)[1]
groups = c(1,1,2,2,3,3)
beta_true = c(-2,1.5,0,0,2,-1.5)
eta = crossprod(t(X), beta_true)
Y = rbinom(n, size=1, prob=1/(1+exp(-eta)))

## K-fold cross-validation. Set parallelize=TRUE for potential speed-up
## If n_cores is not specified, then the function will automatically choose 
# the minimum of either K or the number of available cores minus one.

ssgl_logistic_mods = SSGL_cv(Y, X, groups, family="binomial", parallelize=TRUE, n_cores=2)

## Plot cross-validation curve
plot(ssgl_logistic_mods$lambda0, ssgl_logistic_mods$cve, type="l", xlab="lambda0", ylab="CVE")

## lambda which minimizes mean CVE
ssgl_logistic_mods$lambda0_cve_min
ssgl_logistic_mods$min_cve_index
}
}