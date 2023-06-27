\name{cv_SSGL}
\alias{cv_SSGL}

\title{Cross-Validation for Spike-and-Slab Group Lasso in Group-Regularized Generalized Linear Models (GLMs)}

\description{This function implements \eqn{K}-fold cross-validation for group-regularized GLMs with the spike-and-slab group lasso (SSGL) penalty of Bai et al. (2022) and Bai (2023). The identity link function is used for Gaussian regression, the logit link is used for binomial regression, and the log link is used for Poisson, negative binomial, and gamma regression. 

Although one can choose \code{lambda0} from cross-validation with this function, it can be very time-consuming to do so if the number of groups \eqn{G} and/or the number of total covariantes \eqn{p} is moderate to large. It is \emph{strongly} recommended that the user simply run  the \code{SSGL} function on the training dataset and select the final model according to the \code{lambda0} that minimizes the generalized information criterion (GIC). See description of the \code{SSGL} function for more details.
}

\usage{
cv_SSGL(Y, X, groups, 
        family=c("gaussian","binomial","poisson","negativebinomial","gamma"), 
        nb_size=1, gamma_shape=1, group_weights, n_folds=5, n_lambda0=25,
        lambda0, lambda1=1, a=1, b=dim(X)[2], 
        max_iter=100, tol=1e-6, print_fold=TRUE) 
}

\arguments{
  \item{Y}{\eqn{n \times 1} vector of responses for training data.}
  \item{X}{\eqn{n \times p} design matrix for training data, where the \eqn{j}th column corresponds to the \eqn{j}th overall feature.}
  \item{groups}{\eqn{p}-dimensional vector of group labels. The \eqn{j}th entry in \code{groups} should contain either the group number \emph{or} the factor level name that the feature in the \eqn{j}th column of \code{X} belongs to. \code{groups} must be either a vector of integers or factors.}
  \item{family}{exponential dispersion family of the response variables. Allows for \code{"gaussian"}, \code{"binomial"}, \code{"poisson"}, \code{"negativebinomial"}, and \code{"gamma"}. Note that for \code{"negativebinomial"}, the size parameter must be specified in advance, while for \code{"gamma"}, the shape parameter must be specified in advance.}
    \item{nb_size}{known size parameter \eqn{\alpha} in \eqn{NB(\alpha,\mu_i)} distribution for the responses if the user specifies \code{family="negativebinomial"}. Default is \code{nb_size=1}. Ignored if \code{family} is not \code{"negativebinomial"}.}
  \item{gamma_shape}{known shape parameter \eqn{\nu} in \eqn{G(\mu_i,\nu)} distribution for the responses if the user specifies \code{family="gamma"}. Default is \code{gamma_shape=1}. Ignored if \code{family} is not \code{"gamma"}.}
  \item{group_weights}{group-specific, nonnegative weights for the penalty. Default is to use the square roots of the group sizes.}
    \item{n_folds}{number of folds \eqn{K} to use in \eqn{K}-fold cross-validation. Default is \code{n_folds=5}.}
  \item{n_lambda0}{number of spike hyperparameters \eqn{L}. Default is \code{n_lambda0=25}.}
  \item{lambda0}{grid of \eqn{L} spike hyperparameters \eqn{\lambda_0}. The user may specify either a scalar or a vector. If the user does not provide this, the program chooses the grid automatically.}
  \item{lambda1}{slab hyperparameter \eqn{\lambda_1} in the SSGL prior. Default is \code{lambda1=1}.}
  \item{a}{shape hyperparameter for the \eqn{Beta(a,b)} prior on the mixing proportion in the SSGL prior. Default is \code{a=1}.}
  \item{b}{shape hyperparameter for the \eqn{Beta(a,b)} prior on the mixing proportion in the SSGL prior. Default is \code{b=dim(X)[2]}.}
  \item{max_iter}{maximum number of iterations in the algorithm. Default is \code{max_iter=100}.}
  \item{tol}{convergence threshold for algorithm. Default is \code{tol=1e-6}.}
  \item{print_fold}{Boolean variable for whether or not to print the current fold in the algorithm. Default is \code{print_fold=TRUE}.}
}

\value{The function returns a list containing the following components:
    \item{lambda0}{\eqn{L \times 1} vector of spike hyperparameters \code{lambda0} used to fit the model. \code{lambda0} is displayed in descending order.}
  \item{cve}{\eqn{L \times 1} vector of mean cross-validation error across all \eqn{K} folds. The \eqn{k}th entry in \code{cve} corresponds to the \eqn{k}th spike hyperparameter parameter in \code{lambda0}.}
  \item{cvse}{\eqn{L \times 1} vector of standard errors for cross-validation error across all \eqn{K} folds. The \eqn{k}th entry in \code{cvse} corresponds to the \eqn{k}th spike hyperparameter parameter in \code{lambda0}.}
  \item{lambda0_min}{The value in \code{lambda0} that minimizes mean cross-validation error \code{cve}.}
  \item{min_index}{The index of \code{lambda0_min} in \code{lambda0}.}
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
ssgl_mods = cv_SSGL(Y, X, groups, family="gaussian", lambda0=seq(from=16,to=4,by=-4))

## Plot cross-validation curve
plot(ssgl_mods$lambda0, ssgl_mods$cve, type="l", xlab="lambda0", ylab="CVE")
## lambda which minimizes mean CVE
ssgl_mods$lambda0_min
ssgl_mods$min_index

\donttest{
## Example with Poisson regression

## Generate count responses
eta = crossprod(t(X), beta_true)
Y = rpois(n,exp(eta))

## K-fold cross-validation 
## NOTE: If you do not specify lambda0, the program will automatically choose a suitable grid.
ssgl_poisson_mods = cv_SSGL(Y, X, groups, family="poisson", lambda0=seq(from=20,to=2,by=-4))

## Plot cross-validation curve
plot(ssgl_poisson_mods$lambda0, ssgl_poisson_mods$cve, type="l", xlab="lambda0", ylab="CVE")
## lambda which minimizes mean CVE
ssgl_poisson_mods$lambda0_min
ssgl_poisson_mods$min_index
}
}