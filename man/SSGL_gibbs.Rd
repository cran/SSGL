\name{SSGL_gibbs}
\alias{SSGL_gibbs}

\title{Gibbs sampling for Spike-and-Slab Group Lasso in Group-Regularized Generalized Linear Models (GLMs)}

\description{The \code{SSGL_gibbs} function implements Gibbs sampling for group-regularized GLMs with the spike-and-slab group lasso (SSGL) prior of Bai et al. (2022) and Bai (2023). The identity link function is used for Gaussian regression, the logit link is used for binomial regression, and the log link is used for Poisson regression. 

For binomial and Poisson regression, Polya-gamma data augmentation (Polson et al., 2013) is used to draw MCMC samples. The details are described in Bai (2023).

Note that the \code{SSGL_gibbs} function only returns the posterior mean, the 95 percent posterior credible intervals, and the posterior samples for the elements of the model parameter \eqn{\beta} and the predicted mean response \eqn{\mu_{test} = E(Y_{test})}. This function does not perform variable selection.

It is recommended that you use the \code{SSGL} function to perform variable selection and MAP estimation. If uncertainty quantification is \emph{also} desired, then this \code{SSGL_gibbs} function can be used.
}

\usage{
SSGL_gibbs(Y, X, groups, family=c("gaussian","binomial","poisson"), 
           X_test, group_weights, lambda0=5, lambda1=1, 
           a=1, b=length(unique(groups)),
           burn=1000, n_mcmc=2000, save_samples=TRUE) 
}

\arguments{
  \item{Y}{\eqn{n \times 1} vector of responses for training data.}
  \item{X}{\eqn{n \times p} design matrix for training data, where the \eqn{j}th column corresponds to the \eqn{j}th overall feature.}
  \item{groups}{\eqn{p}-dimensional vector of group labels. The \eqn{j}th entry in \code{groups} should contain either the group number \emph{or} the factor level name that the feature in the \eqn{j}th column of \code{X} belongs to. \code{groups} must be either a vector of integers or factors.}
  \item{family}{exponential dispersion family of the response variables. Allows for \code{"gaussian"}, \code{"binomial"}, and \code{"poisson"}.}
  \item{X_test}{\eqn{n_{test} \times p} design matrix for test data to calculate predictions. \code{X_test} must have the \emph{same} number of columns as X, but not necessarily the same number of rows. If \emph{no} test data is provided or if in-sample predictions are desired, then the function automatically sets \code{X_test=X} in order to calculate \emph{in-sample} predictions.}
  \item{group_weights}{group-specific, nonnegative weights for the penalty. Default is to use the square roots of the group sizes.}
  \item{lambda0}{spike hyperparameter \eqn{\lambda_0} in the SSGL prior. Default is \code{lambda0=5}.}
  \item{lambda1}{slab hyperparameter \eqn{\lambda_1} in the SSGL prior. Default is \code{lambda1=1}.}
  \item{a}{shape hyperparameter for the \eqn{Beta(a,b)} prior on the mixing proportion in the SSGL prior. Default is \code{a=1}.}
  \item{b}{shape hyperparameter for the \eqn{Beta(a,b)} prior on the mixing proportion in the SSGL prior. Default is \code{b=length(unique(groups))}, i.e. the number of groups.}
  \item{burn}{Number of warm-up MCMC samples to discard as burn-in. Default is \code{burn=1000}.}
  \item{n_mcmc}{Number of MCMC samples to save for posterior inference. Default is \code{n_mcmc=2000}.}
  \item{save_samples}{Boolean variable for whether or not to save the MCMC samples for \eqn{\beta} and predicted mean response \eqn{\mu_{test} = E(Y_{text})}. Default is \code{save_samples=TRUE}.}
}

\value{The function returns a list containing the following components:
  \item{beta_hat}{estimated posterior mean of \eqn{p \times 1} regression coefficient vector \eqn{\beta}.}
  \item{Y_pred_hat}{estimated posterior mean of \eqn{n_{test} \times 1} vector of predicted mean response values \eqn{\mu_{test} = E(Y_{test})} based on the test data in \code{X_test} (or training data \code{X} if no argument was specified for \code{X_test}).}
  \item{beta_lower}{\eqn{p \times 1} vector of lower endpoints of the 95 percent posterior credible intervals for \eqn{\beta}.}
  \item{beta_upper}{\eqn{p \times 1} vector of upper endpoints of the 95 percent posterior credible intervals for \eqn{\beta}.}
  \item{Y_pred_lower}{\eqn{n_{test} \times 1} vector of lower endpoints of the 95 percent posterior credible intervals for \eqn{\mu_{test} = E(Y_{test})}.}
  \item{Y_pred_upper}{\eqn{n_{test} \times 1} vector of upper endpoints of the 95 percent posterior credible intervals for \eqn{\mu_{test} = E(Y_{test})}.}
  \item{beta_samples}{\eqn{p \times} \code{n_mcmc} matrix of saved posterior samples for \eqn{\beta}. The \eqn{j}th row of \code{beta_samples} consists of the posterior samples for the \eqn{j}th regression coefficient in \eqn{\beta}. This is not returned if \code{save_samples=FALSE}.}
  \item{Y_pred_samples}{\eqn{n_{test} \times} \code{n_mcmc} matrix of saved posterior samples for \eqn{\beta}. The \eqn{i}th row of \code{Y_pred_samples} consists of the posterior samples of the predicted mean response \eqn{\mu_{i,test} = E(Y_{i,test})} for the \eqn{i}th test point. This is not returned if \code{save_samples=FALSE}.}
}

\references{
Bai, R. (2023). "Bayesian group regularization in generalized linear models with a continuous spike-and-slab prior." \emph{arXiv pre-print arXiv:2007.07021}.

Polson, N. G., Scott, J. G., and Windle, J. (2013). "Bayesian inference for logistic models using Polya-gamma latent variables." \emph{Journal of the American Statistical Association}, \bold{108}: 1339-1349.
}

\examples{
## Generate data
set.seed(1)
X = matrix(runif(200*17), nrow=200)
X_test = matrix(runif(20*17), nrow=20)

n = dim(X)[1]
n_test = dim(X_test)[1]

groups = c(1,1,1,2,2,2,2,3,3,3,4,4,5,5,6,6,6)
true_beta = c(-2,2,2,0,0,0,0,0,0,0,0,0,2.5,-2.5,0,0,0)
Y = crossprod(t(X), true_beta) + rnorm(n)

## Fit SSGL model. You should use the default burn=1000 and n_mcmc=2000

SSGL_mod = SSGL_gibbs(Y, X, groups, family="gaussian", X_test, burn=500, n_mcmc=1000)

## Evaluate results
cbind("True Beta" = true_beta, 
      "Posterior Mean" = SSGL_mod$beta_hat, 
      "95 CI lower" = SSGL_mod$beta_lower, 
      "95 CI upper"= SSGL_mod$beta_upper)

## Predictions on test data
cbind("Predicted E(Y)" = SSGL_mod$Y_pred_hat, 
      "95 CI lower" = SSGL_mod$Y_pred_lower, 
      "95 CI upper" = SSGL_mod$Y_pred_upper)

\donttest{
## Example with binary logistic regression

## Generate data
set.seed(123)
X = matrix(runif(200*16), nrow=200)
X_test = matrix(runif(50*16), nrow=50)
n = dim(X)[1]
n_test = dim(X)[2]
groups = c(1,1,1,1,2,2,2,2,3,4,4,5,5,6,6,6)
true_beta = c(-2,2,2,-2,0,0,0,0,0,0,0,2.5,-2.5,0,0,0)

## Generate binary responses
eta = crossprod(t(X), true_beta)
Y = rbinom(n, 1, 1/(1+exp(-eta)))

## Fit  SSGL logistic model
SSGL_logistic_mod = SSGL_gibbs(Y, X, groups, family="binomial", X_test)

## Evaluate results
cbind("True Beta" = true_beta, 
      "Posterior Mean" = SSGL_logistic_mod$beta_hat, 
      "95 CI lower" = SSGL_logistic_mod$beta_lower, 
      "95 CI upper"= SSGL_logistic_mod$beta_upper)

## Predictions on test data
cbind("Predicted E(Y)" = SSGL_logistic_mod$Y_pred_hat, 
      "95 CI lower" = SSGL_logistic_mod$Y_pred_lower, 
      "95 CI upper" = SSGL_logistic_mod$Y_pred_upper)
}
}