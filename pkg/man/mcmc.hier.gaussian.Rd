\name{mcmc.hier.gaussian}
\alias{mcmc.hier.gaussian}
\title{MCMC for hierarchical normal model}
\usage{
  mcmc.hier.gaussian(ys, method, theta = NULL, niter = 100,
    burnin = 10,
    priors = list(mu = list(mu = 0, sigma = 1), sigma = list(alpha = 2, beta = 0.1)),
    verbose = TRUE, ...)
}
\arguments{
  \item{ys}{list of vectors}

  \item{method}{mcmc method to use}

  \item{theta}{}

  \item{niter}{}

  \item{burnin}{}

  \item{priors}{}

  \item{verbose}{}

  \item{...}{pass extra parameters to the MCMC routine}

  \item{beta}{current}
}
\value{
  list with parameters at each iteration and a vector of
  log posteriors at each iteration
}
\description{
  MCMC for hierarchical normal model
}

