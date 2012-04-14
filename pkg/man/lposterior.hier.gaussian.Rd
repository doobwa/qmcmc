\name{lposterior.hier.gaussian}
\alias{lposterior.hier.gaussian}
\title{Log posterior for hierarchical gaussian model}
\usage{
  lposterior.hier.gaussian(y, theta, mu, sigma,
    priors = list(theta = list(mu = 0, sigma = 1), sigma = list(alpha = 2, beta = 0.1)),
    grad = FALSE, collapse.sigma = TRUE)
}
\arguments{
  \item{y}{}

  \item{theta}{}

  \item{mu}{}

  \item{sigma}{}

  \item{priors}{}

  \item{sigma}{}

  \item{grad}{}

  \item{collapse.sigma}{}
}
\value{
  log posterior density
}
\description{
  Log posterior for hierarchical gaussian model
}

