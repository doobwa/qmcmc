\name{mh}
\alias{mh}
\title{General purpose Metropolis-Hastings sampler}
\usage{
  mh(current, lposterior, olp = NULL, sd = 0.1)
}
\arguments{
  \item{current}{parameter vector}

  \item{lposterior}{function returning an (unnormalized)
  density}

  \item{olp}{log density of current}

  \item{sd}{standard deviation for Gaussian proposal
  density}
}
\value{
  new state and its log density as an attribute
}
\description{
  General purpose Metropolis-Hastings sampler
}

