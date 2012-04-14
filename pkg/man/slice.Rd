\name{slice}
\alias{slice}
\title{General purpose slice sampler (Radford Neal's code)}
\usage{
  slice(current, lposterior, olp = NULL, m = 20)
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
  General purpose slice sampler (Radford Neal's code)
}

