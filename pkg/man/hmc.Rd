\name{hmc}
\alias{hmc}
\title{General purpose HMC (Radford Neal's code)}
\usage{
  hmc(current_q, lposterior, epsilon, L)
}
\arguments{
  \item{current}{parameter vector}

  \item{lposterior}{function returning the log posterior}

  \item{epsilon}{step size}

  \item{L}{number of steps}
}
\value{
  new parameter vector
}
\description{
  General purpose HMC (Radford Neal's code)
}

