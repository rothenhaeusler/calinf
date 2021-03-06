% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distributional_uncertainty.R
\name{drchisq}
\alias{drchisq}
\title{Distributional random number generation (equivalent to rchisq)}
\usage{
drchisq(d_seed, df)
}
\arguments{
\item{d_seed}{A distributional seed as returned by distributional_seed()}

\item{df}{degrees of freedom}
}
\description{
Distributional random number generation (equivalent to rchisq)
}
\examples{
d_seed <- distributional_seed(1000,10)
drchisq(d_seed)
}
