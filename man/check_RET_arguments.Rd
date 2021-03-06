% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{check_RET_arguments}
\alias{check_RET_arguments}
\title{check_RET_arguments}
\usage{
check_RET_arguments(sig.level, power, Delta, n, allocation)
}
\arguments{
\item{sig.level}{A numeric value specifying the significance level (type I error probability)}

\item{power}{A numeric value specifying the target power (1 - type II error probability)}

\item{Delta}{A numeric value specifying the non-inferiority or superiority margin.
Is between 0 and 1 in case of non-inferiority and larger than 1 in case of superiority.}

\item{n}{The total sample size. Needs to be at least 7.}

\item{allocation}{A (non-empty) vector specifying the sample size allocation (nExp/n, nRef/n, nPla/n)}
}
\description{
Check arguments for their respective condition
}
