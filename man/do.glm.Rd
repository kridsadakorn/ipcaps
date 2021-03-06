% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parallelization.R
\name{do.glm}
\alias{do.glm}
\title{(Internal function) Perform regression models, internally used for
parallelization}
\usage{
do.glm(X, PC, method = "linear")
}
\arguments{
\item{X}{A vector of data}

\item{PC}{A matrix of principal components}

\item{method}{Specify a method to be used for regression model, which can be
"linear", "poisson", and "negative.binomial". Default = "linear"}
}
\value{
A vector of residuals computed from regression model
}
\description{
(Internal function) Perform regression models, internally used for
parallelization
}
