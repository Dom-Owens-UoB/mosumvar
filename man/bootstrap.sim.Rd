% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bootstrap.R
\name{bootstrap.sim}
\alias{bootstrap.sim}
\title{Simulate bootstrapped threshold}
\usage{
bootstrap.sim(
  x,
  order = 1,
  G = NULL,
  uni = FALSE,
  method = c("Wald", "Score"),
  estim = c("C", "H"),
  var.estim = c("Local", "Global"),
  rm.cross.terms = TRUE,
  global.resids = TRUE,
  n.bootstrap = 1000
)
}
\description{
Simulate bootstrapped threshold
}
\keyword{internal}
