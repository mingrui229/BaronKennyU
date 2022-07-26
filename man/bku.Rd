% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sensitivity.R
\name{bku}
\alias{bku}
\title{Sensitivity analysis for the Baron-Kenny approach with a single mediator}
\usage{
bku(
  y,
  m,
  a,
  c = NULL,
  r2.max = 0.5,
  randomized.a = FALSE,
  s1 = c(1, -1),
  s2 = c(1, -1),
  s3 = c(1, -1),
  s4 = c(1, -1),
  B = 500,
  shortcut = FALSE,
  progress.display = TRUE
)
}
\arguments{
\item{y}{a vector of the outcome variable}

\item{m}{a vector of the mediator variable}

\item{a}{a vector of the treatment variable}

\item{c}{a vector or matrix of the covariates}

\item{r2.max}{the maximum r^2 to be considered for the grid search}

\item{randomized.a}{if the treatment is randomized}

\item{s1}{(optional) the pre-specified sign S1; see Zhang and Ding (2022)}

\item{s2}{(optional) the pre-specified sign S2; see Zhang and Ding (2022)}

\item{s3}{(optional) the pre-specified sign S3; see Zhang and Ding (2022)}

\item{s4}{(optional) the pre-specified sign S4; see Zhang and Ding (2022)}

\item{B}{the number of bootstrap resampling}

\item{shortcut}{if use the shortcut approach}

\item{progress.display}{if display the progress}
}
\value{
a summary table, and the point estimates, standard errors, t statistics
of the direct and indirect effects with sensitivity parameters generated by
the grid search.
}
\description{
Sensitivity analysis for the Baron-Kenny approach with a single mediator
}
\examples{
set.seed(1234)
a = rnorm(1000)
m = 0.2 * a + rnorm(1000)
u = (-sqrt(1/3.2) - 0.1) * a + 0.5 * m + rnorm(1000)
y = 0.5 * u + rnorm(1000)
cor(a, u)^2
cor(lm(m~a)$residual, lm(u~a)$residual)^2
cor(lm(y~m+a)$residual, lm(u~m+a)$residual)^2
result = bku(y, m, a)
result$coefficients
}
