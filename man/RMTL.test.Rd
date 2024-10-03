\name{RMTL.test}
\alias{RMTL.test}
\title{
 Function to perform multiple RMTL-based tests
}
\description{
 Linear hypotheses of the restricted mean time losts (RMTLs) of k different groups in a competing risks setup can be tested simultaneously.
 Therefore, the multivariate distribution of local Wald-type test statistics is approximated by (1) estimating the covariance between the test statistics (\code{method = "asymptotic"}) or (3) a permutation approach with Bonferroni-correction (\code{method = "permutation"}), respectively. Hence, adjusted p-values can be obtained.
}
\usage{
RMTL.test(
  time = NULL,
  status = NULL,
  group = NULL,
  formula = NULL,
  event = NULL,
  data = NULL,
  hyp_mat,
  hyp_vec = NULL,
  M = NULL,
  tau,
  method = c("permutation", "asymptotic"),
  stepwise = FALSE,
  alpha = 0.05,
  Nres = 4999,
  seed = 1
)
}
\arguments{
  \item{time}{
  A vector containing the observed event times. Default option is \code{NULL}. Either \code{time, status, group} or \code{formula, event, data} needs to be specified.
}
\item{status}{
  A vector of the same length as \code{time} containing the corresponding event indicator with values 0 = censored and 1,...,M for the M different competing events. Default option is \code{NULL}. Either \code{time, status, group} or \code{formula, event, data} needs to be specified.
}
\item{group}{
  A vector of the same length as \code{time} containing the corresponding group labels. Default option is \code{NULL}. Either \code{time, status, group} or \code{formula, event, data} needs to be specified.
}
\item{formula}{
  A model \code{formula} object. The left hand side contains the time variable and the
right hand side contains the factor variables of interest. Default option is \code{NULL}. Either \code{time, status, group} or \code{formula, event, data} needs to be specified.
}
\item{event}{
  The name of censoring status indicator with values 0 = censored and 1,...,M for the M different competing events. Default option is \code{NULL}. Either \code{time, status, group} or \code{formula, event, data} needs to be specified.
}
\item{data}{
  A data.frame or list containing the variables in formula and the censoring status indicator. Default option is \code{NULL}. Either \code{time, status, group} or \code{formula, event, data} needs to be specified.
}
\item{hyp_mat}{
  A list containing all the hypothesis matrices H for the multiple tests or one of the options \code{"Tukey", "Dunnett"} and \code{"center"} for Tukey's or Dunnett's contrasts or the centering matrix, respectively, or \code{"2by2", "2by2 cause-wisely"} for tests on main and interaction effects in a 2-by-2 design without or with cause-wise results, respectively, or a matrix if only one hypothesis is of interest. For the permutation test, all matrices need to fulfill the contrast property as described in Munko et al. (2024).
}
\item{hyp_vec}{
  A list containing all the hypothesis vectors c for the multiple tests or a vector if only one hypothesis is of interest. By default (\code{NULL}), all hypothesis vectors are set to zero vectors of suitable length.
}
\item{M}{
  An integer specifying the number of competing risks. By default (\code{NULL}), the maximum of the values in \code{status} or \code{event} is chosen.
}
\item{tau}{
  A numeric value specifying the end of the relevant time window for the analysis. Default option is \code{NULL}.
}
\item{method}{
 One of the methods \code{"groupwise", "permutation"} and \code{"asymptotic"} that should be used for calculating the critical values. Default option is \code{"groupwise"}.
}
\item{stepwise}{
 A logical vector indicating whether a stepwise extension of the test should be performed. If \code{TRUE}, no confidence intervals can be computed for the linear combinations but it may be that more tests can reject. Default option is \code{FALSE}.
}
\item{alpha}{
  A numeric value specifying the global level of significance. Default option is \code{0.05}.
}
\item{Nres}{
  The number of random variables to approximate the limiting distribution. This is only used if at least one hypothesis matrix is not a row vector. The default option is \code{4999}.
}
\item{seed}{
  A single value, interpreted as an integer, for providing reproducibility of the results or \code{NULL} if reproducibility is not wanted. Default option is \code{1}.
}
}
\details{
The restricted mean time lost (RMTL) of group \eqn{i} and event \eqn{m} is defined as
\deqn{ \mu_{im} := \int\limits_0^{\tau} F_{im}(t) \mathrm{d}t  }
for all \eqn{i\in\{1,...,k\}, m\in\{1,...,M\}}, where \eqn{F_{im}} denotes the cumulative incidence function of group \eqn{i} and event \eqn{m}.
Let
\deqn{\boldsymbol{\mu} := (\mu_{11},\mu_{12},...,\mu_{kM})^{\prime}}
be the vector of the RMTLs and \deqn{\widehat{\boldsymbol{\mu}} := (\widehat{\mu}_{11},\widehat{\mu}_{12},...,\widehat{\mu}_{kM})^{\prime}} be the vector of their estimators.
Let \eqn{\mathbf{H}_{\ell} \in {\rm I\!R}^{r_{\ell}\times kM}} with \eqn{\mathrm{rank}(\mathbf{H}_{\ell}) >0} and \eqn{\mathbf{c}_{\ell} \in {\rm I\!R}^{r_{\ell}}} for all \eqn{\ell\in\{1,...,L\}}.
We are considering the multiple testing problem with null and alternative hypotheses
\deqn{
\mathcal{H}_{0,\ell}: \mathbf{H}_{\ell} \boldsymbol\mu = \mathbf{c}_{\ell} \quad \mathrm{vs.} \quad \mathcal{H}_{1,\ell}: \mathbf{H}_{\ell} \boldsymbol\mu \neq \mathbf{c}_{\ell}, \qquad \mathrm{for }\: \ell\in\{1,...,L\}.}
For the permutation test, the matrices additionally need to fulfill \eqn{\mathbf{H}_{\ell}(\mathbf{1}_k\otimes \mathbf{e}_m) = \mathbf{0}_{r_{\ell}}} for all \eqn{m, \ell}, where \eqn{\mathbf{1}_k\in {\rm I\!R}^{k}, \mathbf{e}_m \in {\rm I\!R}^{M}, \mathbf{0}_{r_{\ell}} \in {\rm I\!R}^{r_{\ell}}} denote the vector of ones, the \eqn{m}th unit vector and the vector of zeros, respectively.

For the argument \code{hyp_mat}, the options \code{"Dunnett"} and \code{"Tukey"} create the cause-wise many-to-one and all-pairs comparisons, respectively, as described in Example 3 of Munko et al. (2024) and the options \code{"2by2"} and \code{"2by2 cause-wisely"} create the hypothesis matrices as described in Example 4 of Munko et al. (2024).

If \code{stepwise = TRUE}, the closed testing procedure is applied. In this case, no confidence intervals can be computed for the linear combinations but potentially more tests can reject.

Note that the calculations for the permutation test may take a while.
}
\value{
A list of class \code{GFDrmst} containing the following components:
\item{method}{A character containing the method which has been used.}
\item{test_stat}{A numeric vector containing the calculated Wald-type test statistics for the local hypotheses.}
\item{p.value}{A numeric vector containing the adjusted p-values for the local hypotheses.}
\item{res}{A list containing the results of the multiple tests including the hypothesis matrices, estimators of the linear combinations of RMTLs, potentially confidence intervals for the linear combinations (if all matrices are row vectors and stepwise = \code{FALSE}), Wald-type test statistics, critical values and the test decisions.}
\item{alpha}{A numeric value containing the global level of significance.}
}
\references{
Munko, M., Dobler, D., Ditzhaus, M. (2024). Multiple tests for restricted mean time lost with competing risks data, arXiv preprint (arXiv:2409.07917). \doi{10.48550/arXiv.2409.07917}
}
\examples{
library(mstate)
data("ebmt2")

# multiple asymptotic tests
out <- RMTL.test(time = ebmt2$time,
                 status = ebmt2$status,
                 group = ebmt2$match,
                 hyp_mat = "Dunnett",
                 tau = 120,
                 method = "asymptotic")
summary(out)
plot(out)

## or, equivalently,
out <- RMTL.test(formula = "time ~ match",
                 event = "status",
                 data = ebmt2,
                 hyp_mat = "Dunnett",
                 tau = 120,
                 method = "asymptotic")
summary(out)
plot(out)

\donttest{
# multiple permutation tests
# this may take a few minutes
out_perm <- RMTL.test(formula = "time ~ match",
                      event = "status",
                      data = ebmt2,
                      hyp_mat = "Dunnett",
                      tau = 120,
                      method = "permutation",
                      Nres = 999)
summary(out_perm)
plot(out_perm)
}
}
