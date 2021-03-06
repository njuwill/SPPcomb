% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DA_FDN2M1M2.R
\name{DA_FDN2M1M2}
\alias{DA_FDN2M1M2}
\title{Data Analysis for Combining (N2,M1) + (N2,M2)}
\usage{
DA_FDN2M1M2(realdata_covariates, realdata_alpha, subset_2, subset_4, p, beta0)
}
\arguments{
\item{realdata_covariates}{a list contains the following data matrics:
CASEZ_2, CASEZhat_2, CASEZhat_22, CONTZ_1, CONTZhat_1,
CONTZhat_2, CONTZhat_22}

\item{realdata_alpha}{a list contains the following data matrics:
prob_case_22, prob_cont_1, prob_cont_2, pwt_cont_2}

\item{subset_2}{A vector of 1:(p-2), which is the subset of \eqn{\hat{Z}_{21}}, i.e.
\eqn{\hat{Z}_{21}^\star} in equation (10) of Huang(2014). \eqn{\hat{Z}_l} may be highly correlated with Z_d, so it is
   removed in the estimation.
 For the view of including more information, you can use the whole dataset.}

\item{subset_4}{A vector of 1:(p-2), which is the subset of \eqn{\hat{Z}_{22}}.}

\item{p}{number of parameters.}

\item{beta0}{an initial parameter for solver "nleqslv".}
}
\value{
A list of estimator and its standard deviation.
}
\description{
The main function to solve the estimating equations constructed by combining pair (N2,M1) and (N2,M2). Since there is just one case data, no selection bias needed.
}
\details{
The function solves estimating equation based on (N2,M1) and (N2,M2), see Huang(2014).
}
\examples{
 #p <- 8
 #subset_2 <- 1:p
 #subset_4 <- 1:p
 #beta0=c(-5.4163,0.7790,-0.1289,0.2773,-0.5510,0.1568,0.4353,-0.6895)
 #DA_FDN2M1M2(realdata_covariates,realdata_alpha,subset_2,subset_4,p=p,beta0=beta0)

}
\references{
Huang, H., Ma, X., Waagepetersen, R., Holford, T.R. , Wang, R., Risch, H., Mueller, L. & Guan, Y. (2014). A New Estimation Approach for Combining Epidemiological Data From Multiple Sources, Journal of the American Statistical Association, 109:505, 11-23.
}

