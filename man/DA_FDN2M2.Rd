% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DA_FDN2M2.R
\name{DA_FDN2M2}
\alias{DA_FDN2M2}
\title{Data Analysis of (N2,M2)}
\usage{
DA_FDN2M2(realdata_covariates, realdata_alpha, p, beta0)
}
\arguments{
\item{realdata_covariates}{a list contains the following data matrics:
CASEZhat_2, CASEZhat_22, CONTZhat_2, CONTZhat_22}

\item{realdata_alpha}{a list contains the following data matrics:
prob_case_22,prob_cont_2, pwt_cont_2}

\item{p}{number of parameters.}

\item{beta0}{an initial parameter for solver "nleqslv".}
}
\value{
A list of estimator and its standard deviation.
}
\description{
The main function to solve the estimating equations constructed by (N2,M2). Since there is just one case data, no selection bias needed.
}
\details{
The function solves estimating equation based on (N2,M2), see Huang(2014).
}
\examples{
 #p <- 8
 #beta0=c(-5.4163,0.7790,-0.1289,0.2773,-0.5510,0.1568,0.4353,-0.6895)
 #DA_FDN2M2(realdata_covariates,realdata_alpha,p=p,beta0=beta0)

}
\references{
Huang, H., Ma, X., Waagepetersen, R., Holford, T.R. , Wang, R., Risch, H., Mueller, L. & Guan, Y. (2014). A New Estimation Approach for Combining Epidemiological Data From Multiple Sources, Journal of the American Statistical Association, 109:505, 11-23.
}

