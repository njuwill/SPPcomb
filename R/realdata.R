#' A data list of matrices containing covariates of cases and controls.
#'
#'  The list includes 10 matrices of covariates of cases and controls from different sources.
#'  Some of them need to impute the missing data, some of them need to estimate the variables
#'  even not missing to make sure the consistent format of input.
#'
#' @format A list of 10 matrices
#'
#'
#' @return
#' in the first case data which has complete cases:
#' \enumerate{
#' \item CASEZ_1=\eqn{[1, Z_d, Z_t, Z_l]}
#' \item  CASEZhat_1=\eqn{[1, Z_d, \hat{Z}_t, Z_l]}
#' }
#'
#' in the second case data(CTR) which has missing lifestyle covariates:
#' \enumerate{
#' \item CASEZ_2=\eqn{[1, Z_d, Z_t]}
#' \item CASEZhat_2=\eqn{[1, Z_d, Z_t, \hat{Z}_l]}
#' \item CASEZhat_22=\eqn{[1, Z_d, \hat{Z}_t, \hat{Z}_l]}
#' }
#'
#' in the 1st control data which has complete controls:
#' \enumerate{
#' \item CONTZ_1=\eqn{[1, Z_d, Z_t, Z_l]}
#' \item CONTZhat_1=\eqn{[1, Z_d, Z_t, \hat{Z}_l]}
#'}
#'
#' in the 2nd control data which has missing traffic covariates(BRFSS):
#' \enumerate{
#' \item CONTZ_2=\eqn{[1, Z_d, Z_l]}
#' \item CONTZhat_2=\eqn{[1, Z_d, \hat{Z}_t, Z_l]}
#' \item CONTZhat_22=\eqn{[1, Z_d, \hat{Z}_t, \hat{Z}_l]}
#'}
#'
#' @examples
#' # For example of each matrix, type the command in R: attributes(realdata_covariates)
#' # to obtain names of 10 bulit-in matrices in the list:
#' # "CASEZ_1", "CASEZhat_1", "CASEZ_2", "CASEZhat_2", "CASEZhat_22", "CONTZ_1",
#' # "CONTZhat_1", "CONTZ_2", "CONTZhat_2", "CONTZhat_22".
#'
 "realdata_covariates"

#' A list of matrices containing value of alpha at each location.
#'
#'
#'
#' @return A list of 8 matrices of calculated value of alpha for case and control points.
#' \describe{
#' \item{counts_agebysex_state}{Age-by-sex stratification in Connecticut, which is a
#' matrix with 18 rows and 2 variables: Male, Female. In this dataset the age-by-sex distribution based on the Census for
#' the following ten age groups: 35-40, 41-45, 46-50, 51-55, 56-60,
#' 61-65, 66-70, 71-75, 76-80, and 81-83.}
#'
#' \item{prob_case_1}{Value of alpha for cases in case-control study,
#'      this \eqn{\alpha_1(s)} is matched to controls' age-by-sex proportion}
#'
#'\item{prob_case_11}{Value of alpha for cases in case-control study,
#'     this \eqn{\alpha_2(s)} is matched to BRFSS age-by-sex proportion}
#'
#' \item{prob_case_2}{Value of alpha for cases in CTR, which is a matrix with 1929 rows and 1 variables: \eqn{\alpha_1(s)}.
#' This \eqn{\alpha_1(s)} in CTR is matched to controls' age-by-sex proportion in case-control study}
#'
#' \item{prob_case_22}{Value of alpha for cases in CTR, which is a matrix with 1929 rows and 1 variables: \eqn{\alpha_2(s)}.
#' This \eqn{\alpha_2(s)} in CTR is matched to BRFSS age-by-sex proportion}
#'
#' \item{prob_cont_1}{Value of alpha for controls in case-control study, which is a matrix with 690 rows and 1 variables:
#' \eqn{\alpha_1(s)}. A dataset of controls' \eqn{\alpha_1(s)} of its own}
#'
#' \item{prob_cont_2}{Another Value of alpha for controls in BRFSS data, which is a matrix with 4459 rows and 1 variables:
#' \eqn{\alpha_2(s)}. A dataset of controls' \eqn{\alpha_2(s)} in BRFSS data of its own}
#'
#' \item{pwt_cont_2}{Value of weights(sampling probability) for controls in BRFSS data, which is \eqn{\alpha_2^{\star}}
#' in equation(14) of Huang(2014), a matrix with 4459 rows and 1 variables}
#' }
#'
#' @examples
#' # For example of each matrix, type the command in R: attributes(realdata_alpha)
#' # to obtain names of 8 matrices in the list:
#' #"counts_agebysex_state", "prob_case_1", "prob_case_11", "prob_case_2",
#' #"prob_case_22", "prob_cont_1", "prob_cont_2", "pwt_cont_2".
#'
"realdata_alpha"

