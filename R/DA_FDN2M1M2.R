#' Data Analysis for Combining (N2,M1) + (N2,M2)
#'
#' The main function to solve the estimating equations constructed by combining pair (N2,M1) and (N2,M2). Since there is just one case data, no selection bias needed.
#'
#' @param realdata_covariates a list contains the following data matrics:
#'        CASEZ_2, CASEZhat_2, CASEZhat_22, CONTZ_1, CONTZhat_1,
#'        CONTZhat_2, CONTZhat_22
#' @param realdata_alpha a list contains the following data matrics:
#'        prob_case_22, prob_cont_1, prob_cont_2, pwt_cont_2
#' @param beta0 an initial parameter for solver "nleqslv".
#' @param subset_2 A vector of 1:(p-2), which is the subset of \eqn{\hat{Z}_{21}}, i.e.
#'      \eqn{\hat{Z}_{21}^\star} in equation (10) of Huang(2014). \eqn{\hat{Z}_l} may be highly correlated with Z_d, so it is
#'         removed in the estimation.
#'       For the view of including more information, you can use the whole dataset.
#' @param subset_4 A vector of 1:(p-2), which is the subset of \eqn{\hat{Z}_{22}}.
#' @param p number of parameters.
#' @return A list of estimator and its standard deviation.
#'
#' @export
#'
#' @details The function solves estimating equation based on (N2,M1) and (N2,M2), see Huang(2014).
#'
#' @examples
#'  #p <- 8
#'  #subset_2 <- 1:p
#'  #subset_4 <- 1:p
#'  #beta0=c(-5.4163,0.7790,-0.1289,0.2773,-0.5510,0.1568,0.4353,-0.6895)
#'  #DA_FDN2M1M2(realdata_covariates,realdata_alpha,subset_2,subset_4,p=p,beta0=beta0)
#'
#' @references
#' Huang, H., Ma, X., Waagepetersen, R., Holford, T.R. , Wang, R., Risch, H., Mueller, L. & Guan, Y. (2014). A New Estimation Approach for Combining Epidemiological Data From Multiple Sources, Journal of the American Statistical Association, 109:505, 11-23.


DA_FDN2M1M2 <- function(realdata_covariates,realdata_alpha,subset_2,subset_4,p,beta0) {
  #use initial parameter beta1 for more convergence
  # beta0=c(-5.556564, 0.9106478, -0.05683087, 0.318503 ,-0.4461375 ,0.2002597 ,0.4306426, -0.6645185);
  # beta0=c(-5.4163,0.7790,-0.1289,0.2773,-0.5510,0.1568,0.4353,-0.6895);

  ## set the sampling probabilities
  prob_cont_1 <- realdata_alpha$prob_cont_1;
  prob_cont_2 <- realdata_alpha$prob_cont_2;
  prob_case_2 <- realdata_alpha$prob_case_2;
  prob_case_22 <- realdata_alpha$prob_case_22;
  pwt_cont_2 <- realdata_alpha$pwt_cont_2;

  #Extract covariates
  CONTZ_1 <- realdata_covariates$CONTZ_1;
  CONTZhat_1 <- realdata_covariates$CONTZhat_1;

  CASEZ_2 <- realdata_covariates$CASEZ_2;
  CASEZhat_2 <- realdata_covariates$CASEZhat_2;
  CASEZhat_22 <- realdata_covariates$CASEZhat_22;

  CONTZ_2 <- realdata_covariates$CONTZ_2;
  CONTZhat_2 <- realdata_covariates$CONTZhat_2;
  CONTZhat_22 <- realdata_covariates$CONTZhat_22;



  dim(beta0) <- c(1,p);

  #DA_FDM2 uses objfun_1011_u() in AnalysisR folder
  f <- function(beta) {
    Y <- DA_FDN2M1M2_inside(beta,CASEZ_2,CASEZhat_2,CASEZhat_22,CONTZ_1,CONTZhat_1,CONTZhat_2,CONTZhat_22,prob_case_2,
                        prob_case_22,prob_cont_1,prob_cont_2,p,c(),c(),subset_2,subset_4,pwt_cont_2);
    Y[[1]]
  }
  reg <- nleqslv::nleqslv(beta0,f,method="Newton");
  betahat <- reg[[1]];
  betahat

  fval <- reg[[2]];
  flag <- reg[[3]];
  if (flag !=1 ){
    reg_tmp <- DA_FDN2M1M2_inside(beta0,CASEZ_2,CASEZhat_2,CASEZhat_22,CONTZ_1,CONTZhat_1,CONTZhat_2,CONTZhat_22,prob_case_2,
                              prob_case_22,prob_cont_1,prob_cont_2,p,c(),c(),subset_2,subset_4,pwt_cont_2);
    f <- reg_tmp[[1]];
    J_u <- reg_tmp[[2]];
    V_u <- reg_tmp[[3]];
    f2 <- function(beta){
      Y2 <- DA_FDN2M1M2_inside(beta,CASEZ_2,CASEZhat_2,CASEZhat_22,CONTZ_1,CONTZhat_1,CONTZhat_2,CONTZhat_22,prob_case_2,
                           prob_case_22,prob_cont_1,prob_cont_2,p,J_u,V_u,subset_2,subset_4,pwt_cont_2);
      Y2[[1]]
    }
    reg_update <- nleqslv::nleqslv(beta0,f2,method="Newton");
    betahat <- reg_update[[1]];
    fval_update <- reg_update[[2]];
    flag_update <- reg_update[[3]];
  }
  beta1 <- betahat;
  dim(beta1) <- c(1,p);
  plugin <- DA_FDN2M1M2_inside(beta1,CASEZ_2,CASEZhat_2,CASEZhat_22,CONTZ_1,CONTZhat_1,CONTZhat_2,CONTZhat_22,prob_case_2,
                           prob_case_22,prob_cont_1,prob_cont_2,p,c(),c(),subset_2,subset_4,pwt_cont_2);
  f <- plugin[[1]];
  J <- plugin[[2]];
  V <- plugin[[3]];
  std <- sqrt(t(abs(diag(solve(t(J)%*%solve(V)%*%J)))));

  return(list(est=beta1, std=std))
}

