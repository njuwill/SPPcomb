#' Data Analysis for Combining (N1,M1) + (N1,M2) + (N2,M1) + (N2,M2)
#'
#' The main function is to solve the estimating equations constructed by combining all pairs (N1,M1), (N1,M2),
#' (N2,M1) and (N2,M2) with selection bias probility \eqn{\pi(s,\eta)}
#' included.
#'
#' @param realdata_covariates a list contains the following data matrics:
#'        CASEZ_1, CASEZ_2, CASEZhat_1, CASEZhat_2, CASEZhat_22,
#'        CONTZ_1, CONTZ_2, CONTZhat_1, CONTZhat_2, CONTZhat_22. For details
#'        please see definition in the help of realdata_covariates. Please be noted
#'        that all the variables have to use the same name as listed above.
#' @param realdata_alpha a list contains the following data matrics:
#'        prob_case_1, prob_case_11, prob_case_2, prob_case_22, prob_cont_1,
#'        prob_cont_2, pwt_cont_2. details please see definition in the
#'        help of realdata_alpha. Please be noted that all the variables have to use the same name as listed above.
#' @param beta0 We need an initial parameter for solver "nleqslv". Default value is beta0=c(-5.4163,0.7790,-0.1289,0.2773,-0.5510,0.1568,0.4353,-0.6895)
#'
#' @return A list of estimator and its standard deviation.
#'
#' @export
#'
#' @details The function solves GMM combined estimating equation
#' with handling selection bias, see Huang(2014).
#'
#' @examples
#' #you can use glm to get the estimate as the initial value of beta0
#' #beta0=c(-5.4163,0.7790,-0.1289,0.2773,-0.5510,0.1568,0.4353,-0.6895)
#' #DA_AllEE5(realdata_covariates,realdata_alpha,beta0=beta0)
#'
#' @references
#' Huang, H., Ma, X., Waagepetersen, R., Holford, T.R. , Wang, R., Risch, H., Mueller, L. & Guan, Y. (2014). A New Estimation Approach for Combining Epidemiological Data From Multiple Sources, Journal of the American Statistical Association, 109:505, 11-23.


DA_AllEE5 <- function(realdata_covariates,realdata_alpha,beta0) {
  #use initial parameter beta0
  #beta0=c(-5.4163,0.7790,-0.1289,0.2773,-0.5510,0.1568,0.4353,-0.6895);
  p <- 8
  subset_2 <- 1:(p-2)
  subset_3 <- 1:p
  subset_4 <- 1:(p-2)
  ## set the sampling probabilities
  prob_cont_1 <- realdata_alpha$prob_cont_1;
  prob_cont_2 <- realdata_alpha$prob_cont_2;
  prob_case_1 <- realdata_alpha$prob_case_1;
  prob_case_11 <- realdata_alpha$prob_case_11;
  prob_case_2 <- realdata_alpha$prob_case_2;
  prob_case_22 <- realdata_alpha$prob_case_22;
  pwt_cont_2 <- realdata_alpha$pwt_cont_2;

  #Extract covariates
  CASEZ_1 <- realdata_covariates$CASEZ_1;
  CASEZhat_1 <- realdata_covariates$CASEZhat_1;

  CONTZ_1 <- realdata_covariates$CONTZ_1;
  CONTZhat_1 <- realdata_covariates$CONTZhat_1;

  CASEZ_2 <- realdata_covariates$CASEZ_2;
  CASEZhat_2 <- realdata_covariates$CASEZhat_2;
  CASEZhat_22 <- realdata_covariates$CASEZhat_22;

  CONTZ_2 <- realdata_covariates$CONTZ_2;
  CONTZhat_2 <- realdata_covariates$CONTZhat_2;
  CONTZhat_22 <- realdata_covariates$CONTZhat_22;

  #Extract additional data
  Z_case_pi_1 <-as.matrix(CASEZ_1[,1:5])
  Z_case_pi_2 <- CASEZ_2[,1:5]
  Z_cont_pi_1 <-CONTZ_1[,1:5]
  Z_cont_pi_2 <-CONTZhat_2[,1:5]

  Z_case_pi_1_t <-as.matrix(CASEZ_1[,1:6])
  Z_case_pi_2_t <- CASEZ_2[,1:6]
  Z_cont_pi_1_t <-CONTZ_1[,1:6]
  Z_cont_pi_2_t <-CONTZhat_2[,1:6]

  # estimate eta
  Y=c(rep(1,dim(CASEZ_1)[1]), rep(0,dim(CASEZ_2)[1]))
  X=rbind(CASEZ_1[,2:6],CASEZ_2[,2:6])
  fit1 <- stats::glm(Y~X[,1]+X[,2]+X[,3]+X[,4],family=stats::binomial())
  etahat1=fit1$coefficients
  pi_case_1=exp(Z_case_pi_1_t[,1:5]%*%etahat1)/(1+exp(Z_case_pi_1_t[,1:5]%*%etahat1))
  pi_case_2=exp(Z_case_pi_2_t[,1:5]%*%etahat1)/(1+exp(Z_case_pi_2_t[,1:5]%*%etahat1))
  pi_cont_1=exp(Z_cont_pi_1_t[,1:5]%*%etahat1)/(1+exp(Z_cont_pi_1_t[,1:5]%*%etahat1))
  pi_cont_2=exp(Z_cont_pi_2_t[,1:5]%*%etahat1)/(1+exp(Z_cont_pi_2_t[,1:5]%*%etahat1))

  fit2 <- stats::glm(Y~X[,1]+X[,2]+X[,3]+X[,4]+X[,5],family=stats::binomial())
  etahat2=fit2$coefficients
  pi_case_1_t=exp(Z_case_pi_1_t%*%etahat2)/(1+exp(Z_case_pi_1_t%*%etahat2))
  pi_case_2_t=exp(Z_case_pi_2_t%*%etahat2)/(1+exp(Z_case_pi_2_t%*%etahat2))
  pi_cont_1_t=exp(Z_cont_pi_1_t%*%etahat2)/(1+exp(Z_cont_pi_1_t%*%etahat2))
  pi_cont_2_t=exp(Z_cont_pi_2_t%*%etahat2)/(1+exp(Z_cont_pi_2_t%*%etahat2))



  dim(beta0) <- c(1,p);

  f <- function(beta){
    Y <- DA_AllEE5_inside(beta,CASEZ_1,CASEZ_2,CASEZhat_1,CASEZhat_2,CASEZhat_22,
                          CONTZ_1,CONTZ_2,CONTZhat_1,CONTZhat_2,CONTZhat_22,
                          prob_case_1,prob_case_11,prob_case_2,prob_case_22,prob_cont_1,
                          prob_cont_2,p,pi_case_1,pi_case_1_t,pi_case_2,pi_case_2_t,
                          pi_cont_1,pi_cont_1_t,pi_cont_2,Z_case_pi_1,Z_case_pi_1_t,
                          Z_case_pi_2,Z_case_pi_2_t,Z_cont_pi_1,Z_cont_pi_1_t,Z_cont_pi_2,
                          c(),c(),pwt_cont_2,subset_2,subset_3,subset_4);
    Y[[1]]
  }
  reg <- nleqslv::nleqslv(beta0,f,method="Newton");
  betahat <- reg[[1]];
  fval <- reg[[2]];
  flag <- reg[[3]];
  if (flag !=1){
    reg_tmp <- DA_AllEE5_inside(beta0,CASEZ_1,CASEZ_2,CASEZhat_1,CASEZhat_2,CASEZhat_22,
                                CONTZ_1,CONTZ_2,CONTZhat_1,CONTZhat_2,CONTZhat_22, prob_case_1,
                                prob_case_11,prob_case_2,prob_case_22,prob_cont_1,prob_cont_2,
                                p,pi_case_1,pi_case_1_t,pi_case_2,pi_case_2_t,pi_cont_1,pi_cont_1_t,
                                pi_cont_2,Z_case_pi_1,Z_case_pi_1_t,Z_case_pi_2,Z_case_pi_2_t,
                                Z_cont_pi_1,Z_cont_pi_1_t,Z_cont_pi_2,c(),c(),pwt_cont_2,subset_2,subset_3,subset_4);
    f_tmp <- reg_tmp[[1]];
    J_u_tmp <- reg_tmp[[2]];
    V_u_tmp <- reg_tmp[[3]];

    f2 <- function(beta){
      Y2 <- DA_AllEE5_inside(beta,CASEZ_1,CASEZ_2,CASEZhat_1,CASEZhat_2,CASEZhat_22,CONTZ_1,
                          CONTZ_2,CONTZhat_1,CONTZhat_2,CONTZhat_22, prob_case_1,prob_case_11,
                          prob_case_2,prob_case_22,prob_cont_1,prob_cont_2,p,pi_case_1,
                          pi_case_1_t,pi_case_2,pi_case_2_t,pi_cont_1,pi_cont_1_t,pi_cont_2,
                          Z_case_pi_1,Z_case_pi_1_t,Z_case_pi_2,Z_case_pi_2_t,Z_cont_pi_1,
                          Z_cont_pi_1_t,Z_cont_pi_2,J_u_tmp,V_u_tmp,pwt_cont_2,subset_2,subset_3,subset_4);
      Y2[[1]]
    }
    reg_update <- nleqslv::nleqslv(beta0,f2,method="Newton");
    betahat <- reg_update[[1]];
    fval <- reg_update[[2]];
    flag <- reg_update[[3]];
  }
  beta1 <- betahat;
  dim(beta1)=c(1,p);
  plugin <- DA_AllEE5_inside(beta1,CASEZ_1,CASEZ_2,CASEZhat_1,CASEZhat_2,CASEZhat_22,
                             CONTZ_1,CONTZ_2,CONTZhat_1,CONTZhat_2,CONTZhat_22, prob_case_1,
                             prob_case_11,prob_case_2,prob_case_22,prob_cont_1,prob_cont_2,p,
                             pi_case_1,pi_case_1_t,pi_case_2,pi_case_2_t,pi_cont_1,pi_cont_1_t,
                             pi_cont_2,Z_case_pi_1,Z_case_pi_1_t,Z_case_pi_2,Z_case_pi_2_t,
                             Z_cont_pi_1,Z_cont_pi_1_t,Z_cont_pi_2,c(),c(),pwt_cont_2,subset_2,subset_3,subset_4);
  f <- plugin[[1]];
  J <- plugin[[2]];
  V <- plugin[[3]];

  std <- sqrt(t(abs(diag(solve(t(J)%*%solve(V)%*%J)))));


  return(list(est=beta1, std=std))
}

