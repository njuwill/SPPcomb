#' Data Analysis of (N2,M2)
#'
#' The main function to solve the estimating equations constructed by (N2,M2). Since there is just one case data, no selection bias needed.
#'
#' @param realdata_covariates a list contains the following data matrics:
#'        CASEZhat_2, CASEZhat_22, CONTZhat_2, CONTZhat_22
#' @param realdata_alpha a list contains the following data matrics:
#'        prob_case_22,prob_cont_2, pwt_cont_2
#' @param beta0 an initial parameter for solver "nleqslv".
#' @param p number of parameters.
#' @return A list of estimator and its standard deviation.
#'
#' @export
#'
#' @details The function solves estimating equation based on (N2,M2), see Huang(2014).
#'
#' @examples
#'  #p <- 8
#'  #beta0=c(-5.4163,0.7790,-0.1289,0.2773,-0.5510,0.1568,0.4353,-0.6895)
#'  #DA_FDN2M2(realdata_covariates,realdata_alpha,p=p,beta0=beta0)
#'
#' @references
#' Huang, H., Ma, X., Waagepetersen, R., Holford, T.R. , Wang, R., Risch, H., Mueller, L. & Guan, Y. (2014). A New Estimation Approach for Combining Epidemiological Data From Multiple Sources, Journal of the American Statistical Association, 109:505, 11-23.



DA_FDN2M2<- function(realdata_covariates,realdata_alpha,p,beta0) {
    #use initial parameter beta1 for more convergence
    # beta0=c(-5.556564, 0.9106478, -0.05683087, 0.318503 ,-0.4461375 ,0.2002597 ,0.4306426, -0.6645185);
    # beta0=c(-5.4163,0.7790,-0.1289,0.2773,-0.5510,0.1568,0.4353,-0.6895);
    # subset_4 <- 1:(p-2)
    ## set the sampling probabilities

  DA_FDN2M2_inside<- function(beta,CASEZhat_2,CASEZhat_22,CONTZhat_2,CONTZhat_22,
                              prob_case_22,prob_cont_2,p,J,V,pwt_cont_2){


    W_case_4=prob_case_22/(prob_case_22+exp(as.matrix(CASEZhat_22)%*%t(beta)))*exp((as.matrix(CASEZhat_22)-as.matrix(CASEZhat_2))%*%t(beta));
    W_cont_4=exp(as.matrix(CONTZhat_2)%*%t(beta))/(prob_cont_2+exp(as.matrix(CONTZhat_22)%*%t(beta)))*prob_cont_2/pwt_cont_2;

    Q_4_cont=colSums(as.matrix(CONTZhat_22)*kronecker(as.matrix(W_cont_4),t(as.matrix(rep(1,p)))));
    Q_4_case=colSums(as.matrix(CASEZhat_22)*kronecker(as.matrix(W_case_4),t(as.matrix(rep(1,p)))));
    Q_4=Q_4_case-Q_4_cont;
    dim(Q_4)=c(p,1);
    Q=Q_4

    if (length(J)+length(V)==0){
      J_4=matrix(0,p,p);
      for (i in 1:p){
        for (j in 1:p){
          ##      J_4[i,j]=sum(CASEZhat_22[,subset_4[i]]*CASEZhat_2[,j]*W_case_4);
          J_4[i,j]=-sum(CASEZhat_22[,i]*(CASEZhat_2[,j]-CASEZhat_22[,j])*W_case_4)-sum(CONTZhat_22[,i]*(CONTZhat_2[,j])*W_cont_4);
        }
      }

      VQQ_44=matrix(0,p,p);
      for (i in 1:p){
        for (j in 1:p){
          VQQ_44[i,j]=sum(CASEZhat_22[,i]*CASEZhat_22[,j]*W_case_4^2)+sum(CONTZhat_22[,i]*CONTZhat_22[,j]*W_cont_4^2);
        }
      }

      V=VQQ_44

      J=J_4
    }
    f=t(J)%*%solve(V)%*%Q;
    return(list(f,J,V))
  }


    prob_cont_2 <- realdata_alpha$prob_cont_2;
    prob_case_22 <- realdata_alpha$prob_case_22;
    pwt_cont_2 <- realdata_alpha$pwt_cont_2;

    #Extract covariates
    CASEZhat_2 <- realdata_covariates$CASEZhat_2;
    CASEZhat_22 <- realdata_covariates$CASEZhat_22;

    CONTZhat_2 <- realdata_covariates$CONTZhat_2;
    CONTZhat_22 <- realdata_covariates$CONTZhat_22;

    dim(beta0) <- c(1,p);

    #DA_FDM2 uses objfun_1011_u() in AnalysisR folder
    f <- function(beta) {
      Y <- DA_FDN2M2_inside(beta,CASEZhat_2,CASEZhat_22,CONTZhat_2,CONTZhat_22,
                              prob_case_22,prob_cont_2,p,c(),c(),pwt_cont_2);
      Y[[1]]
    }
    reg <- nleqslv::nleqslv(beta0,f,method="Newton");
    betahat <- reg[[1]];
    betahat

    fval <- reg[[2]];
    flag <- reg[[3]];
    if (flag !=1 ){
      reg_tmp <- DA_FDN2M2_inside(beta0,CASEZhat_2,CASEZhat_22,CONTZhat_2,CONTZhat_22,
                                    prob_case_22,prob_cont_2,p,c(),c(),pwt_cont_2);
      f <- reg_tmp[[1]];
      J_u <- reg_tmp[[2]];
      V_u <- reg_tmp[[3]];
      f2 <- function(beta){
        Y2 <- DA_FDN2M2_inside(beta,CASEZhat_2,CASEZhat_22,CONTZhat_2,CONTZhat_22,
                                 prob_case_22,prob_cont_2,p,J_u,V_u,pwt_cont_2);
        Y2[[1]]
      }
      reg_update <- nleqslv::nleqslv(beta0,f2,method="Newton");
      betahat <- reg_update[[1]];
      fval_update <- reg_update[[2]];
      flag_update <- reg_update[[3]];
    }
    beta1 <- betahat;
    dim(beta1) <- c(1,p);
    plugin <- DA_FDN2M2_inside(beta1,CASEZhat_2,CASEZhat_22,CONTZhat_2,CONTZhat_22,
                                 prob_case_22,prob_cont_2,p,c(),c(),pwt_cont_2);
    J <- plugin[[2]];
    V <- plugin[[3]];
    # std <- sqrt(t(abs(diag(solve(t(J)%*%solve(V)%*%J)))));
    std <- sqrt(t(abs(diag((solve(J)%*%V%*%solve(t(J)))))));

    return(list(est=beta1, std=std))
  }

