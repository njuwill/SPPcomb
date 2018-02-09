#' Internal calculation of estimating equation for DA_FDN2M1M2
#'
#' The internal function to solve the estimating equations constructed by combining pair (N2,M1) and (N2,M2). Since there is just one case data, no selection bias needed.
#' Since it's a internal function for function DA_FDN2M1M2, thus it's not a necessary
#' or important function.
#'
#' @param beta Parameter \eqn{\beta}.
#' @param CONTZ_1,CONTZhat_1  control data(M1) from case-control study, details please
#'        see definition in the help of realdata_covariates.
#' @param CONTZhat_2,CONTZhat_22 BRFSS data(M2), details please see definition in the help of realdata_covariates.
#' @param CASEZ_2,CASEZhat_2,CASEZhat_22 CTR data(N2), details please see definition in the help of realdata_covariates.
#' @param prob_cont_1,prob_cont_2,prob_case_2,prob_case_22,pwt_cont_2   please see definition in the help of realdata_alpha.
#' @param subset_2 A vector of 1:(p-2).
#' @param subset_4 A vector of 1:(p-2).
#' @param J The derivative of the estimating equation.
#' @param V The variance of the estimating equation.
#' @param p Number of parameters, a constant value of 8.
#'
#' @return A list of (f,J,V)
#' \enumerate{
#' \item f The final form of the estimating equation after adjusting eta.
#' \item J The derivative of the estimating equation.
#' \item V The variance of the estimating equation.
#' }
#'
#' @export
#'
#' @details The function solves estimating equation based on (N2,M1) and (N2, M2) with handling selection bias.
#' It also accounts for the uncertainty due to the estimated value of eta. The function will
#' output the estimating equation at current input value beta. Hence it can be used in "nleqslv" to solve for
#' \eqn{\beta}. Because the function also outputs J and V, the asymptotic variance of \eqn{\beta}
#' can be calculated in a straightforward way. \eqn{\hat{Z}_l} may be highly correlated with Z_d, so it is
#' removed in the estimation.
#'


DA_FDN2M1M2_inside <- function(beta,CASEZ_2,CASEZhat_2,CASEZhat_22,CONTZ_1,CONTZhat_1,CONTZhat_2,CONTZhat_22,prob_case_2,
                           prob_case_22,prob_cont_1,prob_cont_2,p,J,V,subset_2,subset_4,pwt_cont_2){


  lthsub_2=length(subset_2);
  # lthsub_3=length(subset_3);
  lthsub_4=length(subset_4);

  ##  get the derivatives

  W_case_2=prob_case_2/(prob_case_2+exp(as.matrix(CASEZhat_2)%*%t(beta)));
  W_cont_2=exp(as.matrix(CONTZ_1)%*%t(beta))/(prob_cont_1+exp(as.matrix(CONTZhat_1)%*%t(beta)));

  W_case_4=prob_case_22/(prob_case_22+exp(as.matrix(CASEZhat_22)%*%t(beta)))*exp((as.matrix(CASEZhat_22)-as.matrix(CASEZhat_2))%*%t(beta));
  W_cont_4=exp(as.matrix(CONTZhat_2)%*%t(beta))/(prob_cont_2+exp(as.matrix(CONTZhat_22)%*%t(beta)))*prob_cont_2/pwt_cont_2;

  ## get the derivative
  Q_2_cont=colSums(as.matrix(CONTZhat_1)[,subset_2]*kronecker(as.matrix(W_cont_2),t(as.matrix(rep(1,lthsub_2)))));
  Q_2_case=colSums(as.matrix(CASEZhat_2)[,subset_2]*kronecker(as.matrix(W_case_2),t(as.matrix(rep(1,lthsub_2)))));
  Q_2=Q_2_case-Q_2_cont;
  dim(Q_2)=c(lthsub_2,1);

  Q_4_cont=colSums(as.matrix(CONTZhat_22)[,subset_4]*kronecker(as.matrix(W_cont_4),t(as.matrix(rep(1,lthsub_4)))));
  Q_4_case=colSums(as.matrix(CASEZhat_22)[,subset_4]*kronecker(as.matrix(W_case_4),t(as.matrix(rep(1,lthsub_4)))));
  Q_4=Q_4_case-Q_4_cont;
  dim(Q_4)=c(lthsub_4,1);

  Q=rbind(Q_2,Q_4);

  if (length(J)+length(V)==0){
    J_2=matrix(0,length(subset_2),p);
    for (i in 1:length(subset_2)){
      for (j in 1:p){
        J_2[i,j]=-sum(CONTZhat_1[,subset_2[i]]*CONTZ_1[,j]*W_cont_2);
      }
    }



    J_4=matrix(0,length(subset_4),p);
    for (i in 1:length(subset_4)){
      for (j in 1:p){
        ##      J_4[i,j]=sum(CASEZhat_22[,subset_4[i]]*CASEZhat_2[,j]*W_case_4);
        J_4[i,j]=-sum(CASEZhat_22[,subset_4[i]]*(CASEZhat_2[,j]-CASEZhat_22[,j])*W_case_4)-sum(CONTZhat_22[,subset_4[i]]*(CONTZhat_2[,j])*W_cont_4);
      }
    }


    VQQ_22=matrix(0,length(subset_2),length(subset_2));
    for (i in 1:length(subset_2)){
      for (j in 1:length(subset_2)){
        VQQ_22[i,j]=sum(CASEZhat_2[,subset_2[i]]*CASEZhat_2[,subset_2[j]]*W_case_2^2)+sum(CONTZhat_1[,subset_2[i]]*CONTZhat_1[,subset_2[j]]*W_cont_2^2);
      }
    }

    VQQ_24=matrix(0,length(subset_2),length(subset_4));
    for (i in 1:length(subset_2)){
      for (j in 1:length(subset_4)){
        VQQ_24[i,j]=sum(CASEZhat_2[,subset_2[i]]*CASEZhat_22[,subset_4[j]]*W_case_2*W_case_4);
      }
    }

    VQQ_44=matrix(0,length(subset_4),length(subset_4));
    for (i in 1:length(subset_4)){
      for (j in 1:length(subset_4)){
        VQQ_44[i,j]=sum(CASEZhat_22[,subset_4[i]]*CASEZhat_22[,subset_4[j]]*W_case_4^2)+sum(CONTZhat_22[,subset_4[i]]*CONTZhat_22[,subset_4[j]]*W_cont_4^2);
      }
    }

    VQQ=rbind(cbind(VQQ_22,VQQ_24),cbind(t(VQQ_24),VQQ_44));


    V=VQQ;
    J=rbind(J_2,J_4);
  }

  f=t(J)%*%solve(V)%*%Q;
  return(list(f,J,V))
}
