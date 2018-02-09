#' Internal calculation of estimating equation for DA_AllEE5
#'
#' This is the internal function to solve the estimating equation constructed by
#' pair (N1,M1), (N1,M2), (N2,M1) and (N2,M2) with selection bias probility \eqn{\pi(s,\eta)}
#' included. Since it's a internal function for function DA_AllEE5, thus it's not a necessary
#' or important function.
#'
#' @param beta Parameter \eqn{\beta}.
#' @param CONTZ_1,CONTZhat_1  control data(M1) from case-control study, details please
#'        see definition in the help of realdata_covariates.
#' @param CONTZ_2,CONTZhat_2,CONTZhat_22 BRFSS data(M2), details please see definition in the help of realdata_covariates.
#' @param CASEZ_1,CASEZhat_1 case data(N1) from case-control study,
#'        details please see definition in the help of realdata_covariates.
#' @param CASEZ_2,CASEZhat_2,CASEZhat_22 CTR data(N2), details please see definition in the help of realdata_covariates.
#' @param prob_cont_1,prob_cont_2,prob_case_1,prob_case_11,prob_case_2,prob_case_22,pwt_cont_2   please see definition in the help of realdata_alpha.
#' @param pi_case_1,pi_case_1_t,pi_case_2,pi_case_2_t,pi_cont_1,pi_cont_1_t,pi_cont_2 selection bias
#' @param Z_case_pi_1,Z_case_pi_1_t,Z_case_pi_2,Z_case_pi_2_t,Z_cont_pi_1,Z_cont_pi_1_t,Z_cont_pi_2 part
#'       of variables from covariates, used for the estiamtion of variance.
#' @param J_step3 Derivative of the estimating equation.
#' @param V_step3 Variance of the estimating equation.
#' @param subset_2 A vector of 1:(p-2).
#' @param subset_3 A vector of 1:p.
#' @param subset_4 A vector of 1:(p-2).
#' @param p Number of parameters, a constant value of 8.
#'
#' @return A list of (f,J,V)
#' \enumerate{
#' \item f The final form of the estimating equation after adjusting eta.
#' \item J_step3 The derivative of the estimating equation.
#' \item V_step3 The variance of the estimating equation.
#' }
#'
#' @export
#'
#' @details The function solves estimating equation based on GMM combined estimating equations
#' with handling selection bias.
#' It also accounts for the uncertainty due to the estimated value of eta. The function will
#' output the estimating equation at current input value beta. Hence it can be used in "nleqslv" to solve for
#' \eqn{\beta}. Because the function also outputs J and V, the asymptotic variance of \eqn{\beta}
#' can be calculated in a straightforward way. \eqn{\hat{Z}_l} may be highly correlated with Z_d, so it is
#' removed in the estimation. And it has to be careful in constructing f, J and V.
#'

DA_AllEE5_inside <- function(beta,CASEZ_1,CASEZ_2,CASEZhat_1,CASEZhat_2,CASEZhat_22,
                             CONTZ_1,CONTZ_2,CONTZhat_1,CONTZhat_2,CONTZhat_22, prob_case_1,
                             prob_case_11,prob_case_2,prob_case_22,prob_cont_1,prob_cont_2,p,
                             pi_case_1,pi_case_1_t,pi_case_2,pi_case_2_t,pi_cont_1,pi_cont_1_t,
                             pi_cont_2,Z_case_pi_1,Z_case_pi_1_t,Z_case_pi_2,Z_case_pi_2_t,
                             Z_cont_pi_1,Z_cont_pi_1_t,Z_cont_pi_2, J_step3,V_step3,pwt_cont_2,
                             subset_2,subset_3,subset_4){


  lthsub_2=length(subset_2);
  lthsub_3=length(subset_3);
  lthsub_4=length(subset_4);

  # p_pi=dim(Z_case_pi_1)[2];   ## number of parameter in eta_d (no traffic)
  p_pi_t=dim(Z_case_pi_1_t)[2];   ## number of parameter in eta (including traffic)


  ##get the Q ( Q_1 is U_11, Q_2 is U_21, Q_3 is U_12, Q_4 is U_22 in the paper)
  W_case_1=prob_case_1/(prob_case_1+exp(as.matrix(CASEZ_1)%*%t(beta))*pi_case_1_t);
  W_cont_1=pi_cont_1_t*exp(as.matrix(CONTZ_1)%*%t(beta))/(prob_cont_1+exp(as.matrix(CONTZ_1)%*%t(beta))*pi_cont_1_t);

  W_case_2=prob_case_2/(prob_case_2+exp(as.matrix(CASEZhat_2)%*%t(beta))*(1-pi_case_2_t));
  W_cont_2=(1-pi_cont_1_t)*exp(as.matrix(CONTZ_1)%*%t(beta))/(prob_cont_1+exp(as.matrix(CONTZhat_1)%*%t(beta))*(1-pi_cont_1_t));

  W_case_3=prob_case_11/(prob_case_11+exp(as.matrix(CASEZhat_1)%*%t(beta))*pi_case_1)*exp((as.matrix(CASEZhat_1)-as.matrix(CASEZ_1))%*%t(beta))*pi_case_1/pi_case_1_t;
  W_cont_3=pi_cont_2*exp(as.matrix(CONTZhat_2)%*%t(beta))/(prob_cont_2+exp(as.matrix(CONTZhat_2)%*%t(beta))*pi_cont_2)*prob_cont_2/pwt_cont_2;

  W_case_4=prob_case_22/(prob_case_22+exp(as.matrix(CASEZhat_22)%*%t(beta))*(1-pi_case_2))*exp((as.matrix(CASEZhat_22)-as.matrix(CASEZhat_2))%*%t(beta))*(1-pi_case_2)/(1-pi_case_2_t);
  W_cont_4=(1-pi_cont_2)*exp(as.matrix(CONTZhat_2)%*%t(beta))/(prob_cont_2+exp(as.matrix(CONTZhat_22)%*%t(beta))*(1-pi_cont_2))*prob_cont_2/pwt_cont_2;

  Q_1_cont=colSums(as.matrix(CONTZ_1)*kronecker(as.matrix(W_cont_1),t(as.matrix(rep(1,p)))));
  Q_1_case=colSums(as.matrix(CASEZ_1)*kronecker(as.matrix(W_case_1),t(as.matrix(rep(1,p)))));
  Q_1=Q_1_case-Q_1_cont;
  dim(Q_1)=c(p,1);
  Q_2_cont=colSums(as.matrix(CONTZhat_1)[,subset_2]*kronecker(as.matrix(W_cont_2),t(as.matrix(rep(1,lthsub_2)))));
  Q_2_case=colSums(as.matrix(CASEZhat_2)[,subset_2]*kronecker(as.matrix(W_case_2),t(as.matrix(rep(1,lthsub_2)))));
  Q_2=Q_2_case-Q_2_cont;
  dim(Q_2)=c(lthsub_2,1);
  Q_3_cont=colSums(as.matrix(CONTZhat_2)[,subset_3]*kronecker(as.matrix(W_cont_3),t(as.matrix(rep(1,lthsub_3)))));
  Q_3_case=colSums(as.matrix(CASEZhat_1)[,subset_3]*kronecker(as.matrix(W_case_3),t(as.matrix(rep(1,lthsub_3)))));
  Q_3=Q_3_case-Q_3_cont;
  dim(Q_3)=c(lthsub_3,1);
  Q_4_cont=colSums(as.matrix(CONTZhat_22)[,subset_4]*kronecker(as.matrix(W_cont_4),t(as.matrix(rep(1,lthsub_4)))));
  Q_4_case=colSums(as.matrix(CASEZhat_22)[,subset_4]*kronecker(as.matrix(W_case_4),t(as.matrix(rep(1,lthsub_4)))));
  Q_4=Q_4_case-Q_4_cont;
  dim(Q_4)=c(lthsub_4,1);




  ## get the derivative of Q in terms of beta, J is D(\beta) in the paper

  J_1=matrix(0,p,p);
  for (i in 1:p){
    for (j in 1:p){
      J_1[i,j]=-sum(CONTZ_1[,i]*CONTZ_1[,j]*W_cont_1);
    }
  }

  J_2=matrix(0,length(subset_2),p);
  for (i in 1:length(subset_2)){
    for (j in 1:p){
      J_2[i,j]=-sum(CONTZhat_1[,subset_2[i]]*CONTZ_1[,j]*W_cont_2);
    }
  }

  J_3=matrix(0,length(subset_3),p);
  for (i in 1:length(subset_3)){
    for (j in 1:p){
      J_3[i,j]=-sum(CASEZhat_1[,subset_3[i]]*CASEZ_1[,j]*W_case_3);
    }
  }

  J_4=matrix(0,length(subset_4),p);
  for (i in 1:length(subset_4)){
    for (j in 1:p){
      ##      J_4[i,j]=sum(CASEZhat_22[,subset_4[i]]*CASEZhat_2[,j]*W_case_4);
      J_4[i,j]=-sum(CASEZhat_22[,subset_4[i]]*(CASEZhat_2[,j]-CASEZhat_22[,j])*W_case_4)-sum(CONTZhat_22[,subset_4[i]]*(CONTZhat_2[,j])*W_cont_4);
    }
  }






  ##get the naive estimate of Var(Q), VQQ is V(\beta, \eta_0) in the paper

  VQQ_11=matrix(0,p,p);
  for (i in 1:p){
    for (j in 1:p){
      VQQ_11[i,j]=sum(CASEZ_1[,i]*CASEZ_1[,j]*W_case_1^2)+sum(CONTZ_1[,i]*CONTZ_1[,j]*W_cont_1^2);
    }
  }

  VQQ_22=matrix(0,length(subset_2),length(subset_2));
  for (i in 1:length(subset_2)){
    for (j in 1:length(subset_2)){
      VQQ_22[i,j]=sum(CASEZhat_2[,subset_2[i]]*CASEZ_2[,subset_2[j]]*W_case_2^2)+sum(CONTZhat_1[,subset_2[i]]*CONTZhat_1[,subset_2[j]]*W_cont_2^2);
    }
  }

  VQQ_33=matrix(0,length(subset_3),length(subset_3));
  for (i in 1:length(subset_3)){
    for (j in 1:length(subset_3)){
      VQQ_33[i,j]=sum(CASEZhat_1[,subset_3[i]]*CASEZhat_1[,subset_3[j]]*W_case_3^2)+sum(CONTZhat_2[,subset_3[i]]*CONTZhat_2[,subset_3[j]]*W_cont_3^2);
    }
  }

  VQQ_12=matrix(0,p,length(subset_2));
  for (i in 1:p){
    for (j in 1:length(subset_2)){
      VQQ_12[i,j]=sum(CONTZ_1[,i]*CONTZhat_1[,subset_2[j]]*W_cont_1*W_cont_2);
    }
  }

  VQQ_13=matrix(0,p,length(subset_3));
  for (i in 1:p){
    for (j in 1:length(subset_3)){
      VQQ_13[i,j]=sum(CASEZ_1[,i]*CASEZhat_1[,subset_3[j]]*W_case_1*W_case_3);
    }
  }

  VQQ_24=matrix(0,length(subset_2),length(subset_4));
  for (i in 1:length(subset_2)){
    for (j in 1:length(subset_4)){
      VQQ_24[i,j]=sum(CASEZhat_2[,subset_2[i]]*CASEZhat_22[,subset_4[j]]*W_case_2*W_case_4);
    }
  }

  VQQ_34=matrix(0,length(subset_3),length(subset_4));
  for (i in 1:length(subset_3)){
    for (j in 1:length(subset_4)){
      VQQ_34[i,j]=sum(CONTZhat_2[,subset_3[i]]*CONTZhat_22[,subset_4[j]]*W_cont_3*W_cont_4);
    }
  }

  VQQ_44=matrix(0,length(subset_4),length(subset_4));
  for (i in 1:length(subset_4)){
    for (j in 1:length(subset_4)){
      VQQ_44[i,j]=sum(CASEZhat_22[,subset_4[i]]*CASEZhat_22[,subset_4[j]]*W_case_4^2)+sum(CONTZhat_22[,subset_4[i]]*CONTZhat_22[,subset_4[j]]*W_cont_4^2);
    }
  }


  VQQ_23=matrix(0,length(subset_2),length(subset_3));
  VQQ_14=matrix(0,p,length(subset_4));
  VQQ=rbind(cbind(VQQ_11,VQQ_13,VQQ_12,VQQ_14),cbind(t(VQQ_13),VQQ_33,t(VQQ_23),VQQ_34),cbind(t(VQQ_12),VQQ_23,VQQ_22,VQQ_24),cbind(t(VQQ_14),t(VQQ_34),t(VQQ_24),VQQ_44));


  ## Get partial derivative of Q in terms of eta



  H_cont_1=pi_cont_1_t*(1-pi_cont_1_t)*exp(as.matrix(CONTZ_1)%*%t(beta))/(prob_cont_1+exp(as.matrix(CONTZ_1)%*%t(beta))*pi_cont_1_t);
  H_cont_2=(1-pi_cont_1_t)*pi_cont_1_t*exp(as.matrix(CONTZ_1)%*%t(beta))/(prob_cont_1+exp(as.matrix(CONTZhat_1)%*%t(beta))*(1-pi_cont_1_t));
  ## H_cont_3=(1-pi_cont_2)*pi_cont_2*exp(as.matrix(CONTZhat_2)%*%t(beta))/(prob_cont_2+exp(as.matrix(CONTZhat_2)%*%t(beta))*pi_cont_2)*prob_cont_2/pwt_cont_2;
  ## H_cont_4=(1-pi_cont_2)*pi_cont_2*exp(as.matrix(CONTZhat_2)%*%t(beta))/(prob_cont_2+exp(as.matrix(CONTZhat_22)%*%t(beta))*(1-pi_cont_2))*prob_cont_2/pwt_cont_2;
  QE_3=prob_case_11*(1-pi_case_1_t)/(prob_case_11+exp(as.matrix(CASEZhat_1)%*%t(beta))*pi_case_1)*exp((as.matrix(CASEZhat_1)-as.matrix(CASEZ_1))%*%t(beta));
  QE_4=prob_case_22*pi_case_2_t/(prob_case_22+exp(as.matrix(CASEZhat_22)%*%t(beta))*(1-pi_case_2))*exp((as.matrix(CASEZhat_22)-as.matrix(CASEZhat_2))%*%t(beta));

  PQ_1=matrix(0,p,p_pi_t);
  PQ_2=matrix(0,lthsub_2,p_pi_t);
  PQ_3=matrix(0,lthsub_3,p_pi_t);
  PQ_4=matrix(0,lthsub_4,p_pi_t);
  for (i in 1:p_pi_t){
    for (j in 1:p){
      PQ_1[j,i]=-sum(CONTZ_1[,j]*Z_cont_pi_1_t[,i]*H_cont_1);
    }
    for (k in 1:lthsub_2){
      PQ_2[k,i]=sum(CONTZhat_1[,subset_2[k]]*Z_cont_pi_1_t[,i]*H_cont_2);
    }
    for (k in 1:lthsub_3){
      PQ_3[k,i]=-sum(CASEZhat_1[,subset_3[k]]*Z_case_pi_1_t[,i]*QE_3);
    }
    for (l in 1:lthsub_4){
      PQ_4[l,i]=sum(CASEZhat_22[,subset_4[l]]*Z_case_pi_2_t[,i]*QE_4);
    }
  }
  PQ=rbind(PQ_1,PQ_3,PQ_2,PQ_4);


  ## determine var(eta)

  V_eta=matrix(0,p_pi_t,p_pi_t);
  B_eta=matrix(0,p_pi_t,p_pi_t);
  for (i in 1:p_pi_t){
    for (j in 1:p_pi_t){
      V_eta[i,j]=sum(Z_case_pi_1_t[,i]*Z_case_pi_1_t[,j]*(1-pi_case_1_t)^2)+sum(Z_case_pi_2_t[,i]*Z_case_pi_2_t[,j]*pi_case_2_t^2);
      B_eta[i,j]=-sum(Z_case_pi_2_t[,i]*Z_case_pi_2_t[,j]*pi_case_2_t);
    }
  }

  B_eta=solve(B_eta);
  V_eta=B_eta%*%V_eta%*%t(B_eta);

  VPP=PQ%*%V_eta%*%t(PQ);




  ##Determine cov(Q,eta)

  QE_1=prob_case_1*(1-pi_case_1_t)/(prob_case_1+exp(as.matrix(CASEZ_1)%*%t(beta))*pi_case_1_t);
  QE_2=prob_case_2*pi_case_2_t/(prob_case_2+exp(as.matrix(CASEZhat_2)%*%t(beta))*(1-pi_case_2_t));
  QE_3=prob_case_11*(1-pi_case_1_t)/(prob_case_11+exp(as.matrix(CASEZhat_1)%*%t(beta))*pi_case_1)*exp((as.matrix(CASEZhat_1)-as.matrix(CASEZ_1))%*%t(beta));
  QE_4=prob_case_22*pi_case_2_t/(prob_case_22+exp(as.matrix(CASEZhat_22)%*%t(beta))*(1-pi_case_2))*exp((as.matrix(CASEZhat_22)-as.matrix(CASEZhat_2))%*%t(beta));

  CQP_1=matrix(0,p,p_pi_t);
  CQP_2=matrix(0,lthsub_2,p_pi_t);
  CQP_3=matrix(0,lthsub_3,p_pi_t);
  CQP_4=matrix(0,lthsub_4,p_pi_t);
  for (i in 1:p_pi_t){
    for (j in 1:p){
      CQP_1[j,i]=sum(CASEZ_1[,j]*Z_case_pi_1_t[,i]*QE_1);
    }
    for (k in 1:lthsub_2){
      CQP_2[k,i]=-sum(CASEZhat_2[,subset_2[k]]*Z_case_pi_2_t[,i]*QE_2);
    }
    for (l in 1:lthsub_3){
      CQP_3[l,i]=sum(CASEZhat_1[,subset_3[l]]*Z_case_pi_1_t[,i]*QE_3);
    }
    for (l in 1:lthsub_4){
      CQP_4[l,i]=-sum(CASEZhat_22[,subset_4[l]]*Z_case_pi_2_t[,i]*QE_4);
    }
  }

  CQP=rbind(CQP_1,CQP_3,CQP_2,CQP_4);
  CQP=CQP%*%B_eta;
  VQP=CQP%*%t(PQ);
  VPQ=PQ%*%t(CQP);

  V=VQQ+VPP-VQP-VPQ;


  tmp1=length(Q_1)+length(Q_3);
  tmp2=length(Q_1)+length(Q_2)+length(Q_3);

  V_step1=V[1:tmp1,1:tmp1];
  J_step1=rbind(J_1,J_3);
  Q_step1=t(J_step1)%*%solve(V_step1)%*%rbind(Q_1,Q_3);

  COV_12=t(J_step1)%*%solve(V_step1)%*%V[1:tmp1,(tmp1+1):tmp2];
  V2=V[(tmp1+1):tmp2,(tmp1+1):tmp2];
  V_step2=rbind(cbind(t(J_step1)%*%solve(V_step1)%*%J_step1,COV_12),cbind(t(COV_12),V2));
  J_step2=rbind(t(J_step1)%*%solve(V_step1)%*%J_step1,J_2);
  Q_step2=t(J_step2)%*%solve(V_step2)%*%rbind(Q_step1,Q_2);
  if (length(J_step3)+length(V_step3)==0){
    COV_23=t(J_step2)%*%solve(V_step2)%*%rbind(t(J_step1)%*%solve(V_step1)%*%V[1:tmp1,(tmp2+1):length(V[1,])],V[(tmp1+1):tmp2,(tmp2+1):length(V[1,])]);
    V3=V[(tmp2+1):length(V[,1]),(tmp2+1):length(V[1,])];
    V_step3=rbind(cbind(t(J_step2)%*%solve(V_step2)%*%J_step2,COV_23),cbind(t(COV_23),V3));
    J_step3=rbind(t(J_step2)%*%solve(V_step2)%*%J_step2,J_4);
  }
  f=t(J_step3)%*%solve(V_step3)%*%rbind(Q_step2,Q_4);
  return(list(f,J_step3,V_step3));
}
