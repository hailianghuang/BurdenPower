
NCPgen <- function(q,R,K,CaseN,ConN) {
  # q is the AF
  # R is the RR
  
  p <- 1-q;
  
  p_d_aa <- K / (p*p + 2*p*q*R + q*q*R*R);
  p_d_ab <- p_d_aa * R;
  p_d_bb <- p_d_aa * R * R;
  
  p_aa_d <- p_d_aa * p * p/ K ;
  p_ab_d <- p_d_ab * p * q * 2 / K ;
  p_bb_d <- p_d_bb * q * q/ K ;
  
  p_aa_u <- (p*p - K * p_aa_d)  / (1-K) ;
  p_ab_u <- (2*p*q - K * p_ab_d) / (1-K) ; 
  p_bb_u <- (q*q - K * p_bb_d)  / (1-K)	;
  
  p_case <- p_aa_d + 0.5*p_ab_d;
  p_con <- p_aa_u + 0.5*p_ab_u;
  q_case <- 1 - p_case ;
  q_con <- 1 - p_con ;
  
#  NCP <- chisq.test(matrix(c(2*p_case*CaseN,2*p_con*ConN,2*q_case*CaseN,2*q_con*ConN),2,2),correct=F);
  
  summary(as.table(matrix(c(2*p_case*CaseN,2*p_con*ConN,2*q_case*CaseN,2*q_con*ConN),2,2)))$statistic;
  
}

denovoGen <- function(R, q,f,CaseN,ConN) {

    # q is number of de novo mutations 
    # R is the relative risk of having a variant
    # CaseN is the number of cases - N*r/(1+r)
    # CaseN is the number of controls - N/(1+r)
    MutCases = rpois(CaseN, q * 1* (1-f) )  + rpois(CaseN, q * R * f )
    MutCon = rpois(ConN, q * 1)
    
    m <- glm(c(MutCases, MutCon)~ c(rep(1, CaseN), rep(0, ConN)), family="poisson")
    -(m$deviance-m$null.deviance)

}

denovoGenPower <- function(R, q,f,CaseN,ConN, p_cut_denovo) {
  
   pchisq(qchisq(p_cut_denovo, 1, low=F), 1, ncp=denovoGen(R, q,f,CaseN,ConN), low=F)
  
}

getSingleVarPower <- function(f, AF, R, K, N, r, P_cut_single){
  chi <- unlist(lapply(AF[1:f], function(i){ NCPgen(i,R,K,N*r/(1+r),N/(1+r)) })) 
  power_snp <- pchisq(qchisq(P_cut_single, 1, low=F), 1, ncp=chi, low=F)
  ret <- 1-prod(1-power_snp) 
  ret
}

getBurdenPower <- function(f, AF, sum_var, R, K, N, r, p_cut_burden){
  
  ratio <- sum(AF[1:f] * (1-AF[1:f]) ) / sum_var
  chi <- NCPgen(sum(AF) , exp(log(R)*ratio), K, N*r/(1+r),N/(1+r))
  pchisq(qchisq(p_cut_burden, 1, low=F), 1, ncp=chi , low=F) 
}

getDenovoPower_parametric <- function(rr, q,f, N, r,  p_cut_denovo){

  MutCases <- N*r/(1+r)
  MutCon <- N/(1+r)
  x1 <- q * 1* (1-f) +  q * rr * f
  x0 <- q * 1
  chi <- (x1-x0)^2/(x1/MutCases+x0/MutCon)
  power <- pchisq(qchisq(p_cut_denovo, 1, low=F), 1, ncp=chi, low=F)
  power
}
