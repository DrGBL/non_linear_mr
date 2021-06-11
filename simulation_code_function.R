library(tidyverse)
library(TwoSampleMR)

umx_make_MR_data_multi <- function(nSubjects = 1000, 
                                   b_linear_qtl=c(0.02, 0.02),
                                   b_quad_qtl=c(0.02,0.02),
                                   bXY=0.1,
                                   bXY2=0.1,
                                   p_linear_qtl=c(0.5, 0.5),
                                   p_quad_qtl=c(0.5, 0.5),
                                   bUX = 0.5,
                                   bUY = 0.5, 
                                   varX=1,
                                   varY=1,
                                   seed=42) {	
  # nSubjects  = # Individuals
  # bXY  =  Linear component of causal effect of X on Y
  # bXY2 = Quadratic component of causal effect of X on Y
  # bUX  = Confounding effect of U on X
  # bUY  = Confounding effect of U on Y
  # pQTL = Decreaser allele frequency
  # b_linear_qtl = beta of snps associated linearly with exposure
  # b_quad_qtl = beta of snps associated quadratically with exposure
  # p_linear_qtl = allele frequencies of linear snps
  # p_quad_qtl = allele frequencies of quadratic snps
  # varX = variance of error term of X
  # varY = variance of error term of Y
  
  if(length(b_linear_qtl)!=length(p_linear_qtl)){
    print("Unequal linear qtl vectors.")
    return(NA)
  }
  if(length(b_quad_qtl)!=length(p_quad_qtl)){
    print("Unequal quadratic qtl vectors.")
    return(NA)
  }
  
  set.seed(seed)
  
  q_linear_qtl = 1 - p_linear_qtl
  
  q_quad_qtl = 1 - p_quad_qtl 
  
  # Simulate individual genotypic values
  linear_qtl<-as.data.frame(matrix(rbinom(length(q_linear_qtl)*nSubjects,
                                          2,
                                          q_linear_qtl),
                                   ncol=length(q_linear_qtl)))
  colnames(linear_qtl)<-paste0("linear_qtl",c(1:length(q_linear_qtl)))
  
  quad_qtl<-as.data.frame(matrix(rbinom(length(q_quad_qtl)*nSubjects,
                                        2,
                                        q_quad_qtl),
                                 ncol=length(q_quad_qtl)))
  colnames(quad_qtl)<-paste0("quad_qtl",c(1:length(q_quad_qtl)))
  
  #simulate the different relationships
  U <- rnorm(nSubjects, 0, 1) #Confounding variables
  Xa<-rowSums(b_linear_qtl*linear_qtl) #linear QTLs
  
  one_mult<-matrix(sample(x=c(-1,1), size=nSubjects*length(b_quad_qtl), replace=TRUE), ncol=length(b_quad_qtl), byrow = TRUE)
  Xb<-rowSums(sqrt(b_quad_qtl*quad_qtl)*one_mult) #quadratic QTLs
  
  X <- Xa + Xb + bUX * U + rnorm(nSubjects, 0, sqrt(varX)) # X variable
  Y <- bXY * X + bXY2 * X^2 + bUY * U + rnorm(nSubjects, 0, sqrt(varY)) # Y variable
  
  #MR_data
  MR_data = data.frame(X = X, Y = Y, U = U)
  MR_data = cbind(MR_data, linear_qtl, quad_qtl)
  
  #set up vectors for the twosamplemr calls
  #exposure
  SNP_exp<-c(rep(colnames(linear_qtl), each=2),
         rep(colnames(quad_qtl), each=2))
  exposure<-rep(c("linear", "quadratic"), length(b_linear_qtl)+length(b_quad_qtl))
  effect_allele_exp<-rep("A", 2*(length(b_linear_qtl)+length(b_quad_qtl)))
  other_allele_exp<-rep("C", 2*(length(b_linear_qtl)+length(b_quad_qtl)))
  beta_exp<-rep(NA, 2*(length(b_linear_qtl)+length(b_quad_qtl)))
  se_exp<-rep(NA, 2*(length(b_linear_qtl)+length(b_quad_qtl)))
  pval_exp<-rep(NA, 2*(length(b_linear_qtl)+length(b_quad_qtl)))
  eaf_exp<-rep(c(p_linear_qtl, p_quad_qtl), each=2)
  
  #outcome
  SNP_out<-c(colnames(linear_qtl),
             colnames(quad_qtl))
  beta_out<-rep(NA, length(b_linear_qtl)+length(b_quad_qtl))
  se_out<-rep(NA, length(b_linear_qtl)+length(b_quad_qtl))
  pval_out<-rep(NA, length(b_linear_qtl)+length(b_quad_qtl))
  effect_allele_out<-rep("A", length(b_linear_qtl)+length(b_quad_qtl))
  other_allele_out<-rep("C", length(b_linear_qtl)+length(b_quad_qtl))
  outcome<-rep("Y", length(b_linear_qtl)+length(b_quad_qtl))
  id_outcome<-rep("Y", length(b_linear_qtl)+length(b_quad_qtl))
  mr_keep_outcome<-rep("TRUE",length(b_linear_qtl)+length(b_quad_qtl))
  eaf_outcome<-c(p_linear_qtl, p_quad_qtl)
  
  #start with linear SNPs
  for(i in 1:length(b_linear_qtl)){
    #exposure
    mod_lin<-lm(X~linear_qtl[,i]-1)
    mod_lin2<-lm(I(X^2)~linear_qtl[,i]-1)
    
    beta_exp[(i-1)*2+1]<-coefficients(mod_lin)[1]
    beta_exp[(i)*2]<-coefficients(mod_lin2)[1]
    
    se_exp[(i-1)*2+1]<-coef(summary(mod_lin))[2]
    se_exp[(i)*2]<-coef(summary(mod_lin2))[2]
    
    pval_exp[(i-1)*2+1]<-coef(summary(mod_lin))[4]
    pval_exp[(i)*2]<-coef(summary(mod_lin2))[4]
    
    #outcome
    mod_out_lin<-lm(Y~linear_qtl[,i]-1)

    beta_out[i]<-coefficients(mod_out_lin)[1]
    
    se_out[i]<-coef(summary(mod_out_lin))[2]
    
    pval_out[i]<-coef(summary(mod_out_lin))[4]
    
  }
  
  #then with quad SNPs
  for(i in 1:length(b_quad_qtl)){
    #exposure
    mod_quad<-lm(X~quad_qtl[,i]-1)
    mod_quad2<-lm(I(X^2)~quad_qtl[,i]-1)
    
    beta_exp[2*length(linear_qtl)+(i-1)*2+1]<-coefficients(mod_quad)[1]
    beta_exp[2*length(linear_qtl)+(i)*2]<-coefficients(mod_quad2)[1]
    
    se_exp[2*length(linear_qtl)+(i-1)*2+1]<-coef(summary(mod_quad))[2]
    se_exp[2*length(linear_qtl)+(i)*2]<-coef(summary(mod_quad2))[2]
    
    pval_exp[2*length(linear_qtl)+(i-1)*2+1]<-coef(summary(mod_quad))[4]
    pval_exp[2*length(linear_qtl)+(i)*2]<-coef(summary(mod_quad2))[4]
    
    #outcome
    mod_out_quad<-lm(Y~quad_qtl[,i]-1)
    
    beta_out[length(linear_qtl)+i]<-coefficients(mod_out_quad)[1]
    
    se_out[length(linear_qtl)+i]<-coef(summary(mod_out_quad))[2]
    
    pval_out[length(linear_qtl)+i]<-coef(summary(mod_out_quad))[4]
  }
  
  
  #build the exposure and outcome dataframes
  exp_dat<-data.frame(SNP=SNP_exp,
                      exposure=exposure,
                      id.exposure=exposure,
                      effect_allele.exposure=effect_allele_exp, 
                      other_allele.exposure=other_allele_exp,
                      eaf.exposure=eaf_exp,
                      beta.exposure=beta_exp,
                      se.exposure=se_exp,
                      pval.exposure=pval_exp)
  
  out_dat<-data.frame(SNP=SNP_out,
                      beta.outcome=beta_out,
                      se.outcome=se_out,
                      pval.outcome=pval_out,
                      effect_allele.outcome=effect_allele_out, 
                      other_allele.outcome=other_allele_out,
                      outcome=outcome,
                      id.outcome=id_outcome,
                      mr_keep.outcome=mr_keep_outcome,
                      eaf.outcome=eaf_outcome)
  
  
  harmon_data<-mv_harmonise_data(exp_dat, out_dat)
  res<-mv_multiple(harmon_data)
  
  return(list(MR_data=MR_data, exp_dat=exp_dat, out_dat=out_dat, harmon_data=harmon_data, mvmr=res))
}

test<-umx_make_MR_data_multi(nSubjects = 100000, 
                             b_linear_qtl=c(0.5, 0.3, 0.4),
                             b_quad_qtl=c(0.2,0.4, 0.7),
                             bXY=0.8,
                             bXY2=0.2,
                             p_linear_qtl=c(0.3, 0.2, 0.5),
                             p_quad_qtl=c(0.6, 0.4, 0.4),
                             bUX = 0.3,
                             bUY = 0.3, 
                             varX=0.5,
                             varY=0.5,
                             seed=42)



