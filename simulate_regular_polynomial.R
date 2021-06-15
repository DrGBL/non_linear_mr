library(tidyverse)
library(TwoSampleMR)
library(orthopolynom)

setwd("~/Microbiologie/RichardsLab/non_linear_mr")

#see this https://pubmed.ncbi.nlm.nih.gov/25166881/

reg_MR<-function(qtl, MR_pol, MR_y, n=6){
  beta_out<-rep(NA_real_,n)
  se_out<-rep(NA_real_,n)
  pval_out<-rep(NA_real_,n)
  
  #exposure
  for(i in 1:n){
    mod<-lm(MR_pol[,i]~qtl-1)
    beta_out[i]<-coefficients(mod)[1]
    se_out[i]<-coef(summary(mod))[2]
    pval_out[i]<-coef(summary(mod))[4]
  }
  
  #outcome
  mod_y<-lm(MR_y~qtl-1)
  beta_y<-coefficients(mod_y)[1]
  se_y<-coef(summary(mod_y))[2]
  pval_y<-coef(summary(mod_y))[4]
  
  return(list(beta_exp=beta_out, 
              se_exp=se_out, 
              pval_exp=pval_out,
              beta_y=beta_y,
              se_y=se_y,
              pval_y=pval_y))
  
}

umx_make_MR_data_regular_poly <- function(nSubjects = 1000, 
                                          b_linear_qtl=c(0.02, 0.02),
                                          b_quad_qtl=c(0.02,0.02),
                                          b_cub_qtl=c(0.02,0.02),
                                          #bXY=0.1,
                                          #bXY2=0.1,
                                          p_linear_qtl=c(0.5, 0.5),
                                          p_quad_qtl=c(0.5, 0.5),
                                          p_cub_qtl=c(0.5, 0.5),
                                          bUX = 0.5,
                                          bUY = 0.5, 
                                          varX=1,
                                          varY=1,
                                          seed=42,
                                          interval_95=TRUE,
                                          remove_lin_plot=TRUE,
                                          recenter_X=TRUE,
                                          pol_order=6,
                                          FUN,
                                          ...) {	
  # nSubjects  = # Individuals
  # bXY  =  Linear component of causal effect of X on Y
  # bXY2 = Quadratic component of causal effect of X on Y
  # bUX  = Confounding effect of U on X
  # bUY  = Confounding effect of U on Y
  # pQTL = Decreaser allele frequency
  # b_linear_qtl = beta of snps associated linearly with exposure
  # b_quad_qtl = beta of snps associated quadratically with exposure
  # p_linear_qtl = minor allele frequencies of linear snps
  # p_quad_qtl = minor allele frequencies of quadratic snps
  # varX = variance of error term of X
  # varY = variance of error term of Y
  # seed = random seed number
  # interval_95 = whether to plot the final result limited to middle 95% quantiles
  # remove_lin_lot = whether to remove the linear approximation from the plot
  # pol_order = degree of maximal polynomial used
  # FUN = function between X and Y
  
  if(length(b_linear_qtl)!=length(p_linear_qtl)){
    print("Unequal linear qtl vectors.")
    return(NA)
  }
  if(length(b_quad_qtl)!=length(p_quad_qtl)){
    print("Unequal quadratic qtl vectors.")
    return(NA)
  }
  if(length(b_cub_qtl)!=length(p_cub_qtl)){
    print("Unequal cubic qtl vectors.")
    return(NA)
  }
  
  #here the simulation starts
  set.seed(seed)
  
  q_linear_qtl = 1 - p_linear_qtl
  
  q_quad_qtl = 1 - p_quad_qtl 
  
  q_cub_qtl = 1 - p_cub_qtl 
  
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
  
  cub_qtl<-as.data.frame(matrix(rbinom(length(q_cub_qtl)*nSubjects,
                                       2,
                                       q_cub_qtl),
                                ncol=length(q_cub_qtl)))
  colnames(cub_qtl)<-paste0("cub_qtl",c(1:length(q_cub_qtl)))
  
  #simulate the different relationships
  U <- rnorm(nSubjects, 0, 1) #Confounding variables
  Xa<-rowSums(b_linear_qtl*linear_qtl) #linear QTLs
  
  one_mult_quad<-matrix(sample(x=c(-1,1), size=nSubjects*length(b_quad_qtl), replace=TRUE), ncol=length(b_quad_qtl), byrow = TRUE)
  Xb<-rowSums(sqrt(b_quad_qtl*quad_qtl)*one_mult_quad) #quadratic QTLs
  
  one_mult_cub<-matrix(sample(x=c(-1,1), size=nSubjects*length(b_cub_qtl), replace=TRUE), ncol=length(b_cub_qtl), byrow = TRUE)
  Xc<-rowSums(((b_quad_qtl*cub_qtl)^(1/3))*one_mult_cub) #cubic QTLs
  
  #assume X is a sum of all the SNPs effect above
  X <- Xa + Xb + Xc + bUX * U + rnorm(nSubjects, 0, sqrt(varX)) # X variable
  
  if(recenter_X==TRUE){
    X<-X-mean(X)
  }
  
  #Y variable
  Y <- FUN(X,...) + bUY * U + rnorm(nSubjects, 0, sqrt(varY)) # Y variable
  
  #MR_data for the linear regressions below
  MR_data = data.frame(X = X, Y = Y, U = U)
  MR_data = cbind(MR_data, linear_qtl, quad_qtl, cub_qtl)
  
  #find bounds of X to determine bounds of reparameterized Legendre polynomials
  min_x<-min(MR_data$X)
  max_x<-max(MR_data$X)
  
  #add regular polynomials as columns in MR_data
  poly_matrix<-as.data.frame(matrix(data=NA, 
                                    nrow=nSubjects, 
                                    ncol=pol_order+1))
  colnames(poly_matrix)<-paste0("leg",c(0:pol_order))
  
  for(i in 1:(ncol(poly_matrix))){
    poly_matrix[,i]<-X^(i-1)
  }
  
  MR_data<-bind_cols(MR_data,poly_matrix)

  #set up vectors for the twosamplemr calls
  #exposure
  SNP_exp<-c(rep(colnames(linear_qtl), each=pol_order),
             rep(colnames(quad_qtl), each=pol_order),
             rep(colnames(cub_qtl), each=pol_order))
  exposure<-rep(paste0("leg", c(1:pol_order)), length(b_linear_qtl)+length(b_quad_qtl)+length(b_cub_qtl))
  effect_allele_exp<-rep("A", pol_order*(length(b_linear_qtl)+length(b_quad_qtl)+length(b_cub_qtl)))
  other_allele_exp<-rep("C", pol_order*(length(b_linear_qtl)+length(b_quad_qtl)+length(b_cub_qtl)))
  beta_exp<-rep(NA, pol_order*(length(b_linear_qtl)+length(b_quad_qtl)+length(b_cub_qtl)))
  se_exp<-rep(NA, pol_order*(length(b_linear_qtl)+length(b_quad_qtl)+length(b_cub_qtl)))
  pval_exp<-rep(NA, pol_order*(length(b_linear_qtl)+length(b_quad_qtl)+length(b_cub_qtl)))
  eaf_exp<-rep(c(p_linear_qtl, p_quad_qtl, p_cub_qtl), each=pol_order)
  
  #outcome
  SNP_out<-c(colnames(linear_qtl),
             colnames(quad_qtl),
             colnames(cub_qtl))
  beta_out<-rep(NA, length(b_linear_qtl)+length(b_quad_qtl)+length(b_cub_qtl))
  se_out<-rep(NA, length(b_linear_qtl)+length(b_quad_qtl)+length(b_cub_qtl))
  pval_out<-rep(NA, length(b_linear_qtl)+length(b_quad_qtl)+length(b_cub_qtl))
  effect_allele_out<-rep("A", length(b_linear_qtl)+length(b_quad_qtl)+length(b_cub_qtl))
  other_allele_out<-rep("C", length(b_linear_qtl)+length(b_quad_qtl)+length(b_cub_qtl))
  outcome<-rep("Y", length(b_linear_qtl)+length(b_quad_qtl)+length(b_cub_qtl))
  id_outcome<-rep("Y", length(b_linear_qtl)+length(b_quad_qtl)+length(b_cub_qtl))
  mr_keep_outcome<-rep("TRUE",length(b_linear_qtl)+length(b_quad_qtl)+length(b_cub_qtl))
  eaf_outcome<-c(p_linear_qtl, p_quad_qtl, p_cub_qtl)
  
  #start with linear SNPs
  for(i in 1:length(b_linear_qtl)){
    tmp_reg<-reg_MR(qtl=linear_qtl[,i],MR_y=MR_data$Y, MR_pol=MR_data[,(ncol(MR_data)-pol_order+1):ncol(MR_data)], n=pol_order)
    
    #exposure
    beta_exp[((i-1)*pol_order+1):((i-1)*pol_order+pol_order)]<-tmp_reg[["beta_exp"]]
    se_exp[((i-1)*pol_order+1):((i-1)*pol_order+pol_order)]<-tmp_reg[["se_exp"]]
    pval_exp[((i-1)*pol_order+1):((i-1)*pol_order+pol_order)]<-tmp_reg[["pval_exp"]]
    
    #outcome
    beta_out[i]<-tmp_reg[["beta_y"]]
    se_out[i]<-tmp_reg[["se_y"]]
    pval_out[i]<-tmp_reg[["pval_y"]]
    
  }
  
  #then with quad SNPs
  for(i in 1:length(b_quad_qtl)){
    tmp_reg<-reg_MR(qtl=quad_qtl[,i],MR_y=MR_data$Y, MR_pol=MR_data[,(ncol(MR_data)-pol_order+1):ncol(MR_data)], n=pol_order)
    
    #exposure
    beta_exp[(pol_order*length(linear_qtl)+(i-1)*pol_order+1):(pol_order*length(linear_qtl)+(i-1)*pol_order+pol_order)]<-tmp_reg[["beta_exp"]]
    se_exp[(pol_order*length(linear_qtl)+(i-1)*pol_order+1):(pol_order*length(linear_qtl)+(i-1)*pol_order+pol_order)]<-tmp_reg[["se_exp"]]
    pval_exp[(pol_order*length(linear_qtl)+(i-1)*pol_order+1):(pol_order*length(linear_qtl)+(i-1)*pol_order+pol_order)]<-tmp_reg[["pval_exp"]]
    
    #outcome
    beta_out[length(linear_qtl)+i]<-tmp_reg[["beta_y"]]
    se_out[length(linear_qtl)+i]<-tmp_reg[["se_y"]]
    pval_out[length(linear_qtl)+i]<-tmp_reg[["pval_y"]]
  }
  
  #then with cubic SNPs
  for(i in 1:length(b_cub_qtl)){
    tmp_reg<-reg_MR(qtl=cub_qtl[,i],MR_y=MR_data$Y, MR_pol=MR_data[,(ncol(MR_data)-pol_order+1):ncol(MR_data)], n=pol_order)
    
    #exposure
    beta_exp[(pol_order*length(linear_qtl)+pol_order*length(quad_qtl)+(i-1)*pol_order+1):(pol_order*length(linear_qtl)+pol_order*length(quad_qtl)+(i-1)*pol_order+pol_order)]<-tmp_reg[["beta_exp"]]
    se_exp[(pol_order*length(linear_qtl)+pol_order*length(quad_qtl)+(i-1)*pol_order+1):(pol_order*length(linear_qtl)+pol_order*length(quad_qtl)+(i-1)*pol_order+pol_order)]<-tmp_reg[["se_exp"]]
    pval_exp[(pol_order*length(linear_qtl)+pol_order*length(quad_qtl)+(i-1)*pol_order+1):(pol_order*length(linear_qtl)+pol_order*length(quad_qtl)+(i-1)*pol_order+pol_order)]<-tmp_reg[["pval_exp"]]
    
    #outcome
    beta_out[length(linear_qtl)+length(quad_qtl)+i]<-tmp_reg[["beta_y"]]
    se_out[length(linear_qtl)+length(quad_qtl)+i]<-tmp_reg[["se_y"]]
    pval_out[length(linear_qtl)+length(quad_qtl)+i]<-tmp_reg[["pval_y"]]
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
  
  #leg1 MR
  harmon_data_leg1<-harmonise_data(subset(exp_dat, exposure %in% c("leg1")), 
                                   out_dat)
  res_l1<-mr(harmon_data_leg1, method_list=c("mr_ivw"))
  
  #predict
  MR_data$pred_l1<-res_l1$b[1]*MR_data$leg1
  
  for(i in 2:pol_order){
    pred<-data.frame(tmp=rep(0, nSubjects))
    colnames(pred)<-paste0("pred_l",i)
    harmon_data<-mv_harmonise_data(subset(exp_dat, exposure %in% paste0("leg",c(1:i))), 
                                   out_dat)
    res_mvmr<-mv_multiple(harmon_data)
    for(j in 1:i){
      pred<-pred+res_mvmr$result$b[j]*MR_data[,paste0("leg",j)]
    }
    MR_data<-MR_data %>% bind_cols(.,pred)
  }

  quad_fun(bXY=bXY, bXY2=bXY2, X=X)
  MR_plot<-data.frame(X=rep(MR_data$X, pol_order+1),
                      Y=c(FUN(MR_data$X,...),
                          as.vector(unlist(MR_data[,(ncol(MR_data)-pol_order+1):ncol(MR_data)]))),
                      Legendre=rep(c("True",paste0("L", c(1:pol_order))), each=nSubjects))
  
  if(interval_95==TRUE){
    int_95<-quantile(MR_data$X,c(0.025, 0.975))
    MR_plot<-MR_plot %>% filter(X >int_95[1] & X<int_95[2])
  }
  
  if(remove_lin_plot==TRUE){
    MR_plot<-MR_plot %>% filter(Legendre!="L1")
  } 
  
  #View(MR_plot)
  
  plot<-MR_plot %>% ggplot(aes(x=X, y=Y))+geom_line(aes(colour=Legendre, group=Legendre))+xlim(int_95)
  
  return(plot)
  
  #return(list(plot=plot,MR_l1=res_l1, MR_l2=res_l2, MR_l3=res_l3, MR_l4=res_l4, MR_l5=res_l5, MR_l6=res_l6))
}

quad_fun<-function(bXY,bXY2,X){
  return(bXY*X+bXY2*X^2)
}


test<-umx_make_MR_data_regular_poly(nSubjects = 100000, 
                                    b_linear_qtl=rep(c(0.5, 0.3, 0.4),2),
                                    b_quad_qtl=rep(c(0.2,0.4, 0.7),2),
                                    p_linear_qtl=rep(c(0.3, 0.2, 0.5),2),
                                    p_quad_qtl=rep(c(0.6, 0.4, 0.4),2),
                                    bUX = 0.3,
                                    bUY = 0.3, 
                                    varX=0.5,
                                    varY=0.5,
                                    seed=42,
                                    FUN=quad_fun,
                                    bXY=0.8,
                                    bXY2=0.2,
                                    pol_order = 8)

lin_fun<-function(bXY,X){
  return(bXY*X)
}

test<-umx_make_MR_data_regular_poly(nSubjects = 10000, 
                                    b_linear_qtl=c(0.5, 0.3, 0.4),
                                    b_quad_qtl=c(0.2,0.4, 0.7),
                                    p_linear_qtl=c(0.3, 0.2, 0.5),
                                    p_quad_qtl=c(0.6, 0.4, 0.4),
                                    bUX = 0.3,
                                    bUY = 0.3, 
                                    varX=0.5,
                                    varY=0.5,
                                    seed=42,
                                    FUN=lin_fun,
                                    bXY=0.8)

cub_fun<-function(bXY, bXY2, bXY3,X){
  return(bXY*X+bXY2*X^2+bXY3*X^3)
}

test<-umx_make_MR_data_regular_poly(nSubjects = 100000, 
                                    b_linear_qtl=rep(c(0.5, 0.3, 0.4),3),
                                    b_quad_qtl=rep(c(0.2,0.4, 0.7),3),
                                    b_cub_qtl=rep(c(0.2,0.4, 0.7),3),
                                    p_linear_qtl=rep(c(0.3, 0.2, 0.5),3),
                                    p_quad_qtl=rep(c(0.6, 0.4, 0.4),3),
                                    p_cub_qtl=rep(c(0.6, 0.4, 0.4),3),
                                    bUX = 0.3,
                                    bUY = 0.3, 
                                    varX=0.5,
                                    varY=0.5,
                                    seed=42,
                                    FUN=cub_fun,
                                    bXY=0.8,
                                    bXY2=0.4,
                                    bXY3=0.2)

#this is where shit starts to break down for real, and a cubic term might be needed
sig_fun<-function(bXY,X){
  return(1/(1+exp(-bXY*(X-10))))
}

test<-umx_make_MR_data_regular_poly(nSubjects = 400000, 
                                    b_linear_qtl=rep(c(0.5, 0.3, 0.4,0.5, 0.3, 0.4),4),
                                    b_quad_qtl=rep(c(0.2,0.4, 0.7,0.2,0.4, 0.7),4),
                                    b_cub_qtl=c(0.2,0.4, 0.7),
                                    p_linear_qtl=rep(c(0.3, 0.2, 0.5,0.3, 0.2, 0.5),4),
                                    p_quad_qtl=rep(c(0.6, 0.4, 0.4,0.6, 0.4, 0.4),4),
                                    p_cub_qtl=c(0.6, 0.4, 0.4),
                                    bUX = 0.3,
                                    bUY = 0.3, 
                                    varX=5,
                                    varY=2,
                                    seed=42,
                                    FUN=sig_fun,
                                    bXY=1,
                                    pol_order = 9,
                                    remove_lin_plot = FALSE,
                                    recenter_X = FALSE)
#burgess' functions
j_u_shape<-function(bXY=0.2,X){
  return(bXY*(X-1)^2)
}

threshold_func<-function(bXY=1,Shift=2,X){
  return(bXY*ifelse(X-Shift<0,0,X-Shift))
}

test<-umx_make_MR_data_regular_poly(nSubjects = 50000, 
                                    b_linear_qtl=rep(c(0.5, 0.3, 0.4,0.5, 0.3, 0.4),1),
                                    b_quad_qtl=rep(c(0.2,0.4, 0.7,0.2,0.4, 0.7),1),
                                    b_cub_qtl=c(0.2,0.4, 0.7),
                                    p_linear_qtl=rep(c(0.3, 0.2, 0.5,0.3, 0.2, 0.5),1),
                                    p_quad_qtl=rep(c(0.6, 0.4, 0.4,0.6, 0.4, 0.4),1),
                                    p_cub_qtl=c(0.6, 0.4, 0.4),
                                    bUX = 0.3,
                                    bUY = 0.3, 
                                    varX=5,
                                    varY=2,
                                    seed=42,
                                    FUN=j_u_shape,
                                    bXY=0.2,
                                    #Shift=2,
                                    recenter_X = FALSE,
                                    pol_order = 9,
                                    remove_lin_plot=FALSE)

