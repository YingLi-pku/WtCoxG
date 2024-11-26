####calculate batch effect p-value----------------------------------------------

cal.batcheffect.pval = function(n0,         # number of controls
                                n1,         # number of cases
                                n2,         # number of external dataset
                                p0,         # estimated MAF in controls
                                p1,         # estimated MAF in cases
                                p2,         # estimated MAF in external dataset
                                sigma2 = 0,
                                event.rate,
                                er = 0.5,
                                var.ratio=1)
  #GRM = diag(1, n0+n1),
  #vec_weight) # event rate in population
{
  er = n1/(n1+n0)
  w0 = (1-event.rate)/event.rate/((1-er)/er)
  w1 = 1
  
  p01 = sum(p0*w0*n0 + p1*w1*n1) / sum(w0*n0 + w1*n1)
  p012 = sum(p0*w0*n0 + p1*w1*n1 + p2*n2*w0) / sum(n1*w1+n0*w0+ n2*w0)
  
  # if(identical(GRM , diag(1, n0+n1))){
  v =(n1*w1^2+n0*w0^2)/(2*(n1*w1+n0*w0)^2) * p012 * (1-p012) + 1/(2*n2) * p012 * (1-p012) + sigma2
  # }else{
  #   vec_weight = vec_weight/sum(vec_weight)/2
  #   v = t(vec_weight) %*% GRM %*% vec_weight * 2*p012 * (1-p012) + 1/(2*n2) * p012 * (1-p012) + sigma2
  # }
  z = (p01 - p2) / sqrt(v) /var.ratio
  p = 2*pnorm(-abs(z), lower.tail=TRUE)
  
  chisq = z^2
  return(p)
}


#####estimate TPR and sigma2
est_param = function(vec_p_bat, 
                     vec_var_Sbat,
                     vec_cutoff=seq(0.01,0.4,0.1)
){
  
  ########  step1: the propportion of p_bat>cutoff
  vec_p_deno=lapply(vec_cutoff,function(p_cut){
    p_deno=mean(na.omit(vec_p_bat>p_cut))
  })%>%unlist()
  
  ######## optimization function
  opti_fun = function( var_Sbat,vec_p_deno,par ){
    diff=lapply(1: length(vec_cutoff),function(j){
      p_cut=vec_cutoff[j]
      lb = -qnorm(1-p_cut/2) * sqrt(var_Sbat)
      ub = qnorm(1-p_cut/2) *sqrt(var_Sbat)
      p_deno=vec_p_deno[j]
      
      c = pnorm(ub, 0, sqrt(var_Sbat+par[2]), log.p = T)
      d = pnorm(lb, 0, sqrt(var_Sbat+par[2]), log.p = T)
      
      pro.cut = par[1]*( exp(d) * (exp(c-d) - 1) )+(1-par[1])*(1-p_cut)
      t = ((p_deno - pro.cut)/p_deno)^2
    })%>%do.call("sum",.)
    
    return(diff)
    
  }
  
  #######estimate TPR and sigma2 for each SNP
  var.diff = lapply(1:length(vec_var_Sbat), function(i){
    
    if(i%%100==0)cat(i,"\n")
    
    obj = optim(par = c(0.01, 0.01)
                #,method = "SANN"
                #, lower = 0, upper = 1
                , fn = opti_fun,
                vec_p_deno=vec_p_deno,
                var_Sbat=vec_var_Sbat[i])
    TPR = min(1,max(0,obj$par[1]))
    sigma2 = min(1,max(0,obj$par[2]))
    return(cbind(TPR, sigma2))
    
  })%>%
    do.call("rbind",.)%>%as_tibble()
  
  return(var.diff)
} 

###### optimal weight for mu_ext ------------------------------------------
S.root = function(par, prev, R, y, mu1, w , mu, N, n.ext, sigma2, TPR){
  b=par[1]
  
  p.fun= function(b,prev, R, y, mu1, mu, w, N, n.ext, sigma2, TPR){
    
    meanR= mean(R)
    sumR = sum(R)
    
    mu0 = mu
    mu.pop = mu1*prev+mu0*(1-prev)
    
    mu.i = ifelse(y==1, 2*mu1, 2*mu0)
    
    S = sum((R-(1-b)*meanR)*mu.i)-sumR*2*b*mu.pop
    
    w1 = w/(2*sum(w))
    mu = mean(mu.i)/2
    
    var_mu_ext = mu*(1-mu)/(2*n.ext)
    var_Sbat = sum(w1^2)*2*mu*(1-mu) + var_mu_ext
    
    p_cut=0.1
    lb = -qnorm(1-p_cut/2) * sqrt(var_Sbat)
    ub = qnorm(1-p_cut/2) *sqrt(var_Sbat)
    c = pnorm(ub, 0, sqrt(var_Sbat+sigma2),log.p = T)
    d = pnorm(lb, 0, sqrt(var_Sbat+sigma2),log.p = T)
    p_deno = TPR*( exp(d) * (exp(c-d) - 1) )+(1-TPR)*(1-p_cut)
    
    ##sigma2=0
    var.int = sum((R-(1-b)*meanR)^2)*2*mu*(1-mu)
    var_S = var.int +  4*b^2 *sumR^2* var_mu_ext
    cov_Sbat_S = sum(w1*(R-(1-b)*meanR))*2*mu*(1-mu)+2*b*sumR*var_mu_ext
    VAR = matrix(c(var_S, cov_Sbat_S, cov_Sbat_S, var_Sbat),nrow=2)
    p0 = max(0,pmvnorm(lower=c( -Inf, lb),upper=c( -abs(S), ub), mean=c(0,0), sigma=VAR))
    
    ##sigma2!=0
    var_S1 = var.int +  4*b^2 *sumR^2* (var_mu_ext+sigma2)
    cov_Sbat_S1 = sum(w1*(R-(1-b)*meanR))*2*mu*(1-mu)+2*b*sumR*(var_mu_ext+sigma2)
    var_Sbat1 = var_Sbat+sigma2
    VAR1 = matrix(c(var_S1, cov_Sbat_S1, cov_Sbat_S1, var_Sbat1),nrow=2)
    p1=max(0,pmvnorm(lower=c( -Inf, lb),upper=c( -abs(S), ub), mean=c(0,0), sigma=VAR1))
    
    p.con = 2*(TPR*p1+(1-TPR)*p0)/p_deno
    #diff = -log10(p.con)+log10(5e-8)
    diff = -log10(p.con/5e-8)
    
    return(diff)
    
  }
  
  mu1=uniroot(p.fun, lower=mu,upper=1
              
              , b=b, prev=prev,mu=mu,
              R=R, y=y, w=w,N=N, n.ext=n.ext,sigma2=sigma2,TPR=TPR)$root
  
  return(mu1)
  
}

###############################################################################
###SPA method----------------------------------------------------------------------
library(mvtnorm)
M_G0 = function(t, MAF){
  re = (1 - MAF + MAF * exp(t))^2
  return(re)
}
# The first derivative of the MGF of G (genotype)
M_G1 = function(t, MAF){
  re = 2*(MAF * exp(t))*(1 - MAF + MAF * exp(t))
  return(re)                           
}

# The second derivative of the MGF of G (genotype)
M_G2 = function(t, MAF){
  re = 2*(MAF * exp(t))^2 + 2*(MAF * exp(t))*(1 - MAF + MAF * exp(t))
  return(re)
}

# The CGF of G (genotype)
K_G0 = function(t, MAF){
  re = log(M_G0(t, MAF))
  return(re)
}

K_G1 = function(t, MAF){
  re = M_G1(t, MAF)/M_G0(t, MAF)
  return(re)
}

# The second derivative of the CGF of G (genotype)
K_G2 = function(t, MAF){
  re = M_G0(t, MAF)/M_G0(t, MAF) * M_G2(t, MAF)/M_G0(t, MAF) -(M_G1(t, MAF)/M_G0(t, MAF))^2
  return(re)
}

# The CGF of score test statistic 
H_org = function(t, R, MAF,N,n.ext,N.all,sumR,var_mu_ext, g.var.est,meanR,b){
  
  n.t = length(t)
  out = rep(0,n.t)
  R_hat  = sumR/N.all
  
  mu.adj = -2*b*sumR*MAF
  var.adj = 4*b^2 *sumR^2* var_mu_ext
  for(i in 1:n.t){
    t1 = t[i]
    out[i] = sum(K_G0(t1*(R-(1-b)*meanR), MAF)) + mu.adj * t1 + var.adj/2 * t1^2  
    
  }
  
  
  
  
  return(out)
}

# The first derivative of the CGF of score test statistic 
H1_adj = function(t, R, s, MAF, N, n.ext, N.all, sumR,var_mu_ext, g.var.est,meanR,b)
{
  n.t = length(t)
  out = rep(0,n.t)
  R_hat  = sumR/N.all
  
  mu.adj = -2*b*sumR*MAF
  var.adj = 4*b^2 *sumR^2* var_mu_ext
  
  for(i in 1:n.t){
    t1 = t[i]
    out[i] = sum((R-(1-b)*meanR) *K_G1(t1*((R-(1-b)*meanR)), MAF)) + mu.adj +var.adj*t1 - s
  }
  return(out)
}

# The second derivative of the CGF of score test statistic 
H2 = function(t, R, MAF, N, n.ext, N.all, sumR,var_mu_ext, g.var.est,meanR,b)
{
  n.t = length(t)
  out = rep(0,n.t)
  R_hat  = sumR/N.all
  
  mu.adj = -n.ext*R_hat*2*MAF
  var.adj = n.ext*R_hat^2*2*MAF*(1-MAF)
  
  
  for(i in 1:n.t){
    t1 = t[i]
    out[i] = sum((R-(1-b)*meanR)^2*K_G2(t1 * (R-(1-b)*meanR) , MAF)) + var.adj
  }
  return(out)
}

GetProb_SPA_G = function(MAF, R, s, N, n.ext, N.all, sumR,var_mu_ext, g.var.est,meanR,b, lower.tail){
  
  out = uniroot(H1_adj, c(-1,1), extendInt = "yes",
                R=R, s=s, MAF=MAF, N = N, n.ext=n.ext, N.all = N.all, sumR = sumR,
                var_mu_ext = var_mu_ext, g.var.est=g.var.est, meanR =meanR, b=b)
  zeta = out$root
  
  k1 = H_org(zeta, R=R, MAF=MAF, N = N, n.ext=n.ext, N.all = N.all, sumR = sumR,
             var_mu_ext = var_mu_ext, g.var.est=g.var.est, meanR =meanR,b=b)
  k2 = H2(zeta, R=R, MAF=MAF, N = N, n.ext=n.ext, N.all = N.all, sumR = sumR,
          var_mu_ext = var_mu_ext, g.var.est=g.var.est, meanR =meanR, b=b)
  
  temp1 = zeta * s - k1
  
  w = sign(zeta) * (2 *temp1)^{1/2}
  v = zeta * (k2)^{1/2}
  
  pval = pnorm(w + 1/w * log(v/w), lower.tail = lower.tail)
  # pval.norm = pnorm(q2, lower.tail = lower.tail)
  re = pval
  return(re)
}

SPA_G.one.SNP_homo = function(g,
                              R,
                              mu.ext=NA,
                              n.ext=NA,
                              b=0,
                              sigma2=NA,
                              var.ratio=1,
                              Cutoff = 2,
                              impute.method = "fixed",
                              missing.cutoff = 0.15,
                              min.mac = 10,          # update on 2022-08-16 : replace 0.0001 by 0.000001
                              G.model = "Add")
{
  ## calculate MAF and update genotype vector
  
  ##imputation missing SNP
  missing.rate = mean(is.na(g))
  pos.na = which(is.na(g))
  if(missing.rate != 0){
    g[pos.na] = mean(na.omit(g))
  }
  
  
  if(is.na(mu.ext)){
    mu.ext =0
    n.ext=0
    
  }
  
  
  if(sum(g)<min.mac|sum(2-g)<min.mac|missing.rate>missing.cutoff){
    
    MAF= mean(na.omit(g))/2
    
    return(cbind(MAF, missing.rate, NA, NA, NA,NA,NA))      
  }
  
  
  ######################################################################
  
  ## Score statistic
  N=length(g)
  mu.int = mean(g)/2
  MAF = (1-b) * mu.int + b * mu.ext 
  sumR = sum(R)
  N.all = N + n.ext
  S = sum(R *(g-2*MAF))
  S = S / var.ratio
  
  ## estimated variance without adjusting for covariates
  g.var.est = 2 * MAF * (1 - MAF)
  var_mu_ext = ifelse(n.ext==0,0, MAF*(1-MAF)/(2*n.ext)+sigma2) 
  
  
  #  S.var = sum(R^2 * g.var.est + sum(R)^2 * g.var.est/n.ext)
  meanR  = mean(R)
  S.var = sum((R-(1-b)*meanR)^2)*g.var.est +  4*b^2 *sumR^2* var_mu_ext
  
  z = S/sqrt(S.var)
  
  if(abs(z) < Cutoff){
    pval.norm = pnorm(abs(z), lower.tail = FALSE)*2
    return(c(MAF, missing.rate, pval.norm, pval.norm, S, S.var, z))  # update on 2022-10-05 : MAF.est.negative.num 
  }else{
    pval1 = GetProb_SPA_G(MAF, R = R, abs(S), N=N, n.ext=n.ext, N.all=N.all,
                          var_mu_ext = var_mu_ext, g.var.est=g.var.est,meanR =meanR,
                          sumR=sumR, b=b, lower.tail = FALSE) # EmpSPA-G p value 
    pval2 = GetProb_SPA_G(MAF, R = R, -abs(S), N=N, n.ext=n.ext, N.all=N.all,
                          var_mu_ext = var_mu_ext, g.var.est=g.var.est,meanR =meanR,
                          sumR=sumR, b=b, lower.tail = TRUE) # EmpSPA-G p value 
    
    pval3 = pnorm(abs(z), lower.tail = FALSE) # Normal
    pval4 = pnorm(-abs(z), lower.tail = TRUE) # Normal
    
    pval.spa.G = pval1 + pval2
    pval.norm = pval3 + pval4
    
    # if(abs(z) < Cutoff){
    #   pval.spa.G = pval.norm
    # }
    
    pval = c(pval.spa.G, pval.norm) # 2 elements: element 1 is from EmpSPA-G, element 2 is from Normal
    return(c(MAF, missing.rate, pval, S, S.var, z))
    
    
  }
  
}




#####calculate p-value --------------------------------------------------

cal.cutoff.p = function(g,
                        R,
                        w,
                        TPR,
                        p_bat,
                        var.ratio = 1,
                        var.ratio.w0 = 1,
                        var.ratio.w1 = 1,
                        var.ratio0 = 1,
                        var.ratio1 = 1,
                        mu.ext=NA,
                        n.ext=NA,
                        b=0,
                        sigma2=NA,
                        p_cut=0.1
){
  ##imputation missing SNP
  missing.rate = mean(is.na(g))
  pos.na = which(is.na(g))
  if(missing.rate != 0){
    g[pos.na] = mean(na.omit(g))
  }
  
  if(is.na(mu.ext)){
    n.ext=0
    mu.ext=0
    b=0
    cov_Sbat_S=cov_Sbat_S1=NA
    meanR= mean(R)
    sumR = sum(R)
    
    mu.int = mean(g)/2
    S= sum(R*(g-2*mu.int))
    mu=mu.int
    p.con = SPA_G.one.SNP_homo(g=g, R=R, mu.ext=NA, n.ext= 0,sigma2=0,var.ratio = var.ratio)[3]
    p_deno =NA
    return(cbind(p.con, sigma2, S, p_deno, cov_Sbat_S, cov_Sbat_S1) )
  }
  
  
  if(p_bat<p_cut|is.na(p_bat)|sum(g)<10|sum(2-g)<10 ){
    
    p.con=b=S=p_deno=cov_Sbat_S=cov_Sbat_S1=NA
    return(cbind(p.con, b, S, p_deno, cov_Sbat_S, cov_Sbat_S1) )
  }else{
    
    
    meanR= mean(R)
    sumR = sum(R)
    mu.int = mean(g)/2
    N = length(g)
    
    
    mu.com = (1-b)*mu.int + b*mu.ext
    S= sum(R*(g-2*mu.com))
    # S=sum((R-(1-b)*meanR)*g)-sumR*2*b*mu.ext
    
    mu=mu.com
    w1 = w/(2*sum(w))
    
    var_mu_ext = mu*(1-mu)/(2*n.ext)
    var_Sbat = sum(w1^2)*2*mu*(1-mu) + var_mu_ext
    
    
    lb = -qnorm(1-p_cut/2) * sqrt(var_Sbat)*var.ratio.w0
    ub = qnorm(1-p_cut/2) *sqrt(var_Sbat)*var.ratio.w0
    c = pnorm(ub/var.ratio.w1, 0, sqrt(var_Sbat+sigma2), log.p = T)
    d = pnorm(lb/var.ratio.w1, 0, sqrt(var_Sbat+sigma2), log.p = T)
    p_deno = TPR*(  exp(d) * (exp(c-d) - 1) )+(1-TPR)*(1-p_cut)
    
    ##sigma2=0
    p_spa_s0 = SPA_G.one.SNP_homo(g=g, R=R, b=b ,mu.ext=mu.ext, n.ext= n.ext, sigma2=0,
                                  var.ratio = var.ratio0)[3]
    var_S = (S/var.ratio0)^2/qchisq(p_spa_s0, 1, ncp = 0, lower.tail = F)
    
    
    var.int = sum((R-(1-b)*meanR)^2)*2*mu*(1-mu)
    # var_S = var.int +  4*b^2 *sumR^2* var_mu_ext
    cov_Sbat_S = sum(w1*(R-(1-b)*meanR))*2*mu*(1-mu)+2*b*sumR*var_mu_ext
    cov_Sbat_S = cov_Sbat_S*sqrt(var_S/(var.int +  4*b^2 *sumR^2* var_mu_ext))
    VAR = matrix(c(var_S, cov_Sbat_S, cov_Sbat_S, var_Sbat),nrow=2)
    p0 =max(0,min(1, pmvnorm(lower=c( -Inf, lb/var.ratio.w0),upper=c( -abs(S/var.ratio0), ub/var.ratio.w0), mean=c(0,0), sigma=VAR)))
    
    ##sigma2!=0
    p_spa_s1 = SPA_G.one.SNP_homo(g=g, R=R, b=b, mu.ext=mu.ext, n.ext= n.ext, sigma2=sigma2,
                                  var.ratio = var.ratio1)[3]
    var_S1 = (S/var.ratio1)^2/qchisq(p_spa_s1, 1, ncp = 0, lower.tail = F)
    
    #var_S1 = var.int +  4*b^2 *sumR^2* (var_mu_ext+sigma2)
    cov_Sbat_S1 = sum(w1*(R-(1-b)*meanR))*2*mu*(1-mu)+2*b*sumR*(var_mu_ext+sigma2)
    cov_Sbat_S1 = cov_Sbat_S1*sqrt(var_S1/(var.int +  4*b^2 *sumR^2* (var_mu_ext+sigma2)))
    var_Sbat1 = var_Sbat+sigma2
    VAR1 = matrix(c(var_S1, cov_Sbat_S1, cov_Sbat_S1, var_Sbat1),nrow=2)
    p1= max(0,min(1,pmvnorm(lower=c( -Inf, lb/var.ratio.w1),upper=c( -abs(S/var.ratio1), ub/var.ratio.w1), mean=c(0,0), sigma=VAR1)))
    
    p.con = 2*(TPR*p1+(1-TPR)*p0)/p_deno
    
    
    
    return(cbind(p.con, b,  p_deno, cov_Sbat_S, cov_Sbat_S1) )
    
    
  }
  
  
}
