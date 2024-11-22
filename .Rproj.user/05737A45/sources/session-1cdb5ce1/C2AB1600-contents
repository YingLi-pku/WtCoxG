# cd /lustre/home/2211110003/code
# sbatch -p C064M0256G --qos=low --ntasks-per-node=1 -J set3 -t 1-0:0 --array=1-1080 -o log/%A_%a.log --wrap='Rscript SNPinfo-unif.R $SLURM_ARRAY_TASK_ID'
args=commandArgs(TRUE)
print(args)
print(sessionInfo())
n.cpu=as.numeric(args)
##########################################

library(survival)
library(ggplot2)
library(patchwork)
library(reshape2)
library(data.table)
library(tidyr)
library(dplyr)
library(msm)
library(LTFHPlus)
library(truncnorm)
library(KScorrect)

source("./function-1021.R")

##################################################################################################
params = expand.grid(prev = c(0.01,0.05,0.1)                     # disease prevalence
                     ,h2 = c(0.01, seq(0.012,0.0155,0.0005))
                     ,n1 = c(2000)                              # ncase=ncontrol
                     ,bc = c(0, 0.03, 0.12, 0.5)                 # proportion of batch effect
                     ,theta = c( 0.5 )                     # level of batch effect
                     ,group=1:10
)

h2 = params[n.cpu, "h2"]
prev = params[n.cpu, "prev"]
bc= params[n.cpu, "bc"]
n1= params[n.cpu, "n1"]
theta=params[n.cpu, "theta"]
group = params[n.cpu, "group"]



for(v in ((group-1)*10 + 1):(group*10)){
  cat(v,"\n")

  outputfile = paste0("./SNPinfo_unif/prev",prev,
                      "_her",h2,
                      "_bc",bc,
                      "_theta",theta,
                      "_n",n1,
                      "_v",v,
                      ".RData")

  if(!file.exists(outputfile)){

    ##load pheno
    cat("load phenotype -----------------------------------\n")
    load(paste0("./Pheno/prev",prev,
                "_her",h2,
                "_n",n1,
                "_v",v,
                ".RData"))

    ##update SNP info--------------------------------------------------------------------------
    cat("load causal SNP  -----------------------------------\n")
    load(paste0("./Geno/causal_prev",prev,
                "_n",n1,
                ".RData"))



    info.SNP.causal = info.SNP.causal%>%
      mutate(mu0 = apply(G.ds.causal[Phen.ds$event==0,], 2, function(x){mean(x)/2}),  # MAF in control
             mu1 = apply(G.ds.causal[Phen.ds$event==1,], 2, function(x){mean(x)/2}),  # MAF in case
             mu.int = (mu0 + mu1)/2,                                                  # MAF in internal samples
             mu.w = prev*mu1 + (1-prev)*mu0)%>%                                       # weighted MAF in internal samples
      mutate(batcheffect = rbinom(n(),1, bc),                                         # assign batch effect=0/1 for each SNPs
             trueMAF_ref = ifelse(batcheffect==0, trueMAF,                            # if batch effect=0 then mu.ext = mu
                                  runif(n(), 0.1*trueMAF, 4*trueMAF )))%>% # if batch effect=1 then mu.ext ~N(mu,sigma2)
      mutate(causal = "yes", trueMAF_ref = pmax(0, pmin(1, trueMAF_ref)))             # mu.ext is truncated within [0,1]


    cat("load non causal SNP  -----------------------------------\n")

    load(paste0("./Geno/null_n",n1,
                ".RData"))

    # mu0 = apply(G.ds.null[Phen.ds$event==0,], 2, function(x){mean(x)/2})
    # mu1 = apply(G.ds.null[Phen.ds$event==1,], 2, function(x){mean(x)/2})




    cat("calculate MAF  -----------------------------------\n")
    mu_data = lapply(1:ncol(G.ds.null), function(i){

      if(i%%100000==0)cat(i,"\n")

      mu0 = mean(G.ds.null[Phen.ds$event==0,i])/2
      mu1 = mean(G.ds.null[Phen.ds$event==1,i])/2


      return(cbind(mu0, mu1))
    })%>%
      do.call("rbind",.)%>%as_tibble()

    rm(G.ds.null);gc()
    info.SNP = info.SNP %>% mutate(mu0 = mu_data$mu0,
                                   mu1 = mu_data$mu1,
                                   mu.int = (mu0 + mu1)/2,
                                   mu.w = prev * mu1 + (1-prev) * mu0)%>%
      mutate(batcheffect = rbinom(n(),1, bc),
             trueMAF_ref = ifelse(batcheffect==0, trueMAF,
                                  runif(n(), 0.1*trueMAF, 4*trueMAF)))%>%
      mutate(causal = "no", trueMAF_ref = pmax(0, pmin(1, trueMAF_ref)))

    rm(mu_data);gc()
    ####combine all SNPs-------------------------------------------------------------------
    cat("combine all SNPs-------------------------------\n")
    info.SNP.all = rbind(info.SNP, info.SNP.causal)%>%
      mutate( sigma2.true = theta*trueMAF*(1-trueMAF))

    rm(info.SNP, info.SNP.causal);gc()

    n.ref = c("4e3", "1e4", "5e4", "2e5") #size of reference population
    refMAFs = lapply(n.ref, function(x){
      x = as.numeric(x)
      mu.ref = rbinom(nrow(info.SNP.all), 2*x, info.SNP.all$trueMAF_ref)/(2*x) # estimated MAF from reference population
    }) %>% do.call("cbind",.) %>%
      as_tibble() %>% setNames(paste0(rep("maf",4),n.ref))

    info.SNP.all=cbind(info.SNP.all, refMAFs)


    ####calculate batch effect p-value-------------------------------------------------------
    cat("calculate batch effect p-value-------------------------------\n")
    p_bat = lapply(1:nrow(info.SNP.all), function(i){

      if(i%%10000==0)cat(i,"\n")
      temp = lapply(n.ref , function(x){

        n.ext = as.numeric(x)
        p = cal.batcheffect.pval(n0 = n1,
                                 n1 = n1,
                                 n2 = n.ext,
                                 p0 = info.SNP.all$mu0[i],
                                 p1 = info.SNP.all$mu1[i],
                                 p2 = info.SNP.all[i,paste0("maf",x)] ,
                                 event.rate=prev)
        return(p)

      })%>%
        do.call("cbind",.)

      return(temp)
    })%>%do.call("rbind",.)%>%as_tibble()%>%
      setNames(paste0(rep("p",4),  n.ref))
    info.SNP.all = cbind(info.SNP.all, p_bat)


    #####estimate TPR and sigma2---------------------------------------------------------------
    cat("estimate TPR and sigma2--------------\n")
    t1 = Sys.time()
    maf.group = c(seq(0, 0.4, 0.01),max(info.SNP.all$mu.int))
    info.SNP.all =lapply(1:(length(maf.group)-1), function(i){
      cat(i,"\n")
      data.ref = info.SNP.all %>% filter(mu.int >=max(maf.group[i]-0.1,0) & mu.int <min(1,maf.group[i+1]+0.1) )
      #data.ref=info.SNP.all
      data = info.SNP.all %>% filter(mu.int >maf.group[i] & mu.int <=maf.group[i+1] )

      mu = (maf.group[i]+maf.group[i+1])/2
      w1 = Phen.ds$weight/(2*sum(Phen.ds$weight))


      temp = lapply(n.ref , function(x){

        n.ext = as.numeric(x)
        var_mu_ext = mu*(1-mu)/(2*n.ext)
        var_Sbat = sum(w1^2)*2*mu*(1-mu) + var_mu_ext


        obj=est_param(vec_p_bat=data.ref[,paste0("p",x)],
                      vec_var_Sbat=var_Sbat)


        b = optim(par = 0.5, method = "L-BFGS-B"
                  , lower = 0, upper = 1
                  , fn = S.root
                  , prev = prev
                  , y = Phen.ds$event
                  , R = Phen.ds$Resid.coxph.w
                  , w = Phen.ds$weight
                  , mu = mu
                  , N = 2*n1
                  , n.ext = n.ext
                  , sigma2 = obj$sigma2
                  , TPR = obj$TPR
        )$par[1]


        return(cbind(obj$TPR, obj$sigma2, b=b))

      })%>%
        do.call("cbind",.) %>%as_tibble()%>%
        setNames(paste0(rep(c("TPR_","sigma2_","b_"),length(n.ref)),
                        rep(n.ref,each=3) ))



      data=data%>%cbind(.,temp)

    })%>%
      do.call("rbind",.)%>%as_tibble()%>%arrange(id)
    t2=Sys.time()
    print(t2-t1)

    info.SNP = info.SNP.all %>%filter(causal=="no")%>%arrange(id)
    info.SNP.causal = info.SNP.all %>%filter(causal=="yes")%>%arrange(id)
    rm(info.SNP.all);gc()

    save(info.SNP.causal, info.SNP, file = outputfile )

  }

}
