# cd /gdata01/user/liying/simulation-2023-04-12/code
# sbatch -J SimuPheno --partition=bi1 --exclude=node03,node02 --mem=4000M -t 1-0:0 --array=1-75 -o log/%A_%a.log --wrap='Rscript Generate-Pheno-LTM.R $SLURM_ARRAY_TASK_ID'
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
library(LTFHPlus)

####################################################################################################
params = expand.grid(prev = c(0.01,0.05,0.1)
                     ,n1 = c( 2000 )
                     ,h2 = seq(0.002, 0.01, 0.002)
                     ,group=1:5)

prev = params[n.cpu, "prev"]
n1 =params[n.cpu, "n1"]
h2 = params[n.cpu, "h2"]
group=params[n.cpu, "group"]

###load genotypes
load(paste0("./Geno/causal_prev",prev,
            "_n",n1,
            ".RData"))

###standardized genotype matrix
G.stad = lapply(1: ncol(G.pop), function(i){
  mu = info.SNP.causal$trueMAF[i]
  g = (G.pop[,i] - 2*mu)/sqrt(2*mu*(1-mu))
}) %>%
  do.call("cbind",.) %>% as.matrix()
###sum genotypes up for each individuals
nSNP = ncol(G.stad)
combG = apply(G.stad, 1, sum)
rm(G.stad);gc()
N = nrow(G.pop)
n = 2*n1

lapply(((group-1)* 20+ 1):(group*20), function(v){
  cat("v:",v,"-----------------------------------------------------------------------------","\n")

  outputfile = paste0("./Pheno_LTM/prev",prev,
                      "_her",h2,
                      "_n",n1,
                      "_v",v,
                      ".RData")
  gamma1 =  sqrt( h2)

  if(!file.exists(outputfile)){

    ##generate population under LTM
    Phen.pop = data.frame(ID = paste0("IID-",1:N),
                          IDnum = 1:N
    )%>%
      mutate(true_gen_lia = gamma1* combG,
             true_full_lia = true_gen_lia + rnorm(n(),0,sqrt(1 - nSNP*h2)),
             event = (true_full_lia > qnorm(prev,lower.tail = F)) + 0L)

    aoo = lapply(1:nrow(Phen.pop), function(i){
      t=convert_liability_to_aoo(liability = Phen.pop$true_full_lia[i],pop_prev = prev)
    })%>%unlist()

    Phen.pop = Phen.pop %>%
      mutate(aoo=aoo,
             time = ifelse(event==1, aoo, runif(n(), 10, 100)),
             fam_id = IDnum,
             role = rep("o",N) )%>%
      arrange(time)%>% mutate(prop  = pmin((cumsum(event) + 1)/n(), prev))


    ## case oversampling--------------------------------------------------------

    Phen.ds = rbind(Phen.pop[Phen.pop$event == 1,],
                    Phen.pop[Phen.pop$event == 0,] %>% slice_sample(n=n-sum(Phen.pop$event))) %>%
      mutate(er = mean(event),
             weight = ifelse(event==1, 1, (1-prev)/prev/((1-er)/er)))%>%
      mutate(est   = qnorm(prop, lower.tail = F),
             age   = LTFHPlus::convert_cir_to_age(cir = prop, pop_prev = prev, slope = 0.2, mid_point = 50),
             lower = ifelse(event, LTFHPlus::convert_age_to_thresh(age = age, pop_prev = prev, mid_point = 50, slope = 0.2), -Inf),
             upper = ifelse(event, LTFHPlus::convert_age_to_thresh(age = age, pop_prev = prev, mid_point = 50, slope = 0.2), #alternatively, set upper bound to Inf for cases
                            LTFHPlus::convert_age_to_thresh(age = age, pop_prev = prev, mid_point = 50, slope = 0.2)),
             status = event + 0L)

    G.ds.causal = G.pop[Phen.ds$ID,]

    ## estimate liability for oversampled individuals---------------------------
    t1=Sys.time()
    gen_lia = estimate_liability(.tbl = Phen.ds[,c("fam_id", "ID", "role", "lower", "upper")],
                                 h2 = h2*nSNP, pid = "ID", fam_id = "fam_id")
    t2=Sys.time()
    print(t2-t1)

    Phen.ds = Phen.ds %>% mutate(gen_est = as.numeric(gen_lia$genetic_est)) %>%
      select(ID, IDnum, event, time, er, weight, age, true_gen_lia, true_full_lia, gen_est )

    ## fit null Weighted CoxPH model using oversampled individuals--------------
    obj.coxph.w = coxph(Surv(time, event) ~ 1 , Phen.ds,
                        weight = weight,
                        robust = T) # weighted coxph model
    Phen.ds = Phen.ds %>% mutate(Resid.coxph.w = obj.coxph.w$residuals)

    save(Phen.ds, G.ds.causal, file = outputfile )

  }

})









