# sbatch -J SimuPheno --mem=4000M -t 1-0:0 --array=1-1500 -o log/%A_%a.log --wrap='Rscript Simulate-Pheno-Cox.R $SLURM_ARRAY_TASK_ID'
args=commandArgs(TRUE)
print(args)
print(sessionInfo())
n.cpu=as.numeric(args)
##########################################
library(data.table)
library(tidyr)
library(dplyr)

#### lower function to estimate beta0 given an event rate. Will be used in data.simu.surv().
f.surv = function(N,                  # Sample size
                  event.rate,         # Event rate
                  lamda,              # Scale parameter
                  gamma,              # Genetic effect
                  g,                  # Genotype vector
                  seed,               # random seed
                  bVec = 0)           # Additional effect, could be random effect. If bVec = 0 (default), then no additional effect is included.
{
  scale0 = lamda
  shape0 = 2
  set.seed(seed)
  X1 = rnorm(N)            # Covariate vector 1
  X2 = rbinom(N, 1, 0.5)   # Covariate vector 2
  betas = c(0.5, 0.5)      # Coefficient vector of fixed effects
  cens = rweibull(N, shape = 1, scale = 0.15)
  eps <- runif(N, 0, 1)
  time = (-log(eps) * exp(- betas[1] * X1 - betas[2] * X2 - gamma * g - bVec)) ^ (1/shape0) * scale0
  surv.time = pmin(time, cens)
  event = ifelse(time < cens, 1, 0)
  re = mean(event) - event.rate
  return(re)
}

#### Main function to simulate time-to-event phenotype
data.simu.surv = function(N,              # Sample size
                          event.rate,     # Event rate
                          gamma,          # Genetic effect
                          g,              # Genotype vector
                          bVec = 0)       # Additional effect, could be random effect. If bVec = 0 (default), then no additional effect is included.
{
  seed = sample(1e6, 1);
  cat("random number seed is", seed)
  scale0 = uniroot(f.surv, c(-1000000, 1000000),N = N, event.rate = event.rate, gamma = gamma, g = g, seed = seed, bVec = bVec, extendInt = "yes")
  scale0 = scale0$root
  shape0 = 2
  set.seed(seed)
  X1 = rnorm(N)          # Covariate vector 1
  X2 = rbinom(N,1,0.5)   # Covariate vector 2
  betas = c(0.5, 0.5)    # Coefficient vector of fixed effects
  cens = rweibull(N, shape = 1, scale = 0.15)
  eps <- runif(N, 0, 1)
  time = (-log(eps) * exp(- betas[1] * X1 - betas[2] * X2 - gamma * g - bVec)) ^ (1/shape0) * scale0
  surv.time = pmin(time, cens)
  event = ifelse(time < cens, 1, 0)
  out = data.frame(surv.time, event, X1, X2, bVec)
  return(out)
}


####################################################################################################
params = expand.grid(prev = c(0.01,0.05,0.1)
                     ,n1 = c( 2000 )
                     ,h2 = c(0.01, seq(0.012,0.0155,0.0005))   # heritability
                     ,group=1:50)

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
combG = apply(G.stad, 1, sum)
rm(G.stad);gc()
N = nrow(G.pop)
n = 2*n1

lapply(((group-1)*10 + 1):(group*10), function(v){
  cat("v:",v,"-----------------------------------------------------------------------------","\n")

  outputfile = paste0("./Pheno/prev",prev,
                      "_her",h2,
                      "_n",n1,
                      "_v",v,
                      ".txt")
  print(params[n.cpu,])

  if(!file.exists(outputfile)){
    nSNP = ncol(G.pop)
    gamma1 =  sqrt(  (5/16/(1-h2*nSNP)-5/16 )/ nSNP)   ##effect size
    ## generate population phenotype--------------------------------------------
    data = data.simu.surv(N = N, event.rate = prev, gamma = gamma1, g = combG)
    Phen.pop = data.frame(ID = paste0("IID-",1:N),
                          IDnum = 1:N,
                          event = data$event,
                          time = data$surv.time,
                          Cov1 = data$X1,
                          Cov2 = data$X2
    ) %>%
      as_tibble %>%
      mutate(fam_id = IDnum, role = rep("o",N))%>%
      arrange(time)%>%
      mutate(prop  = pmin((cumsum(event) + 1)/n(), prev))

    rm(data);gc()

    ## case oversampling--------------------------------------------------------

    Phen.ds = rbind(Phen.pop[Phen.pop$event == 1,],
                    Phen.pop[Phen.pop$event == 0,] %>% slice_sample(n=n-sum(Phen.pop$event)))%>%
      select(ID, event, time, Cov1, Cov2)

    G.ds.causal = G.pop[Phen.ds$ID,]


    save(G.ds.causal, info.SNP.causal, file = paste0("./Geno/prev",prev,
                                     "_her",h2,
                                     "_n",n1,
                                     "_v",v,
                                     ".RData") )
    fwrite(Phen.ds, file = outputfile, sep="\t", col.names=T)

  }

})




