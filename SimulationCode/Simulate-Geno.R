# source activate myEnvR
# cd /gdata01/user/liying/simulation-2023-04-12/code
# sbatch -J geno --partition=bi1 --mem=20000M -t 1-0:0 --array=1-1 -o log/%A_%a.log --wrap='Rscript Generate-Geno-0620.R $SLURM_ARRAY_TASK_ID'
args=commandArgs(TRUE)
print(args)
print(sessionInfo())
n.cpu=as.numeric(args)
##########################################--------------------------------
library(tidyr)
library(dplyr)
library(data.table)
set.seed(n.cpu)

params = expand.grid(n1 = 2000,
                     prev = c(0.01, 0.05, 0.1))
prev  =params[n.cpu, "prev"]                          ### sample size after oversampling
n1 = params[n.cpu, "n1"]                              ### sample size after oversampling

# ###generate population genotypes for power----------------------------------------
nSNP.causal = 60
N = n1 / prev
set.seed(n.cpu)

mu.causal = c(runif(nSNP.causal/3, 0.01, 0.05),
              runif(nSNP.causal/3, 0.05, 0.1),
              runif(nSNP.causal/3, 0.1, 0.5))

G.pop= lapply(mu.causal, function(m){return(rbinom(N,2,m))}) %>%
  do.call("cbind",.) %>% as.matrix()
rownames(G.pop) = paste0("IID-",1 : N)
colnames(G.pop) = paste0("causalSNP-",1 : nSNP.causal)

info.SNP.causal = data.frame(id = colnames(G.pop), trueMAF = mu.causal)
save( G.pop
     ,info.SNP.causal
     ,file = paste0("./Geno/causal_prev",prev,
                    "_n",n1,
                    ".RData"))
rm(G.pop)
rm(info.SNP.causal)
gc()

cat("Generate genotypes for power \n")

##generate oversampled genotypes for type 1 error--------------------
nSNP.null = 3e5
n = 2*n1
for(part in 1:10){
  cat(part,"\n")
  set.seed(part)
  mu.null = c(runif(nSNP.null/3, 0.01, 0.05),
              runif(nSNP.null/3, 0.05, 0.1),
              runif(nSNP.null/3, 0.1, 0.5))


  G.ds.null= lapply(mu.null, function(m){return(rbinom(n,2,m))}) %>%
    do.call("cbind",.) %>% as.matrix()
  rownames(G.ds.null) = paste0("IID-",1 : n)
  colnames(G.ds.null) = paste0("nullSNP-",1 : nSNP.null)

  info.SNP = data.frame(id = colnames(G.ds.null), trueMAF = mu.null)
  save( G.ds.null,info.SNP,
        file = paste0("./Geno/null_n",n1,
                      "_part",part,
                      ".RData"))
  cat("Generate genotypes for type1error \n")

}





