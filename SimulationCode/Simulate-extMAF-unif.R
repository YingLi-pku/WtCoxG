# sbatch -J unif -t 1-0:0 --array=1-1080 -o log/%A_%a.log --wrap='Rscript Simulate-extMAF-unif.R $SLURM_ARRAY_TASK_ID'
args=commandArgs(TRUE)
print(args)
print(sessionInfo())
n.cpu=as.numeric(args)
##########################################

library(data.table)
library(tidyr)
library(dplyr)

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
                      ".txt")

  if(!file.exists(outputfile)){


    ##update SNP info--------------------------------------------------------------------------
    cat("load causal SNP  -----------------------------------\n")
    load( paste0("./Geno/prev",prev,
                 "_her",h2,
                 "_n",n1,
                 "_v",v,
                 ".RData") )
    #colnames(G.ds.causal) = paste0("causalSNP-", 1:ncol(G.ds.causal))
    cat("load non causal SNP  -----------------------------------\n")

    load(paste0("./Geno/null_n",n1,
                ".RData"))
    #colnames(G.ds.null) = paste0("nullSNP-", 1:ncol(G.ds.null))


    G.all = cbind(G.ds.null, G.ds.causal)
    rm(G.ds.null, G.ds.causal)
    info.SNP.all = rbind(info.SNP, info.SNP.causal) %>%
      mutate(CHR=1, POS=1:n(),
             REF= "A", ALT = "T",
             ID = colnames(G.all),
             batcheffect = rbinom(n(),1, bc),
             trueMAF_ref = ifelse(batcheffect==0, trueMAF,                            # if batch effect=0 then mu.ext = mu
                                  runif(n(), 0.1*trueMAF, 4*trueMAF )),
             trueMAF_ref = pmax(0, pmin(1, trueMAF_ref)),
             AN_ref = 8e3,
             AF_ref = rbinom(n(), AN_ref, trueMAF_ref)/AN_ref
             )


    fwrite(info.SNP.all, file =outputfile, col.names=T, sep="\t")

  }

}
