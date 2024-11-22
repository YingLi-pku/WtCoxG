# cd /gdata01/user/liying/simulation-2023-04-12/code
# sbatch -J type1 --partition=bi2 -q low --exclude=node04 --mem=20000M -t 1-0:0 --array=1-1440 -o log/%A_%a.log --wrap='Rscript type1error0620.R $SLURM_ARRAY_TASK_ID'
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

###################################################################################################
params = expand.grid(prev = c(0.01,0.05,0.1)
                     #,h2 = c(0.01, seq(0.012,0.0155,0.0005))###2000
                     ,h2 =  c(0.012,seq(0.014,0.016,0.0005))##1000
                     #,h2 = seq(0.004,0.014,0.002)
                     #,h2 = seq(0.01,0.015,0.001)
                     ,n1 = c(1000)
                     ,bc = c(0, 0.03, 0.12, 0.5)
                     ,theta = c( 0.5 )
                     ,group=1:20
)

h2 = params[n.cpu, "h2"]
prev = params[n.cpu, "prev"]
bc= params[n.cpu, "bc"]
n1= params[n.cpu, "n1"]
theta=params[n.cpu, "theta"]
group=params[n.cpu, "group"]




for(part in 1:1){
  load( paste0("./Geno/null_n",n1,
               ".RData"))


  lapply(((group-1)*1+ 1):(group*1), function(v){
    cat("part",part,"v:",v,"-----------------------------------------------------------------------------","\n")


    outputfile = paste0("./type1error_n",n1,
                        "/prev",prev,
                        "_her",h2,
                        "_batcheffect",bc,
                        "_theta",theta,
                        "_v",v,
                        "_part",part,
                        ".csv")


    if(!file.exists(outputfile)){


      ##load pheno
      cat("load phenotype and SNPinfo -----------------------------------\n")
      load(paste0("./Pheno/prev",prev,
                  "_her",h2,
                  "_n",n1,
                  "_v",v,
                  ".RData"))

      #load(paste0("/gdata01/user/liying/simulation-2023-04-12/simu0620/SNPinfo1_3M/prev",prev,
      load(paste0("./SNPinfo_norm/prev",prev,
                  "_her",h2,
                  "_bc",bc,
                  "_theta",theta,
                  "_n",n1,
                  "_v",v,
                  #"_part",part,
                  ".RData"))
      rm(info.SNP.causal)

      G = G.ds.null
      rm(G.ds.null)
      Phen.data = Phen.ds
      # n.ref = c("4e3", "1e4", "5e4", "2e5") #size of reference population
      n.ref = c("4e3","2e5") #size of reference population

      t1 = Sys.time()
      GWAS = lapply(1:ncol(G), function(i){
        if(i%%100000==0)cat(i,"\n")
        g = G[,i]
        y = Phen.data$event
        R = Phen.data$Resid.coxph.w

        #WtCoxG0k---------------------------------------
        WtCoxG0k= cal.cutoff.p(g=g,
                               R=R,
                               w=Phen.data$weight,
                               TPR=0,
                               p_bat=1)[1]

        #WtCoxG.ref---------------------------------------
        temp = lapply(n.ref , function(x){
          TPR = as.numeric(info.SNP[i,paste0("TPR_",x)])
          p_bat = as.numeric(info.SNP[i,paste0("p",x)])
          mu.ext = as.numeric(info.SNP[i,paste0("maf",x)])
          sigma2 = as.numeric(info.SNP[i,paste0("sigma2_",x)])
          # sigma2 =info.SNP$sigma2.true[i]

          b = as.numeric(info.SNP[i,paste0("b_",x)])


          WtCoxG= cal.cutoff.p(g=g,
                               R=R,
                               w=Phen.data$weight,
                               TPR=TPR,
                               p_bat=p_bat,
                               mu.ext = mu.ext,
                               n.ext = as.numeric(x),
                               sigma2 =sigma2,
                               b= b,
                               #b=b,
                               p_cut =0.1)[1]

          return(WtCoxG)

        })%>%
          do.call("cbind",.)


        return( cbind( WtCoxG0k, temp ))

      })%>%do.call("rbind",.)%>%as_tibble()%>%
        setNames(c("WtCoxG0k",paste0("WtCoxG",n.ref )))

      t2=Sys.time()
      print(t2-t1)
      #apply(GWAS,2,function(x){mean(na.omit(x<5e-5))})

      GWAS = cbind(GWAS, info.SNP)


      fwrite(GWAS, outputfile )

    }



  })
}

