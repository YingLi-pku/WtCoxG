# sbatch -J Power --mem=12000M -t 1-0:0 --array=1-60 -o log/%A_%a.log --wrap='Rscript power.R $SLURM_ARRAY_TASK_ID'
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
library(SPACox)
source("/gdata01/user/liying/simulation-2023-04-12/code/function-1021.R")
###################################################################################################
params = expand.grid(prev = c(0.01,0.05,0.1)
                     ,h2 = seq(0.002, 0.01, 0.002)
                     ,h2 = c(0.01, seq(0.012,0.0155,0.0005))
                     ,n1 = c( 2000)
                     ,bc = c(0, 0.03, 0.12, 0.5)
                     ,theta = c( 0.5 )

)

h2 = params[n.cpu, "h2"]
prev = params[n.cpu, "prev"]
bc= params[n.cpu, "bc"]
n1= params[n.cpu, "n1"]
theta=params[n.cpu, "theta"]

outputfile = paste0("./power_n",n1,
                    "/prev",prev,
                    "_her",h2,
                    "_batcheffect",bc,
                    "_theta",theta,
                    ".csv")
if(!file.exists(outputfile)){

  results = lapply(1:100, function(v){
    cat("v:",v,"-----------------------------------------------------------------------------","\n")

    print(params[n.cpu,])

    ##load pheno
    cat("load phenotype and geno and SNPinfo -----------------------------------\n")
    load(paste0("./Pheno/prev",prev,
                "_her",h2,
                "_n",n1,
                "_v",v,
                ".RData"))

    load(paste0("./SNPinfo_norm/prev",prev,
                "_her",h2,
                "_bc",bc,
                "_theta",theta,
                "_n",n1,
                "_v",v,
                ".RData"))

    rm(info.SNP)

    obj.null = SPACox_Null_Model(Surv(time,event)~Cov1 + Cov2, data=Phen.ds,
                                 pIDs=Phen.ds$ID, gIDs=rownames(G.ds.causal))

    G = G.ds.causal
    Phen.data = Phen.ds
    n.ref = c("4e3", "1e4", "5e4", "2e5") #size of reference population



    t1 = Sys.time()
    GWAS = lapply(1:ncol(G), function(i){
      if(i%%100==0)cat(i,"\n")
      g = G[,i]
      y = Phen.data$event
      R = Phen.data$Resid.coxph.w
      trueMAF=info.SNP.causal$trueMAF[i]
      N = length(g)

      #ADuLT
      if(sum(g)<10|sum(2-g)<10){
        WtCoxG.true=ADuLT=WtCoxG0k = SPACox=NA
      }else{
        #adult
        obj.lm = lm(gen_est ~ g + Cov1 + Cov2, Phen.data)
        t.score = summary(obj.lm)$coefficients["g","t value"]
        ADuLT = summary(obj.lm)$coefficients["g","Pr(>|t|)"]

        #true MAF
        S = sum(R*(g - 2*trueMAF))
        S.var = sum(R^2)*2*trueMAF*(1-trueMAF)
        z=S/sqrt(S.var)
        WtCoxG.true = pnorm(abs(z), lower.tail =F)*2

        #. WtCoxG0k
        WtCoxG0k= cal.cutoff.p(g=g,
                               R=R,
                               w=Phen.data$weight,
                               TPR=0,
                               p_bat=1)[1]

        #. SPACox
        g= as.matrix(g)
        colnames(g) = colnames(G)[i]
        rownames(g) = rownames(G )

        scox = SPACox(obj.null, g) %>% as_tibble()
        SPACox = scox$p.value.spa

      }



      temp = lapply(n.ref , function(x){
        TPR = as.numeric(info.SNP.causal[i,paste0("TPR_",x)])
        p_bat = as.numeric(info.SNP.causal[i,paste0("p",x)])
        mu.ext = as.numeric(info.SNP.causal[i,paste0("maf",x)])
        sigma2 = as.numeric(info.SNP.causal[i,paste0("sigma2_",x)])
        #  sigma2 =info.SNP.causal$sigma2.true[i]

        b = as.numeric(info.SNP.causal[i,paste0("b_",x)])


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
        do.call("cbind",.)%>%as_tibble()%>%
        setNames(paste0("WtCoxG",n.ref ))


      return(cbind(temp, WtCoxG.true, ADuLT, WtCoxG0k,SPACox))

    })%>%do.call("rbind",.)%>%as_tibble()

    t2=Sys.time()
    print(t2-t1)
    #apply(na.omit(GWAS),2,function(x){mean(x<5e-8)})

    GWAS = cbind(GWAS, info.SNP.causal)

    return(GWAS)

  })%>%do.call("rbind",.)%>%as_tibble()
  fwrite(results, outputfile )

}

