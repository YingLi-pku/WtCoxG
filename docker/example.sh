# step1: fit null model and calculate batch effect p-values
docker run -v /docker/inst/:/test \
  wtcoxg Rscript WtCoxG_step1.R \
  --GenoFile=/test/GenoMat2.bed \
  --GenoMarkerIndexFile=/test/GenoMat2.bim \
  --GenoSampleFile=/test/GenoMat2.fam \
  --PhenoFile=/test/simuPHENO_WtSPAG.txt \
  --IndicatorColumn=SurvEvent \
  --SurvTimeColumn=SurvTime \
  --SampleIDColumn=IID \
  --RefAfFile=/test/RefMAFs1.txt \
  --RefPrevalence=0.1 \
  --formula='Surv(SurvTime,Indicator) ~ Cov1 + Cov2' \
  --OutputFile=/test/dockertest1 

# step2: gwas
docker run -v /docker/inst/:/test \
  wtcoxg Rscript WtCoxG_step2.R \
  --GenoFile=/test/GenoMat2.bed \
  --objWtCoxGFile=/test/dockertest1.RData \
  --RefPrevalence=0.1 \
  --OutputFile=/test/dockertest1

  
  
