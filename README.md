# WtCoxG
WtCoxG is an accurate, powerful, and computationally efficient Cox-based approach to perform genome-wide time-to-event data analyses in study cohorts with case ascertainment.  

## Table of contents
  * [How to install WtCoxG](#How-to-install-WtCoxG)
  * [Introduction of WtCoxG](#Introduction-of-WtCoxG)
  * [Step-by-step Workflow](#Step-by-step-Workflow)
  * [Example of slurm scripts](#Example-of-slurm-scripts)
  * [Reproducibility](#Reproducibility)

## How to install WtCoxG
```
library(devtools)  # author version: 0.0.9
install_github("YingLi-pku/WtCoxG",dependencies = TRUE)
library(WtCoxG)
?WtCoxG  # manual of WtCoxG package
example(WtCoxG)  # run example
```
WtCoxG supports both BGEN file and PLINK file.  
Please do not hesitate to contact me (yingli@stu.pku.edu.cn) if you have any problem or suggestion!  

## Introduction of WtCoxG
**WtCoxG is a Cox-regression based method designed for time-to-event GWAS, which accounts for case ascertainment and gains power by utilizing external allele frequencies (AFs) from publicly available datasets.** WtCoxG consists of three main steps:  
- **Step 0:** If the external MAFs are available, we first test for batch effect between internal and external MAFs and calculate a batch effect p-value for each genetic variant. And parameters of batch effect are estimated according to the genome-wide batch effect p values.   
- **Step 1:** In the context of case ascertainment, we fit an Inverse Probability Weighting (IPW) Cox PH null model, in which subjects were assigned with different sampling weights. The sampling weight can be calculated according to disease prevalence in the population. The covariates includes but not limited to age, gender, principle components (PCs).  
- **Step 2:** we calculate score statistics for each genetic variant. To further boost statistical power, we incorporate external MAFs with batch effect p value > 0.1 into the score statistics. Then we approximate the distribution of the score statistics by Saddle approximation.  
![plot](https://github.com/YingLi-pku/WtCoxG/blob/main/Figure/Pipeline.png)  

### We support Dense GRM to adjust for sample relatedness  
To account for sample relatedness, we follow the strategy of GATE (Dey et al., 2022, Nature communications), which calculates the ratio of the variance of the score statistic with and without GRM. Therefore, when performing the genome-wide scan in Step 2, the score statistic is calibrated using the variance ratio.  

## Step-by-step Workflow

## Example of slurm scripts

## Reproducibility
Scripts to reproduce the experiments performed for the manuscript: **Applying weighted Cox regression to boost powers for time-to-event genome-wide association studies**  available in the **SimulationCode** folder.  
### Simulation studies
We conducted extensive numeric simulations to evaluate WtCoxG in terms of type I error rates and powers in the context of case-ascertainment. And We varied the settings of the external sample sizes and batch effect scenarios, and evaluated the power improvement of WtCoxG brought by external MAFs. And we compare WtCoxG with ADuLT, SPACox, GATE and gwasurvivr.  
### 1. Simulate Genotype
We simulated 60 causal variants and 300,000 null variants following a binomial distribution Binom(2, MAFs), in which MAFs were simulated following a uniform distribution. For both causal and null variants, 1/3 are rare variants simulated following U(0.01,0.05), one-third are low-frequency variants following U(0.05,0.1), and the remaining one-third are common variants following U(0.1,0.5).
```
Simulate-Geno.R  # R script to simulate genotypes
```
### 2. Simulate Phenotype
To fairly compare WtCoxG (Cox PH model-based) and ADuLT (liability threshold model-based, LTM), we simulated time-to-event phenotypes under Cox PH model and LTM model, respectively. To mimic the case oversampling process, we randomly selected 2000 cases and 2000 controls from the population. 
```
Simulate-Pheno-Cox.R  # R script to simulate time-to-event phenotypes under Cox PH model
Simulate-Pheno-LTM.R  # R script to simulate time-to-event phenotypes under LTM model
```
### 3. Simulate external MAFs and test for batch effect
Leveraging external MAFs could boost powers, but potential batch effect between internal and external data is inevitable. Here we simulated external MAFs $\mu_{ext}$ following a mixed distribution. A proportion of genetic variants are of high quality and free from batch effect, that is, the external MAFs are identical to the internal ones ($\mu_{ext}=\mu$). Meanwhile, the other genetic variants suffer from batch effect and the external MAFs are different from the internal ones. We increase the proportion of the variants with batch effect from 0%, 3%, 12% to 50%. For the variants with batch effect, we considered three settings to simulate external MAFs:  
- **Setting 1:** a normal distribution $N(\mu,\theta∙\mu(1-\mu))$ with a large variance $\theta=0.5$ (Setting 1A), a normal distribution with with a small variance θ=0.01 (Setting 1B)  
- **Setting 2:** a truncated normal distribution $TN(\mu,\theta∙\mu(1-\mu),0,1)$  
- **Setting 3:** a uniform distribution $U(0.1\mu,4\mu)$.  
We set the external sample sizes to 4,000, 10,000, 50,000, and 200,000. For each setting, we test for batch effect between internal and external MAFs across genome, and estimate the batch effect parameters.
```
Simulate-extMAF-norm.R        # R script to simulate external MAFs following Setting1
Simulate-extMAF-truncnorm.R   # R script to simulate external MAFs following Setting2
Simulate-extMAF-unif.R        # R script to simulate external MAFs following Setting3
```
### Type I error
In each scenario, we conduct 1e9 type I error simualtions
```
type1error.R 
```
### Power
In each scenario, we conduct 1e3 power simualtions
```
power.R
```
## References
Dey, R., Zhou, W., Kiiskinen, T. et al. Efficient and accurate frailty model approach for genome-wide survival association analysis in large-scale biobanks. Nat Commun 13, 5437 (2022). https://doi.org/10.1038/s41467-022-32885-x



