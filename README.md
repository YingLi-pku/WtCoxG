# WtCoxG
WtCoxG is an accurate, powerful, and computationally efficient Cox-based approach to perform genome-wide time-to-event data analyses in study cohorts with case ascertainment.  

## How to install and load this package
```
library(devtools)  # author version: 0.0.9
install_github("YingLi-pku/WtCoxG")
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
To account for sample relatedness, we following the strategy of GATE (Rounak Dey et al., 2022, Nature communications), which calculates the ratio of the variance of the score statistic with and without GRM. Therefore, when performing the genome-wide scan in Step 2, the score statistic is calibrated using the variance ratio.  
## Reproducibility
Scripts to reproduce the experiments performed for the manuscript: **Applying weighted Cox regression to boost powers for time-to-event genome-wide association studies**  
### Simulation studies








