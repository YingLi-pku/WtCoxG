#!/usr/bin/env Rscript

options(stringsAsFactors=F)

## load R libraries
library(WtCoxG)
#library(SAIGE, lib.loc="/net/hunt/zhowei/project/imbalancedCaseCtrlMixedModel/Rpackage_SPAGMMAT/installSAIGEFolder/GATE")
require(optparse) #install.packages("optparse")

print(sessionInfo())

## set list of cmd line arguments
option_list <- list(
  make_option("--GenoFile", type="character",default=NULL,
              help="Path to plink file/BGEN file, bed file or bgen file"),
  make_option("--GenoMarkerIndexFile", type="character", default=NULL,
              help="marker index file corresponding to GenoFile, bim file or bgi file."),
  make_option("--GenoSampleFile", type="character", default=NULL,
              help="sample file corresponding to GenoFile, fam file or sample file."),
              
  make_option("--GenomtxFile", type="character",default=NULL,
              help="Path to genotype matrix RData file, which includes the genotype matrix nameed Genomtx"),
  make_option("--AlleleOrder", type="character",default="ref-first",
              help="ref-first or alt-first"),
  make_option("--AllMarkers",type="logical",default=TRUE,
              help="if analysis all markers, TRUE or FALSE, the default value is TRUE, if FALSE, then IDtoIncludeFile must be specified"),
  make_option("--IDtoIncludeFile", type="character",default=NULL,
              help="Path to the IDtoIncludeFile, a file to store rsID of analysed SNPs"),

  make_option("--PhenoFile", type="character", default=NULL,
              help="Path to the phenotype file. The phenotype file must have at least three columns: the column of personal identifiers for all individuals, the column of whether the event occurred (0 or 1 or NA), the column of the time of occurrence."),
  make_option("--IndicatorColumn", type="character",default=NULL,
              help="Column name of the event indicator"),
  make_option("--SampleIDColumn", type="character",default=NULL,
              help="Column name of the sample ID"),
  make_option("--SurvTimeColumn", type="character",default=NULL,
              help="Column name of the survival time"),

  make_option("--RefAfFile", type="character", default=NULL,
              help="Path to the refrence AF file."),
  make_option("--RefPrevalence", type="double",default = 0.1,
              help="The value of population prevalence"),
  make_option("--formula", type="character",default = 'Surv(SurvTime , Indicator)~.',
              help="The formula for fitting null model, e.g., 'Surv(SurvTime , Indicator)~Cov1+Cov2'"),
  make_option("--sparseGRMFile", type="character",default = NULL,
              help="Path to the sparseGRM .txt file,  a three-column sparse GRM file with the first column as ID1,the second column as ID2, and the last column as Value (i.e., two times of kinship coefficient) "),
  make_option("--OutputFile", type="character",default="~/",
              help="Path to the output file (Prefix)")
  
  )
  

## list of options
parser <- OptionParser(usage="%prog [options]", option_list=option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)


if(!is.null(opt$GenomtxFile)){
  load(opt$GenomtxFile)
}else{Genomtx = NULL}


if(!is.null(opt$sparseGRMFile)){
  sparseGRM = fread("sparseGRMFile")
}else{sparseGRM = NULL}

GenoFileIndex = c(opt$GenoMarkerIndexFile, opt$GenoSampleFile)


obj.WtCoxG = TestforBatchEffect(GenoFile = opt$GenoFile,
                             GenoFileIndex = opt$GenoFileIndex,
                             OutputFile = paste0(opt$OutputFile,".markerInfo.txt"),
                             control=list(AlleleOrder = opt$AlleleOrder,
                                          AllMarkers = opt$AllMarkers,
                                          IDtoIncludeFile =  opt$IDtoIncludeFile,
                                          IndicatorColumn = opt$IndicatorColumn, 
                                          SampleIDColumn = opt$SampleIDColumn, 
                                          SurvTimeColumn = opt$SurvTimeColumn), # specify the column names of sampleID, event, and time
                             PhenoFile = opt$PhenoFile,
                             RefAfFile = opt$RefAfFile,
                             RefPrevalence = opt$RefPrevalence,
                             formula = as.formula(opt$formula),
                             sparseGRM = sparseGRM,
                             SNPnum=1e4)

fwrite (obj.WtCoxG$PhenoData, 
        file = paste0(opt$OutputFile,".resid.txt"),
        row.names = F, col.names = T, quote = F, sep = "\t")

save (obj.WtCoxG, 
        file = paste0(opt$OutputFile,".RData"))

