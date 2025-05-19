#!/usr/bin/env Rscript

options(stringsAsFactors=F)

## load R libraries
library(WtCoxG)
library(dplyr)
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
              help="if analysis all markers"),
  make_option("--IDtoIncludeFile", type="character",default=NULL,
              help="Path to the IDtoIncludeFile"),
  
  make_option("--objWtCoxGFile", type="character",default=NULL,
              help="Path to the output .RData file of WtCoxG_step1"),
  make_option("--OutputFile", type="character",default="~/",
              help="Path to the output file"),
  make_option("--RefPrevalence", type="double",default = 0.1,
              help="The value of population prevalence"),

  make_option("--PhenoFile", type="character", default=NULL,
              help="Path to the phenotype file. The phenotype file must have at least three columns: the column of personal identifiers for all individuals, the column of whether the event occurred (0 or 1 or NA), the column of the time of occurrence."),
  make_option("--mergeGenoInfoFile", type="character", default=NULL,
              help="Path to the refrence AF file.")
  
)


## list of options
parser <- OptionParser(usage="%prog [options]", option_list=option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)


if(!is.null(opt$GenomtxFile)){
  load(opt$GenomtxFile)
}else{Genomtx = NULL}


if(!is.null(opt$objWtCoxGFile)){
  load(opt$objWtCoxGFile)
}else{obj.WtCoxG = NULL}

GenoFileIndex = c(opt$GenoMarkerIndexFile, opt$GenoSampleFile)

gwas = WtCoxG(GenoFile = opt$GenoFile,
              GenoFileIndex = opt$GenoFileIndex,
              Geno.mtx = Genomtx ,
              control=list(AlleleOrder = opt$AlleleOrder,
                           AllMarkers = opt$AllMarkers,
                           IDtoIncludeFile = opt$IDtoIncludeFile
                           ),
              obj.WtCoxG = obj.WtCoxG,
              PhenoFile = opt$PhenoFile,
              mergeGenoInfoFile = opt$mergeGenoInfoFile,
              RefPrevalence = opt$RefPrevalence,
              OutputFile = paste0(opt$OutputFile,".gwas.txt"))

