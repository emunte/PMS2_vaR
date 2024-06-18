# Description: Run the PMS2_vaR pipeline
# USAGE: Rscript run_PMS2_vaR.R [-t tools_file] [-b bamTxt_file] [-r reference_bam_file] [-g genome_assembly] [-v vardictJava_params] [-o outputdir] [-n folder_name]

#libs
library(yaml)
library("optparse")
library(dplyr)
library(stringr)
library(vcfR)
library(GenomicRanges)
################################################################################

library(jsonlite)
library(httr)
library(VariantAnnotation)
library(stringr)
library(purrr)
library(Rsamtools)
# Check run_PMS2_vaR.R is being called from PMS2_vaR folder

if (length(list.files(pattern = "run_PMS2_vaR.R"))== 0){
  cat("Sorry, run_PMS2_vaR.R .R should be called from PMS2_vaR.R  folder\n")
  quit()
}


# Build options list
 option_list <- list(
   make_option(c("-t", "--tools"), type="character", default="params/tools.yaml",
               help="Path to tools file (yaml)", metavar="character"),
   make_option(c("-b", "--bamTxt"), type="character", default="",
               help="Path to bams TXT routes. Copy all full paths of your bam files to a txt-file", metavar="character"),
   make_option(c("-r", "--reference"), type="character", default="",
               help="Full path to reference modified file (fa), with no PMS2CL sequence.", metavar="character"),
   make_option(c("-g", "--genome"), type="character", default="hg19",
               help="Assembly used hg19 or hg38", metavar="character"),
   make_option(c("-v", "--vardictjava"), type="character", default="params/vardictjavaParams.yaml",
               help="Full path to vardict java params", metavar="character"),
   make_option(c("-o", "--outputdir"), type="character", default="",
               help="Output directory to store results", metavar="character"),
   make_option(c("-n", "--samplesname"), type="character", default="",
               help="A name to store the results of all the samples", metavar="character")
 );

opt_parser <- OptionParser(option_list=option_list);

# Load params
args <- parse_args(opt_parser);
tools <- yaml.load_file(args$tools)
vardict <- yaml.load_file(args$vardictjava)

#Checks
assertthat::assert_that(file.exists(args$reference), msg="Please enter a reference file")
assertthat::assert_that(stringr::str_detect(args$reference, ".fa"), msg="Please enter a valid reference file (fasta format)")
assertthat::assert_that((args$genome=="hg19" | args$genome=="hg38"), msg="Assembly can only be hg19 or hg38")


#create logs/output folder if not exists
#dir.create("logs", showWarnings = F)
#dir.create("output", showWarnings = F)


#Load functions and variables
source(file.path(getwd(), "functions/PMS2_vaR_functions.R"))
source(file.path(getwd(), "utils/ROIS.R"))
source(file.path(getwd(), "utils/classification.R"))

#Execute the pipeline
output.samples.dir <- file.path(args$outputdir, paste(args$samplesname, Sys.Date(), sep="_"))
print(output.samples.dir)
message("Re-alignment")
pms2_realignment(args, tools, vardict, output.samples.dir)
message("zip tabix function")
zip_tabix(output.samples.dir, "*VARDICT.vcf$")
message("convert to txt file")
convert_vcf_txt(resultsDir = output.samples.dir, vcf.pattern = "*VARDICT.vcf.gz$", rng, args)
message("merge pipelines")
merge_pipelines(resultsDir = output.samples.dir, vars.paralogous=vars.paralogous, classification=classification, args=args)














