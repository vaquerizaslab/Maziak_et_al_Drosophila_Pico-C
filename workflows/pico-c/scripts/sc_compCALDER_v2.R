#!/usr/bin/env Rscript

# R script for CALDER subcompartments with arguments

## Load libraries
library(optparse)
library(CALDER)
# library.path <- .libPaths()
# library(CALDER, lib.loc = '/home/slloren/R/x86_64-pc-linux-gnu-library/4.1')

## Set up arguments
option_list = list(
  make_option(c('-f', '--contact_file_folder'), 
              type='character', 
              metavar='character', 
              default=NULL, 
              help='folder with contact file for each chromosome in the format pos_x pos_y contact_value'),
  make_option(c('-t', '--featuretrack'),
              type='character',
              metavar='character',
              default=NULL,
              help='feature track file in bigWig format (see CALDER v2.0 docs)'),
  make_option(c('-c', '--chrom'), 
              type='character', 
              metavar='character',
              default=NULL, 
              help='chromosome of the contact file, without the prefix "chr"'),
  make_option(c('-b', '--binsize'), 
              type='character',
              metavar='character',
              default=NULL,
              help='binsize of the contact file'),
  make_option(c('-o', '--outfolder'),
              type='character', 
              metavar='character',
              default=NULL,
              help='output folder name'))

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

## Checking all inputs are provided
# if (is.null(opt$contact_file) | is.null(opt$chrom) | is.null(opt$binsize) | is.null(opt$outfolder)) {
#   print_help(opt_parser)
#   stop('Check inputs: contact file, chromosome, binsize and output folder need to be provided')}

if (is.null(opt$contact_file) | is.null(opt$binsize) | is.null(opt$outfolder)) {
  print_help(opt_parser)
  stop('Check inputs: contact file, binsize and output folder need to be provided')}



### Accessing to arguments
# opt$contact_file
# opt$chrom
# opt$binsize
# opt$outfolder


## Check binsize and convert to bp if needed
bp_magnitudes <- c('kb'=1e3, 
                   'mb'=1e6)

#binsize_num <- regmatches(opt$binsize, regexpr('\\d+', opt$binsize))
#binsize_magnitud <- regmatches(opt$binsize, regexpr('\\D+', opt$binsize))
#
#if (length(binsize_magnitud) == 0) {
#  binsize_parsed = as.numeric(opt$binsize)
#} else if (length(binsize_magnitud) == 1) {
#  binsize_parsed = as.numeric(binsize_num)*bp_magnitudes[binsize_magnitud]
#}
#binsize_parsed

binsize_num <- as.numeric(gsub('[^0-9.]', '', opt$binsize))
binsize_magnitud <- gsub('[0-9.]', '', opt$binsize)

if (binsize_magnitud %in% names(bp_magnitudes)) {
  binsize_parsed <- binsize_num * bp_magnitudes[[binsize_magnitud]]
} else {
  binsize_parsed <- binsize_num
}

## prepare feature_track
library(rtracklayer)
feature_track_raw  = import(opt$featuretrack)
feature_track = data.table::as.data.table(feature_track_raw)[, c(1:3, 6)]

## prepare contact files

chrs = c("2L", "2R", "3L", "3R", "X")
contact_file_dump_list = as.list(list.files(opt$contact_file_folder, pattern = paste0(".*", opt$binsize, ".*_sparseMatrixReg_slim.txt"), full.names = T))
names(contact_file_dump_list) = c("chr2L", "chr2R", "chr3L", "chr3R", "chrX")
contact_file_dump_list

## CALDER function 

CALDER(contact_file_dump = contact_file_dump_list,
        chr=chrs,
        bin_size=binsize_parsed,
        genome = 'others',
        save_dir=opt$outfolder,
        feature_track=feature_track,
        sub_domains=FALSE,
        save_intermediate_data=TRUE,
        n_cores=1)


## Run like: Rscript --vanilla sc_subcompCALDERarg.R -f ${CONTACT_FILE} -c ${CHROM} -b ${BINSIZE} -o ${OUTFOLDER}
## e.g.,
## Rscript --vanilla sc_subcompCALDERarg.R -f scripts/20220420_EvaXenlae_HiC_subcompartmentsCalder/Endoderm_q42_500kb_downsample_filt_chr7S_MatrixRegions_slim.txt -c 7S -b 500kb -o tmp/20220420_CALDERtest/xlae10_chr7S_Rarg


message(date())
sessionInfo()

