# adapted from: https://bitbucket.org/schwarzlab/refphase/src/master/
library(ASCAT)
library(optparse)

option_list = list(
  make_option(
    c("-i", "--input_directory"),
    type = "character",
    default = NULL,
    help = "Path to ASCAT input",
    metavar = "path"
  ),
  make_option(
    c("-o", "--output_directory"),
    type = "character",
    default = NULL,
    help = "Path to Refphase output directory",
    metavar = "path"
  ),
  make_option(
    c("-n", "--normal_sample_name"),
    type = "character",
    default = NULL,
    help = "Name of the normal sample",
    metavar = "string"
  ),
  make_option(
    c("-t", "--samples"),
    type = "character",
    default = NULL,
    help = "Array of tumour sample names",
    metavar = "string"
  ),
  make_option(
    c("-c", "--case_id"),
    type = "character",
    default = NULL,
    help = "case id",
    metavar = "string"
  )
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

ascat_dir = opt$input_directory
refphase_output = opt$output_directory
normal_sample_name = opt$normal_sample_name
all_samples = unlist(strsplit(opt$samples, " "))
tumor_samples = setdiff(all_samples, normal_sample_name)
case_id = opt$case_id

input_dir = paste0(ascat_dir, "/input")
output_dir = paste0(ascat_dir, "/output")

ascat_input_name_rds = paste0(output_dir,"/","ascat_input.rds")
ascat_output_name_rds = paste0(output_dir,"/","ascat_output.rds")

# change to the input directory:
setwd(input_dir)


# Do some extra work to create a table that will be needed when we run refphase in the next step. This is not needed by ASCAT itself.
refphase_sample_data <- data.frame(sample_id = tumor_samples, segmentation = paste0(tumor_samples, "_refphase_segs.tsv"), snps = paste0(tumor_samples, "_refphase_snps.tsv"), purity = NA, ploidy = NA, row.names = tumor_samples)

# ascat_input and ascat_output are required for refphase later on
# if ascat outptu already exists, load it:
if (file.exists(ascat_input_name_rds) && file.exists(ascat_output_name_rds)) {
  print('loading ascat input and output from RDS files')
  ascat_input <- readRDS(ascat_input_name_rds)
  ascat_output <- readRDS(ascat_output_name_rds)
} else {
  print('running ascat')
  ascat_input <- list()
  ascat_output <- list()

  for (tumor_sample in tumor_samples) {
    cur_ascat_input <- ascat.loadData(Tumor_LogR_file = paste0(tumor_sample, "_logr_tumor.tsv"),
                              Tumor_BAF_file = paste0(tumor_sample, "_baf_tumor.tsv"),
                              Germline_LogR_file = paste0(tumor_sample, "_logr_normal.tsv"),
                              Germline_BAF_file = paste0(tumor_sample, "_baf_normal.tsv"))
    ascat.plotRawData(cur_ascat_input)
    cur_ascat_input <- ascat.aspcf(cur_ascat_input)
    ascat.plotSegmentedData(cur_ascat_input)
    cur_ascat_output <- ascat.runAscat(cur_ascat_input, gamma = 1.0) #
    cur_ascat_input$SNPpos$Chromosome = cur_ascat_input$SNPpos$chrs
    cur_ascat_input$SNPpos$Position = cur_ascat_input$SNPpos$pos
    ascat_input[[tumor_sample]] <- cur_ascat_input
    ascat_output[[tumor_sample]] <- cur_ascat_output

    # Save the segmentation in the default ASCAT format
    write.table(cur_ascat_output$segments, file = paste0(output_dir,"/",tumor_sample, "_segments.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

    refphase_sample_data[tumor_sample, "purity"] <- cur_ascat_output$aberrantcellfraction[[1]]
    refphase_sample_data[tumor_sample, "ploidy"] <- cur_ascat_output$ploidy[[1]]
  }
  # save ascat_input and ascat_output lists as RDS files:
  saveRDS(ascat_input, file = ascat_input_name_rds)
  saveRDS(ascat_output, file = ascat_output_name_rds)
  out_name = paste0(output_dir,"/","refphase-sample-data.tsv")
  write.table(refphase_sample_data, file = out_name, sep = "\t", quote = FALSE, row.names = FALSE)
}


# adjust data:
for (tumor_sample in tumor_samples){
  print(tumor_sample)
  ascat_input[[tumor_sample]]$SNPpos$Chromosome = ascat_input[[tumor_sample]]$SNPpos$chrs
  ascat_input[[tumor_sample]]$SNPpos$Position = ascat_input[[tumor_sample]]$SNPpos$pos
}


# run refphase:
print('running refphase')
library(refphase)
library(GenomicRanges)
library(IRanges)
# Load the ascat input and output
refphase_input <- refphase_load(
    data_format = "ascat3", samples = tumor_samples,
    ascat_input = ascat_input, ascat_output = ascat_output
)

# (Optional) If your data shows reference bias, which presents itself as BAFs
# in regions with a balanced copy number that are systematically shifted away
# from 0.5 (usually something like 0.47), this can try to correct for that.
# refphase_input <- center_baf(refphase_input)

# (Optional) Fit SNP logr data to improve copy number re-estimation in refphase,
# when using the default ASCAT formula-based method fo re-estimating copy
# numbers
# refphase_input <- fit_logr_to_ascat(refphase_input)

# Run refphase on the experimental data
results <- refphase(refphase_input)

# go to output dir:
setwd(refphase_output)

#store entire results object as RDS:
#saveRDS(results, file = paste0(case_id, "-refphase-results.rds"))
#print ('results saved as RDS')
# for downstream compatibility, save also as RData:
save(results, file = paste0(case_id, "-refphase-results.RData"))


write_segs(results$phased_segs, file = paste0(case_id, "-refphase-segmentation.tsv"))

# (optional) output the SNPs, including phasing information
write_snps(results$phased_snps, file = paste0(case_id, "-refphase-phased-snps.tsv.gz"))

# Ploidy might have changed, if we have updated any copy numbers
write.table(results$sample_data, file = paste0(case_id, "-refphase-sample-data-updated.tsv"), sep = "\t", row.names = FALSE)
