# adapted from: https://bitbucket.org/schwarzlab/refphase/src/master/

library(readr)
library(gtools)
library(optparse)

option_list = list(
  make_option(
    c("-i", "--input_directory"),
    type = "character",
    default = NULL,
    help = "Path to mpileup files",
    metavar = "path"
  ),
  make_option(
    c("-o", "--output_directory"),
    type = "character",
    default = NULL,
    help = "Path to the ASCAT input directory (where output of this script will be placed)",
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
  )
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

input_directory = opt$input_directory
output_directory = opt$output_directory
normal_sample_name = opt$normal_sample_name
all_samples = unlist(strsplit(opt$samples, " "))
tumor_samples = setdiff(all_samples, normal_sample_name)

# print arguments:
print('Input directory:')
print(input_directory)
print('Output directory:')
print(output_directory)
print('Normal sample name:')
print(normal_sample_name)
print('Tumor samples:')
print(tumor_samples)
print('starting...')

chrom2integer = function(x) {
  print(unique(x))
  as.integer(gsub("M", "25", gsub("X", "23", gsub("Y", "24", gsub("chr", "", x)))))
}

replace_dots = function(df) {
df$alt <- replace(df$alt, df$alt == ".", 0)
df$ref <- replace(df$ref, df$ref == ".", 0)
df$alt = as.integer(df$alt)
df$ref = as.integer(df$ref)
  return(df)
}

column_names = c("chrom", "pos", "ref", "alt")
normal = as.data.frame(read_tsv(paste0(input_directory,'/',normal_sample_name,".mpileup.gz"), col_names = column_names, progress = FALSE))
normal = replace_dots(normal)

for (tumor_sample in tumor_samples) {
  tumor = as.data.frame(read_tsv(paste0(input_directory,'/',tumor_sample, ".mpileup.gz"), col_names = column_names, progress = FALSE))
  tumor = replace_dots(tumor)
  # Merge tumor and normal samples to only keep the positions that are shared
  # between both samples
  merged = merge(normal, tumor, by = c("chrom", "pos"), suffixes = c("_normal", "_tumor"))

  merged$pos = as.integer(merged$pos)
  merged = merged[order(chrom2integer(merged$chrom), merged$pos), ]
  # replace dots with zeros:
  
  total_reads_tumor = sum(merged$alt_tumor, na.rm = TRUE) + sum(merged$ref_tumor, na.rm = TRUE)
  total_reads_normal = sum(merged$alt_normal, na.rm = TRUE) + sum(merged$ref_normal, na.rm = TRUE)

  merged$baf_tumor = merged$alt_tumor / (merged$alt_tumor + merged$ref_tumor)
  merged$logr_tumor = log2(((merged$alt_tumor + merged$ref_tumor) / total_reads_tumor) / ((merged$alt_normal + merged$ref_normal) / total_reads_normal))

  merged$baf_normal = merged$alt_normal / (merged$alt_normal + merged$ref_normal)

  # logR of normal should always be 0
  merged$logr_normal = merged$logr_tumor
  merged$logr_normal[] = 0

  # Remove NAs and +/-Inf
  keep_numeric = apply(merged[, c("logr_tumor", "logr_normal", "baf_normal", "baf_tumor")], 1, function(x) all(unlist(Map(is.finite, x))))
  merged = merged[keep_numeric, ]

  # Write the specific tsv format that we can later use with refphase
  # Set positions with germline BAF < 0.1 or BAF > 0.9 as germline homozygous,
  # and the rest as heterozygous
  refphase_snps=data.frame(chrom = merged$chrom, pos = merged$pos, baf = merged$baf_tumor, logr = merged$logr_tumor, germline_zygosity = ifelse(abs(merged$baf_normal - 0.5) > 0.4, "hom", "het"))

  hom = refphase_snps[refphase_snps$germline_zygosity == "hom",]
  het = refphase_snps[refphase_snps$germline_zygosity == "het",]
  write.table(refphase_snps, file = paste0(output_directory,'/',tumor_sample, "_refphase_snps.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  # ASCAT only takes integer chromosome names
  merged$chrom = gsub("Y", "24",
                  gsub("X", "23",
                  gsub("chr", "", merged$chrom)))

  # Write the specific input format files that ASCAT expects
  rownames(merged) = paste0("SNP", seq_len(nrow(merged)))
  write.table(data.frame(chrs = merged$chrom, pos = merged$pos, sample = merged$baf_tumor, row.names = rownames(merged)), file = paste0(output_directory,'/',tumor_sample, "_baf_tumor.tsv"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

  # Tumor
  dat = data.frame(chrs = merged$chrom, pos = merged$pos, sample = merged$baf_tumor, row.names = rownames(merged))
  colnames(dat)[[3]] = tumor_sample

  write.table(dat, file = paste0(output_directory,'/',tumor_sample, "_baf_tumor.tsv"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

  dat[[tumor_sample]] = merged$logr_tumor
  write.table(dat, file = paste0(output_directory,'/',tumor_sample, "_logr_tumor.tsv"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

  # Normal
  dat[[tumor_sample]] = merged$baf_normal
  write.table(dat, file = paste0(output_directory,'/',tumor_sample, "_baf_normal.tsv"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

  dat[[tumor_sample]] = merged$logr_normal
  write.table(dat, file = paste0(output_directory,'/',tumor_sample, "_logr_normal.tsv"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

}
