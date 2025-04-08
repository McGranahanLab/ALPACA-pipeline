#!/bin/bash

# ==============================
# Setup
# ==============================


# To run this pipeline, you need to have Singularity v.3.11.3 available
# and a Gurobi license. The pipeline is designed to run in a Singularity container
# via Singularity shell.
# Start by building the Singularity image using the provided definition file:
# singularity build singularity/pipeline.sif singularity/pipeline.def

# ALPACA requires Gurobi license  - it is free for academic use and can be obtained via 
# https://www.gurobi.com/academia/academic-program-and-licenses/

# The Singularity image does not contain the license, make sure that your license is available
# from within the container if you run this pipeline with Singulairty:
image_path="singularity/pipeline.sif"
license_path=$GRB_LICENSE_FILE
work_dir=$(pwd)
singularity shell -B "$work_dir" -B "$license_path" $image_path 

# initialise conda:
source /opt/miniforge3/etc/profile.d/conda.sh
# or use `conda init` if you are not using singularity


# ==============================
# Assets
# ==============================


# define base path:
base_dir=$(pwd)
assets="${base_dir}/_assets"

# NB these files are not provided, please use your own reference genome, germline variants and panel of normals
reference="${assets}/reference/ucsc.hg19.fasta"
gnomad_af="${assets}/gnomad/af-only-gnomad.raw.sites.b37.vcf.gz"
pon="${assets}/pon/gatk4_mutect2_4136_pon_hg19_liftover.vcf.gz"
bams_path="${assets}/bams/"


# ==============================
# Parameters
# ==============================


# Specify the normal sample name - we assume there is only one normal sample.
# We assume that all the BAM files are named according to the following pattern:
# <sample_name>.bam
normal_sample_name="LTX0028_BS_GL--164c6d2dc495"
bam_normal_path="${bams_path}/${normal_sample_name}.bam"
# to get create an array of sample names:
SAMPLES=()
for bam_file in "$bams_path"/*.bam; do
    SAMPLE=$(basename "$bam_file" .bam)
    SAMPLES+=("$SAMPLE")
done
echo "Sample names: ${SAMPLES[@]}"

# conipher requires case identifier composed of a 'prefix' and a case number:
case_id=LTX0028
prefix=LTX

# retain only germline SNPS with a minimum depth of X
germline_snp_min_depth_normal=20
germline_snp_min_depth_tumour=20


# mutect memory
mutect2_memory=8
# mutect cpus:
mutect2_cpus=24

# CONIPHER input: filter variants based on these thresholds
depth_threrehold=30 # Threshold for total depth of mutation in each sample
varcount_threshold=10 # Threshold for minimum variant count of mutation in sample
germline_variant_read_count_threshold=5 # Reject variant if we find more than X reads for this variant in germline
germline_vaf_threshold=0.01 # Reject variant if gemline VAF is above X
variant_vaf_threshold=0.05  # Reject variant if VAF is below X


# ==============================================
# Define file paths and create directories
# ==============================================


# germline variants:
germline_dir="${base_dir}/germline_variants"
germline_variants="${germline_dir}/germline-variants.bcf"
germline_variants_pos="${germline_dir}/germline-variants.pos.gz"
mpileup=${base_dir}/mpileup

# somatic variants:
somatic_dir="${base_dir}/somatic"

# ascat dirs:
ascat_dir="${base_dir}/ascat"
ascat_input_dir="${ascat_dir}/input"
ascat_output_dir="${ascat_dir}/output"

# refphase dirs:
refphase_dir="${base_dir}/refphase"
refphase_dir_output="${refphase_dir}/output"

# conipher dirs:
conipher_dir="${base_dir}/conipher"
conipher_dir_input="${conipher_dir}/input"
conipher_dir_output="${conipher_dir}/output"

# alpaca dirs:
alpaca_directory="${base_dir}/alpaca"
alpaca_input_dir="${alpaca_directory}/input/${case_id}"
alpaca_output_dir="${alpaca_directory}/output/${case_id}"

# gnomad:
gnomad_common_biallelic="${assets}/gnomad/gnomad.common_biallelic.vcf.gz"

mkdir -p $mpileup
mkdir -p $germline_dir
mkdir -p $ascat_dir
mkdir -p $ascat_input_dir
mkdir -p $ascat_output_dir
mkdir -p $refphase_dir
mkdir -p $refphase_dir_output
mkdir -p $somatic_dir
mkdir -p $conipher_dir_input
mkdir -p $conipher_dir_output
mkdir -p $alpaca_input_dir
mkdir -p $alpaca_output_dir


# ==============================
# FUNCTIONS
# ==============================


get_chromosomes_style() {
    local fai_file="$1"
    local first_chr
    first_chr=$(head -n 1 "$fai_file" | cut -f1)
    if [[ "$first_chr" == chr* ]]; then
        echo "chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 \
              chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 \
              chr21 chr22 chrX chrY"
    else
        echo "1 2 3 4 5 6 7 8 9 10 \
              11 12 13 14 15 16 17 18 19 20 \
              21 22 X Y"
    fi

}


run_mutect2() {
    local interval=$1
    local output_vcf="${temp_dir}/mutect2_${interval}.vcf"
    local f1r2_output="${temp_dir}/f1r2_${interval}.tar.gz"
    cmd="/opt/GATK/gatk-4.6.1.0/gatk --java-options \"-Xmx${mutect2_memory}g\" \
      Mutect2 \
      -R ${reference} \
      $(echo ${input_args}) \
      -normal ${normal_sample_name} \
      --germline-resource ${gnomad_af} \
      --panel-of-normals ${pon} \
      -L ${interval} \
      -O ${output_vcf} \
      --f1r2-tar-gz ${f1r2_output}"
    eval "$cmd"
    echo "${output_vcf}"
}





# ==============================
# GERMLINE VARIANTS
# ==============================


# find germline variants in normal sample and keep only those with a minimum depth of germline_snp_min_depth_normal
bcftools mpileup -Ou -f $reference $bam_normal_path | \
bcftools call -mv -Ou | \
bcftools view -v snps | \
bcftools filter -e "INFO/DP[*] < $germline_snp_min_depth_normal" -Ob -o $germline_variants

# throw error if $germline_variants is empty or does not exist:
if [ ! -s $germline_variants ]; then
    echo "Error: $germline_variants is empty or does not exist."
    exit 1
fi


bcftools query --format '%CHROM\t%POS\n' $germline_variants | bgzip -c > $germline_variants_pos
bgzip --reindex $germline_variants_pos
echo "germline variants pos done"

# create a text file with bam paths:
printf "%s\n" "${SAMPLES[@]}" | xargs -I{} echo "${bams_path}/{}.bam" > "${mpileup}/bam_list.txt"

bcftools mpileup --regions-file "${germline_variants_pos}" \
                 --fasta-ref "${reference}" \
                 --annotate FORMAT/AD,FORMAT/DP \
                 --bam-list "${mpileup}/bam_list.txt" \
                 -Ou | \
bcftools call -m -Ou | \
bcftools view -v snps | \
bcftools filter -e "FORMAT/DP[*] < $germline_snp_min_depth_tumour" -Ob -o "${mpileup}/all_samples.bcf"

# Index the output file
bcftools index "${mpileup}/all_samples.bcf"

# throw error if "${mpileup}/all_samples.bcf" is empty or does not exist:
if [ ! -s "${mpileup}/all_samples.bcf" ]; then
    echo "Error: "${mpileup}/all_samples.bcf" is empty or does not exist."
    exit 1
fi

# Split into individual sample files and extract required information
for sample in "${SAMPLES[@]}"; do
    echo "Processing $sample"
    bcftools view -s "$sample" "${mpileup}/all_samples.bcf" -Ob -o "${mpileup}/${sample}.bcf"
    bcftools query --format "%CHROM\t%POS\t[%AD{0}]\t[%AD{1}]\n" "${mpileup}/${sample}.bcf" | \
    gzip -c > "${mpileup}/${sample}.mpileup.gz"
done


# ==============================
# ASCAT and REFPHASE
# ==============================


conda activate /opt/miniforge3/envs/refphase # or just conda activate refphase if running without singularity

echo "running ascat preprocessing"
Rscript ${base_dir}/scripts/ascat_preproc.R \
    --input_directory "${mpileup}" \
    --output_directory "${ascat_input_dir}" \
    --normal_sample_name "${normal_sample_name}" \
    --samples "${SAMPLES[*]}"

echo "running ascat"

# run ASCAT followed by refphase:
Rscript ${base_dir}/scripts/ascat_refphase_run.R  \
    --input_directory $ascat_dir \
    --output_directory $refphase_dir_output \
    --normal_sample_name "${normal_sample_name}" \
    --samples "${SAMPLES[*]}" \
    --case_id $case_id 

conda deactivate


# ==============================
# SOMATIC VARIANTS
# ==============================


# activate gatk conda:
conda activate /opt/miniforge3/envs/gatk # or just conda activate gatk if running without singularity

# create dictionary:
/opt/GATK/gatk-4.6.1.0/gatk CreateSequenceDictionary \
    -R ${reference}

# index germline variants:
/opt/GATK/gatk-4.6.1.0/gatk IndexFeatureFile \
    -I ${gnomad_af}

# index panel of normals:
/opt/GATK/gatk-4.6.1.0/gatk IndexFeatureFile \
    -I ${pon}

# create gnomad_common_biallelic file:
/opt/GATK/gatk-4.6.1.0/gatk SelectVariants \
    -V ${gnomad_af} \
    --restrict-alleles-to BIALLELIC \
    --select-type-to-include SNP \
    --max-nocall-fraction 0.1 \
    --select "AF > 0.05" \
    -O "${gnomad_common_biallelic}"

# create temporary directory
temp_dir="${somatic_dir}/temp"
mkdir -p $temp_dir

# read bams files from the file:
input_args=""
while read bam_path; do
  input_args+=" -I ${bam_path}"
done < "${mpileup}/bam_list.txt"

chromosomes=($(get_chromosomes_style ${reference}.fai))

export -f run_mutect2
export reference input_args normal_sample_name gnomad_af pon temp_dir mutect2_memory

interval_files=""
for chrom in "${chromosomes[@]}"; do
    interval_files="${interval_files} ${temp_dir}/mutect2_${chrom}.vcf"
done

total_jobs=${#chromosomes[@]}
current_job=0

while [ $current_job -lt $total_jobs ]; do
    # Start a new batch of jobs
    active_jobs=0
    
    while [ $active_jobs -lt $mutect2_cpus ] && [ $current_job -lt $total_jobs ]; do
        chrom=${chromosomes[$current_job]}
        run_mutect2 "$chrom" &
        
        current_job=$((current_job + 1))
        active_jobs=$((active_jobs + 1))
    done
    
    # Wait for this batch to complete before starting the next batch
    wait
done

output_vcf="${somatic_dir}/unfiltered.vcf"

# Merge VCF files
echo "Merging VCF files..."
/opt/GATK/gatk-4.6.1.0/gatk --java-options "-Xmx${mutect2_memory}g" \
  MergeVcfs \
  $(for vcf in ${temp_dir}/mutect2_*.vcf; do echo "-I $vcf"; done) \
  -O "${output_vcf}"

# Merge f1r2 files
echo "Merging f1r2 tar.gz files..."
/opt/GATK/gatk-4.6.1.0/gatk --java-options "-Xmx${mutect2_memory}g" \
  LearnReadOrientationModel \
  $(for f1r2 in ${temp_dir}/f1r2_*.tar.gz; do echo "-I $f1r2"; done) \
  -O "${somatic_dir}/read-orientation-model.tar.gz"

# Merge stats files
echo "Merging stats  files..."
stats_inputs=""
for chrom in "${chromosomes[@]}"; do
    stats_inputs="${stats_inputs} --stats ${temp_dir}/mutect2_${chrom}.vcf.stats"
done

/opt/GATK/gatk-4.6.1.0/gatk --java-options "-Xmx${mutect2_memory}g" \
  MergeMutectStats \
  ${stats_inputs} \
  -O "${output_vcf}.stats"


# Clean up temporary files if successful
if [ -f "${output_vcf}" ]; then
  echo "Cleaning up temporary files..."
  rm -rf "${temp_dir}"
fi

# get summaries for each sample:
for sample in "${SAMPLES[@]}"; do
    /opt/GATK/gatk-4.6.1.0/gatk GetPileupSummaries \
        -I "${bams_path}/${sample}.bam" \
        -V ${gnomad_common_biallelic} \
        -L ${gnomad_common_biallelic} \
        -O "${somatic_dir}/${sample}.pileups.table"

    # Calculate contamination
    # skip normal sample:
    if [ "$sample" == "$normal_sample_name" ]; then
        continue
    fi
    /opt/GATK/gatk-4.6.1.0/gatk CalculateContamination \
        -I "${somatic_dir}/${sample}.pileups.table" \
        -matched "${somatic_dir}/${normal_sample_name}.pileups.table" \
        -O "${somatic_dir}/${sample}.contamination.table"
done

input_args_cont=""
for sample in "${SAMPLES[@]}"; do
    # skip normal sample:
    if [ "$sample" == "$normal_sample_name" ]; then
        continue
    fi
    path_to_cont_table="${somatic_dir}/${sample}.contamination.table"
    input_args_cont+=" --contamination-table  ${path_to_cont_table}"
done

/opt/GATK/gatk-4.6.1.0/gatk FilterMutectCalls \
    -V "${somatic_dir}/_unfiltered.vcf" \
    -R ${reference} \
    $(echo ${input_args_cont}) \
    --ob-priors "${somatic_dir}/_read-orientation-model.tar.gz" \
    -O "${somatic_dir}/filtered.vcf"


# apply filters and create allele counts tables for each tumour sample:
# filters:
python ${base_dir}/scripts/vcf_to_allele_calls.py \
    --input_vcf "${somatic_dir}/filtered.vcf" \
    --output_file "${somatic_dir}/allele_counts.tsv" \
    --normal_sample $normal_sample_name \
    --vaf_threshold $variant_vaf_threshold \
    --germline_vaf_threshold $germline_vaf_threshold \
    --germline_variant_read_count_threshold $germline_variant_read_count_threshold


# combine allele_counts with refphase output to create conipher input
# this script aplies filetering based on total depth across samples and variant count per sample - modify it to suit your needs

conda deactivate
# use alpaca env as it contains the required packages:
conda activate /opt/miniforge3/envs/alpaca/
python ${base_dir}/scripts/create_conipher_input.py \
  -c "${refphase_dir_output}/${case_id}-refphase-segmentation.tsv" \
  -a "${somatic_dir}/allele_counts.tsv" \
  -p "${refphase_dir_output}/${case_id}-refphase-sample-data-updated.tsv" \
  -o "${conipher_dir_input}" \
  -t "${case_id}" \
  -b "${bams_path}" \
  --depth-thr $depth_threrehold \
  --var-threshold $varcount_threshold

conda deactivate


# ==============================
# CONIPHER
# ==============================


# activate CONIPHER env:

conda activate /opt/miniforge3/envs/conipher
echo "running conipher"
input_tsv_loc="${conipher_dir_input}/clustering_input.tsv"
Rscript ${base_dir}/scripts/run_conipher.R $case_id $prefix "$input_tsv_loc" "$conipher_dir_output"
conda deactivate


# ==============================
# IMPORTANT
# ==============================
# Before running ALPACA, make sure the solutions provided by CONIPHER and Refphase
# are satisfactory. For example, check if purity and ploidy estimates are reasonable
# and ensure that the clonal cluster in CONIPHER output tree has CCF values of around 100 in each sample
# CCF values can be found in conipher/output/Trees/pytree_and_bar.pdf plot on left hand side.


# ==============================
# ALPACA
# ==============================


conda activate /opt/miniforge3/envs/alpaca/
# convert refphase and conipher outputs to alpaca input
CONIPHER_tree_object="${conipher_dir_output}/Trees/${case_id}.tree.RDS"
refphase_rData="${refphase_dir_output}/${case_id}-refphase-results.RData"
input_conversion \
 --tumour_id $case_id \
 --refphase_rData $refphase_rData \
 --CONIPHER_tree_object $CONIPHER_tree_object \
 --output_dir $alpaca_input_dir


# run alpaca:
alpaca \
    --input_tumour_directory "${alpaca_input_dir}" \
    --output_directory "${alpaca_output_dir}" > "${alpaca_output_dir}/alpaca.log"

conda deactivate

