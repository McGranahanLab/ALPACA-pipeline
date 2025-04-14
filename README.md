# ALPACA Pipeline

This pipeline allows the user to run ALPACA starting from BAM files. It requires multiple intermediate steps, including running [CONIPHER](https://github.com/McGranahanLab/CONIPHER/blob/main/README.md) to create phylogenetic tree and [Refphase](https://bitbucket.org/schwarzlab/refphase/src/master/) to obtain fractional copy numbers from a multi-sample input.

## Requirements

ALPACA is designed to work with a multi-sample data, i.e. multiple samples of the same tumour acquired via 2 or more separate samplings. Therefore, it requires at least two separate tumour BAM files. This tutorial below assumes that on top of at least two BAM files, the matched normal BAM is also available.

### Files

Additionally, to call the  germline and somatic variants, the user should have:

- Reference genome
- Panel of normals
- Germline variants

These files can be acquired from [GATK resource bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle) and must match genome build used to align the samples (i.e. the same build which was used for creation of the BAM files.)

### Software

This tutorial has been tested using Singularity (v. 3.11.3) image alongside three separate conda/mamba virtual environments. If you want to run the code without using Singularity, ensure that your base environment has all the necessary libraries (see Singularity definition file `singularity/pipeline.def`) as well as `samtools` `tabix` and `bcftools` installed.

Importantly, ALPACA requires Gurobi solver to work. The user can obtain free academic license at https://www.gurobi.com/academia/academic-program-and-licenses.

After obtaining the license, log in to your Gurobi account, go to "Review your current licenses" and press download icon. You will see detailed instructions to activate the license. Please take note of where the gurobi.lic is saved.

### Hardware

The tutorial is implemented as a simple sequential bash script. Running Mutect2 is the longest component, but the process can be sped up by using multiple CPUs. Since Mutect2 script in this pipeline operates on each chromosome separately, we recommend to use as many CPUs as chromosomes present in the sample. Number of CPUs is controlled by the `mutect2_cpus` parameter.

## Running the pipeline

All the necessary command to run the pipeline can be found in the `pipeline.sh` script. However, instead of running the entire script we advise the user to either create their own pipeline using a pipeline management software such as Nextflow or Snakemake or alternatively execute each block sequentially from a stable environment (e.g. using TMUX: `tmux new -s pipeline`).

## Running via Singularity

In order to guarantee the reproducibility, we recommend using Singularity container. To build the image, please use Singularity 3.11.3 - using other versions might cause issues during the build. If you have sudo privileges on the host machine, use:

`sudo singularity build singularity/pipeline.sif singularity/pipeline.def`

If you don't have sudo available, you can use remote build option (`singularity build --remote pipeline.sif pipeline.def`) see: https://cloud.sylabs.io/builder.

On Windows of Mac singularity must be run via a [Virtual Machine](https://docs.sylabs.io/guides/3.0/user-guide/installation.html)

Since Gurobi license is needed to run ALPACA, find the path to your Gurobi license, it should be present under `$GRB_LICENSE_FILE`. If your are using free academic, store the path to the license file under `GRB_LICENSE_FILE variable`.

Once the image is built, activate the singularity shell with:

```bash
image_path=path/to/your/image.sif
singularity shell -B path/to/host/filesystem -B path/to/gurobi.lic $image_path
```

E.g.

```bash
image_path=singularity/pipeline.sif
singularity shell -B /Users/GitHub/ALPACA-pipeline -B $GRB_LICENSE_FILE $image_path
```

## Running without Singularity

You can also run the pipeline without the Singularity image. In this case, ensure that all the necessary libraries as well as `samtools` `tabix` and `bcftools` installed in your base environment. You will also need conda or mamba to create three virtual environments.

### GATK environment

For basic bioinformatic processing tools, create a `gatk` conda environment with:

```bash
mkdir GATK && cd GATK
wget https://github.com/broadinstitute/gatk/releases/download/4.6.1.0/gatk-4.6.1.0.zip
unzip gatk-4.6.1.0.zip
cd gatk-4.6.1.0
conda env create -n gatk -f gatkcondaenv.yml
cd ../..
```

### Refphase environment

In order to run ASCAT and Refphase create a `refphase` conda environment with:

```bash
conda create -n refphase -c bioconda -c conda-forge r-base r-tidyverse r-gtools r-jsonlite r-optparse r-devtools bioconductor-genomicranges -y
conda run -n refphase R -e "devtools::install_bitbucket('schwarzlab/refphase')"
conda run -n refphase R -e "devtools::install_bitbucket('schwarzlab/ascat_v3_fork/ASCAT')"
```

### CONIPHER environment

In order to run CONIPHER create a `conipher` conda environment with:

```bash
conda create -n conipher -c conda-forge -c bioconda conipher -y
```

### ALPACA environment

Finally, create `alpaca` environment with:

```bash
# download repositories
model_repo_url='https://github.com/McGranahanLab/ALPACA-model.git'
git clone $model_repo_url
conda env create -f ALPACA-model/environment.yml -n alpaca -y
conda run -n alpaca pip install repo/ALPACA-model/dist/*.whl
conda run -n alpaca pip install pysam
```

## The pipeline

See pipeline.sh for the full pipeline example. Below we present the same code but with additional comments.

### Setup

We recommend running the pipeline via TMUX and Singularity shell.
Start by creating a new TMUX seession:

```bash
tmux new -s pipeline
```

Next, define the paths to Singularity image, Guorbi license and working directory and open Singularity shell

```bash
image_path="singularity/pipeline.sif"
license_path=$GRB_LICENSE_FILE
work_dir=$(pwd)
singularity shell -B $work_dir -B $license_path $image_path 
```

Finally, initialize conda with:

```bash
source /opt/miniforge3/etc/profile.d/conda.sh
```

### Assets

Start by defining file path to all the required assets. In this tutorial we use "_assets" directory which holds reference files, gnomad resource, panel of normals and BAM files:

```bash
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
```

Example _assets structure:

```bash
_assets/
├── bams
│   ├── LTX000_BS_GL--7c096035hwyv.bam
│   ├── LTX000_BS_GL--7c096035hwyv.bam.bai
│   ├── LTX000_SU_T1-R1--d864fcbbhwyv.bam
│   ├── LTX000_SU_T1-R1--d864fcbbhwyv.bam.bai
│   ├── LTX000_SU_T1-R2--c384bdbchwyv.bam
│   ├── LTX000_SU_T1-R2--c384bdbchwyv.bam.bai
├── gnomad
│   ├── af-only-gnomad.raw.sites.b37.vcf.gz
│   ├── af-only-gnomad.raw.sites.b37.vcf.gz.tbi
│   ├── gnomad.common_biallelic.vcf.gz
│   └── gnomad.common_biallelic.vcf.gz.tbi
├── pon
│   ├── gatk4_mutect2_4136_pon_hg19_liftover.vcf.gz
│   └── gatk4_mutect2_4136_pon_hg19_liftover.vcf.gz.tbi
└── reference
    ├── ucsc.hg19.dict
    ├── ucsc.hg19.fasta
    └── ucsc.hg19.fasta.fai
```

### Parameters

Next, we define parameters used in the pipeline:

```bash

# ==============================
# Parameters
# ==============================


# Specify the normal sample name - we assume there is only one normal sample.
# We assume that all the BAM files are named according to the following pattern:
# <sample_name>.bam
normal_sample_name="LTX000-GL--164c6d2dc495"
bam_normal_path="${bams_path}/${normal_sample_name}.bam"
# to get create an array of sample names:
SAMPLES=()
for bam_file in "$bams_path"/*.bam; do
    SAMPLE=$(basename "$bam_file" .bam)
    SAMPLES+=("$SAMPLE")
done
echo "Sample names: ${SAMPLES[@]}"

# CONIPHER requires case identifier composed of a 'prefix' and a case number:
case_id=LTX0000
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
germline_vaf_threshold=0.01 # Reject variant if germline VAF is above X
variant_vaf_threshold=0.05  # Reject variant if VAF is below X
```

### Paths

Define all the necessary paths and initialize empty directories:

```bash
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
```

### Functions

We define these functions to facilitate processing:

```bash
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
```

### Input preparation

Make sure that your BAMs and reference files are properly indexed. If not, you can create the index by running `samtools index` and `samtools faidx`.
E.g.:

```bash
# ==============================
# BAM indexing
# ==============================


for bam_file in "$bams_path"/*.bam; do
    samtools index $bam_file
done

samtools faidx $reference
```

### Germline variants

Next, we call germline variants using bcftools:

```bash
# ==============================
# GERMLINE VARIANTS
# ==============================


# find germline variants in normal sample and keep only those with a minimum depth of germline_snp_min_depth_normal
bcftools mpileup -Ou -f $reference $bam_normal_path | \
bcftools call -mv -Ou | \
bcftools view -v snps | \
bcftools filter -e "FORMAT/DP[*]" < $germline_snp_min_depth_normal -Ob -o $germline_variants

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

# Split into individual sample files and extract required information
for sample in "${SAMPLES[@]}"; do
    echo "Processing $sample"
    bcftools view -s "$sample" "${mpileup}/all_samples.bcf" -Ob -o "${mpileup}/${sample}.bcf"
    bcftools query --format "%CHROM\t%POS\t[%AD{0}]\t[%AD{1}]\n" "${mpileup}/${sample}.bcf" | \
    gzip -c > "${mpileup}/${sample}.mpileup.gz"
done


```

### ASCAT and Refphase

We start by activaing the refphase conda and then we run a preprocessing script, followed by ASCAT and Refphase.

```bash
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
```

Troubleshooitng:
If ASCAT fails to find optimal purity and ploidy, you migh need to adujust SNP filtering thresholds (germline_snp_min_depth_normal).

### Somatic variants

In order to run CONIPHER we need to identify somatic mutations - this step takes a long time if done on a single CPU. The script below can accomodate both single and multiple cpus. For maximum efficient, use separate cpu for each chromosome

Before running this commands, remember to deactivate previous conda environment and activate `gatk`

```bash
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
  # rm -rf "${temp_dir}"
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
    -V "${somatic_dir}/unfiltered.vcf" \
    -R ${reference} \
    $(echo ${input_args_cont}) \
    --ob-priors "${somatic_dir}/read-orientation-model.tar.gz" \
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
```

### CONIPHER

In this setup we transfrom somatic variants into CONIPHER input table and run the CONIPHER tool to obtain phylogenetic trees and clone proportions. Before proceeding, remember to deactivate previous environment and activate `conipher`.

```bash
# ==============================
# CONIPHER
# ==============================


# activate CONIPHER env:

conda activate /opt/miniforge3/envs/conipher
echo "running conipher"
input_tsv_loc="${conipher_dir_input}/clustering_input.tsv"
Rscript ${base_dir}/scripts/run_conipher.R $case_id $prefix "$input_tsv_loc" "$conipher_dir_output"
conda deactivate
```

### ALPACA

Finally, we will use Refphase and CONIPHER outputs to run Alpaca. Before proceeding, please remember to deactivate `conipher` environment and activate `alpaca`. Before running ALPACA, make sure the solutions provided by CONIPHER and Refphase are satisfactory. For example, check if purity and ploidy estimates are reasonable and ensure that the clonal cluster in CONIPHER output tree has CCF values of around 100 in each sample CCF values can be found in `conipher/output/Trees/pytree_and_bar.pdf` plot on left hand side.

```bash
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
```
