import argparse
import pandas as pd
import pysam
from typing import Tuple

def parse_arguments():
    parser = argparse.ArgumentParser(description="Ingest input files for copy number and allele counts.")
    
    parser.add_argument(
        "-c", "--copy-number", 
        nargs='+',
        required=True,
        help="Paths to TSV files containing copy number data."
    )
    
    parser.add_argument(
        "-a", "--allele-counts", 
        required=True,
        help="Paths to TSV file containing allele count data for mutations."
    )
    
    parser.add_argument(
        "-p", "--purity-ploidy", 
        required=True,
        help="Paths to TSV files containing purity and ploidy per sample."
    )
    
    parser.add_argument(
        "-o", "--output-directory", 
        required=True,
        help="Path to output directory."
    )
    
    parser.add_argument(
        "-t", "--case-id", 
        required=True,
        help="Tumour identifier."
    )
    
    parser.add_argument(
        "-b", "--bam-files", 
        required=True,
        help="Path to directory containing BAM files."
    )
    
    parser.add_argument(
        "-v", "--var-threshold", 
        default=10,
        type=int,
        help="Threshold for minimum variant count of mutation in sample."
    )
    
    parser.add_argument(
        "-d", "--depth-thr", 
        default=30,
        type=int,
        help="Threshold for total depth of mutation in all samples."
    )
    
    parser.add_argument(
        "-e", "--test", 
        default=False,
        action='store_true',
        help="If True, random mutations will be removed from each sample. Use this to speedup the process for troubleshooting."
    )
    
    
    
    args = parser.parse_args()
    return args


def force_call_mutations(allele_count_df:pd.DataFrame, bam_files:str) -> pd.DataFrame:
    # for each sample, find sites not mutated in this sample, but mutated in other samples:
    all_mutations = set(allele_count_df['MUT_ID'])
    force_called_mutations = {'SAMPLE':[],'CHR':[],'POS':[],'REF':[],'ALT':[],'REF_COUNT':[],'VAR_COUNT':[],'DEPTH':[]}
    for sample_name,df in allele_count_df.groupby('SAMPLE'):
        mutations_in_sample = set(df['MUT_ID'])
        mutations_to_force_call = all_mutations - mutations_in_sample
        sample_bam_path = f'{bam_files}/{sample_name}.bam'
        bamfile = pysam.AlignmentFile(sample_bam_path, "rb")
        contigs_names = [contig['SN'] for contig in bamfile.header['SQ']]
        for mut in mutations_to_force_call:
            mut_chr, mut_pos, mut_ref, mut_alt = mut.split('_')
            mut_pos = int(mut_pos)
            if mut_chr not in contigs_names:
                mut_chr = 'chr' + mut_chr
                assert mut_chr in contigs_names
            reference_count, variant_count = get_variant_and_reference_counts(bamfile, mut_chr, mut_pos, mut_ref, mut_alt)
            force_called_mutations['SAMPLE'].append(sample_name)
            force_called_mutations['CHR'].append(mut_chr.replace('chr',''))
            force_called_mutations['POS'].append(mut_pos)
            force_called_mutations['REF'].append(mut_ref)
            force_called_mutations['ALT'].append(mut_alt)
            force_called_mutations['REF_COUNT'].append(reference_count)
            force_called_mutations['VAR_COUNT'].append(variant_count)
            force_called_mutations['DEPTH'].append(reference_count + variant_count)  
        bamfile.close()
    force_called_mutations_df = pd.DataFrame(force_called_mutations)
    force_called_mutations_df['CHR'] = force_called_mutations_df['CHR'].astype(str)
    force_called_mutations_df['MUT_ID'] = force_called_mutations_df['CHR'].astype(str) + '_' + force_called_mutations_df['POS'].astype(str) + '_' + force_called_mutations_df['REF'] + '_' + force_called_mutations_df['ALT']
    return force_called_mutations_df


def get_variant_and_reference_counts(bamfile: pysam.AlignmentFile, mut_chr: str, mut_pos: int, mut_ref: str, mut_alt: str) -> Tuple[int, int]:
    variant_count = 0
    reference_count = 0
    for pileupcolumn in bamfile.pileup(mut_chr, mut_pos - 1, mut_pos):
        if pileupcolumn.pos == mut_pos - 1:
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del and not pileupread.is_refskip:
                    base = pileupread.alignment.query_sequence[pileupread.query_position]
                    if base == mut_ref:
                        reference_count += 1
                    elif base == mut_alt:
                        variant_count += 1
    return reference_count, variant_count


def remove_random_mutations(allele_count_df: pd.DataFrame,frac: float=0.01) -> pd.DataFrame:
    dfs = []
    for sample_name,df in allele_count_df.groupby('SAMPLE'):
        random_mutations = df.sample(frac=frac)[['MUT_ID']]
        # remove random mutations from clustering_input df:
        df = df[~df['MUT_ID'].isin(random_mutations['MUT_ID'])]
        dfs.append(df)
    return pd.concat(dfs)



def main():
    args = parse_arguments()
    copy_number_files = args.copy_number
    allele_counts = args.allele_counts
    purity_ploidy = args.purity_ploidy
    output_directory = args.output_directory
    case_id = args.case_id
    bam_files = args.bam_files
    var_threshold = args.var_threshold
    depth_thr = args.depth_thr
    test = args.test
    if test:
        print('-------------------------------------------------------------')
        print('Test mode: random mutations will be removed from each sample.')
        print('-------------------------------------------------------------')
    print('Creating conipher input files...')
    # conipher expects this columns in this order
    clustering_columns=[
    'CASE_ID',
    'SAMPLE',
    'CHR',
    'POS',
    'REF',
    'ALT',
    'REF_COUNT',
    'VAR_COUNT',
    'DEPTH',
    'COPY_NUMBER_A',
    'COPY_NUMBER_B',
    'ACF',
    'PLOIDY',
    'MUT_TYPE']

    copy_number_df = pd.concat([pd.read_csv(f, sep='\t') for f in copy_number_files])
    copy_number_df.rename(columns={'chrom':'CHR',
                                'sample_id':'SAMPLE',
                                'cn_a':'COPY_NUMBER_A',
                                'cn_b':'COPY_NUMBER_B',
                                }, inplace=True)
    # ensure type consistency:
    copy_number_df['CHR'] = copy_number_df['CHR'].astype(str)
    purity_ploidy_df = pd.read_csv(purity_ploidy, sep='\t')
    purity_ploidy_df.rename(columns={'sample_id':'SAMPLE',
                                    'purity':'ACF',
                                    'ploidy':'PLOIDY',
                                }, inplace=True)
    print('Reading allele count files...')

    allele_count_df = pd.read_table(allele_counts)


    allele_count_df.rename(columns={'CHROM':'CHR', 'TOTAL_COUNT':'DEPTH','ALT_COUNT':'VAR_COUNT'}, inplace=True)
    allele_count_df['CHR'] = allele_count_df['CHR'].str.replace('chr','')
    allele_count_df = allele_count_df[allele_count_df['CHR'].isin(copy_number_df['CHR'].astype(str).unique())]
    allele_count_df['CHR'] = allele_count_df['CHR'].astype(str)
    allele_count_df['MUT_ID'] = allele_count_df['CHR'].astype(str) + '_' + allele_count_df['POS'].astype(str) + '_' + allele_count_df['REF'] + '_' + allele_count_df['ALT']

    print('Removing low confidence mutations and forcing call mutations...')
    # remove low confidence mutations:
    count_agg = allele_count_df.groupby('MUT_ID').agg({'VAR_COUNT':'max','DEPTH':'min'}).reset_index().rename(columns={'VAR_COUNT':'max_VAR_COUNT','DEPTH':'min_DEPTH'})
    loci_COV_threshold_pass = count_agg[count_agg['min_DEPTH'] >= depth_thr]
    loci_MUTCOV_threshold_pass = count_agg[count_agg['max_VAR_COUNT'] >= var_threshold]
    loci_pass = set(loci_COV_threshold_pass['MUT_ID']).intersection(set(loci_MUTCOV_threshold_pass['MUT_ID']))
    allele_count_df = allele_count_df[allele_count_df['MUT_ID'].isin(loci_pass)]

    #allele_count_df = allele_count_df[allele_count_df['VAR_COUNT']<var_threshold]
    #total_depth_per_mutation = allele_count_df.groupby('MUT_ID')['DEPTH'].sum()
    # keep only samples with at least depth_thr reads in all the samples:
    # allele_count_df = allele_count_df[allele_count_df['MUT_ID'].isin(total_depth_per_mutation[total_depth_per_mutation>=depth_thr].index)]
    

    if test:
        allele_count_df = remove_random_mutations(allele_count_df)
    

    # force call mutations:
    force_called_muts_df = force_call_mutations(allele_count_df, bam_files)
    allele_count_df = pd.concat([allele_count_df, force_called_muts_df])
    allele_count_df = allele_count_df.drop_duplicates()
    allele_count_df['CASE_ID'] = case_id
    allele_count_df['MUT_TYPE'] = 'SNV'
    assert (allele_count_df['MUT_ID'].value_counts()==allele_count_df['SAMPLE'].nunique()).all()
    copy_number_df_pp = copy_number_df.merge(purity_ploidy_df, on='SAMPLE', how='left')
    assert copy_number_df_pp['CHR'].dtype == allele_count_df['CHR'].dtype
    copy_number_df_pp_mut = copy_number_df_pp.merge(allele_count_df, on=['CHR','SAMPLE'], how='inner')
    # ensure that POS, start and end are integers:
    copy_number_df_pp_mut['POS'] = copy_number_df_pp_mut['POS'].astype(int)
    copy_number_df_pp_mut['start'] = copy_number_df_pp_mut['start'].astype(int)
    copy_number_df_pp_mut['end'] = copy_number_df_pp_mut['end'].astype(int)
    copy_number_df_pp_mut = copy_number_df_pp_mut[(copy_number_df_pp_mut['POS']>=copy_number_df_pp_mut['start'])& (copy_number_df_pp_mut['POS']<=copy_number_df_pp_mut['end'])]
    print('Saving output files...')
    # save main and test dfs:
    if test:
        random_mutations = copy_number_df_pp_mut.sample(frac=0.05)['MUT_ID'].unique()
        copy_number_df_pp_mut_test = copy_number_df_pp_mut[copy_number_df_pp_mut['MUT_ID'].isin(random_mutations)]
        copy_number_df_pp_mut_test = copy_number_df_pp_mut_test[clustering_columns].drop_duplicates()
        copy_number_df_pp_mut_test.to_csv(f'{output_directory}/clustering_input_test.tsv', sep='\t', index=False)
    clustering_input = copy_number_df_pp_mut[clustering_columns].drop_duplicates()
    clustering_input.to_csv(f'{output_directory}/clustering_input.tsv', sep='\t', index=False)
    unique_mutations = clustering_input[['CHR','POS']].drop_duplicates().shape[0]
    print(f'{unique_mutations} unique mutations in CONIPHER input.')
    print('Done!')
    
if __name__ == "__main__":
    main()