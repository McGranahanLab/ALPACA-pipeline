import csv
import argparse

def parse_vcf_to_allele_counts(vcf_file, output_file, normal_sample, vaf_threshold, germline_vaf_threshold, germline_variant_read_count_threshold):

    with open(vcf_file, "r") as infile, open(output_file, "w", newline="") as outfile:
        writer = csv.writer(outfile, delimiter="\t")
        writer.writerow([
            "SAMPLE",
            "CHROM",
            "POS",
            "REF",
            "ALT",
            "REF_COUNT",
            "ALT_COUNT",
            "TOTAL_COUNT",
            "ALT_FREQ"
        ])
        for line in infile:
             # Skip header lines
            if line.startswith("##"):
                continue
            # Extract column names from the header line
            if line.startswith("#CHROM"):
                header = line.strip().split("\t")
                continue
            parse_line(line, writer, header,normal_sample, vaf_threshold, germline_vaf_threshold, germline_variant_read_count_threshold)


def parse_line(line, writer, header, normal_sample, vaf_threshold, germline_vaf_threshold, germline_variant_read_count_threshold):
    format_index = header.index("FORMAT")
    # Process data lines
    fields = line.strip().split("\t")
    chrom = fields[0]
    pos = fields[1]
    ref = fields[3]
    alt = fields[4]
    filter = fields[6]
    if filter != "PASS":
        return
    # keep only single base mutations, with only one variant
    if (len(ref) > 1) or (len(alt) > 1):
        return
    # Get FORMAT field and tumor sample data
    normal_index = header.index(normal_sample)
    normal_data = fields[normal_index].split(":")
    format_fields = fields[format_index].split(":")
    # Create a dictionary to map format fields to tumor data
    normal_dict = dict(zip(format_fields, normal_data))
    # Extract allele depth (AD), reference count is first, alt count is second
    germline_ref_count, germline_alt_count = map(
        int, normal_dict["AD"].split(",")
    )
    germline_total_count = germline_ref_count + germline_alt_count
    # Extract allele frequency (AF)
    germline_alt_freq = (
        float(normal_dict["AF"])
        if "AF" in normal_dict
        else (
            germline_alt_count / germline_total_count
            if germline_total_count > 0
            else 0
        )
    )

    if germline_alt_freq > germline_vaf_threshold:
        return
    if germline_alt_count >= germline_variant_read_count_threshold:
        return
    for sample_name in header[9:]:
        if sample_name == normal_sample:
            continue
        sample_index = header.index(sample_name)
        tumor_data = fields[sample_index].split(":")
        format_fields = fields[format_index].split(":")
        # Create a dictionary to map format fields to tumor data
        tumor_dict = dict(zip(format_fields, tumor_data))
        # Extract allele depth (AD), reference count is first, alt count is second
        ref_count, alt_count = map(int, tumor_dict["AD"].split(","))
        total_count = ref_count + alt_count
        # Extract allele frequency (AF)
        alt_freq = float(tumor_dict["AF"])
        if alt_freq < vaf_threshold:
            continue
        # Write to output file
        new_line = [sample_name, chrom, pos, ref, alt, ref_count, alt_count, total_count, f"{alt_freq:.4f}"]
        writer.writerow(new_line)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_vcf", help="Input VCF file")
    parser.add_argument("--output_file", help="Output file")
    parser.add_argument("--normal_sample", help="Normal sample name")
    parser.add_argument("--vaf_threshold", help="VAF threshold", type=float, default=0.05)
    parser.add_argument("--germline_vaf_threshold", help="Germline VAF threshold", type=float, default=0.01)
    parser.add_argument("--germline_variant_read_count_threshold", help="Germline variant read count threshold", type=int, default=5)
    args = parser.parse_args()
    parse_vcf_to_allele_counts(args.input_vcf, args.output_file, args.normal_sample, args.vaf_threshold, args.germline_vaf_threshold, args.germline_variant_read_count_threshold)