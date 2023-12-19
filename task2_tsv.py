import subprocess
import vcf

def annotate_vcf(input_vcf, output_annotated_vcf):
    # Annotate VCF using VEP
    vep_command = f"vep --cache --offline --dir_cache /home/roshny/AbiomixTestFiles " \
                  f"--input_file {input_vcf} --output_file {output_annotated_vcf} " \
                  "--format vcf --vcf --force_overwrite"
    subprocess.run(vep_command, shell=True)

def is_exonic(vep_annotations):
    for annotation_entry in vep_annotations:
        for annotation in annotation_entry.split(','):
            consequences = annotation.split('&')
            if any('exon' in c for c in consequences):
                return True
    return False

def get_variant_class(record):
    if len(record.REF) < len(record.ALT[0]):
        return 'INS'
    elif len(record.REF) > len(record.ALT[0]):
        return 'DEL'
    else:
        return 'SNP'

def parse_vcf(vcf_file, output_tsv):
    #Extract Relevant Fields and Generate TSV
    vcf_reader = vcf.Reader(filename=vcf_file)
    
    with open(output_tsv, 'w') as output_file:
        # Print headers
        output_file.write("CHROM\tPOS\tREF\tALT\tGT\tDP\tAF\tISExonic\tVariantClass\n")

        # Iterate through VCF records
        for record in vcf_reader:
            is_exonic_variant = is_exonic(record.INFO.get('CSQ', []))
            variant_class = get_variant_class(record)

            # Write data to the output file
            output_file.write(f"{record.CHROM}\t{record.POS}\t{record.REF}\t{record.ALT[0]}\t"
                              f"{record.samples[0]['GT']}\t{record.samples[0]['DP']}\t{record.samples[0]['AF']}\t"
                              f"{is_exonic_variant}\t{variant_class}\n")

def main():
    input_vcf = 'testsample.vcf.gz'
    output_annotated_vcf = 'annotated_testsample'
    output_tsv = 'task2_tsv_output.tsv'
    
    # Annotate VCF
    annotate_vcf(input_vcf, output_annotated_vcf)
    
    #Parse VCF and Generate TSV
    parse_vcf(output_annotated_vcf, output_tsv)

if __name__ == "__main__":
    main()

