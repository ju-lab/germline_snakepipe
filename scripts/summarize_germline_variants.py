'''Script to summarize germline variant calls from Freebayes
Originally modified from variant_fishing_no_inheritance_model.py
Input file is VCF file that has been left-normalized and decomposed with vt and annotated with VEP. 
2018.06.21 CJY

Different from PDP analysis in that population frequency filter has been loosened to 0.1 from 0.01 
since no variant was found with 0.01 cutoff.

Modified 2018.07.01. Fixing variant_type and vaf calculation errors. 
Modiifed 2018.07.02. Output all varians regardless of max_af_vaf and variant classification to include MT variants and just to make the table very inclusive. 
Modified 2018.07.03 Input any vcf with any number of samples, without having to manually load their names.

'''

import re
import os
import cyvcf2
from collections import Counter
import sys
import numpy as np 
import argparse

def argument_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_vcf', required=True, help='Input vcf or vcf.gz file to be summarized into a table')
    parser.add_argument('-o', '--output_dir', required=False, default=os.getcwd(), help='Output directory')
    parser.add_argument('-a', '--af_threshold', required=False, type=float, default=0.01, help='Allele Fraction cutoff for writing to the table')
    parser.add_argument('-p', '--inclusion_pattern', required=False, default='missense|stop|splice', help='Includion criteria for filtering variant consequences')
    args = vars(parser.parse_args())
    
    return args['input_vcf'], args['output_dir'], args['af_threshold'], args['inclusion_pattern']

def variant_type(genotype):
    '''convert genotype from cyvcf2 into a string for variant classification
    if a variant is not genotyped, then index 0 is -1. need to account for this
    '''
    if genotype[0] == -1:
        return 'NA'
    else:
        if genotype == [0, 0]:
            return 'homo_ref'
        elif genotype == [0, 1] or genotype == [1, 0]:
            return 'het'
        elif genotype != [0, 0] and genotype[0] == genotype[1]:
            return 'homo_alt'
        elif genotype != [0, 1] and genotype != [0, 1] and genotype[0] != genotype[1]:
            return 'het'
        else:
            print('Not prepared for this type of variant genotypes')
            raise ValueError

def variants_genotyping(variant_genotypeINFO):
    '''calls variant_type function for a given list of variant samples'''
    variants_genotype = []
    for genotypeinfo in variant_genotypeINFO:
        geno = variant_type(genotypeinfo[0:2])
        variants_genotype.append(geno)

    return variants_genotype

def find_canonical_annotation(vep_annotation_string): 

    """VEP annotates with many alternative transcripts as well as canonical transcript
    this function finds the canonical transcript within vep_annotation_string.
    If there is no canonical transcript, which is usually the case fore intergenic,
    will just report the first annotation.
    """
    annotations = vep_annotation_string.split(',') 
    return_status = 0
    for annotation in annotations: 
        CANONICAL = annotation.split('|')[26] # CANONICAL
        if CANONICAL == 'YES': 
            return_status = 1
            return annotation 

        if return_status == 0: 
            return vep_annotation_string.split(',')[0]

def calculate_vaf(alt_depth, total_depth):
    '''for a given depth of both total and alt from a variant info, will calulate the vaf'''
    if total_depth != 0:
        return round(float(alt_depth/total_depth), 3)
    else:
        return 'NA'

def sample_vafs(variant):
    '''calculate vafs of each variants'''
    variant_vafs = []

    for index in range(0, len(variant.gt_depths)):
        vaf = calculate_vaf(variant.gt_alt_depths[index], variant.gt_depths[index])
        variant_vafs.append(vaf)

    return variant_vafs

def main():
    input_vcf, output_dir, af_threshold, inclusion_pattern= argument_parser()
    pattern_compile = re.compile(inclusion_pattern)

    output = os.path.join(output_dir, os.path.basename(input_vcf) + '.germline_summary.tsv')
    with open(output, 'w') as g:
        # write header
        vcfHandle = cyvcf2.VCF(input_vcf)
        sampleNames = vcfHandle.samples
        nSamples = len(sampleNames)
        sample_geno_vaf_header_string = ''
        for sample in sampleNames:
            sample_geno_vaf_header_string = sample_geno_vaf_header_string + sample + '_geno\t' + sample + '_vaf\t'

        print(sample_geno_vaf_header_string)
        g.write(f'CHROM\tPOS\tREF\tALT\tmax_af_gnomad\tgene\tconsequence\tprotein_change\tsift\tpolyphen\t{sample_geno_vaf_header_string}avgdepth')

        for variant in vcfHandle:
            variant_genotypeINFO = variant.genotypes
            variant_genotypes = variants_genotyping(variant_genotypeINFO)
            variant_vafs = sample_vafs(variant)
            variant_geno_vaf_string = ''
            for genotype, vaf in (zip(variant_genotypes, variant_vafs)):
                variant_geno_vaf_string = variant_geno_vaf_string + genotype + '\t'
                variant_geno_vaf_string = variant_geno_vaf_string + str(vaf)+ '\t'

            csq = variant.INFO.get('CSQ')
            canonical = find_canonical_annotation(csq)
            consequence = canonical.split('|')[1]
            gene = canonical.split('|')[3]
            protein_change = canonical.split('|')[11]
            sift = canonical.split('|')[36]
            polyphen = canonical.split('|')[37]
            domains = canonical.split('|')[38]

            max_af = find_canonical_annotation(csq).split('|')[57]
            gnomad_af = find_canonical_annotation(csq).split('|')[48]
            if max_af.strip() == '':
                max_af = 0

            if gnomad_af.strip() == '':
                gnomad_af = 0


            max_af_gnomad = max(float(max_af), float(gnomad_af))
            totaldepth = variant.INFO['DP']
            avgdepth = round(float(np.sum(np.array(variant.format('DP'))) / 9), 3)
            # only report rare (10% threshold) that are not intronic/regulatory/synonymous... etc and 
            # not in very high depth region of the genome. 

            if max_af_gnomad < af_threshold and re.search(pattern_compile, consequence):
                g.write(f'\n{variant.CHROM}\t{variant.POS}\t{variant.REF}\t{variant.ALT[0]}\t{max_af_gnomad}\t{gene}\t{consequence}\t{protein_change}\t{sift}\t{polyphen}\t{variant_geno_vaf_string}{avgdepth}')

                

if __name__=='__main__':
    main()
