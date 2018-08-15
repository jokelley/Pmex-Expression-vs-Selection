################################################################################################
#
#           3. Calculate FST between adjacent populations                      
#
#
################################################################################################

# Calculate FST between populations using vcftools
# Need to calculate FST between adjacent population pairs
# puya_nonsulfidic_samples.txt contains a list of sample IDs for Puyacatengo non-sulfidic individuals
# puya_sulfidic_samples.txt contains a list of sample IDs for Puyacatengo sulfidic individuals
# taco_nonsulfidic_samples.txt contains a list of sample IDs for Tacotalpa non-sulfidic individuals
# taco_sulfidic_samples.txt contains a list of sample IDs for Tacotalpa sulfidic individuals

vcftools --vcf merged_snps_minDP8_maxmiss0.9.ann.vcf --weir-fst-pop puya_nonsulfidic_samples.txt --weir-fst-pop puya_sulfidic_samples.txt --out puya_fst
vcftools --vcf merged_snps_minDP8_maxmiss0.9.ann.vcf --weir-fst-pop taco_nonsulfidic_samples.txt --weir-fst-pop taco_sulfidic_samples.txt --out taco_fst