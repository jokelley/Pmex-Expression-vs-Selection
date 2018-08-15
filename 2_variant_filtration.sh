################################################################################################
#
#           2. Calling/Filtering SNPs                      
#
#
################################################################################################

# Use GATK to call SNPs on a per population basis
# Create sequence dictionary for reference prior to calling SNPs

java -jar picard.jar CreateSequenceDictionary \ 
      R=<path to reference> \
      O=<path to reference without .fa extension>.dict

# Use UnifiedGenotyper in GATK to call SNPs per population
# Use EMIT_ALL_SITES to show all sites in the resulting vcf files (facilitates merging of vcf files later)
# Call SNPs for the Puyacatengo non-sulfidic population
	  
java -Xmx2048m -jar GenomeAnalysisTK.jar \
-T UnifiedGenotyper \
-R <path to reference> \
-o puya_nonsulfidic.vcf.gz \
--output_mode EMIT_ALL_SITES \
-I MX50.bam \
-I MX51.bam \
-I MX52.bam \
-I MX53.bam \
-I MX54.F.bam \
-I MX55.F.bam

# Call SNPs for the Puyacatengo sulfidic population
	
java -Xmx2048m -jar GenomeAnalysisTK.jar \
-T UnifiedGenotyper \
-R <path to reference> \
-o puya_sulfidic.vcf.gz \
--output_mode EMIT_ALL_SITES \
-I MX44.bam \
-I MX45.bam \
-I MX46.bam \
-I MX48.bam \
-I MX49.F.bam \

# Call SNPs for the Tacotalpa non-sulfidic population
java -Xmx2048m -jar GenomeAnalysisTK.jar \
-T UnifiedGenotyper \
-R <path to reference> \
-o taco_nonsulfidic.vcf.gz \
--output_mode EMIT_ALL_SITES \
-I MX57.bam \
-I MX59.bam \
-I MX60.bam \
-I MX61.bam \
-I MX62.F.bam \
-I MX63.F.bam

# Call SNPs for the Tacotalpa sulfidic population

java -Xmx2048m -jar GenomeAnalysisTK.jar \
-T UnifiedGenotyper \
-R <path to reference> \
-o taco_sulfidic.vcf.gz \
--output_mode EMIT_ALL_SITES \
-I MX71.bam \
-I MX73.bam \
-I MX74.bam \
-I MX75.bam \
-I MX76.F.bam \
-I MX77.F.bam

# Use vcf-merge (Perl module from vcftools) to merge population vcf files into one vcf
# Move into directory that contains vcftools Perl modules
# Use bgzip to zip the merged vcf file, use tabix to index the merged vcf file

./vcf-merge puya_nonsulfidic.vcf.gz puya_sulfidic.vcf.gz taco_nonsulfidic.vcf.gz taco_sulfidic.vcf.gz | bgzip -c > merged.vcf.gz
tabix merged.vcf.gz

# Filter merged vcf file with GATK
# Use SelectVariants to pull out biallelic SNPs

java -jar GenomeAnalysisTK.jar \
    -T SelectVariants \
    -R <path to reference> \
    -V merged.vcf.gz \
    -selectType SNP \
    -restrictAllelesTo BIALLELIC \
    -o merged_just_biallelic_snps.vcf

# Use VariantFiltration to apply standard GATK hard filters

java -jar GenomeAnalysisTK.jar \
    -T VariantFiltration \
    -R <path to reference> \
    -V merged_just_biallelic_snps.vcf \
    --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
    --filterName "my_snp_filter" \
    -o filter_applied_biallelic_snps.vcf

# Use SelectVariants to exclude sites that didn't pass the hard filters

java -jar GenomeAnalysisTK.jar \
    -T SelectVariants \
    -R <path to reference> \
    -V filter_applied_biallelic_snps.vcf \
    --excludeFiltered \
    -o filter_excluded_biallelic_snps.vcf

# Use vcftools to further filter the SNPs. --minDP 8 only keeps genotypes supported by at least 8x coverage. --max-missing 0.9 only keeps sites where 90% of individuals have an inferred genotype.
# The combination of min-alleles 2 and max-alleles 2 ensures that only biallelic sites are kept
# Vcftools adds a suffix to the end of the out file automatically (".recode.vcf")

vcftools --vcf filter_excluded_biallelic_snps.vcf --min-alleles 2 --max-alleles 2 --minDP 8 --max-missing 0.9 --recode --recode-INFO-all --out merged_snps_minDP8_maxmiss0.9

# Annotate VCF using SNPEff
java -Xmx4g -jar snpEff.jar Xipmac4.4.2.75 merged_snps_minDP8_maxmiss0.9.recode.vcf > merged_snps_minDP8_maxmiss0.9.ann.vcf
