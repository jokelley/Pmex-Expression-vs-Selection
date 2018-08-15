################################################################################################
#
#           4. Preparing to perform demographic history analyses                      
#
#
################################################################################################

# These analyses were performed with pooled sequencing data from the European Nucleotide Archive, accession PRJEB8912
# This is an example with just the Puyacatengo populations
# Follow the same pipeline as in the 3 previous scripts, until you get to SNP calling, use samtools instead (recommended by the PoPoolation2 manual)
# Prerequisites: samtools, PoPoolation2
# Call SNPs with samtools mpileup (NS = non-sulfidic bam, S = sulfidic bam)

samtools mpileup -B Puy_NS.sort.rg.rm_dup.bam Puy_S.sort.rg.rm_dup.bam > Puy.mpileup

# Pull out sites with at least 20x coverage in each population

awk '$4 >= 20 && $7 >= 20' Puy.mpileup > Puy_20.mpileup

# Convert to sync format for PoPoolation2

java -ea -Xmx7g -jar mpileup2sync.jar --input Puy_20.mpileup --output Puy_20java.sync --fastq-type sanger --min-qual 20 --threads 2

# Separate allele counts into separate columns in the sync file

tr -s ':' $'\t' <Puy_20java.sync >Puy_20_columns_java.sync

# Add new columns that are the proportions of allele counts
# Columns 16-19 contain the proportions for each nucleotide (one per column) at each site for the NS population, 20-23 contain the proportions for the S population
# These counts are rounded to the nearest whole number

awk -v OFS='\t' '{$16 = ($4 != 0) ? sprintf("%.3f", $4 / ($4 + $5 + $6 + $7)) : "0"}1' Puy_20_columns_java.sync > Puy_20_columns_java16.sync

awk -v OFS='\t' '{$17 = ($5 != 0) ? sprintf("%.3f", $5 / ($4 + $5 + $6 + $7)) : "0"}1' Puy_20_columns_java16.sync > Puy_20_columns_java17.sync

rm Puy_20_columns_java16.sync

awk -v OFS='\t' '{$18 = ($6 != 0) ? sprintf("%.3f", $6 / ($4 + $5 + $6 + $7)) : "0"}1' Puy_20_columns_java17.sync > Puy_20_columns_java18.sync

rm Puy_20_columns_java17.sync

awk -v OFS='\t' '{$19 = ($7 != 0) ? sprintf("%.3f", $7 / ($4 + $5 + $6 + $7)) : "0"}1' Puy_20_columns_java18.sync > Puy_20_columns_java19.sync

rm Puy_20_columns_java18.sync

awk -v OFS='\t' '{$20 = ($10 != 0) ? sprintf("%.3f", $10 / ($10 + $11 + $12 + $13)) : "0"}1' Puy_20_columns_java19.sync > Puy_20_columns_java20.sync

rm Puy_20_columns_java19.sync

awk -v OFS='\t' '{$21 = ($11 != 0) ? sprintf("%.3f", $11 / ($10 + $11 + $12 + $13)) : "0"}1' Puy_20_columns_java20.sync > Puy_20_columns_java21.sync

rm Puy_20_columns_java20.sync

awk -v OFS='\t' '{$22 = ($12 != 0) ? sprintf("%.3f", $12 / ($10 + $11 + $12 + $13)) : "0"}1' Puy_20_columns_java21.sync > Puy_20_columns_java22.sync

rm Puy_20_columns_java21.sync

awk -v OFS='\t' '{$23 = ($13 != 0) ? sprintf("%.3f", $13 / ($10 + $11 + $12 + $13)) : "0"}1' Puy_20_columns_java22.sync > Puy_20_columns_java23.sync

rm Puy_20_columns_java22.sync

# Keep only the site information (scaffold and position) as well as the proportion of each allele

cut -f3,4,5,6,7,8,9,10,11,12,13,14,15 --complement Puy_20_columns_java23.sync > Puy_20_columns_java_only_proportions.sync

# Multiply the proportions by 20 to get a normalized count for each allele in each population

awk -v OFS='\t' '{$11 = ($3 != 0) ? sprintf("%.0f", $3*20) : "0"}1' Puy_20_columns_java_only_proportions.sync > Puy_20_columns_java_only_proportions11.sync

awk -v OFS='\t' '{$12 = ($4 != 0) ? sprintf("%.0f", $4*20) : "0"}1' Puy_20_columns_java_only_proportions11.sync > Puy_20_columns_java_only_proportions12.sync

rm Puy_20_columns_java_only_proportions11.sync

awk -v OFS='\t' '{$13 = ($5 != 0) ? sprintf("%.0f", $5*20) : "0"}1' Puy_20_columns_java_only_proportions12.sync > Puy_20_columns_java_only_proportions13.sync

rm Puy_20_columns_java_only_proportions12.sync

awk -v OFS='\t' '{$14 = ($6 != 0) ? sprintf("%.0f", $6*20) : "0"}1' Puy_20_columns_java_only_proportions13.sync > Puy_20_columns_java_only_proportions14.sync

rm Puy_20_columns_java_only_proportions13.sync

awk -v OFS='\t' '{$15 = ($7 != 0) ? sprintf("%.0f", $7*20) : "0"}1' Puy_20_columns_java_only_proportions14.sync > Puy_20_columns_java_only_proportions15.sync

rm Puy_20_columns_java_only_proportions14.sync

awk -v OFS='\t' '{$16 = ($8 != 0) ? sprintf("%.0f", $8*20) : "0"}1' Puy_20_columns_java_only_proportions15.sync > Puy_20_columns_java_only_proportions16.sync

rm Puy_20_columns_java_only_proportions15.sync

awk -v OFS='\t' '{$17 = ($9 != 0) ? sprintf("%.0f", $9*20) : "0"}1' Puy_20_columns_java_only_proportions16.sync > Puy_20_columns_java_only_proportions17.sync

rm Puy_20_columns_java_only_proportions16.sync

awk -v OFS='\t' '{$18 = ($10 != 0) ? sprintf("%.0f", $10*20) : "0"}1' Puy_20_columns_java_only_proportions17.sync > Puy_20_columns_java_only_proportions18.sync

rm Puy_20_columns_java_only_proportions17.sync

# Keep only the columns that contain the site information (scaffold and position the normalized counts for each allele for each population

cut -f3,4,5,6,7,8,9,10 --complement Puy_20_columns_java_only_proportions18.sync > Puy_20_columns_java_only_count.sync

# Filter out non-variant sites between populations (any sites where each population has 20 counts of the same allele)
awk -v OFS='\t' '!(($3==20 && $7==20) || ($4==20 && $8==20) || ($5==20 && $9==20) || ($6==20 && $10==20))' Puy_20_columns_java_only_count.sync >\
Puy_20_columns_java_only_count_keepsnps.sync

# Obtain count information for each bin in a folded site frequency spectrum
# The count is saved in a file that includes the minor allele counts for the bin (e.g. Puy_0_1.txt means 0 in the NS, and 1 in the S)
# These counts can be assembled into a shared site frequency spectrum that can then be used for demographic analysis

awk -v OFS='\t' '(($3==20 && $7==19) || \
($4==20 && $8==19) || \
($5==20 && $9==19) || \
($6==20 && $10==19))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_0_1.txt

awk -v OFS='\t' '(($3==20 && $7==18) || \
($4==20 && $8==18) || \
($5==20 && $9==18) || \
($6==20 && $10==18))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_0_2.txt

awk -v OFS='\t' '(($3==20 && $7==17) || \
($4==20 && $8==17) || \
($5==20 && $9==17) || \
($6==20 && $10==17))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_0_3.txt

awk -v OFS='\t' '(($3==20 && $7==16) || \
($4==20 && $8==16) || \
($5==20 && $9==16) || \
($6==20 && $10==16))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_0_4.txt

awk -v OFS='\t' '(($3==20 && $7==15) || \
($4==20 && $8==15) || \
($5==20 && $9==15) || \
($6==20 && $10==15))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_0_5.txt

awk -v OFS='\t' '(($3==20 && $7==14) || \
($4==20 && $8==14) || \
($5==20 && $9==14) || \
($6==20 && $10==14))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_0_6.txt

awk -v OFS='\t' '(($3==20 && $7==13) || \
($4==20 && $8==13) || \
($5==20 && $9==13) || \
($6==20 && $10==13))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_0_7.txt

awk -v OFS='\t' '(($3==20 && $7==12) || \
($4==20 && $8==12) || \
($5==20 && $9==12) || \
($6==20 && $10==12))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_0_8.txt

awk -v OFS='\t' '(($3==20 && $7==11) || \
($4==20 && $8==11) || \
($5==20 && $9==11) || \
($6==20 && $10==11))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_0_9.txt

awk -v OFS='\t' '(($3==20 && $7==10) || \
($4==20 && $8==10) || \
($5==20 && $9==10) || \
($6==20 && $10==10))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_0_10.txt

awk -v OFS='\t' '(($3==20 && $7==9) || \
($4==20 && $8==9) || \
($5==20 && $9==9) || \
($6==20 && $10==9))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_0_11.txt

awk -v OFS='\t' '(($3==20 && $7==8) || \
($4==20 && $8==8) || \
($5==20 && $9==8) || \
($6==20 && $10==8))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_0_12.txt

awk -v OFS='\t' '(($3==20 && $7==7) || \
($4==20 && $8==7) || \
($5==20 && $9==7) || \
($6==20 && $10==7))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_0_13.txt

awk -v OFS='\t' '(($3==20 && $7==6) || \
($4==20 && $8==6) || \
($5==20 && $9==6) || \
($6==20 && $10==6))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_0_14.txt

awk -v OFS='\t' '(($3==20 && $7==5) || \
($4==20 && $8==5) || \
($5==20 && $9==5) || \
($6==20 && $10==5))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_0_15.txt

awk -v OFS='\t' '(($3==20 && $7==4) || \
($4==20 && $8==4) || \
($5==20 && $9==4) || \
($6==20 && $10==4))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_0_16.txt

awk -v OFS='\t' '(($3==20 && $7==3) || \
($4==20 && $8==3) || \
($5==20 && $9==3) || \
($6==20 && $10==3))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_0_17.txt

awk -v OFS='\t' '(($3==20 && $7==2) || \
($4==20 && $8==2) || \
($5==20 && $9==2) || \
($6==20 && $10==2))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_0_18.txt

awk -v OFS='\t' '(($3==20 && $7==1) || \
($4==20 && $8==1) || \
($5==20 && $9==1) || \
($6==20 && $10==1))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_0_19.txt

awk -v OFS='\t' '(($3==19 && $7==19) || \
($4==19 && $8==19) || \
($5==19 && $9==19) || \
($6==19 && $10==19))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_1_1.txt

awk -v OFS='\t' '(($3==19 && $7==18) || \
($4==19 && $8==18) || \
($5==19 && $9==18) || \
($6==19 && $10==18))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_1_2.txt

awk -v OFS='\t' '(($3==19 && $7==17) || \
($4==19 && $8==17) || \
($5==19 && $9==17) || \
($6==19 && $10==17))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_1_3.txt

awk -v OFS='\t' '(($3==19 && $7==16) || \
($4==19 && $8==16) || \
($5==19 && $9==16) || \
($6==19 && $10==16))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_1_4.txt

awk -v OFS='\t' '(($3==19 && $7==15) || \
($4==19 && $8==15) || \
($5==19 && $9==15) || \
($6==19 && $10==15))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_1_5.txt

awk -v OFS='\t' '(($3==19 && $7==14) || \
($4==19 && $8==14) || \
($5==19 && $9==14) || \
($6==19 && $10==14))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_1_6.txt

awk -v OFS='\t' '(($3==19 && $7==13) || \
($4==19 && $8==13) || \
($5==19 && $9==13) || \
($6==19 && $10==13))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_1_7.txt

awk -v OFS='\t' '(($3==19 && $7==12) || \
($4==19 && $8==12) || \
($5==19 && $9==12) || \
($6==19 && $10==12))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_1_8.txt

awk -v OFS='\t' '(($3==19 && $7==11) || \
($4==19 && $8==11) || \
($5==19 && $9==11) || \
($6==19 && $10==11))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_1_9.txt

awk -v OFS='\t' '(($3==19 && $7==10) || \
($4==19 && $8==10) || \
($5==19 && $9==10) || \
($6==19 && $10==10))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_1_10.txt

awk -v OFS='\t' '(($3==19 && $7==9) || \
($4==19 && $8==9) || \
($5==19 && $9==9) || \
($6==19 && $10==9))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_1_11.txt

awk -v OFS='\t' '(($3==19 && $7==8) || \
($4==19 && $8==8) || \
($5==19 && $9==8) || \
($6==19 && $10==8))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_1_12.txt

awk -v OFS='\t' '(($3==19 && $7==7) || \
($4==19 && $8==7) || \
($5==19 && $9==7) || \
($6==19 && $10==7))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_1_13.txt

awk -v OFS='\t' '(($3==19 && $7==6) || \
($4==19 && $8==6) || \
($5==19 && $9==6) || \
($6==19 && $10==6))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_1_14.txt

awk -v OFS='\t' '(($3==19 && $7==5) || \
($4==19 && $8==5) || \
($5==19 && $9==5) || \
($6==19 && $10==5))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_1_15.txt

awk -v OFS='\t' '(($3==19 && $7==4) || \
($4==19 && $8==4) || \
($5==19 && $9==4) || \
($6==19 && $10==4))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_1_16.txt

awk -v OFS='\t' '(($3==19 && $7==3) || \
($4==19 && $8==3) || \
($5==19 && $9==3) || \
($6==19 && $10==3))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_1_17.txt

awk -v OFS='\t' '(($3==19 && $7==2) || \
($4==19 && $8==2) || \
($5==19 && $9==2) || \
($6==19 && $10==2))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_1_18.txt

awk -v OFS='\t' '(($3==19 && $7==1) || \
($4==19 && $8==1) || \
($5==19 && $9==1) || \
($6==19 && $10==1))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_1_19.txt

awk -v OFS='\t' '(($3==18 && $7==19) || \
($4==18 && $8==19) || \
($5==18 && $9==19) || \
($6==18 && $10==19))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_2_1.txt

awk -v OFS='\t' '(($3==18 && $7==18) || \
($4==18 && $8==18) || \
($5==18 && $9==18) || \
($6==18 && $10==18))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_2_2.txt

awk -v OFS='\t' '(($3==18 && $7==17) || \
($4==18 && $8==17) || \
($5==18 && $9==17) || \
($6==18 && $10==17))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_2_3.txt

awk -v OFS='\t' '(($3==18 && $7==16) || \
($4==18 && $8==16) || \
($5==18 && $9==16) || \
($6==18 && $10==16))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_2_4.txt

awk -v OFS='\t' '(($3==18 && $7==15) || \
($4==18 && $8==15) || \
($5==18 && $9==15) || \
($6==18 && $10==15))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_2_5.txt

awk -v OFS='\t' '(($3==18 && $7==14) || \
($4==18 && $8==14) || \
($5==18 && $9==14) || \
($6==18 && $10==14))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_2_6.txt

awk -v OFS='\t' '(($3==18 && $7==13) || \
($4==18 && $8==13) || \
($5==18 && $9==13) || \
($6==18 && $10==13))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_2_7.txt

awk -v OFS='\t' '(($3==18 && $7==12) || \
($4==18 && $8==12) || \
($5==18 && $9==12) || \
($6==18 && $10==12))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_2_8.txt

awk -v OFS='\t' '(($3==18 && $7==11) || \
($4==18 && $8==11) || \
($5==18 && $9==11) || \
($6==18 && $10==11))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_2_9.txt

awk -v OFS='\t' '(($3==18 && $7==10) || \
($4==18 && $8==10) || \
($5==18 && $9==10) || \
($6==18 && $10==10))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_2_10.txt

awk -v OFS='\t' '(($3==18 && $7==9) || \
($4==18 && $8==9) || \
($5==18 && $9==9) || \
($6==18 && $10==9))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_2_11.txt

awk -v OFS='\t' '(($3==18 && $7==8) || \
($4==18 && $8==8) || \
($5==18 && $9==8) || \
($6==18 && $10==8))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_2_12.txt

awk -v OFS='\t' '(($3==18 && $7==7) || \
($4==18 && $8==7) || \
($5==18 && $9==7) || \
($6==18 && $10==7))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_2_13.txt

awk -v OFS='\t' '(($3==18 && $7==6) || \
($4==18 && $8==6) || \
($5==18 && $9==6) || \
($6==18 && $10==6))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_2_14.txt

awk -v OFS='\t' '(($3==18 && $7==5) || \
($4==18 && $8==5) || \
($5==18 && $9==5) || \
($6==18 && $10==5))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_2_15.txt

awk -v OFS='\t' '(($3==18 && $7==4) || \
($4==18 && $8==4) || \
($5==18 && $9==4) || \
($6==18 && $10==4))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_2_16.txt

awk -v OFS='\t' '(($3==18 && $7==3) || \
($4==18 && $8==3) || \
($5==18 && $9==3) || \
($6==18 && $10==3))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_2_17.txt

awk -v OFS='\t' '(($3==18 && $7==2) || \
($4==18 && $8==2) || \
($5==18 && $9==2) || \
($6==18 && $10==2))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_2_18.txt

awk -v OFS='\t' '(($3==18 && $7==1) || \
($4==18 && $8==1) || \
($5==18 && $9==1) || \
($6==18 && $10==1))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_18_1.txt

awk -v OFS='\t' '(($3==17 && $7==19) || \
($4==17 && $8==19) || \
($5==17 && $9==19) || \
($6==17 && $10==19))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_3_1.txt

awk -v OFS='\t' '(($3==17 && $7==18) || \
($4==17 && $8==18) || \
($5==17 && $9==18) || \
($6==17 && $10==18))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_3_2.txt

awk -v OFS='\t' '(($3==17 && $7==17) || \
($4==17 && $8==17) || \
($5==17 && $9==17) || \
($6==17 && $10==17))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_3_3.txt

awk -v OFS='\t' '(($3==17 && $7==16) || \
($4==17 && $8==16) || \
($5==17 && $9==16) || \
($6==17 && $10==16))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_3_4.txt

awk -v OFS='\t' '(($3==17 && $7==15) || \
($4==17 && $8==15) || \
($5==17 && $9==15) || \
($6==17 && $10==15))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_3_5.txt

awk -v OFS='\t' '(($3==17 && $7==14) || \
($4==17 && $8==14) || \
($5==17 && $9==14) || \
($6==17 && $10==14))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_3_6.txt

awk -v OFS='\t' '(($3==17 && $7==13) || \
($4==17 && $8==13) || \
($5==17 && $9==13) || \
($6==17 && $10==13))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_3_7.txt

awk -v OFS='\t' '(($3==17 && $7==12) || \
($4==17 && $8==12) || \
($5==17 && $9==12) || \
($6==17 && $10==12))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_3_8.txt

awk -v OFS='\t' '(($3==17 && $7==11) || \
($4==17 && $8==11) || \
($5==17 && $9==11) || \
($6==17 && $10==11))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_3_9.txt

awk -v OFS='\t' '(($3==17 && $7==10) || \
($4==17 && $8==10) || \
($5==17 && $9==10) || \
($6==17 && $10==10))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_3_10.txt

awk -v OFS='\t' '(($3==17 && $7==9) || \
($4==17 && $8==9) || \
($5==17 && $9==9) || \
($6==17 && $10==9))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_3_11.txt

awk -v OFS='\t' '(($3==17 && $7==8) || \
($4==17 && $8==8) || \
($5==17 && $9==8) || \
($6==17 && $10==8))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_3_12.txt

awk -v OFS='\t' '(($3==17 && $7==7) || \
($4==17 && $8==7) || \
($5==17 && $9==7) || \
($6==17 && $10==7))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_3_13.txt

awk -v OFS='\t' '(($3==17 && $7==6) || \
($4==17 && $8==6) || \
($5==17 && $9==6) || \
($6==17 && $10==6))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_3_14.txt

awk -v OFS='\t' '(($3==17 && $7==5) || \
($4==17 && $8==5) || \
($5==17 && $9==5) || \
($6==17 && $10==5))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_3_15.txt

awk -v OFS='\t' '(($3==17 && $7==4) || \
($4==17 && $8==4) || \
($5==17 && $9==4) || \
($6==17 && $10==4))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_3_16.txt

awk -v OFS='\t' '(($3==17 && $7==3) || \
($4==17 && $8==3) || \
($5==17 && $9==3) || \
($6==17 && $10==3))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_3_17.txt

awk -v OFS='\t' '(($3==17 && $7==2) || \
($4==17 && $8==2) || \
($5==17 && $9==2) || \
($6==17 && $10==2))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_17_2.txt

awk -v OFS='\t' '(($3==17 && $7==1) || \
($4==17 && $8==1) || \
($5==17 && $9==1) || \
($6==17 && $10==1))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_17_1.txt

awk -v OFS='\t' '(($3==16 && $7==19) || \
($4==16 && $8==19) || \
($5==16 && $9==19) || \
($6==16 && $10==19))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_4_1.txt

awk -v OFS='\t' '(($3==16 && $7==18) || \
($4==16 && $8==18) || \
($5==16 && $9==18) || \
($6==16 && $10==18))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_4_2.txt

awk -v OFS='\t' '(($3==16 && $7==17) || \
($4==16 && $8==17) || \
($5==16 && $9==17) || \
($6==16 && $10==17))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_4_3.txt

awk -v OFS='\t' '(($3==16 && $7==16) || \
($4==16 && $8==16) || \
($5==16 && $9==16) || \
($6==16 && $10==16))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_4_4.txt

awk -v OFS='\t' '(($3==16 && $7==15) || \
($4==16 && $8==15) || \
($5==16 && $9==15) || \
($6==16 && $10==15))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_4_5.txt

awk -v OFS='\t' '(($3==16 && $7==14) || \
($4==16 && $8==14) || \
($5==16 && $9==14) || \
($6==16 && $10==14))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_4_6.txt

awk -v OFS='\t' '(($3==16 && $7==13) || \
($4==16 && $8==13) || \
($5==16 && $9==13) || \
($6==16 && $10==13))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_4_7.txt

awk -v OFS='\t' '(($3==16 && $7==12) || \
($4==16 && $8==12) || \
($5==16 && $9==12) || \
($6==16 && $10==12))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_4_8.txt

awk -v OFS='\t' '(($3==16 && $7==11) || \
($4==16 && $8==11) || \
($5==16 && $9==11) || \
($6==16 && $10==11))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_4_9.txt

awk -v OFS='\t' '(($3==16 && $7==10) || \
($4==16 && $8==10) || \
($5==16 && $9==10) || \
($6==16 && $10==10))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_4_10.txt

awk -v OFS='\t' '(($3==16 && $7==9) || \
($4==16 && $8==9) || \
($5==16 && $9==9) || \
($6==16 && $10==9))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_4_11.txt

awk -v OFS='\t' '(($3==16 && $7==8) || \
($4==16 && $8==8) || \
($5==16 && $9==8) || \
($6==16 && $10==8))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_4_12.txt

awk -v OFS='\t' '(($3==16 && $7==7) || \
($4==16 && $8==7) || \
($5==16 && $9==7) || \
($6==16 && $10==7))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_4_13.txt

awk -v OFS='\t' '(($3==16 && $7==6) || \
($4==16 && $8==6) || \
($5==16 && $9==6) || \
($6==16 && $10==6))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_4_14.txt

awk -v OFS='\t' '(($3==16 && $7==5) || \
($4==16 && $8==5) || \
($5==16 && $9==5) || \
($6==16 && $10==5))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_4_15.txt

awk -v OFS='\t' '(($3==16 && $7==4) || \
($4==16 && $8==4) || \
($5==16 && $9==4) || \
($6==16 && $10==4))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_4_16.txt

awk -v OFS='\t' '(($3==16 && $7==3) || \
($4==16 && $8==3) || \
($5==16 && $9==3) || \
($6==16 && $10==3))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_16_3.txt

awk -v OFS='\t' '(($3==16 && $7==2) || \
($4==16 && $8==2) || \
($5==16 && $9==2) || \
($6==16 && $10==2))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_16_2.txt

awk -v OFS='\t' '(($3==16 && $7==1) || \
($4==16 && $8==1) || \
($5==16 && $9==1) || \
($6==16 && $10==1))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_16_1.txt

awk -v OFS='\t' '(($3==15 && $7==19) || \
($4==15 && $8==19) || \
($5==15 && $9==19) || \
($6==15 && $10==19))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_5_1.txt

awk -v OFS='\t' '(($3==15 && $7==18) || \
($4==15 && $8==18) || \
($5==15 && $9==18) || \
($6==15 && $10==18))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_5_2.txt

awk -v OFS='\t' '(($3==15 && $7==17) || \
($4==15 && $8==17) || \
($5==15 && $9==17) || \
($6==15 && $10==17))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_5_3.txt

awk -v OFS='\t' '(($3==15 && $7==16) || \
($4==15 && $8==16) || \
($5==15 && $9==16) || \
($6==15 && $10==16))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_5_4.txt

awk -v OFS='\t' '(($3==15 && $7==15) || \
($4==15 && $8==15) || \
($5==15 && $9==15) || \
($6==15 && $10==15))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_5_5.txt

awk -v OFS='\t' '(($3==15 && $7==14) || \
($4==15 && $8==14) || \
($5==15 && $9==14) || \
($6==15 && $10==14))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_5_6.txt

awk -v OFS='\t' '(($3==15 && $7==13) || \
($4==15 && $8==13) || \
($5==15 && $9==13) || \
($6==15 && $10==13))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_5_7.txt

awk -v OFS='\t' '(($3==15 && $7==12) || \
($4==15 && $8==12) || \
($5==15 && $9==12) || \
($6==15 && $10==12))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_5_8.txt

awk -v OFS='\t' '(($3==15 && $7==11) || \
($4==15 && $8==11) || \
($5==15 && $9==11) || \
($6==15 && $10==11))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_5_9.txt

awk -v OFS='\t' '(($3==15 && $7==10) || \
($4==15 && $8==10) || \
($5==15 && $9==10) || \
($6==15 && $10==10))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_5_10.txt

awk -v OFS='\t' '(($3==15 && $7==9) || \
($4==15 && $8==9) || \
($5==15 && $9==9) || \
($6==15 && $10==9))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_5_11.txt

awk -v OFS='\t' '(($3==15 && $7==8) || \
($4==15 && $8==8) || \
($5==15 && $9==8) || \
($6==15 && $10==8))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_5_12.txt

awk -v OFS='\t' '(($3==15 && $7==7) || \
($4==15 && $8==7) || \
($5==15 && $9==7) || \
($6==15 && $10==7))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_5_13.txt

awk -v OFS='\t' '(($3==15 && $7==6) || \
($4==15 && $8==6) || \
($5==15 && $9==6) || \
($6==15 && $10==6))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_5_14.txt

awk -v OFS='\t' '(($3==15 && $7==5) || \
($4==15 && $8==5) || \
($5==15 && $9==5) || \
($6==15 && $10==5))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_5_15.txt

awk -v OFS='\t' '(($3==15 && $7==4) || \
($4==15 && $8==4) || \
($5==15 && $9==4) || \
($6==15 && $10==4))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_15_4.txt

awk -v OFS='\t' '(($3==15 && $7==3) || \
($4==15 && $8==3) || \
($5==15 && $9==3) || \
($6==15 && $10==3))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_15_3.txt

awk -v OFS='\t' '(($3==15 && $7==2) || \
($4==15 && $8==2) || \
($5==15 && $9==2) || \
($6==15 && $10==2))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_15_2.txt

awk -v OFS='\t' '(($3==15 && $7==1) || \
($4==15 && $8==1) || \
($5==15 && $9==1) || \
($6==15 && $10==1))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_15_1.txt

awk -v OFS='\t' '(($3==14 && $7==19) || \
($4==14 && $8==19) || \
($5==14 && $9==19) || \
($6==14 && $10==19))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_6_1.txt

awk -v OFS='\t' '(($3==14 && $7==18) || \
($4==14 && $8==18) || \
($5==14 && $9==18) || \
($6==14 && $10==18))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_6_2.txt

awk -v OFS='\t' '(($3==14 && $7==17) || \
($4==14 && $8==17) || \
($5==14 && $9==17) || \
($6==14 && $10==17))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_6_3.txt

awk -v OFS='\t' '(($3==14 && $7==16) || \
($4==14 && $8==16) || \
($5==14 && $9==16) || \
($6==14 && $10==16))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_6_4.txt

awk -v OFS='\t' '(($3==14 && $7==15) || \
($4==14 && $8==15) || \
($5==14 && $9==15) || \
($6==14 && $10==15))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_6_5.txt

awk -v OFS='\t' '(($3==14 && $7==14) || \
($4==14 && $8==14) || \
($5==14 && $9==14) || \
($6==14 && $10==14))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_6_6.txt

awk -v OFS='\t' '(($3==14 && $7==13) || \
($4==14 && $8==13) || \
($5==14 && $9==13) || \
($6==14 && $10==13))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_6_7.txt

awk -v OFS='\t' '(($3==14 && $7==12) || \
($4==14 && $8==12) || \
($5==14 && $9==12) || \
($6==14 && $10==12))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_6_8.txt

awk -v OFS='\t' '(($3==14 && $7==11) || \
($4==14 && $8==11) || \
($5==14 && $9==11) || \
($6==14 && $10==11))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_6_9.txt

awk -v OFS='\t' '(($3==14 && $7==10) || \
($4==14 && $8==10) || \
($5==14 && $9==10) || \
($6==14 && $10==10))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_6_10.txt

awk -v OFS='\t' '(($3==14 && $7==9) || \
($4==14 && $8==9) || \
($5==14 && $9==9) || \
($6==14 && $10==9))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_6_11.txt

awk -v OFS='\t' '(($3==14 && $7==8) || \
($4==14 && $8==8) || \
($5==14 && $9==8) || \
($6==14 && $10==8))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_6_12.txt

awk -v OFS='\t' '(($3==14 && $7==7) || \
($4==14 && $8==7) || \
($5==14 && $9==7) || \
($6==14 && $10==7))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_6_13.txt

awk -v OFS='\t' '(($3==14 && $7==6) || \
($4==14 && $8==6) || \
($5==14 && $9==6) || \
($6==14 && $10==6))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_6_14.txt

awk -v OFS='\t' '(($3==14 && $7==5) || \
($4==14 && $8==5) || \
($5==14 && $9==5) || \
($6==14 && $10==5))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_14_5.txt

awk -v OFS='\t' '(($3==14 && $7==4) || \
($4==14 && $8==4) || \
($5==14 && $9==4) || \
($6==14 && $10==4))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_14_4.txt

awk -v OFS='\t' '(($3==14 && $7==3) || \
($4==14 && $8==3) || \
($5==14 && $9==3) || \
($6==14 && $10==3))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_14_3.txt

awk -v OFS='\t' '(($3==14 && $7==2) || \
($4==14 && $8==2) || \
($5==14 && $9==2) || \
($6==14 && $10==2))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_14_2.txt

awk -v OFS='\t' '(($3==14 && $7==1) || \
($4==14 && $8==1) || \
($5==14 && $9==1) || \
($6==14 && $10==1))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_14_1.txt

awk -v OFS='\t' '(($3==13 && $7==19) || \
($4==13 && $8==19) || \
($5==13 && $9==19) || \
($6==13 && $10==19))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_7_1.txt

awk -v OFS='\t' '(($3==13 && $7==18) || \
($4==13 && $8==18) || \
($5==13 && $9==18) || \
($6==13 && $10==18))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_7_2.txt

awk -v OFS='\t' '(($3==13 && $7==17) || \
($4==13 && $8==17) || \
($5==13 && $9==17) || \
($6==13 && $10==17))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_7_3.txt

awk -v OFS='\t' '(($3==13 && $7==16) || \
($4==13 && $8==16) || \
($5==13 && $9==16) || \
($6==13 && $10==16))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_7_4.txt

awk -v OFS='\t' '(($3==13 && $7==15) || \
($4==13 && $8==15) || \
($5==13 && $9==15) || \
($6==13 && $10==15))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_7_5.txt

awk -v OFS='\t' '(($3==13 && $7==14) || \
($4==13 && $8==14) || \
($5==13 && $9==14) || \
($6==13 && $10==14))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_7_6.txt

awk -v OFS='\t' '(($3==13 && $7==13) || \
($4==13 && $8==13) || \
($5==13 && $9==13) || \
($6==13 && $10==13))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_7_7.txt

awk -v OFS='\t' '(($3==13 && $7==12) || \
($4==13 && $8==12) || \
($5==13 && $9==12) || \
($6==13 && $10==12))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_7_8.txt

awk -v OFS='\t' '(($3==13 && $7==11) || \
($4==13 && $8==11) || \
($5==13 && $9==11) || \
($6==13 && $10==11))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_7_9.txt

awk -v OFS='\t' '(($3==13 && $7==10) || \
($4==13 && $8==10) || \
($5==13 && $9==10) || \
($6==13 && $10==10))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_7_10.txt

awk -v OFS='\t' '(($3==13 && $7==9) || \
($4==13 && $8==9) || \
($5==13 && $9==9) || \
($6==13 && $10==9))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_7_11.txt

awk -v OFS='\t' '(($3==13 && $7==8) || \
($4==13 && $8==8) || \
($5==13 && $9==8) || \
($6==13 && $10==8))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_7_12.txt

awk -v OFS='\t' '(($3==13 && $7==7) || \
($4==13 && $8==7) || \
($5==13 && $9==7) || \
($6==13 && $10==7))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_7_13.txt

awk -v OFS='\t' '(($3==13 && $7==6) || \
($4==13 && $8==6) || \
($5==13 && $9==6) || \
($6==13 && $10==6))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_13_6.txt

awk -v OFS='\t' '(($3==13 && $7==5) || \
($4==13 && $8==5) || \
($5==13 && $9==5) || \
($6==13 && $10==5))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_13_5.txt

awk -v OFS='\t' '(($3==13 && $7==4) || \
($4==13 && $8==4) || \
($5==13 && $9==4) || \
($6==13 && $10==4))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_13_4.txt

awk -v OFS='\t' '(($3==13 && $7==3) || \
($4==13 && $8==3) || \
($5==13 && $9==3) || \
($6==13 && $10==3))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_13_3.txt

awk -v OFS='\t' '(($3==13 && $7==2) || \
($4==13 && $8==2) || \
($5==13 && $9==2) || \
($6==13 && $10==2))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_13_2.txt

awk -v OFS='\t' '(($3==13 && $7==1) || \
($4==13 && $8==1) || \
($5==13 && $9==1) || \
($6==13 && $10==1))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_13_1.txt

awk -v OFS='\t' '(($3==12 && $7==19) || \
($4==12 && $8==19) || \
($5==12 && $9==19) || \
($6==12 && $10==19))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_8_1.txt

awk -v OFS='\t' '(($3==12 && $7==18) || \
($4==12 && $8==18) || \
($5==12 && $9==18) || \
($6==12 && $10==18))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_8_2.txt

awk -v OFS='\t' '(($3==12 && $7==17) || \
($4==12 && $8==17) || \
($5==12 && $9==17) || \
($6==12 && $10==17))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_8_3.txt

awk -v OFS='\t' '(($3==12 && $7==16) || \
($4==12 && $8==16) || \
($5==12 && $9==16) || \
($6==12 && $10==16))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_8_4.txt

awk -v OFS='\t' '(($3==12 && $7==15) || \
($4==12 && $8==15) || \
($5==12 && $9==15) || \
($6==12 && $10==15))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_8_5.txt

awk -v OFS='\t' '(($3==12 && $7==14) || \
($4==12 && $8==14) || \
($5==12 && $9==14) || \
($6==12 && $10==14))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_8_6.txt

awk -v OFS='\t' '(($3==12 && $7==13) || \
($4==12 && $8==13) || \
($5==12 && $9==13) || \
($6==12 && $10==13))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_8_7.txt

awk -v OFS='\t' '(($3==12 && $7==12) || \
($4==12 && $8==12) || \
($5==12 && $9==12) || \
($6==12 && $10==12))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_8_8.txt

awk -v OFS='\t' '(($3==12 && $7==11) || \
($4==12 && $8==11) || \
($5==12 && $9==11) || \
($6==12 && $10==11))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_8_9.txt

awk -v OFS='\t' '(($3==12 && $7==10) || \
($4==12 && $8==10) || \
($5==12 && $9==10) || \
($6==12 && $10==10))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_8_10.txt

awk -v OFS='\t' '(($3==12 && $7==9) || \
($4==12 && $8==9) || \
($5==12 && $9==9) || \
($6==12 && $10==9))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_8_11.txt

awk -v OFS='\t' '(($3==12 && $7==8) || \
($4==12 && $8==8) || \
($5==12 && $9==8) || \
($6==12 && $10==8))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_8_12.txt

awk -v OFS='\t' '(($3==12 && $7==7) || \
($4==12 && $8==7) || \
($5==12 && $9==7) || \
($6==12 && $10==7))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_12_7.txt

awk -v OFS='\t' '(($3==12 && $7==6) || \
($4==12 && $8==6) || \
($5==12 && $9==6) || \
($6==12 && $10==6))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_12_6.txt

awk -v OFS='\t' '(($3==12 && $7==5) || \
($4==12 && $8==5) || \
($5==12 && $9==5) || \
($6==12 && $10==5))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_12_5.txt

awk -v OFS='\t' '(($3==12 && $7==4) || \
($4==12 && $8==4) || \
($5==12 && $9==4) || \
($6==12 && $10==4))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_12_4.txt

awk -v OFS='\t' '(($3==12 && $7==3) || \
($4==12 && $8==3) || \
($5==12 && $9==3) || \
($6==12 && $10==3))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_12_3.txt

awk -v OFS='\t' '(($3==12 && $7==2) || \
($4==12 && $8==2) || \
($5==12 && $9==2) || \
($6==12 && $10==2))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_12_2.txt

awk -v OFS='\t' '(($3==12 && $7==1) || \
($4==12 && $8==1) || \
($5==12 && $9==1) || \
($6==12 && $10==1))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_12_1.txt

awk -v OFS='\t' '(($3==11 && $7==19) || \
($4==11 && $8==19) || \
($5==11 && $9==19) || \
($6==11 && $10==19))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_9_1.txt

awk -v OFS='\t' '(($3==11 && $7==18) || \
($4==11 && $8==18) || \
($5==11 && $9==18) || \
($6==11 && $10==18))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_9_2.txt

awk -v OFS='\t' '(($3==11 && $7==17) || \
($4==11 && $8==17) || \
($5==11 && $9==17) || \
($6==11 && $10==17))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_9_3.txt

awk -v OFS='\t' '(($3==11 && $7==16) || \
($4==11 && $8==16) || \
($5==11 && $9==16) || \
($6==11 && $10==16))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_9_4.txt

awk -v OFS='\t' '(($3==11 && $7==15) || \
($4==11 && $8==15) || \
($5==11 && $9==15) || \
($6==11 && $10==15))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_9_5.txt

awk -v OFS='\t' '(($3==11 && $7==14) || \
($4==11 && $8==14) || \
($5==11 && $9==14) || \
($6==11 && $10==14))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_9_6.txt

awk -v OFS='\t' '(($3==11 && $7==13) || \
($4==11 && $8==13) || \
($5==11 && $9==13) || \
($6==11 && $10==13))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_9_7.txt

awk -v OFS='\t' '(($3==11 && $7==12) || \
($4==11 && $8==12) || \
($5==11 && $9==12) || \
($6==11 && $10==12))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_9_8.txt

awk -v OFS='\t' '(($3==11 && $7==11) || \
($4==11 && $8==11) || \
($5==11 && $9==11) || \
($6==11 && $10==11))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_9_9.txt

awk -v OFS='\t' '(($3==11 && $7==10) || \
($4==11 && $8==10) || \
($5==11 && $9==10) || \
($6==11 && $10==10))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_9_10.txt

awk -v OFS='\t' '(($3==11 && $7==9) || \
($4==11 && $8==9) || \
($5==11 && $9==9) || \
($6==11 && $10==9))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_9_11.txt

awk -v OFS='\t' '(($3==11 && $7==8) || \
($4==11 && $8==8) || \
($5==11 && $9==8) || \
($6==11 && $10==8))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_11_8.txt

awk -v OFS='\t' '(($3==11 && $7==7) || \
($4==11 && $8==7) || \
($5==11 && $9==7) || \
($6==11 && $10==7))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_11_7.txt

awk -v OFS='\t' '(($3==11 && $7==6) || \
($4==11 && $8==6) || \
($5==11 && $9==6) || \
($6==11 && $10==6))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_11_6.txt

awk -v OFS='\t' '(($3==11 && $7==5) || \
($4==11 && $8==5) || \
($5==11 && $9==5) || \
($6==11 && $10==5))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_11_5.txt

awk -v OFS='\t' '(($3==11 && $7==4) || \
($4==11 && $8==4) || \
($5==11 && $9==4) || \
($6==11 && $10==4))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_11_4.txt

awk -v OFS='\t' '(($3==11 && $7==3) || \
($4==11 && $8==3) || \
($5==11 && $9==3) || \
($6==11 && $10==3))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_11_3.txt

awk -v OFS='\t' '(($3==11 && $7==2) || \
($4==11 && $8==2) || \
($5==11 && $9==2) || \
($6==11 && $10==2))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_11_2.txt

awk -v OFS='\t' '(($3==11 && $7==1) || \
($4==11 && $8==1) || \
($5==11 && $9==1) || \
($6==11 && $10==1))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_11_1.txt

awk -v OFS='\t' '(($3==10 && $7==19) || \
($4==10 && $8==19) || \
($5==10 && $9==19) || \
($6==10 && $10==19))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_10_1.txt

awk -v OFS='\t' '(($3==10 && $7==18) || \
($4==10 && $8==18) || \
($5==10 && $9==18) || \
($6==10 && $10==18))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_10_2.txt

awk -v OFS='\t' '(($3==10 && $7==17) || \
($4==10 && $8==17) || \
($5==10 && $9==17) || \
($6==10 && $10==17))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_10_3.txt

awk -v OFS='\t' '(($3==10 && $7==16) || \
($4==10 && $8==16) || \
($5==10 && $9==16) || \
($6==10 && $10==16))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_10_4.txt

awk -v OFS='\t' '(($3==10 && $7==15) || \
($4==10 && $8==15) || \
($5==10 && $9==15) || \
($6==10 && $10==15))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_10_5.txt

awk -v OFS='\t' '(($3==10 && $7==14) || \
($4==10 && $8==14) || \
($5==10 && $9==14) || \
($6==10 && $10==14))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_10_6.txt

awk -v OFS='\t' '(($3==10 && $7==13) || \
($4==10 && $8==13) || \
($5==10 && $9==13) || \
($6==10 && $10==13))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_10_7.txt

awk -v OFS='\t' '(($3==10 && $7==12) || \
($4==10 && $8==12) || \
($5==10 && $9==12) || \
($6==10 && $10==12))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_10_8.txt

awk -v OFS='\t' '(($3==10 && $7==11) || \
($4==10 && $8==11) || \
($5==10 && $9==11) || \
($6==10 && $10==11))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_10_9.txt

awk -v OFS='\t' '(($3==10 && $7==10) || \
($4==10 && $8==10) || \
($5==10 && $9==10) || \
($6==10 && $10==10))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_10_10.txt

awk -v OFS='\t' '(($3==20 && ($7==20 || $8==20 || $9==20 || $10==20)) || \
($4==20 && ($7==20 || $8==20 || $9==20 || $10==20)) || \
($5==20 && ($7==20 || $8==20 || $9==20 || $10==20)) || \
($6==20 && ($7==20 || $8==20 || $9==20 || $10==20)))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_0_20.txt

awk -v OFS='\t' '(($3==19 && $7==20) || \
($4==19 && $8==20) || \
($5==19 && $9==20) || \
($6==19 && $10==20))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_1_0.txt

awk -v OFS='\t' '(($3==18 && $7==20) || \
($4==18 && $8==20) || \
($5==18 && $9==20) || \
($6==18 && $10==20))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_2_0.txt

awk -v OFS='\t' '(($3==17 && $7==20) || \
($4==17 && $8==20) || \
($5==17 && $9==20) || \
($6==17 && $10==20))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_3_0.txt

awk -v OFS='\t' '(($3==16 && $7==20) || \
($4==16 && $8==20) || \
($5==16 && $9==20) || \
($6==16 && $10==20))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_4_0.txt

awk -v OFS='\t' '(($3==15 && $7==20) || \
($4==15 && $8==20) || \
($5==15 && $9==20) || \
($6==15 && $10==20))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_5_0.txt

awk -v OFS='\t' '(($3==14 && $7==20) || \
($4==14 && $8==20) || \
($5==14 && $9==20) || \
($6==14 && $10==20))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_6_0.txt

awk -v OFS='\t' '(($3==13 && $7==20) || \
($4==13 && $8==20) || \
($5==13 && $9==20) || \
($6==13 && $10==20))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_7_0.txt

awk -v OFS='\t' '(($3==12 && $7==20) || \
($4==12 && $8==20) || \
($5==12 && $9==20) || \
($6==12 && $10==20))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_8_0.txt

awk -v OFS='\t' '(($3==11 && $7==20) || \
($4==11 && $8==20) || \
($5==11 && $9==20) || \
($6==11 && $10==20))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_9_0.txt

awk -v OFS='\t' '(($3==10 && $7==20) || \
($4==10 && $8==20) || \
($5==10 && $9==20) || \
($6==10 && $10==20))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_10_0.txt

awk -v OFS='\t' '(($3==9 && $7==20) || \
($4==9 && $8==20) || \
($5==9 && $9==20) || \
($6==9 && $10==20))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_11_0.txt

awk -v OFS='\t' '(($3==8 && $7==20) || \
($4==8 && $8==20) || \
($5==8 && $9==20) || \
($6==8 && $10==20))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_12_0.txt

awk -v OFS='\t' '(($3==7 && $7==20) || \
($4==7 && $8==20) || \
($5==7 && $9==20) || \
($6==7 && $10==20))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_13_0.txt

awk -v OFS='\t' '(($3==6 && $7==20) || \
($4==6 && $8==20) || \
($5==6 && $9==20) || \
($6==6 && $10==20))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_14_0.txt

awk -v OFS='\t' '(($3==5 && $7==20) || \
($4==5 && $8==20) || \
($5==5 && $9==20) || \
($6==5 && $10==20))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_15_0.txt

awk -v OFS='\t' '(($3==4 && $7==20) || \
($4==4 && $8==20) || \
($5==4 && $9==20) || \
($6==4 && $10==20))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_16_0.txt

awk -v OFS='\t' '(($3==3 && $7==20) || \
($4==3 && $8==20) || \
($5==3 && $9==20) || \
($6==3 && $10==20))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_17_0.txt

awk -v OFS='\t' '(($3==2 && $7==20) || \
($4==2 && $8==20) || \
($5==2 && $9==20) || \
($6==2 && $10==20))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_18_0.txt

awk -v OFS='\t' '(($3==1 && $7==20) || \
($4==1 && $8==20) || \
($5==1 && $9==20) || \
($6==1 && $10==20))' \
Puy_20_columns_java_only_count_keepsnps.sync | wc -l > Puy_19_0.txt

# Compile the shared site frequency spectrum using the individual files, it will look something like this in dadi format:
# The first line lists the number of chromosomes+1 in each population, then the names of the populations
# The second line lists the counts in each bin
# The third line indicates whether each bin should be masked (1 = masked, 2 = not masked)
# An example of this is provided in the file puy_shared_sfs.txt