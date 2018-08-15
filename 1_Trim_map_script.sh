################################################################################################
#
#           1. Read trimming and mapping                      
#
#
################################################################################################

# The general structure of this script could be used for either RNA-seq or DNA-seq data
# First part of the script locates all fastq files from current directory to all subfolders and then runs FastQC on them. Note, the output directory is hardcoded so you might want to change it
####################################################################################################################################

# Move into directory that contains raw files
cd raw_files

# Find all files ending in .fastq.gz from current directory and all subfolders. Output full path (`pwd` needed otherwise you'd only get the path from the current directory). -type specifies a normal file. -name allows us to look for things based on file name.
# find `pwd` -type f -name *.fastq.gz
# Saves names of directories to an array
fileList=($(find `pwd` -type f -name '*.fastq'))

# Loops through all elements in the array. Note that "i" holds the value for the element in the array during each iteration. In this case "i" holds the path to the fastq file. Path to output directory is fixed. Runs all FastQC on all fastq files located in array.
for i in ${fileList[@]}
do

fastqc --outdir fastqc_dir --noextract --nogroup $i

done

# Runs Trim Galore on the files specified in the "TrimGalore_fastq_Locations.txt"

# The for loop is hard coded to the number to lines in the "TrimGalore_fastq_Locations.txt" file
# "cat" opens the file
# "sed" is grabbing the ith line of the file
# "awk is grabbing the 1st or 2nd column from the line above

for i in {1..23}
 
do

# TrimGalore_fastq_Locations_psulph.txt has two columns
# Column one contains the paths to raw read one fastq files
# Column two contains the paths to raw read two fastq files

read1=$(cat TrimGalore_fastq_Locations.txt | sed -n ''$i'p' | awk '{print $1}')
read2=$(cat TrimGalore_fastq_Locations.txt | sed -n ''$i'p' | awk '{print $2}')
 
trim_galore --quality 20 --fastqc_args "--noextract --nogroup" --stringency 6 --gzip --length 50 --paired $read1 $read2
 
done

# Index reference genome
bwa index <path to reference genome>
samtools faidx <path to reference genome>

# Map to reference genome using BWA-mem
# Use same loop structure seen above

for i in {1..23}

do

# trimmed_fastq_Locations.txt has three columns
# Column one contains the paths to trimmed read one fastq files
# Column two contains the paths to trimmed read two fastq files
# Column three contains the name of each sample

read1=$(cat trimmed_fastq_Locations.txt | sed -n ''$i'p' | awk '{print $1}')
read2=$(cat trimmed_fastq_Locations.txt | sed -n ''$i'p' | awk '{print $2}')
filename=$(cat trimmed_fastq_Locations.txt | sed -n ''$i'p' | awk '{print $3}')

bwa mem <path to reference genome> $read1 $read2 > ${filename}.sam

done


# Process the sam files after mapping to the reference genome

for i in {1..23}

do

# sam_locations.txt has five columns
# Column one contains the paths to sam files
# Column two contains the sample names
# Column three contains the readgroup information
# Column four contains the lanenumber
# Column five contains the population

samfile=$(cat sam_locations.txt | sed -n ''$i'p' | awk '{print $1}')
samplename=$(cat sam_locations.txt | sed -n ''$i'p' | awk '{print $2}')
readgroup=$(cat sam_locations.txt | sed -n ''$i'p' | awk '{print $3}')
population=$(cat sam_locations.txt | sed -n ''$i'p' | awk '{print $4}')

# First, convert to bam file to save space
# -u = Output uncompressed BAM, -b = Output in the BAM format, -h = Include the header in the output, -t = links to .fai index
samtools view -ubht <path to reference>.fai $samfile -o ${samplename}.bam

# CleanSam: Cleans the provided SAM/BAM, soft-clipping beyond-end-of-reference alignments and setting MAPQ to 0 for unmapped reads
java -Xmx6G -jar picard.jar CleanSam INPUT=${samplename}.bam OUTPUT=${samplename}.clean.bam VALIDATION_STRINGENCY=SILENT

# FixMateInformation: Verify mate-pair information between mates and fix if needed
java -Xmx6G -jar picard.jar FixMateInformation INPUT=${samplename}.clean.bam OUTPUT=${samplename}.clean.fix.bam VALIDATION_STRINGENCY=SILENT

# ValidateSamFile: This tool reports on the validity of a SAM or BAM file relative to the SAM format specification
java -Xmx6G -jar picard.jar ValidateSamFile INPUT=${samplename}.clean.fix.bam OUTPUT=${samplename}.validation MODE=SUMMARY VALIDATION_STRINGENCY=LENIENT

# SortSam: Sort the bam file by coordinate
java -Xmx4G -jar picard.jar SortSam INPUT=${samplename}.clean.fix.bam OUTPUT=${samplename}.F.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT

# Index the bam file
samtools index ${samplename}.F.bam

# Check alignment statistics using bamtools
bamtools stats -in ${samplename}.F.bam > ${samplename}.STATS

# AddOrReplaceReadGroups: Add readgroups to the bam files for SNP calling
java -Xmx6G -jar /home/repository_software/picard-tools-1.133/picard.jar AddOrReplaceReadGroups INPUT=${samplename}.F.bam OUTPUT=${samplename}.F.readgroup.bam SORT_ORDER=coordinate RGPL=Illumina RGLB=$samplename RGSM=$readgroup RGPU=$population

done