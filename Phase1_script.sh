# Create the environment (software version not specified as using Mac)
conda create -n vc -c bioconda -c conda-forge fastqc bwa samtools freebayes bcftools

# Activate the environment
conda activate vc

# Create an output directory for the reports
mkdir quality

# Run FastQC on all your FASTQ files
fastqc --outdir quality --noextract --nogroup --format fastq --threads 2 --quiet *.fastq

# Index the reference sequence
bwa index GeneX_reference.fa

# Map reads, convert to BAM, and sort in one pipeline
bwa mem -t 2 -R '@RG\tID:hlzon\tSM:hlzon\tPL:ILLUMINA\tLB:WGS\tPU:HiSeq2500' GeneX_reference.fa mutant_GeneX_R1.fastq mutant_GeneX_R2.fastq | \
samtools view -b -@ 2 | \
samtools sort -@ 2 -o hlzon_mutant_sorted.bam -

# Index the BAM file for later visualisation
samtools index -@ 2 hlzon_mutant_sorted.bam

# Call variants
freebayes --ploidy 1 --min-mapping-quality 30 --min-base-quality 25 \
--fasta-reference GeneX_reference.fa hlzon_mutant_sorted.bam > hlzon_raw_variants.vcf

# Filter the VCF file for high-quality variants
bcftools view -e 'QUAL < 10 || INFO/DP < 10 || AF = 0' -O v \
-o hlzon_filtered_variants.vcf hlzon_raw_variants.vcf

# Compress the filtered VCF file
bgzip hlzon_filtered_variants.vcf

# Index the compressed VCF file
bcftools index hlzon_filtered_variants.vcf.gz

# Construct the consensus sequence
bcftools consensus -s hlzon -f GeneX_reference.fa -o hlzon_consensus.fa hlzon_filtered_variants.vcf.gz

# Update the FASTA header to reflect the mutant genome
sed -i '' 's/>.*/>Mutant_GeneX/' hlzon_consensus.fa

# Without unzipping the file, grep the lines that is not starting with # for mutation
zgrep -v "^#" hlzon_filtered_variants.vcf.gz
