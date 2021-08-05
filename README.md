# CoV Chimeras
Scripts that identify and quantify inter-virus recombination from RNA-seq and Nanopore datasets
Note: This pipeline is unpublished and currently under construction.

## Before you start:

1) RNA-seq datasets should be aligned to viral genomes using ViReMa. You will have to make your own FASTA file containing the genomes of interest.

Example workflow:

```
input="./samples.txt"
while IFS= read -r line
do
	python ./ViReMa_0.21/ViReMa.py  virus_genomes_index ${line}_virema.fastq ${line}_virema.sam --p 32 --Output_Tag ${line}_virema -FuzzEntry --Defuzz 0 --MicroInDel_Length 5 --Output_Dir ${line}_virema -BED
	samtools view -b -@ 32 ${line}_virema.sam > ${line}_virema.bam
	samtools sort -@ 32 -o ${line}_virema.sort.bam ${line}_virema.bam
	samtools index -@ 32 -b ${line}_virema.sort.bam ${line}_virema.sort.bam.bai
	./bbmap/pileup.sh in=${line}_virema.sam basecov=${line}_virema_coverage.txt delcoverage=f 32bit=t -Xmx64g
done < "$input"
```

{sample}_virema_Recombination_Results.txt must be manually parsed to isolate virus1 --> virus2 and virus2 -->virus1 recombination junctions. Future releases will contain a script/functionality for automatic parsing. Save files as {sample}_virus1_to_virus2.txt and {sample}_virus2_to_virus1.txt

2) Nanopore datasets should also be aligned to the genomes but this time using minimap2. Please download FLAIR to utilize its junction detection abilities. See https://github.com/BrooksLabUCSC/flair

Example workflow:

```
touch chimeric_recombination_stats.txt
input="samples.txt"
while IFS= read -r line
do
	cd target_directory/
	python flair-master/flair.py align -m ./minimap2/minimap2 -o ${line}_virus1 -t 32 -v1.3 -r ${line}.fastq -g virus1_genome.fasta
	samtools view -b -@ 32 ${line}_virus1.sam > ${line}_virus1.bam
	samtools sort -@ 32 -o ${line}_virus1.sort.bam ${line}_virus1.bam
	samtools index -@ 32 -b ${line}_virus1.sort.bam ${line}_virus1.sort.bam.bai
	python /home/denison-thelio/flair-master/bin/bam2Bed12.py -i ${line}_virus1.sort.bam > ${line}_virus1.bed
	./bbmap/pileup.sh in=${line}_virus1.sam basecov=${line}_virus1_coverage.txt delcoverage=f 32bit=t -Xmx64g
	echo "Statistics for " ${line} " Virus1 alignment:" >> chimeric_recombination_stats.txt
	samtools idxstats ${line}_virus1.sort.bam >> chimeric_recombination=_stats.txt
	python flair-master/flair.py align -m ./minimap2/minimap2 -o ${line}_virus2 -t 32 -v1.3 -r ${line}.fastq -g virus2_genome.fasta	
  samtools view -b -@ 32 ${line}_virus2.sam > ${line}_virus2.bam
	samtools sort -@ 32 -o ${line}_virus2.sort.bam ${line}_virus2.bam
	samtools index -@ 32 -b ${line}_virus2.sort.bam ${line}_virus2.sort.bam.bai
	python ./flair-master/bin/bam2Bed12.py -i ${line}_virus2.sort.bam > ${line}_virus2.bed
	./bbmap/pileup.sh in=${line}_virus2.sam basecov=${line}_virus2_coverage.txt delcoverage=f 32bit=t -Xmx64g
	echo "Statistics for " ${line} " Virus2 alignment:" >> chimeric_recombination_stats.txt
	samtools idxstats ${line}_virus2.sort.bam >> chimeric_recombination_stats.txt
	cd ~
done < "$input"
```

## R script usage (short-read Illlumina RNA-seq datasets)
**Junction_Pattern_Plots.R**
Requires editing of scripts in RStudio. Future releases will be in Python and command-line executable.
Required packages:
ggplot2

1) edit input data:
```
data <- read.table("./{sample}_MERS_to_SARS2_junctions.txt", header = TRUE)
```
2) edit output names:

```
ggsave(filename = "{sample}_virus1_to_virus2_junctions.tiff", plot = last_plot(), device = "tiff", path = "path/to/Junction_Plots/", scale = 1, width = 6, height = 4, units = "in", dpi = 1200, limitsize = TRUE)
```

## Python script usage (direct RNA Nanopore sequencing datasets)
Section under construction.
