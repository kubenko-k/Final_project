#!/bin/bash

SAMPLE="AI-70_S61"

fastqc --outdir fastqc -t 8 sequencing_files/${SAMPLE}_R1_001.fastq.gz sequencing_files/${SAMPLE}_R2_001.fastq.gz -q

trimmomatic PE -phred33 -threads 8 sequencing_files/AI-70_S61_R1_001.fastq.gz sequencing_files/AI-70_S61_R2_001.fastq.gz trimmed_reads/AI-70_S61_R1_001_paired trimmed_reads/AI-70_S61_R1_001_unpaired trimmed_reads/AI-70_S61_R2_001_paired trimmed_reads/AI-70_S61_R2_001_unpaired LEADING:15 TRAILING:15 SLIDINGWINDOW:10:20 MINLEN:20

bwa index sequencing_files/T5_sequence.fasta

bwa mem -t 8 sequencing_files/T5_sequence.fasta trimmed_reads/${SAMPLE}_R1_001_paired trimmed_reads/${SAMPLE}_R2_001_paired | samtools view -Sb>mapped/T5_${SAMPLE}_aligned.bam


samtools sort -@ 8 mapped/T5_${SAMPLE}_aligned.bam -o mapped_sorted/T5_${SAMPLE}_aligned_sorted.bam

bcftools mpileup -Ou -f Sequencing_files/T5_sequence.fasta mapped_sorted/T5_${SAMPLE}_aligned_sorted.bam | bcftools call -Ou -mv --ploidy 1 | b$


