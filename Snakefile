SAMPLES = ['AI-69_S60', 'AI-70_S61', 'AI-71_S62', 'AI-72_S63', 'AI-73_S64']
REFERENCES= ['T5']

rule all:
    input:
        expand(os.path.join("results", 'sorted_{reference}_{sample}_flagstat.log'),
         reference=REFERENCES,sample=SAMPLES),
        expand(os.path.join("calling", 'variants_filtered_{reference}_{sample}.vcf'),
         reference=REFERENCES,sample=SAMPLES)

rule index_ref:
 input:
    os.path.join("Sequencing_files",'{reference}_sequence.fasta')
 output:
    os.path.join("Sequencing_files",'{reference}_sequence.fasta.amb'),
    os.path.join("Sequencing_files",'{reference}_sequence.fasta.ann'),
    os.path.join("Sequencing_files",'{reference}_sequence.fasta.bwt'),
    os.path.join("Sequencing_files",'{reference}_sequence.fasta.pac'),
    os.path.join("Sequencing_files",'{reference}_sequence.fasta.sa')
 conda:
    os.path.join('envs','bwa_env.yml')
 shell:
    "bwa index {input}"

rule bwa_map:
    input:
        index_amb = os.path.join("Sequencing_files", '{reference}_sequence.fasta.amb'),
        index_ann = os.path.join("Sequencing_files", '{reference}_sequence.fasta.ann'),
        index_bwt = os.path.join("Sequencing_files", '{reference}_sequence.fasta.bwt'),
        index_pac = os.path.join("Sequencing_files", '{reference}_sequence.fasta.pac'),
        index_sa = os.path.join("Sequencing_files", '{reference}_sequence.fasta.sa'),
        fasta = os.path.join("Sequencing_files", '{reference}_sequence.fasta'),
        fastq1 = os.path.join("Sequencing_files", '{sample}_R1_001.fastq.gz'),
        fastq2 = os.path.join("Sequencing_files", '{sample}_R2_001.fastq.gz')
    output:
        bam = os.path.join("mapped", '{reference}_{sample}.bam')
    conda:
        os.path.join('envs', 'bwa_env.yml')
    threads: 8
    shell:
        "bwa mem -t {threads} {input.fasta} {input.fastq1} {input.fastq2} | samtools view -Sb > {output.bam}"

rule get_statistics:
     input:
        os.path.join("mapped_sorted",'sorted_{reference}_{sample}.bam')
     output:
        os.path.join("results",'sorted_{reference}_{sample}_flagstat.log')
     threads: 1
     conda:
        os.path.join('envs','samtools_env.yml')
     shell:
        "samtools flagstat {input} > {output}"

rule sam_sort:
    input:
        os.path.join("mapped", '{reference}_{sample}.bam')
    output:
        os.path.join("mapped_sorted", 'sorted_{reference}_{sample}.bam')
    threads: 8
    conda:
        os.path.join('envs', 'samtools_env.yml')
    shell:
        "samtools sort -@ 8 {input} -o {output}"

rule index:
    input:
     os.path.join("mapped_sorted", "{reference}_{sample}_sorted.bam")
    output:
      os.path.join("mapped_sorted","{reference}_{sample}_sorted.bam.bai")
    threads: 4
    conda:  os.path.join("envs", "samtools_env.yml")
    shell:
        "samtools index -@ {threads} {input} {output}"


rule variant_calling:
    input:
        fasta=os.path.join("Sequencing_files", '{reference}_sequence.fasta'),
        bam=os.path.join("mapped_sorted", 'sorted_{reference}_{sample}.bam'),
        index_bai=os.path.join("mapped_sorted", 'sorted_{reference}_{sample}.bam.bai')
    output:
        os.path.join("calling", 'variants_filtered_{reference}_{sample}.vcf')
    conda:
        os.path.join('envs', 'bcftools_env.yml')
    shell:
        """
        bcftools mpileup -Ou -f {input.fasta} {input.bam} | \
        bcftools call -Ou -mv --ploidy 1 | \
        bcftools filter -s LowQual > {output}
        """
