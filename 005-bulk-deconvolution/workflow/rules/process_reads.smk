# ----------------------------------------------------- #
# EXAMPLE WORKFLOW                                      #
# ----------------------------------------------------- #


# fetch genome sequence from Ensembl
# -----------------------------------------------------
rule get_genome:
    output:
        fasta="results/get_genome/genome.fna",
    conda:
        "../envs/get_genome.yml"
    message:
        """--- Downloading genome sequence."""
    params:
        ensembl_ftp=lookup(within=config, dpath="get_genome/ensembl_ftp"),
    log:
        "results/get_genome/genome.log",
    shell:
        "wget -O results/get_genome/genome.fna.gz {params.ensembl_ftp} > {log} 2>&1 && "
        "gunzip results/get_genome/genome.fna.gz >> {log} 2>&1"


# fetch genome annotations from Ensembl
# -----------------------------------------------------
rule get_genome_gtf:
    output:
        gtf="results/get_genome_gtf/genome.gtf",
    conda:
        "../envs/get_genome.yml"
    message:
        """--- Downloading genome annotations."""
    params:
        ensembl_ftp=lookup(within=config, dpath="get_genome_gtf/ensembl_ftp"),
    log:
        "results/get_genome_gtf/genome.log",
    shell:
        "wget -O results/get_genome_gtf/genome.gtf.gz {params.ensembl_ftp} > {log} 2>&1 && "
        "gunzip results/get_genome_gtf/genome.gtf.gz >> {log} 2>&1"

# index genome sequence with STAR
# source <https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/bio/star/index.html>
# -----------------------------------------------------
rule star_index:
    input:
        fasta="results/get_genome/genome.fna",
        gtf="results/get_genome_gtf/genome.gtf",
    output:
        directory("results/star_index"),
    message:
        """--- Building STAR index."""
    threads: 16
    resources:
        mem=lookup(within=config, dpath="star_index/mem"),
        runtime=lookup(within=config, dpath="star_index/runtime"),
    params:
        sjdbOverhang=lookup(within=config, dpath="star_index/sjdbOverhang"),
        extra="",
    log:
        "logs/star_index/star_index.log",
    wrapper:
        "v7.1.0/bio/star/index"


def get_star_bam_files(wildcards):
    """
    Get all STAR BAM files for the given sample wildcard.
    """
    sample_data=samples[samples['sample'] == wildcards.sample].filter(items=['directory', 'read1', 'read2'])
    fq1=sample_data['read1'].to_list()
    fq2=sample_data['read2'].to_list()
    return {
        'fq1': fq1,
        'fq2': fq2,
    }


# map reads with STAR
# source <https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/bio/star/align.html>
# -----------------------------------------------------
rule star_pe:
    input:
        # use a list for multiple fastq files for one sample
        # usually technical replicates across lanes/flowcells
        #fq1=["reads/{sample}_R1.1.fastq", "reads/{sample}_R1.2.fastq"],
        # paired end reads needs to be ordered so each item in the two lists match
        #fq2=["reads/{sample}_R2.1.fastq", "reads/{sample}_R2.2.fastq"],  #optional
        unpack(get_star_bam_files),
        # path to STAR reference genome index
        idx="results/star_index",
    output:
        # see STAR manual for additional output files
        aln="results/star_pe/{sample}/pe_aligned.bam",
        log="results/star_pe/{sample}/Log.out",
        sj="results/star_pe/{sample}/SJ.out.tab",
        unmapped=["results/star_pe/{sample}/unmapped.1.fastq.gz","results/star_pe/{sample}/unmapped.2.fastq.gz"],
    message:
        """--- Mapping reads using STAR."""
    log:
        "logs/star_pe/{sample}.log",
    params:
        # optional parameters
        extra="--outSAMtype BAM SortedByCoordinate --twopassMode Basic",
    threads: 8
    resources:
        mem=lookup(within=config, dpath="star_pe/mem"),
        runtime=lookup(within=config, dpath="star_pe/runtime"),
    wrapper:
        "v7.1.0/bio/star/align"


# mark duplicates with GATK
# -----------------------------------------------------
rule mark_duplicate:
    input:
        bam="results/star_pe/{sample}/pe_aligned.bam",
    output:
        bam="results/mark_duplicate/{sample}.bam",
        metrics="results/mark_duplicate/{sample}.metrics",
    conda:
        "../envs/gatk.yml"
    message:
        """--- Running GATK MarkDuplicates."""
    log:
        "logs/mark_duplicate/{sample}.log",
    shell:
        "gatk MarkDuplicates"
        " --INPUT {input.bam}"
        " --OUTPUT {output.bam}"
        " --CREATE_INDEX true"
        " --VALIDATION_STRINGENCY SILENT"
        " --METRICS_FILE {output.metrics} > {log} 2>&1 "


# mark duplicates with picard tools
# source <https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/bio/gatk/markduplicatesspark.html>
# -----------------------------------------------------
rule mark_duplicates_spark:
    input:
        "results/star_pe/{sample}/pe_aligned.bam",
    output:
        bam="results/mark_duplicates_spark/{sample}.bam",
        metrics="results/mark_duplicates_spark/{sample}.metrics",
    message:
        """--- Running GATK MarkDuplicates."""
    log:
        "logs/mark_duplicates_spark/{sample}.log",
    params:
        extra="",  # optional
        java_opts="",  # optional
        #spark_runner="",  # optional, local by default
        #spark_v7.1.0="",  # optional
        #spark_extra="", # optional
    resources:
        # Memory needs to be at least 471859200 for Spark, so 589824000 when
        # accounting for default JVM overhead of 20%. We round round to 650M.
        mem_mb=lambda wildcards, input: max([input.size_mb * 0.25, 650]),
    threads: 8
    wrapper:
        "v7.1.0/bio/gatk/markduplicatesspark"

# rule mark_duplicate:
#     input:
#         bam="results/star_pe/{sample}/pe_aligned.bam",
#     output:
#         bam="results/mark_duplicate/{sample}.bam",
#         metrics="results/mark_duplicate/{sample}.metrics",
#     conda:
#         "../envs/gatk.yml"
#     message:
#         """--- Running GATK MarkDuplicates."""
#     log:
#         "logs/mark_duplicate/{sample}.log",
#     shell:
#         "gatk MarkDuplicates"
#         " --INPUT {input.bam}"
#         " --OUTPUT {output.bam}"
#         " --CREATE_INDEX true"
#         " --VALIDATION_STRINGENCY SILENT"
#         " --METRICS_FILE {output.metrics} > {log} 2>&1 "

# validate genome sequence file
# -----------------------------------------------------
# rule validate_genome:
#     input:
#         fasta=rules.get_genome.output.fasta,
#     output:
#         fasta="results/validate_genome/genome.fna",
#     conda:
#         "../envs/validate_genome.yml"
#     message:
#         """--- Validating genome sequence file."""
#     log:
#         "results/validate_genome/genome.log",
#     script:
#         "../scripts/validate_fasta.py"


# simulate read data using DWGSIM
# -----------------------------------------------------
# rule simulate_reads:
#     input:
#         fasta=rules.validate_genome.output.fasta,
#     output:
#         multiext(
#             "results/simulate_reads/{sample}",
#             read1=".bwa.read1.fastq.gz",
#             read2=".bwa.read2.fastq.gz",
#         ),
#     conda:
#         "../envs/simulate_reads.yml"
#     message:
#         """--- Simulating read data with DWGSIM."""
#     params:
#         output_type=1,
#         read_length=lookup(within=config, dpath="simulate_reads/read_length"),
#         read_number=lookup(within=config, dpath="simulate_reads/read_number"),
#     log:
#         "results/simulate_reads/{sample}.log",
#     shell:
#         "output_prefix=`echo {output.read1} | cut -f 1 -d .`;"
#         "dwgsim "
#         " -1 {params.read_length}"
#         " -2 {params.read_length}"
#         " -N {params.read_number}"
#         " -o {params.output_type}"
#         " {input.fasta}"
#         " ${{output_prefix}}"
#         " > {log} 2>&1"


# make QC report
# -----------------------------------------------------
# rule fastqc:
#     input:
#         fastq="results/simulate_reads/{sample}.bwa.{read}.fastq.gz",
#     output:
#         html="results/fastqc/{sample}.bwa.{read}_fastqc.html",
#         zip="results/fastqc/{sample}.bwa.{read}_fastqc.zip",
#     params:
#         extra="--quiet",
#     message:
#         """--- Checking fastq files with FastQC."""
#     log:
#         "results/fastqc/{sample}.bwa.{read}.log",
#     threads: 1
#     wrapper:
#         "v6.0.0/bio/fastqc"


# run multiQC on tool output
# -----------------------------------------------------
# rule multiqc:
#     input:
#         expand(
#             "results/fastqc/{sample}.bwa.{read}_fastqc.{ext}",
#             sample=samples.index,
#             read=["read1", "read2"],
#             ext=["html", "zip"],
#         ),
#     output:
#         report="results/multiqc/multiqc_report.html",
#     params:
#         extra="--verbose --dirs",
#     message:
#         """--- Generating MultiQC report for seq data."""
#     log:
#         "results/multiqc/multiqc.log",
#     wrapper:
#         "v6.0.0/bio/multiqc"
