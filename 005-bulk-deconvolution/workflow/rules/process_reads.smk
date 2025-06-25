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
# -----------------------------------------------------
rule star_index:
    input:
        fasta="results/get_genome/genome.fna",
        gtf="results/get_genome_gtf/genome.gtf",
    output:
        directory("results/star_index"),
    message:
        "Testing STAR index"
    threads: 16
    resources:
        mem=lookup(within=config, dpath="star_index/mem"),
    params:
        sjdbOverhang=lookup(within=config, dpath="star_index/sjdbOverhang"),
        extra="",
    log:
        "logs/star_index/star_index.log",
    wrapper:
        "v3.3.0/bio/star/index"

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
