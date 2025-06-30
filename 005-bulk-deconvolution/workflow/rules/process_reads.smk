# ----------------------------------------------------- #
# EXAMPLE WORKFLOW                                      #
# ----------------------------------------------------- #


# copy genome sequence from cellranger
# -----------------------------------------------------
rule get_genome:
    input:
        fasta=config["genome_fasta"],
        fai=config["genome_fasta"] + ".fai",
    output:
        fasta="results/get_genome/genome.fa",
        fai="results/get_genome/genome.fa.fai",
    message:
        """--- Downloading genome sequence."""
    log:
        "results/get_genome/genome.log",
    shell:
        "cp {input.fasta} {output.fasta} > {log} 2>&1 && "
        "cp {input.fai} {output.fai} >> {log} 2>&1 "


# copy genome annotations from cellranger
# -----------------------------------------------------
rule get_genome_gtf:
    input:
        gtf=config["genome_gtf"],
    output:
        gtf="results/get_genome_gtf/genome.gtf",
    message:
        """--- Downloading genome annotations."""
    log:
        "results/get_genome_gtf/genome.log",
    shell:
        "cp {input.gtf} {output.gtf} > {log} 2>&1 "


# fetch dbSNP VCF file from NCBI
# -----------------------------------------------------
# rule get_dbsnp_vcf:
#     output:
#         vcf="results/get_dbsnp_vcf/genome.vcf.gz",
#     conda:
#         "../envs/get_genome.yml"
#     message:
#         """--- Downloading genome annotations."""
#     params:
#         ncbi_ftp=lookup(within=config, dpath="get_dbsnp_vcf/ncbi_ftp"),
#     log:
#         "results/get_dbsnp_vcf/get_dbsnp_vcf.log",
#     shell:
#         "wget -O {output.vcf} {params.ensembl_ftp} > {log} 2>&1 "


# index genome sequence with STAR
# source <https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/bio/star/index.html>
# Necessary because the cellranger index was generated with an old version of STAR
# -----------------------------------------------------
rule star_index:
    input:
        fasta="results/get_genome/genome.fa",
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


# helper function to get paired fastq files for a sample
# -----------------------------------------------------
def get_paired_fastq_files(wildcards):
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


# make uBAM from FASTQ files
# extra step to enter the GATK Best Practices pipeline
# -----------------------------------------------------
rule paired_fastqs_to_ubam:
    input:
        unpack(get_paired_fastq_files),
    output:
        ubam="results/paired_fastqs_to_ubam/{sample}.bam",
    message:
        """--- Running GATK FastqToSam."""
    threads: 16
    resources:
        mem=lookup(within=config, dpath="paired_fastqs_to_ubam/mem"),
        runtime=lookup(within=config, dpath="paired_fastqs_to_ubam/runtime"),
    conda:
        "../envs/gatk.yml"
    params:
        read_group=lambda wildcards, input: samples['read_group'][wildcards.sample],
        library=lambda wildcards, input: samples['library'][wildcards.sample],
    log:
        "logs/paired_fastqs_to_ubam/{sample}.log",
    shell:
        "gatk FastqToSam"
        " --FASTQ {input.fq1}"
        " --FASTQ2 {input.fq2}"
        " --OUTPUT {output.ubam}"
        " --READ_GROUP_NAME {params.read_group}"
        " --SAMPLE_NAME {wildcards.sample}"
        " --LIBRARY_NAME pooled"
        " --PLATFORM_UNIT {params.read_group}"
        " --PLATFORM illumina > {log} 2>&1"


# fetch VCF file for reference variants
# -----------------------------------------------------
rule get_reference_variants:
    output:
        vcf="results/get_reference_variants/{chr}.vcf.gz",
    message:
        """--- Downloading reference variants."""
    params:
        url=lambda wildcards: ref_vcfs["url"][wildcards.chr],
    log:
        "logs/get_reference_variants/{chr}.log",
    shell:
        "wget -O {output.vcf} {params.url} > {log} 2>&1"


# index VCF files
# necessary for merging VCF files later
# -----------------------------------------------------
rule bcftools_index_reference_variants:
    input:
        "results/get_reference_variants/{chr}.vcf.gz",
    output:
        "results/get_reference_variants/{chr}.vcf.gz.csi",
    log:
        "logs/bcftools_index_reference_variants/{chr}.log",
    params:
        extra="",  # optional parameters for bcftools index
    wrapper:
        "v7.1.0/bio/bcftools/index"


rule filter_gtf_mrna:
    input:
        gtf="results/get_genome_gtf/genome.gtf",
    output:
        bed="results/filter_gtf_mrna/mRNAs.bed",
    message:
        """--- Filter GTF to exons and UTRs."""
    threads: 1
    log:
        "logs/filter_gtf_mrna/filter_gtf_mrna.log",
    shell:
        """awk '$3 == "exon" || $3 == "UTR" {{print $1 "\\t" $4 "\\t" $5}}' {input.gtf} > {output.bed} 2> {log}"""


# filter reference VCF file for common exonic SNVs
# -----------------------------------------------------
rule filter_common_variants:
    input:
        "results/get_reference_variants/{chr}.vcf.gz",
        index="results/get_reference_variants/{chr}.vcf.gz.csi",
        regions="results/filter_gtf_mrna/mRNAs.bed",
    output:
        "results/filter_common_variants/{chr}.vcf.gz",
    message:
        """--- Filtering common variants."""
    log:
        "logs/filter_common_variants/{chr}.log",
    params:
        extra="--exclude 'AF<0.01' --exclude-types indels", # TODO: add filter for exonic variants
    threads: 8,
    resources:
        mem=lookup(within=config, dpath="bcftools_call/mem"),
        runtime=lookup(within=config, dpath="bcftools_call/runtime"),
    wrapper:
        "v3.7.0/bio/bcftools/view"


rule bcftools_concat:
    input:
        calls=expand("results/filter_common_variants/{chr}.vcf.gz", chr=ref_vcfs.index.unique()),
    output:
        "results/bcftools_concat/common_snvs.vcf.gz",
    log:
        "logs/bcftools_concat/bcftools_concat.log",
    params:
        uncompressed_bcf=False,
        extra="",  # optional parameters for bcftools concat (except -o)
    threads: 4
    resources:
        mem_mb=10,
    wrapper:
        "v7.1.0/bio/bcftools/concat"


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
        unpack(get_paired_fastq_files),
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


# create sequence dictionary for genome
# needed by some GATK tools
# -----------------------------------------------------
rule create_sequence_dictionary:
    input:
        "results/get_genome/genome.fa",
    output:
        "results/get_genome/genome.dict",
    conda:
        "../envs/gatk.yml"
    message:
        """--- Creating sequence dictionary for genome."""
    log:
        "logs/create_sequence_dictionary.log",
    threads: 8
    resources:
        mem=lookup(within=config, dpath="create_sequence_dictionary/mem"),
        runtime=lookup(within=config, dpath="create_sequence_dictionary/runtime"),
    shell:
        "gatk CreateSequenceDictionary"
        " -R {input}"
        " -O {output} > {log} 2>&1"


# merge reads mapped using STAR with unmapped BAM files
# -----------------------------------------------------
rule merge_bam_alignment:
    input:
        ubam="results/paired_fastqs_to_ubam/{sample}.bam",
        bam="results/star_pe/{sample}/pe_aligned.bam",
        fasta="results/get_genome/genome.fa",
        dict="results/get_genome/genome.dict",
    output:
        bam="results/merge_bam_alignment/{sample}.bam",
    message:
        """--- Running GATK MergeBamAlignment."""
    threads: 16
    resources:
        mem=lookup(within=config, dpath="merge_bam_alignment/mem"),
        runtime=lookup(within=config, dpath="merge_bam_alignment/runtime"),
    conda:
        "../envs/gatk.yml"
    params:
        read_group=lambda wildcards, input: samples['read_group'][wildcards.sample],
        library=lambda wildcards, input: samples['library'][wildcards.sample],
    log:
        "logs/merge_bam_alignment/{sample}.log",
    shell:
        "gatk MergeBamAlignment"
        " --REFERENCE_SEQUENCE {input.fasta}"
        " --UNMAPPED_BAM {input.ubam}"
        " --ALIGNED_BAM {input.bam}"
        " --OUTPUT {output}"
        " --INCLUDE_SECONDARY_ALIGNMENTS false"
        " --VALIDATION_STRINGENCY SILENT > {log} 2>&1"


# mark duplicates with picard tools
# source <https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/bio/gatk/markduplicatesspark.html>
# -----------------------------------------------------
# rule mark_duplicates_spark:
#     input:
#         "results/merge_bam_alignment/{sample}.bam",
#     output:
#         bam="results/mark_duplicates_spark/{sample}.bam",
#         metrics="results/mark_duplicates_spark/{sample}.metrics",
#     message:
#         """--- Running GATK MarkDuplicates."""
#     log:
#         "logs/mark_duplicates_spark/{sample}.log",
#     params:
#         extra="",  # optional
#         java_opts="",  # optional
#         #spark_runner="",  # optional, local by default
#         #spark_v7.1.0="",  # optional
#         #spark_extra="", # optional
#     resources:
#         # Memory needs to be at least 471859200 for Spark, so 589824000 when
#         # accounting for default JVM overhead of 20%. We round round to 650M.
#         mem=lookup(within=config, dpath="mark_duplicates_spark/mem"),
#         runtime=lookup(within=config, dpath="mark_duplicates_spark/runtime"),
#     threads: 8
#     wrapper:
#         "v7.1.0/bio/gatk/markduplicatesspark"


## not needed because cellranger already indexed the genome
# rule samtools_faidx:
#     input:
#         fa="results/get_genome/genome.fa",
#     output:
#         bai="results/get_genome/genome.fa.fai",
#     message:
#         """--- Running samtools faidx."""
#     log:
#         "logs/samtools_faidx/genome.log",
#     resources:
#         mem=lookup(within=config, dpath="samtools_faidx/mem"),
#         runtime=lookup(within=config, dpath="samtools_faidx/runtime"),
#     conda:
#         "../envs/samtools.yml"
#     threads: 8
#     shell:
#         "samtools faidx {input.fa} > {log} 2>&1"


# split reads with N CIGAR operations
# source <https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/bio/gatk/splitncigarreads.html>
# -----------------------------------------------------
rule splitncigarreads:
    input:
        bam="results/merge_bam_alignment/{sample}.bam",
        ref="results/get_genome/genome.fa",
        dict="results/get_genome/genome.dict",
    output:
        "results/splitncigarreads/{sample}.bam",
    message:
        """--- Running GATK SplitNCigarReads."""
    log:
        "logs/splitncigarreads/{sample}.log",
    params:
        extra="",  # optional
        java_opts="",  # optional
    resources:
        # Memory needs to be at least 471859200 for Spark, so 589824000 when
        # accounting for default JVM overhead of 20%. We round round to 650M.
        mem=lookup(within=config, dpath="splitncigarreads/mem"),
        runtime=lookup(within=config, dpath="splitncigarreads/runtime"),
    wrapper:
        "v7.1.0/bio/gatk/splitncigarreads"


# Run BCFtools mpileup
# source <https://samtools.github.io/bcftools/howtos/variant-calling.html>
# -----------------------------------------------------
rule bcftools_mpileup:
    input:
        alignments=["results/splitncigarreads/{sample}.bam",],
        ref="results/get_genome/genome.fa",
        index="results/get_genome/genome.fa.fai",
    output:
        pileup="results/bcftools_mpileup/{sample}.pileup.bcf",
    params:
        uncompressed_bcf=False,
        extra="--max-depth 500 --min-BQ 20",
    log:
        "logs/bcftools_mpileup/{sample}.log",
    threads: 8,
    message:
        """--- Running BCFtools mpileup."""
    resources:
        mem=lookup(within=config, dpath="bcftools_mpileup/mem"),
        runtime=lookup(within=config, dpath="bcftools_mpileup/runtime"),
    wrapper:
        "v7.1.0/bio/bcftools/mpileup"


# Run BCFtools call
# source <https://samtools.github.io/bcftools/howtos/variant-calling.html>
# -----------------------------------------------------
rule bcftools_call:
    input:
        pileup="results/bcftools_mpileup/{sample}.pileup.bcf",
    output:
        calls="results/bcftools_call/{sample}.calls.bcf",
    params:
        uncompressed_bcf=False,
        caller="-m",  # valid options include -c/--consensus-caller or -m/--multiallelic-caller
        extra="-v -Ob",
    log:
        "logs/bcftools_call/{sample}.log",
    threads: 8,
    resources:
        mem=lookup(within=config, dpath="bcftools_call/mem"),
        runtime=lookup(within=config, dpath="bcftools_call/runtime"),
    wrapper:
        "v7.1.0/bio/bcftools/call"


# filter BCF calls
# source <https://samtools.github.io/bcftools/howtos/variant-calling.html>
# - quality >= 30
# - exclude indels
# -----------------------------------------------------
rule bcftools_view:
    input:
        "results/bcftools_call/{sample}.calls.bcf",
    output:
        "results/bcftools_view/{sample}.calls.filtered.vcf.gz",
    log:
        "logs/bcftools_view/{sample}.log",
    params:
        extra="--include 'QUAL>=30' --exclude-types indels",
    threads: 8,
    resources:
        mem=lookup(within=config, dpath="bcftools_view/mem"),
        runtime=lookup(within=config, dpath="bcftools_view/runtime"),
    wrapper:
        "v7.1.0/bio/bcftools/view"


# index VCF files
# necessary for merging VCF files later
# -----------------------------------------------------
rule bcftools_index_calls:
    input:
        "results/bcftools_view/{sample}.calls.filtered.vcf.gz",
    output:
        "results/bcftools_view/{sample}.calls.filtered.vcf.gz.csi",
    log:
        "logs/bcftools_index_calls/{sample}.log",
    params:
        extra="",  # optional parameters for bcftools index
    wrapper:
        "v7.1.0/bio/bcftools/index"


# merge VCF files
# prepares the file for deconvolution by popscle demuxlet
# -----------------------------------------------------
rule bcftools_merge:
    input:
        calls=expand("results/bcftools_view/{sample}.calls.filtered.vcf.gz", sample=samples['sample'].unique()),
        idx=expand("results/bcftools_view/{sample}.calls.filtered.vcf.gz.csi", sample=samples['sample'].unique()),
    output:
        "results/bcftools_merge/all.vcf.gz",
    log:
        "logs/bcftools_merge/all.log",
    params:
        extra="",  # optional parameters for bcftools concat (except -o)
    threads: 8,
    resources:
        mem=lookup(within=config, dpath="bcftools_merge/mem"),
        runtime=lookup(within=config, dpath="bcftools_merge/runtime"),
    wrapper:
        "v7.1.0/bio/bcftools/merge"


# first step of deconvolution with popscle
# -----------------------------------------------------
rule popscle_dsc:
    input:
        vcf="results/bcftools_concat/common_snvs.vcf.gz",
    output:
        pileup="results/popscle_dsc/{pool}.pileup",
    params:
        bam=lambda wildcards: pooled_bams['bam'][wildcards.pool],
    conda:
        "../envs/popscle.yml"
    message:
        """--- Running popscle dsc-pileup."""
    log:
        "logs/popscle_dsc/{pool}.log",
    threads: 1
    resources:
        mem=lookup(within=config, dpath="popscle_dsc/mem"),
        runtime=lookup(within=config, dpath="popscle_dsc/runtime"),
    shell:
        "popscle dsc-pileup"
        " --sam {params.bam}"
        " --vcf {input.vcf}"
        " --out {output.pileup} > {log} 2>&1 &&"
        " touch {output.pileup}"


# second step of deconvolution with popscle
# -----------------------------------------------------
rule popscle_demuxlet:
    input:
        pileup="results/popscle_dsc/{pool}.pileup",
        vcf="results/bcftools_merge/all.vcf.gz",
    output:
        pileup="results/popscle_demuxlet/{pool}",
    conda:
        "../envs/popscle.yml"
    message:
        """--- Running popscle demuxlet."""
    log:
        "logs/popscle_demuxlet/{pool}.log",
    threads: 1
    resources:
        mem=lookup(within=config, dpath="popscle_demuxlet/mem"),
        runtime=lookup(within=config, dpath="popscle_demuxlet/runtime"),
    shell:
        "popscle demuxlet"
        " --plp {input.pileup}"
        " --vcf {input.vcf}"
        " --field GT"
        " --out {output.pileup} > {log} 2>&1 &&"
        " touch {output.pileup}"

# rule gatk_baserecalibratorspark:
#     input:
#         bam="results/splitncigarreads/{sample}.bam",
#         ref="results/get_genome/genome.fa",
#         dict="results/create_sequence_dictionary/genome.dict",
#         known="results/get_dbsnp_vcf/genome.vcf.gz",
#     output:
#         recal_table="results/gatk_baserecalibratorspark/{sample}.grp",
#     log:
#         "logs/gatk_baserecalibratorspark/{sample}.log",
#     params:
#         extra="--use-original-qualities",  # optional
#         java_opts="",  # optional
#         #spark_runner="",  # optional, local by default
#         #spark_v7.1.0="",  # optional
#         #spark_extra="", # optional
#     resources:
#         mem=lookup(within=config, dpath="gatk_baserecalibratorspark/mem"),
#         runtime=lookup(within=config, dpath="gatk_baserecalibratorspark/runtime"),
#     wrapper:
#         "v7.1.0/bio/gatk/baserecalibratorspark"


# merge reads mapped using STAR with unmapped BAM files
# see <https://github.com/statgen/popscle?tab=readme-ov-file>
# -----------------------------------------------------
# rule dsc_pileup:
#     input:
#         bam=TODO,
#         vcf=TODO
#     output:
#         pileup=TODO
#     message:
#         """--- Running popscle dsc-pileup."""
#     threads: TODO
#     resources:
#         mem=lookup(within=config, dpath="dsc_pileup/mem"),
#         runtime=lookup(within=config, dpath="dsc_pileup/runtime"),
#     conda:
#         "../envs/gatk.yml"
#     params:
#         read_group=lambda wildcards, input: samples['read_group'][wildcards.sample],
#         library=lambda wildcards, input: samples['library'][wildcards.sample],
#     log:
#         "logs/dsc_pileup/{sample}.log",
#     shell:
#         "popscle dsc-pileup"
#         " --sam /data/$bam"
#         " --vcf /data/$ref_vcf"
#         " --out /data/$pileup > {log} 2>&1"


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
