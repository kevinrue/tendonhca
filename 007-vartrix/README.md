## Vartrix

GitHub repository <https://github.com/10XGenomics/vartrix>

## Installation

Documentation: <https://github.com/10XGenomics/vartrix?tab=readme-ov-file#installation>

Steps:

Download [vartrix_linux](https://github.com/10XGenomics/vartrix/releases/download/v1.1.22/vartrix_linux) v1.1.22

Run `chmod u+x vartrix_linux`.

Run `./vartrix_linux --help`.

The help is displayed. Success!

## Usage

Documentation: <https://github.com/10XGenomics/vartrix?tab=readme-ov-file#usage>

```bash
./vartrix_linux -v <path_to_input_vcf> -b <path_to_cellranger_bam> -f <path_to_fasta_file> -c <path_to_cell_barcodes_file> -o <path_for_output_matrix>
```

Example files:

path_to_input_vcf: `/ceph/project/tendonhca/albrecht/resources/agne_poorva_vcf/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf`

path_to_cellranger_bam: `/project/tendonhca/albrecht/004-deconvolution/results/cellranger_count/12G/outs/possorted_genome_bam.bam`

path_to_fasta_file: `/project/tendonhca/albrecht/004-deconvolution/resources/refdata-gex-GRCh38-2020-A/fasta/genome.fa`

path_to_cell_barcodes_file: `/project/tendonhca/albrecht/004-deconvolution/results/cellranger_count/12G/outs/filtered_feature_bc_matrix/barcodes.tsv.gz`

path_for_output_matrix: `results/vartrix/12G.matrix.mtx`

## Snakemake

```bash
snakemake
```