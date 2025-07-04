# import basic packages
import pandas as pd
from snakemake.utils import validate


# read sample sheet
samples = (
    pd.read_csv(config["samplesheet"], sep="\t", dtype={"sample": str})
    .set_index("sample", drop=False)
    .sort_index()
)

# read reference VCFs sheet
ref_vcfs = (
    pd.read_csv(config["ref_vcf_sheet"], sep="\t", dtype={"chr": str})
    .set_index("chr", drop=False)
    .sort_index()
)

# read pooled BAMs sheet
pooled_bams = (
    pd.read_csv(config["poolsheet"], sep="\t", dtype={"id": str})
    .set_index("id", drop=False)
    .sort_index()
)

# validate sample sheet and config file
validate(samples, schema="../../config/schemas/samples.schema.yml")
validate(config, schema="../../config/schemas/config.schema.yml")
