import pandas as pd
import re
from datetime import datetime

datetime_start = None
datetime_end = None

samples = (
    pd.read_csv(snakemake.config["samples"], sep="\t", dtype={"sample_name": str})
    .set_index("sample_name", drop=False)
    .sort_index()
)

sample_names = samples.index.tolist()

with open(snakemake.output[0], "wt") as filestream_out:
    filestream_out.write(f"sample_name\thours\tminutes\tseconds\n")
    for sample_name in sample_names:
        logfile = f"./results/spaceranger_count/{sample_name}/_log"
        with open(logfile, "rt") as filestream_in:
            for line in filestream_in:
                datetime_match = re.match("\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}", line)
                if datetime_match is None:
                    continue
                datetime_str = datetime_match.group(0)
                if datetime_start is None:
                    datetime_start = datetime.strptime(datetime_str, '%Y-%m-%d %H:%M:%S')
            datetime_end = datetime.strptime(datetime_str, '%Y-%m-%d %H:%M:%S')

        duration = datetime_end - datetime_start
        duration_seconds = duration.total_seconds()

        complete_hours = int(divmod(duration_seconds, 3600)[0])
        leftover_seconds = duration_seconds - complete_hours*3600
        complete_minutes = int(divmod(leftover_seconds, 60)[0])
        complete_seconds = int(leftover_seconds - complete_minutes*60)

        filestream_out.write(f"{sample_name}\t{complete_hours}\t{complete_minutes}\t{complete_seconds}\n")
