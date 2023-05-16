import re
from datetime import datetime

datetime_start = None
datetime_end = None

sample_name = "OMB1556_Ach_MB_H"

with open(f"./results/spaceranger_count/{sample_name}/_log", "rt") as logfile:
    for line in logfile:
        datetime_match = re.match("\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}", line)
        if datetime_match is None:
            continue
        datetime_str = datetime_match.group(0)
        if datetime_start is None:
            datetime_start = datetime.strptime(datetime_str, '%Y-%m-%d %H:%M:%S')
    datetime_end = datetime.strptime(datetime_str, '%Y-%m-%d %H:%M:%S')

duration = datetime_end - datetime_start
duration_seconds = duration.total_seconds()

complete_hours = divmod(duration_seconds, 3600)[0]
leftover_seconds = duration_seconds - complete_hours*3600
complete_minutes = divmod(leftover_seconds, 60)[0]
complete_seconds = leftover_seconds - complete_minutes*60

print(f"sample_name\thours\tminutes\tseconds")
print(f"{sample_name}\t{complete_hours}\t{complete_minutes}\t{complete_seconds}")
