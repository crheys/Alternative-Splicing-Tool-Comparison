#!/usr/bin/env python3

import gzip

#Define paths to the counts file and output groups files
counts_file = "/mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/HD/AS_analysis/leafcutter/intron_clustering/leafcutter_perind_numers.counts.gz"
output_file = "/mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/HD/AS_analysis/leafcutter/groups.txt"

# Define conditions for each sample
# Adjust this to match your experiment
sample_to_group = {
    "SRR3306823": "con1",
    "SRR3306824": "con1",
    "SRR3306825": "con1",
    "SRR3306826": "con1",
    "SRR3306827": "con1",
    "SRR3306828": "con1",
    "SRR3306829": "con1",
    "SRR3306830": "con2",
    "SRR3306831": "con2",
    "SRR3306832": "con2",
    "SRR3306833": "con2",
    "SRR3306834": "con2",
    "SRR3306835": "con2",
    "SRR3306836": "con2"
}

with gzip.open(counts_file, 'rt') as f:
    header = f.readline().strip()

samples = header.split()

# Sanity check
missing = [s for s in samples if s not in sample_to_group]
if missing:
    raise ValueError(f"Missing group assignment for: {missing}")

with open(output_file, 'w') as out:
    for sample in samples:
        out.write(f"{sample} {sample_to_group[sample]}\n")

print(f"groups.txt written with {len(samples)} samples.")
