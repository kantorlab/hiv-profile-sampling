import csv
import sys
from collections import defaultdict
from itertools import groupby
from operator import itemgetter

min_coverage = 10

codon_files = sys.argv[1:-1]
out_file    = sys.argv[-1]

codons = defaultdict(list)

for codon_file in codon_files:
    dataset = codon_file.split("/")[-2]
    with open(codon_file) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for i, rows in groupby(reader, key=itemgetter("pos")):
            row = next(rows)
            if int(row["freq"]) > min_coverage:
                codons[dataset].append(row["codon"])

with open(out_file, "w") as f:
    for dataset in sorted(codons):
        print(">{}".format(dataset), file=f)
        print("".join(codons[dataset]), file=f)

# vim: syntax=python expandtab sw=4 ts=4
