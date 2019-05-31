import csv
import numpy as np
import sys
from itertools import groupby
from operator import itemgetter

min_coverage = 10

nsamples    = int(sys.argv[1])
seed        = int(sys.argv[2])
codon_files = sys.argv[3:-1]
out_file    = sys.argv[-1]

np.random.seed(seed)

codons = []

# Load codon frequencies
for codon_file in codon_files:
    with open(codon_file) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for _, rows in groupby(reader, key=itemgetter("pos")):
            position = {"codon": [], "freq": []}
            for row in rows:
                if int(row["freq"]) > min_coverage:
                    position["codon"].append(row["codon"])
                    position["freq"].append(float(row["freq"]))
            if position["freq"]:
                position["freq"] = np.array(position["freq"])
                position["freq"] /= position["freq"].sum()
                assert abs(position["freq"].sum() - 1) < 0.001, position["freq"]
            codons.append(position)

# Sample from codon frequencies
with open(out_file, "w") as f:
    for i in range(nsamples):
        print(">{}".format(i), file=f)
        for position in codons:
            if position["codon"]:
                f.write(np.random.choice(position["codon"], p=position["freq"]))
        f.write("\n")

# vim: syntax=python expandtab sw=4 ts=4
