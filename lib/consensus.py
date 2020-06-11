import pandas as pd
import sys
from Bio import Seq
from collections import defaultdict

min_coverage = 10

codon_files = sys.argv[1:-1]
out_file    = sys.argv[-1]

_ambiguous = dict(("".join(sorted(b)), a) for a, b in Seq.IUPAC.IUPACData.ambiguous_dna_values.items())
consensus = defaultdict(list)

for codon_file in codon_files:

    dataset = codon_file.split("/")[-2]
    codons = pd.read_csv(codon_file, sep="\t", index_col="pos").fillna("")

    # Sum counts by site
    sums = codons.groupby(level=0)["freq"].sum()
    sums = sums.to_frame(name="sum")
    codons = codons.join(sums, how="left")

    # Filter out low-coverage sites
    codons = codons[codons["sum"] >= min_coverage]

    # Normalize counts
    codons["freq"] = codons["freq"] * (1.0 / codons["sum"])

    # 20% consensus
    codons = codons[codons["freq"] >= 0.2]

    for i, rows in codons.groupby(level=0):
        for nt in zip(*rows["codon"].tolist()):
            consensus[dataset].append(_ambiguous["".join(sorted(set(nt)))])

with open(out_file, "w") as f:
    for dataset in sorted(consensus):
        print(">{}".format(dataset), file=f)
        print("".join(consensus[dataset]), file=f)

# vim: syntax=python expandtab sw=4 ts=4
