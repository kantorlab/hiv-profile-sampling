import os
import pandas as pd
import sys

from Bio import Seq

min_coverage = int(sys.argv[1])
min_freq     = float(sys.argv[2])
csvfiles     = sys.argv[3:-2]
ntfile       = sys.argv[-2]
aafile       = sys.argv[-1]

with open(ntfile, "w") as nt, open(aafile, "w") as aa:

    for csvfile in csvfiles:

        dataset = os.path.basename(csvfile).partition(".")[0]
        print(">{}".format(dataset), file=nt)
        print(">{}".format(dataset), file=aa)

        codons = pd.read_csv(csvfile)
        codons[codons < min_coverage] = 0
        codons = codons.div(codons.sum(axis=1), axis=0)
        codons[codons < min_freq] = 0

        seq = []
        for _, row in codons.iterrows():
            if row.sum() > 0:
                codon = max(row.index, key=lambda x: row[x])
                if codon != "del" and codon[:3] != "ins":
                    seq.append(codon)

        print("".join(seq), file=nt)
        print("".join(Seq.translate(codon) for codon in seq), file=aa)

# vim: syntax=python expandtab sw=4 ts=4
