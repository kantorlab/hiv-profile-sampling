import numpy as np
import os
import pandas as pd
import sys

from Bio import Seq
from collections import defaultdict

seed         = int(sys.argv[1])
min_coverage = int(sys.argv[2])
min_freq     = float(sys.argv[3])
csvfiles     = sys.argv[4:-2]
ntfile       = sys.argv[-2]
aafile       = sys.argv[-1]

np.random.seed(seed)

with open(ntfile, "w") as nt, open(aafile, "w") as aa:

    for csvfile in csvfiles:

        dataset = os.path.basename(csvfile).partition(".")[0]
        print(">{}".format(dataset), file=nt)
        print(">{}".format(dataset), file=aa)

        codons = defaultdict(dict)
        for _, row in pd.read_csv(csvfile, sep="\t").fillna("").iterrows():
            if (row.freq >= min_coverage):
                codons[row.pos][row.codon] = row.freq

        seq = []
        for pos in sorted(codons):
            row = pd.Series(codons[pos])
            row = row.divide(row.sum())
            row = row[row > min_freq]
            if row.sum() > 0:
                row = row.divide(row.sum())
                seq.append(np.random.choice(row.index, p=row.values))

        print("".join(seq), file=nt)
        print("".join(Seq.translate(codon) for codon in seq), file=aa)

# vim: syntax=python expandtab sw=4 ts=4
