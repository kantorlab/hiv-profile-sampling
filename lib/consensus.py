import os
import pandas as pd
import sys

from Bio import Seq
from itertools import groupby
from operator import attrgetter

min_coverage = int(sys.argv[1])
csvfiles     = sys.argv[2:-2]
ntfile       = sys.argv[-2]
aafile       = sys.argv[-1]

with open(ntfile, "w") as nt, open(aafile, "w") as aa:

    for csvfile in csvfiles:

        dataset = os.path.basename(csvfile).partition(".")[0]
        print(">{}".format(dataset), file=nt)
        print(">{}".format(dataset), file=aa)

        codons = pd.read_csv(csvfile, sep="\t").fillna("")
        seq = []
        for _, group in groupby(codons.itertuples(), key=attrgetter("pos")):
            row = next(group)
            if (row.freq >= min_coverage):
                seq.append(row.codon)

        print("".join(seq), file=nt)
        print("".join(Seq.translate(codon) for codon in seq), file=aa)

# vim: syntax=python expandtab sw=4 ts=4
