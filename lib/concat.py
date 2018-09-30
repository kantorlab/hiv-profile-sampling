import sys
from Bio import SeqIO
from itertools import chain

infiles = sys.argv[1:-1]
outfile = sys.argv[-1]

alignments = [SeqIO.index(f, "fasta") for f in infiles]

lengths = [max(map(len, alignment.values())) for alignment in alignments]
print(lengths)

datasets = sorted(frozenset(chain(*[list(alignment.keys()) for alignment in alignments])))
print(datasets)

with open(outfile, "w") as f:
    for dataset in datasets:
        print(">{}".format(dataset), file=f)
        for n, alignment in zip(lengths, alignments):
            if dataset in alignment:
                f.write(str(alignment[dataset].seq))
            else:
                f.write("-" * n)
        f.write("\n")
