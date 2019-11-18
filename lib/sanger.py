import sys
from Bio import AlignIO
from Bio import Seq

in_file, datasets, prrt_file, int_file, env_file, wgs_file = sys.argv[1:]

datasets = frozenset("MC" + i for i in datasets.split(","))

regions = {
    "prrt": (prrt_file, 2253, 3869),
    "int": (int_file, 4230, 5096),
    "env": (env_file, 6225, 8795)
}

align = AlignIO.read(in_file, "fasta")
N = align.get_alignment_length()

# find hxb2 alignment
for a in align:
    if a.id == "hxb2":
        hxb2 = a
        break

# map hxb2 coordinates to alignment coordinates
coords = {}
j = 0
for i, nt in enumerate(hxb2.seq):
    if nt != "-":
        coords[j] = i
        j += 1

# non-gap columns
nogap = set()
for i in range(N):
    if sum(1 for a in align if a.seq[i] != "-" and a.id in datasets):
        nogap.add(i)

# trim alignment for each region, excluding hxb2
for region in regions:

    out_file, start, end = regions[region]

    # adjust 1-index inclusive to 0-index exclusive range
    start = coords[int(start) - 1]
    end   = coords[int(end)]

    with open(out_file, "w") as f:
        for a in align:
            if a.id in datasets:
                print(">" + a.description, file=f)
                print("".join(a.seq[i] for i in range(start, end) if i in nogap), file=f)

# wgs alignment includes all coordinates, excluding hxb2
with open(wgs_file, "w") as f:
    for a in align:
        if a.id in datasets:
            print(">" + a.description, file=f)
            print("".join(a.seq[i] for i in range(N) if i in nogap), file=f)

# vim: expandtab sw=4 ts=4
