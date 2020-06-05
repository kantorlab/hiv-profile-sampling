import sys
from Bio import SeqIO
from collections import defaultdict

# Load alignments and convert to columns
ids = set()
columns = defaultdict(dict)
start = 0
for f in sys.argv[1:]:
    for record in SeqIO.parse(f, "fasta"):
        ids.add(record.id)
        for j, nt in enumerate(str(record.seq)):
            if nt != "-":
                columns[start+j][record.id] = nt.upper()
    if columns:
        start = max(columns)

# Remove gapped columns
idx = sorted(j for j in columns if sum(1 for i in columns[j] if i != "outgroup") > 0)

# Write output
for i in sorted(ids):
    print(">{}".format(i))
    for j in idx:
        sys.stdout.write(columns[j].get(i, "-"))
    sys.stdout.write("\n")
