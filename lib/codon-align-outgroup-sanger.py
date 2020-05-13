import sys
from Bio import Seq, SeqIO
from collections import defaultdict
from subprocess import Popen, PIPE

in_file, outgroup_file, out_aa_file, out_nt_file = sys.argv[1:]

# Read nucleotide sequences
if in_file == "-":
    nt_sequences = list(SeqIO.parse(sys.stdin, "fasta"))
else:
    with open(in_file, "r") as f:
        nt_sequences = list(SeqIO.parse(f, "fasta"))

# Append outgroup
with open(outgroup_file, "r") as f:
    nt_sequences += list(SeqIO.parse(f, "fasta"))

# Translate to AA
aa_sequences = "\n".join(">{}\n{}".format(record.id, str(record.seq.translate())) for record in nt_sequences)

# Align with mafft
with open(out_aa_file, "w") as f:
    with Popen(["mafft", "-"], stdin=PIPE, stdout=f) as mafft:
        mafft.communicate(bytes(aa_sequences, "ascii"))

# Convert aligned AA to codons
nt_sequences = dict((record.id, record.seq) for record in nt_sequences)
ids = []
columns = defaultdict(list)
for record in SeqIO.parse(out_aa_file, "fasta"):
    ids.append(record.id)
    aa_seq = str(record.seq)
    nt_seq = str(nt_sequences[record.id])
    assert len(nt_seq) // 3 == sum(1 for aa in aa_seq if aa != "-"), record.id
    i = 0
    for j, aa in enumerate(aa_seq):
        if aa == "-":
            columns[j].append("---")
        else:
            codon = nt_seq[i:i+3]
            assert str(Seq.translate(codon)) == aa.upper(), codon
            columns[j].append(codon)
            i += 3
assert ids[-1] == "outgroup"

# Remove gapped columns
N = len(nt_sequences) - 1
idx = [i for i in sorted(columns) if columns[i].count("---") < N]
with open(out_nt_file, "w") as f:
    for j, id in enumerate(ids):
        print(">{}\n{}".format(id, "".join(columns[i][j] for i in idx)), file=f)
