import sys
from Bio import SeqIO

fafile, outfile = sys.argv[1:]

sequences = list(SeqIO.parse(fafile, "fasta"))

def length(s):
    return sum(c != "-" for c in s)

def distance(s1, s2):
    return sum(c1 != c2 for (c1, c2) in zip(s1, s2) if c1 != "-" and c2 != "-")

with open(outfile, "w") as f:
    print("id1,id2,min_length,distance", file=f)
    for i, seq1 in enumerate(sequences):
        for seq2 in sequences[i:]:
            print(",".join([seq1.id,
                            seq2.id,
                            str(min(length(seq1.seq), length(seq2.seq))),
                            str(distance(str(seq1.seq), str(seq2.seq)))]),
                  file=f)

# vim: syntax=python expandtab sw=4 ts=4
