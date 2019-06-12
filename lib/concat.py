import sys
from Bio import SeqIO

alignments = []

for fa_file in sys.argv[1:]:
    alignments.append(SeqIO.to_dict(SeqIO.parse(fa_file, "fasta")))

ids = sorted(set(sum([list(alignment) for alignment in alignments], [])))
lengths = [max(len(record.seq) for record in alignment.values()) for alignment in alignments]

for i in ids:
    print(">{}".format(i))
    for alignment, length in zip(alignments, lengths):
        if i in alignment:
            sys.stdout.write(str(alignment[i].seq))
        else:
            print("-" * length)
    sys.stdout.write("\n")
