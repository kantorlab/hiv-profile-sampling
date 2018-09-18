from Bio import SeqIO
import sys
n = 0
for record in SeqIO.parse(sys.argv[1], "fasta"):
    seq = str(record.seq).replace("#", "-").replace("?", "-").replace("$", "*")
    if n == 0:
         assert seq[-1] == "*"
         n = len(seq)
    if len(seq) < n:
        d = n - len(seq) - 1
        seq = "".join((seq, "-" * d, "*"))
    print(">{}".format(record.id))
    print(seq)
