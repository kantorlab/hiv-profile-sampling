import Bio.SeqIO
import os
import sys

for fa_file in sys.argv[2:]:
    dataset = os.path.basename(fa_file).split(".")[1]
    print(">{}".format(dataset))
    print(Bio.SeqIO.index(fa_file, "fasta")[sys.argv[1]].seq)

