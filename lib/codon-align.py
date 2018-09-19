import sys

from Bio import SearchIO
from Bio import Seq
from Bio import SeqIO

reads = sys.argv[1]
hmmer = sys.argv[2]

print("indexing reads", file=sys.stderr)
reads = SeqIO.index(reads, "fasta")

print("parsing hmmer", file=sys.stderr)
hmmer = SearchIO.read(hmmer, "hmmer3-text")

nref = max(f.query_end for f in hmmer.fragments)

for hit in hmmer.hits:

    print(">" + hit.id)
    seq = reads[hit.id].seq
    j = 0

    for hsp in hit.hsps:

        sys.stdout.write("---" * (hsp.query_start - j))

        i = 3*hsp.hit_start
        j = hsp.query_start

        for aa in hsp.aln[1].seq:
            if aa.islower():
                i += 3
            elif aa == "-":
                sys.stdout.write("---")
                j += 1
            else:
                sys.stdout.write(str(seq[i:i+3]))
                i += 3
                j += 1
            if j >= nref: break

    sys.stdout.write("---" * (nref - j))
    sys.stdout.write("\n")

# vim: expandtab sw=4 ts=4
