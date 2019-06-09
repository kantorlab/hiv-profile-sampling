import sys
from Bio import Seq, SeqIO
from io import StringIO
from subprocess import Popen, PIPE

in_file, ref_file, out_aa_file, out_nt_file = sys.argv[1:]

# Read nucleotide sequences
if in_file == "-":
    nt_sequences = list(SeqIO.parse(sys.stdin, "fasta"))
else:
    with open(in_file, "r") as f:
        nt_sequences = list(SeqIO.parse(f, "fasta"))

# Translate to AA
aa_sequences = "\n".join(">{}\n{}".format(record.id, record.seq.translate()) for record in nt_sequences)

# Align with hmmer
with Popen(["hmmalign", "--amino", "--outformat", "pfam", "-o", out_aa_file, ref_file, "-"], stdin=PIPE) as hmmalign:
    hmmalign.communicate(bytes(aa_sequences, "ascii"))

# Convert aligned AA to codons
nt_sequences = dict((record.id, record.seq) for record in nt_sequences)
with open(out_nt_file, "w") as f:
    for record in SeqIO.parse(out_aa_file, "stockholm"):
        print(">{}".format(record.id), file=f)
        aa_seq = str(record.seq)
        nt_seq = str(nt_sequences[record.id])
        assert len(nt_seq) == 3*sum(1 for aa in aa_seq if aa != "-")
        i = 0
        for aa in aa_seq:
            if aa == "-":
                f.write("---")
            else:
                codon = nt_seq[i:i+3]
                assert Seq.translate(codon) == aa.upper(), "{} = {}".format(Seq.translate(codon), aa)
                f.write(codon)
                i += 3
        f.write("\n")

