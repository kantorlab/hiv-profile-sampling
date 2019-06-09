import sys
from Bio import AlignIO, Seq, SeqIO
from io import StringIO
from subprocess import Popen, PIPE

in_file, ref_file, out_file = sys.argv[1:]

# Read nucleotide sequences
if in_file == "-":
    nt_sequences = list(SeqIO.parse(sys.stdin, "fasta"))
else:
    with open(in_file, "r") as f:
        nt_sequences = list(SeqIO.parse(f, "fasta"))

# Translate to AA
aa_sequences = [record.seq.translate() for record in nt_sequences]
aa_fasta = "\n".join(">{}\n{}".format(record.id, record.seq) for record in aa_sequences)

# Align with hmmer
with Popen(["hmmalign", "--amino", ref_file, "-"], stdin=PIPE, stdout=PIPE) as hmmalign:
    aa_alignment = AlignIO.read(StringIO(hmmalign.communicate(aa_fasta)), "stockholm")

# Convert aligned AA to codons
nt_sequences = dict((record.id, record.seq) for record in nt_sequences)
with open(out_file, "w") as f:
    for record in aa_alignment:
        print(">{}".format(record.id), file=f)
        aa_seq = record.seq
        nt_seq = nt_sequences[record.id].seq
        assert len(nt_seq) == 3*len(1 for aa in aa_seq if aa != "-")
        i = 0
        for aa in aa_seq:
            if aa == "-":
                f.write("---")
            else:
                codon = nt_seq[i:i+3]
                assert Seq.translate(codon) == aa
                f.write(codon)
                i += 3
        f.write("\n")

