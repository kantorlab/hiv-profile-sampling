import sys
from Bio import SeqIO
from collections import defaultdict

consensusfile = sys.argv[1]
sangerfile    = sys.argv[2]
samplefiles   = sys.argv[3:-1]
outfile       = sys.argv[-1]

N = len(samplefiles)

def read_fa_headers(filename):
    """
    Read cluster structure from the header lines in the fasta files outputted
    by ClusterPicker.
    """
    clusters = defaultdict(list)
    with open(filename) as f:
        for record in SeqIO.parse(f, "fasta"):
            i, mc = record.id.split("_")
            clusters[i].append(mc)
    return [",".join(sorted(cluster)) for cluster in clusters.values()]

# Calculate support as frequency of cluster in samples

support = {}

for clusters in map(read_fa_headers, samplefiles):
    for cluster in clusters:
        support[cluster] = support.get(cluster, 0) + 1

consensus = read_fa_headers(consensusfile)
sanger    = read_fa_headers(sangerfile)

# Add any clusters in the consensus/sanger trees that are missing from the
# samples, at 0 support.
for cluster in consensus:
    support[cluster] = support.get(cluster, 0)
for cluster in sanger:
    support[cluster] = support.get(cluster, 0)

# Prune redundant clusters
clusters = list(support)
while clusters:
    cluster = clusters.pop(0).split(",")
    if len(cluster) > 2:
        for i in range(len(cluster)):
            subcluster = ",".join(cluster[:i] + cluster[i+1:])
            if subcluster in support:
                support.pop(subcluster)
            clusters.append(subcluster)

with open(outfile, "w") as f:
    print("cluster,N,support,consensus,sanger", file=f)
    for cluster in sorted(support, key=support.get, reverse=True):
        print("\"{}\"".format(cluster),
              support[cluster],
              100.0*support[cluster]/N,
              int(cluster in consensus),
              int(cluster in sanger),
              sep=",", file=f)

