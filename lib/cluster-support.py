import sys
from Bio import SeqIO
from collections import defaultdict

consensusfile = sys.argv[1]
samplefiles   = sys.argv[2:-1]
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

def subsupport(cluster):
    """
    Recursive function to add support for all subclusters.
    """
    cluster = cluster.split(",")
    if len(cluster) > 2:
        for i in range(len(cluster)):
            subcluster = ",".join(cluster[:i] + cluster[i+1:])
            support[subcluster] = support.get(subcluster, 0) + 1
            subsupport(subcluster)

for clusters in map(read_fa_headers, samplefiles):
    for cluster in clusters:
        support[cluster] = support.get(cluster, 0) + 1
        subsupport(cluster)

consensus = read_fa_headers(consensusfile)

with open(outfile, "w") as f:
    print("cluster,N,support,consensus", file=f)
    for cluster in sorted(support, key=support.get, reverse=True):
        print("\"{}\"".format(cluster),
              support[cluster],
              100.0*support[cluster]/N,
              int(cluster in consensus),
              sep=",", file=f)

