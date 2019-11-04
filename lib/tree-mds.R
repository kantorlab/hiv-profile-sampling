library(ape)

args <- commandArgs(trailingOnly=TRUE)
distfile <- args[1]
outfile  <- args[2]

load(distfile, verbose=TRUE)

# Compute multi-dimensional scaling.
mds <- cmdscale(distance, k=2, eig=TRUE)

save(genes, consensus, mds, file=outfile)

