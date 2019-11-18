library(distory)
library(stringr)

args <- commandArgs(trailingOnly=TRUE)
filenames <- args[1:(length(args)-1)]
outfile <- args[length(args)]

trees <- c()
genes <- c()
names <- c()

for (filename in filenames) 
{
    # Find gene name in filename path.
    path <- str_split(filename, "/")[[1]]
    gene <- path[length(path)-1]

    # Load and append trees.
    t <- read.tree(filename)
    if (length(trees) == 0) {
        trees <- t
    } else {
        trees <- c(trees, t)
    }

    # Append gene and consensus values.
    if (endsWith(filename, "consensus")) {
        n <- 1
        names <- c(names, "NGS Consensus")
    } else if (endsWith(filename, "sanger")) {
        n <- 1
        names <- c(names, "Sanger Consensus")
    } else {
        n <- length(t)
        names <- c(names, rep(NA, n))
    }
    genes <- c(genes, rep(gene, n))
}
.compressTipLabel(trees)

# Compute geodesic distance.
distance <- dist.multiPhylo(trees, force.multi2di=TRUE)

save(trees, genes, names, distance, file=outfile)

