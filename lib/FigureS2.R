library(ape)
library(tidyverse)
library(ggtree)

args <- commandArgs(trailingOnly=TRUE)
treefiles <- args[1:(length(args)-1)]
outfile <- args[length(args)]

plots <- list()
genes <- c()
order <- c()

for (treefile in treefiles)
{
    # Find gene name in filename path.
    path <- str_split(treefile, "/")[[1]]
    gene <- path[length(path)-1]
    genes <- c(genes, gene)

    # Load tree
    trees <- read.tree(treefile)

    if (length(order) == 0) {
        tree <- ladderize(trees[[1]])
        is_tip <- tree$edge[,2] <= length(tree$tip.label)
        ordered_tips <- tree$edge[is_tip, 2]
        order <- rev(tree$tip.label[ordered_tips])
    }

    # Plot tree
    plots[[gene]] <- ggdensitree(trees, alpha=0.1*(100/length(trees)), tip.order=order) +
                     geom_tiplab(size=1.6)
}

pdf(outfile, width=18, height=6)
multiplot(plotlist=plots, ncol=4, labels=genes)
dev.off()
