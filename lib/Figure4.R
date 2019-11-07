library(ape)
library(tidyverse)
library(ggtree)

args <- commandArgs(trailingOnly=TRUE)
infiles <- args[1:(length(args)-1)]
outfile <- args[length(args)]

plots <- list()
genes <- c()
colors <- list()

kelly_palette <- c(
    red = "#be0032",
    yellow = "#f3c300",
    blue = "#0067a5",
    olivegreen = "#2b3d26",
    yellowgreen = "#8db600",
    purplishpink = "#e68fac",
    orange = "#f38400",
    purple = "#875692",
    reddishbrown = "#882d17",
    green = "#008856",
    buff = "#c2b280",
    lightblue = "#a1caf1",
    yellowishpink = "#f99379",
    gray = "#848482",
    yellowishbrown = "#654522",
    reddishorange = "#e25822",
    purplishred = "#b3446c",
    greenishyellow = "#dcd300",
    orangeyellow = "#f6a600",
    violet = "#604e97"
)

for (i in seq(1, length(infiles), 2))
{
    treefile = infiles[i]
    clusterfile = infiles[i+1]

    # Find gene name in filename path.
    path <- str_split(treefile, "/")[[1]]
    gene <- path[length(path)-1]
    genes <- c(genes, gene)

    # Load tree
    tree <- root(read.tree(treefile), "MC50")

    # Load cluster support
    support <- read_csv(clusterfile) %>% filter(consensus == 1)
    support$tips <- lapply(support$cluster, function(x) strsplit(x, ",")[[1]])
    support$MRCA <- sapply(support$tips, function(x) getMRCA(tree, x))
    print(support)

    # Plot tree
    xscale <- max(node.depth.edgelength(tree))
    g <- ggtree(tree) +
         geom_tiplab(size=2) +
         geom_treescale(y=-1) +
         xlim_tree(1.1*xscale)
    for (i in 1:nrow(support)) {
        cluster <- support$cluster[i]
        if (!(cluster %in% names(colors))) {
            colors[[cluster]] <- kelly_palette[length(colors)+1]
        }
        g <- g + geom_cladelabel(support$MRCA[i], label=round(support$support[i]), offset=0.1*xscale, color=colors[[cluster]])
    }
    plots[[gene]] <- g
}

pdf(outfile, width=18, height=6)
multiplot(plotlist=plots, ncol=4, labels=genes)
dev.off()
