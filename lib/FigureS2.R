library(tidyverse)
library(ggtree)
library(gridExtra)
library(grid)
library(devEMF)
library(ape)
library(treeio)

args <- commandArgs(trailingOnly=TRUE)
infiles <- args[1:(length(args)-1)]
outfile <- args[length(args)]

plots <- list()
genes <- c()

for (i in seq(1, length(infiles), 2))
{
    treefile = infiles[i]
    clusterfile = infiles[i+1]

    # Find tree name in filename path.
    if (endsWith(treefile, "consensus")) {
        name <- "consensus"
    } else if (endsWith(treefile, "sanger")) {
        name <- "sanger"
    }

    # Find gene name in filename path.
    path <- str_split(treefile, "/")[[1]]
    gene <- path[length(path)-1]
    genes <- c(genes, gene)

    # Load tree
    tree <- drop.tip(read.raxml(treefile), "outgroup")

    # Load cluster support
    support <- read_csv(clusterfile)
    support$tips <- lapply(support$cluster, function(x) strsplit(x, ",")[[1]])
    support$level <- sapply(support$cluster, function(x) str_count(x, ","))
    support$MRCA <- sapply(support$tips, function(x) getMRCA(as.phylo(tree), x))
    print(support)

    # Plot tree
    #xscale <- max(node.depth.edgelength(tree))
    xscale <- 0.4
    g <- ggtree(tree) +
         geom_tiplab(size=1.8, hjust=-0.1) +
         geom_treescale(y=-1, width=0.01, fontsize=1.8) +
         xlim_tree(1.1*xscale) +
         ggtitle(gene) +
	 geom_text2(aes(subset=(!isTip & bootstrap >= 70), label=bootstrap), hjust=1.2, nudge_y=0.5, size=1.5) +
         theme(plot.title=element_text(face="bold"))
    support <- filter(support, support >= 1)
    for (i in 1:nrow(support)) {
        g <- g + geom_cladelabel(support$MRCA[i], label=round(support$support[i]), offset=(0.09 + 0.04*support$level[i])*xscale, color="blue", hjust=-0.2, size=1, fontsize=1.8)
    }
    support <- filter(support, get(name) == 1)
    for (i in 1:nrow(support)) {
        cluster <- support$cluster[i]
      	g <- g + geom_cladelabel(support$MRCA[i], label="", offset=(0.08 + 0.04*support$level[i])*xscale, color="red", size=1)
    }
    plots[[gene]] <- g
}

pdf(outfile, width=6.5, height=8)
grid.arrange(grobs=plots, ncol=2, nrow=2)
dev.off()
