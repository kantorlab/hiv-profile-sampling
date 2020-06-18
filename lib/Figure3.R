library(ape)
library(tidyverse)
library(devEMF)

args <- commandArgs(trailingOnly=TRUE)
filenames <- args[1:(length(args)-1)]
outfile <- args[length(args)]

sum_branches <- function(tree) {
    sum(tree$edge.length)
}

sum_tip_branches <- function(tree) {
    sum(tree$edge.length[1:length(tree$tip.label)])
}

Gene <- c()
BranchLengths <- c()
Name <- c()
Tips <- c()

for (filename in filenames)
{
    # Find gene name in filename path.
    path <- str_split(filename, "/")[[1]]
    gene <- path[length(path)-1]

    # Load trees.
    trees <- read.tree(filename)

    # Append values.
    if (endsWith(filename, "consensus")) {
        n <- 1
        BranchLengths <- c(BranchLengths, sum_branches(trees), sum_tip_branches(trees))
        Name <- c(Name, rep("NGS Consensus", 2*n))
    } else if (endsWith(filename, "sanger")) {
        n <- 1
        BranchLengths <- c(BranchLengths, sum_branches(trees), sum_tip_branches(trees))
        Name <- c(Name, rep("Sanger Consensus", 2*n))
    } else {
        n <- length(trees)
        BranchLengths <- c(BranchLengths, sapply(trees, sum_branches), sapply(trees, sum_tip_branches))
        Name <- c(Name, rep(NA, 2*n))
    }
    Gene <- c(Gene, rep(gene, 2*n))
    Tips <- c(Tips, rep(FALSE, n), rep(TRUE, n))
}

data <- tibble(Gene=factor(Gene, levels=c("prrt", "int", "env", "wgs")),
               BranchLengths=BranchLengths,
               Name=factor(Name, levels=c("NGS Consensus", "Sanger Consensus")),
               Tips=factor(Tips, levels=c("FALSE", "TRUE"), labels=c("All Branches", "Tip Branches")))

g <- ggplot(data, aes(x=Gene, y=BranchLengths, color=Gene)) +
     geom_violin(color=NA, fill="gray", show.legend=FALSE) +
     geom_point(aes(shape=Name), data=filter(data, !is.na(Name)), size=3) +
     labs(y="Sum of Branch Lengths in Tree") +
     facet_grid(Tips ~ .) +
     scale_y_continuous(limits=c(0, 6), breaks=seq(0, 6, 1)) +
     scale_color_manual(values=c("#DF8F44", "#00A1D5", "#374E55", "#B24745")) +
     coord_flip() +
     theme_minimal() +
     theme(legend.position="bottom",
           legend.title=element_blank(),
           panel.grid.minor=element_line(size=0.25),
           panel.grid.major=element_line(size=0.25),
           panel.grid.major.y=element_blank(),
           axis.text=element_text(color="black"))

emf(outfile, width=6.5, height=4)
print(g)
dev.off()
