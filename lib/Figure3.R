library(ape)
library(tidyverse)

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
        Name <- c(Name, rep("NGS consensus", 2*n))
    } else if (endsWith(filename, "sanger")) {
        n <- 1
        BranchLengths <- c(BranchLengths, sum_branches(trees), sum_tip_branches(trees))
        Name <- c(Name, rep("Sanger consensus", 2*n))
    } else {
        n <- length(trees)
        BranchLengths <- c(BranchLengths, sapply(trees, sum_branches), sapply(trees, sum_tip_branches))
        Name <- c(Name, rep("NGS profile-sampled", 2*n))
    }
    Gene <- c(Gene, rep(gene, 2*n))
    Tips <- c(Tips, rep(FALSE, n), rep(TRUE, n))
}

data <- tibble(Gene=factor(Gene, levels=c("prrt", "int", "env", "wgs")),
               BranchLengths=BranchLengths,
               Name=Name,
               Tips=factor(Tips, levels=c("FALSE", "TRUE"), labels=c("All Branches", "Tip Branches")))

g <- ggplot(data, aes(x=Gene, y=BranchLengths, fill=Gene, shape=Name)) +
     geom_violin(colour="white") +
     geom_point(data=filter(data, Name!="NGS profile-sampled"), colour="black") +
     expand_limits(y=0) +
     labs(y="Sum of Branch Lengths in Tree") +
     facet_grid(. ~ Tips) +
     theme(legend.position="none")

ggsave(outfile, g, width=6, height=6, units="in")
