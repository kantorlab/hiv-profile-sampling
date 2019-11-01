library(ape)
library(tidyverse)

args <- commandArgs(trailingOnly=TRUE)
ngenes <- (length(args)-1)/2
outfile <- args[length(args)]

print(paste(ngenes, "genes"))

sum_branches <- function(tree) {
    sum(tree$edge.length)
}

sum_tip_branches <- function(tree) {
    sum(tree$edge.length[1:length(tree$tip.label)])
}

sampled_trees <- function(filename) {
    path <- str_split(filename, "/")[[1]]
    gene <- path[length(path)-1]
    trees <- read.tree(filename)
    tibble(Gene=rep(gene, 2*length(trees)),
           BranchLengths=c(sapply(trees, sum_branches), sapply(trees, sum_tip_branches)),
           Consensus=rep(FALSE, 2*length(trees)),
           Tips=c(rep(FALSE, length(trees)), rep(TRUE, length(trees))))
}

consensus_trees <- function(filename) {
    path <- str_split(filename, "/")[[1]]
    gene <- path[length(path)-1]
    tree <- read.tree(filename)
    tibble(Gene=c(gene, gene),
           BranchLengths=c(sum_branches(tree), sum_tip_branches(tree)),
           Consensus=c(TRUE, TRUE),
           Tips=c(FALSE, TRUE))
}

data <- bind_rows(lapply(args[1:ngenes], sampled_trees),
                  lapply(args[(ngenes+1):(2*ngenes)], consensus_trees)) %>%
        mutate(Gene=factor(Gene, levels=c("prrt", "int", "env", "wgs")),
               Tips=factor(Tips, levels=c("FALSE", "TRUE"), labels=c("All Branches", "Tip Branches")))
print(head(data))

cons <- data %>% filter(Consensus==TRUE)

g <- ggplot(data, aes(x=Gene, y=BranchLengths)) +
     geom_violin() +
     geom_point(data=cons, colour="red") +
     expand_limits(y=0) +
     labs(y="Sum of Branch Lengths in Tree") +
     facet_grid(. ~ Tips)

ggsave(outfile, g, width=6, height=6, units="in")
