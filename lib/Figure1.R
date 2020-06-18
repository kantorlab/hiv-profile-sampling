library(tidyverse)
library(devEMF)

args <- commandArgs(trailingOnly=TRUE)
distfiles <- args[1:(length(args)-1)]
outfile  <- args[length(args)]

Gene <- c()
Dataset <- c()
Diversity <- c()

for (distfile in distfiles)
{
    # Find gene name in filename path.
    path <- str_split(distfile, "/")[[1]]
    gene <- path[length(path)-1]

    # Find dataset name in filename path.
    path <- str_split(distfile, "\\.")[[1]]
    dataset <- path[length(path)-1]

    dist <- read_csv(distfile)
    dist$diversity <- dist$distance / dist$min_length

    Gene <- c(Gene, gene)
    Dataset <- c(Dataset, dataset)
    Diversity <- c(Diversity, mean(dist$diversity))
}

data <- tibble(Gene=factor(Gene, levels=c("prrt", "int", "env", "wgs")),
               Dataset=factor(Dataset),
               Diversity=Diversity)

# Print min/max ranges for each region
print(data %>%
      group_by(Gene) %>%
      summarise(Min=min(Diversity), Max=max(Diversity)))

# Order by ascending env diversity
env <- filter(data, Gene=="env") %>% arrange(Diversity)
data <- mutate(data, Dataset=factor(Dataset, levels=env$Dataset))

g <- ggplot(data, aes(x=Dataset, y=Diversity, color=Gene, group=Gene)) +
     geom_point(size=1) +
     facet_grid(Gene ~ .) +
     labs(x="Patient", y="Within-Host Diversity", y2="Gene") +
     scale_y_continuous(limits=c(0, 0.04), breaks=c(0.01, 0.02, 0.03, 0.04), labels=scales::percent_format(accuracy=1)) +
     scale_color_manual(values=c("#DF8F44", "#00A1D5", "#374E55", "#B24745")) +
     theme(legend.position="none",
           panel.grid.minor=element_line(size=0.25),
           panel.grid.major=element_line(size=0.25),
           axis.ticks=element_blank(),
           axis.text=element_text(color="black"),
           axis.text.x=element_text(angle=90, vjust=0.5))

emf(outfile, width=6.5, height=5)
print(g)
dev.off()
