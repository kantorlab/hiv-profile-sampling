library(tidyverse)

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

g <- ggplot(data, aes(x=Dataset, y=Diversity, colour=Gene)) +
     geom_point() +
     labs(x="Patient", y="Intra-patient Genetic Diversity") +
     scale_y_continuous(labels=scales::percent_format(accuracy=1)) +
     theme(legend.position="right",
           axis.text.x=element_text(angle = 90))

ggsave(outfile, g, width=12, height=3.5, units="in")
