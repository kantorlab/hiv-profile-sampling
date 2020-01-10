library(tidyverse)
library(gridExtra)

args <- commandArgs(trailingOnly=TRUE)

samplefile <- args[1]
distfiles  <- args[2:(length(args)-1)]
outfile    <- args[length(args)]

samples <- read_csv(samplefile) %>%
           filter(Included == 1 & Analyte == "PL" & abs(VLToEnrollDays) < 90) %>%
           mutate(Dataset=factor(MC)) %>%
           select("Dataset", "VL")

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
               Diversity=Diversity) %>%
        left_join(samples, by="Dataset")

# Print min/max ranges for each region
print(data %>%
      group_by(Gene) %>%
      summarise(Min=min(Diversity), Max=max(Diversity)))

# Order by ascending env diversity
env <- filter(data, Gene=="env") %>% arrange(Diversity)
data <- mutate(data, Dataset=factor(Dataset, levels=env$Dataset))

p1 <- ggplot(data, aes(x=Dataset, y=VL)) +
      geom_bar(stat="identity") +
      labs(x="Patient", y="Viral Load") +
      theme(axis.text.x=element_text(angle=90),
            axis.text=element_text(size=6),
            axis.title=element_text(size=9))

p2 <- ggplot(data, aes(x=Dataset, y=Diversity, colour=Gene)) +
      geom_point() +
      labs(x="Patient", y="Within-host Genetic Diversity") +
      scale_y_continuous(labels=scales::percent_format(accuracy=1)) +
      theme(legend.position="bottom",
            axis.text.x=element_text(angle=90),
            axis.text=element_text(size=6),
            axis.title=element_text(size=9),
            legend.text=element_text(size=11))

pdf(outfile, width=10, height=5)
grid.arrange(p1, p2, ncol=1, heights=c(1,2))
dev.off()
