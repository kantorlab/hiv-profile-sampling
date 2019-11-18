library(ape)
library(tidyverse)
library(ggtree)
library(gridExtra)

args <- commandArgs(trailingOnly=TRUE)
infiles <- args[1:(length(args)-1)]
outfile <- args[length(args)]

colors <- list()
data   <- data.frame(Cluster=character(),
                     Gene=character(),
                     Name=character(),
                     Detected=integer(),
                     Support=double())

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

for (infile in infiles)
{
    # Find gene name in filename path.
    path <- str_split(infile, "/")[[1]]
    gene <- path[length(path)-1]

    # Load cluster support
    support <- read_csv(infile)
    support$N <- NULL

    # Assign color to cluster
    for (i in 1:nrow(support)) {
        cluster <- support$cluster[i]
        if (!(cluster %in% names(colors))) {
            colors[[cluster]] <- kelly_palette[length(colors)+1]
        }
    }

    # Unpivot and append to dataframe
    support <- gather(support, "consensus", "sanger", key="Name", value="Detected")
    support$Gene <- gene
    colnames(support)[1:2] <- c("Cluster", "Support")
    data <- rbind(data, support)
    print(tail(data))
}

data$Name <- factor(data$Name, levels=c("consensus", "sanger"), labels=c("NGS Consensus", "Sanger Consensus"))
data$Gene <- factor(data$Gene, levels=c("prrt", "int", "env", "wgs"))
data$Cluster <- factor(data$Cluster, levels=names(colors))

data <- filter(data, Detected==1)

g <- ggplot(data, aes(x=Gene, y=Cluster, fill=Cluster, label=Support)) +
     geom_tile(color="black", alpha=0.33) +
     geom_text(size=2) +
     scale_fill_manual(values=unlist(colors, use.names=FALSE)) +
     facet_grid(. ~ Name) + 
     theme(legend.position="none",
           panel.background=element_blank())

ggsave(outfile, g, width=5, height=5, units="in")
