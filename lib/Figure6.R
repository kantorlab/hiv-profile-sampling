library(ape)
library(tidyverse)
library(ggtree)
library(gridExtra)

args <- commandArgs(trailingOnly=TRUE)
infiles <- args[1:(length(args)-1)]
outfile <- args[length(args)]

data <- bind_rows(lapply(infiles, function(infile) {

    # Find gene name in filename path.
    path <- str_split(infile, "/")[[1]]
    gene <- path[length(path)-1]

    # Load cluster support
    read_csv(infile) %>% mutate(Cluster=cluster,
                                Gene=gene,
                                Support=support,
                                Consensus=case_when(consensus == 1 & sanger == 1 ~ "Both",
                                                    consensus == 1               ~ "NGS Consensus",
                                                    sanger == 1                  ~ "Sanger Consensus",
                                                    TRUE                         ~ "Neither")) %>%
                         select(Cluster, Gene, Support, Consensus)
}))

# Adjust subclusters
data <- bind_rows(data,
                  tibble(Gene=c("wgs", "prrt", "wgs", "env", "prrt"),
                         Cluster=c("MC17,MC21", "MC17,MC21", "MC47,MC56", "MC47,MC56", "MC47,MC56"),
                         Support=c(0, 0, 0, 0, 0),
                         Consensus=c("Both", "Neither", "NGS Consensus", "Sanger Consensus", "Neither")))

for (gene in c("wgs", "env", "int", "prrt")) {
    data[which(data$Gene == gene & data$Cluster == "MC17,MC21"), c("Support")] <- data[which(data$Gene == gene & data$Cluster == "MC17,MC21"), c("Support")] +
                                                                                  data[which(data$Gene == gene & data$Cluster == "MC17,MC20,MC21"), c("Support")]
}
for (gene in c("wgs", "env", "prrt")) {
    data[which(data$Gene == gene & data$Cluster == "MC47,MC56"), c("Support")] <- data[which(data$Gene == gene & data$Cluster == "MC47,MC56"), c("Support")] +
                                                                                  data[which(data$Gene == gene & data$Cluster == "MC41,MC47,MC56"), c("Support")]
}
data[which(data$Gene == "wgs" & data$Cluster == "MC47,MC56"), c("Support")] <- data[which(data$Gene == "wgs" & data$Cluster == "MC47,MC56"), c("Support")] +
                                                                               data[which(data$Gene == "wgs" & data$Cluster == "MC37,MC41,MC47,MC53,MC56"), c("Support")]
data[which(data$Gene == "wgs" & data$Cluster == "MC41,MC47,MC56"), c("Support")] <- data[which(data$Gene == "wgs" & data$Cluster == "MC41,MC47,MC56"), c("Support")] +
                                                                                    data[which(data$Gene == "wgs" & data$Cluster == "MC37,MC41,MC47,MC53,MC56"), c("Support")]
data[which(data$Gene == "wgs" & data$Cluster == "MC37,MC53"), c("Support")] <- data[which(data$Gene == "wgs" & data$Cluster == "MC37,MC53"), c("Support")] +
                                                                               data[which(data$Gene == "wgs" & data$Cluster == "MC37,MC41,MC47,MC53,MC56"), c("Support")]
#
# Order genes
data$Gene <- factor(data$Gene, levels=c("prrt", "int", "env", "wgs"))

# Order consensus
data$Consensus <- factor(data$Consensus, levels=c("NGS Consensus", "Sanger Consensus", "Both", "Neither"))

# Percent support
data$Percent <- paste0(formatC(data$Support, format="f", digits=1), "%")
data$Support <- data$Support * 0.01

g <- ggplot(data, aes(x=Gene, y=Cluster, fill=Consensus, label=Percent, alpha=Support)) +
     geom_tile(color=NA) +
     geom_text(size=3.3, color="black") +
     scale_fill_manual(values=c("Both"="#ffffbf", "NGS Consensus"="#fc8d59", "Sanger Consensus"="#99d594", "Neither"="white")) +
     scale_alpha_continuous(labels=scales::percent_format()) +
     theme(legend.position="right",
           panel.background=element_blank(),
           axis.text.x=element_text(angle=90))

ggsave(outfile, g, width=6.5, height=5, units="in")
