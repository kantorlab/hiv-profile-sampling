library(tidyverse)

args <- commandArgs(trailingOnly=TRUE)
mdsfile <- args[1]
outfile <- args[2]

load(mdsfile, verbose=TRUE)

print("Percent of variance explained by MDS axes:")
print(round(mds$eig*100/sum(mds$eig), 1))

data <- tibble(Gene=factor(genes, levels=c("prrt", "int", "env", "wgs")),
               Consensus=consensus,
               X=mds$points[,1],
               Y=mds$points[,2]) %>%
        arrange(Consensus)

g <- ggplot(data, aes(x=X, y=Y, colour=Gene)) +
     geom_point(shape=3) +
     geom_point(data=filter(data, Consensus==TRUE), shape=16, colour="black") +
     labs(x="MDS Axis 1", y="MDS Axis 2") +
     theme(legend.position="bottom")

ggsave(outfile, g, width=12, height=12, units="in")
