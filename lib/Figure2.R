library(tidyverse)

args <- commandArgs(trailingOnly=TRUE)
distfile <- args[1]
mdsfile  <- args[2]
outfile  <- args[3]

load(distfile, verbose=TRUE)
load(mdsfile, verbose=TRUE)

print("Percent of variance explained by MDS axes:")
print(round(mds$eig*100/sum(mds$eig), 1))

data <- tibble(Gene=factor(genes, levels=c("prrt", "int", "env", "wgs")),
               Consensus=factor(consensus),
               X=mds$points[,1],
               Y=mds$points[,2])

g <- ggplot(data, aes(x=X, y=Y, colour=Gene, shape=Consensus)) +
     geom_point(alpha=0.3) +
     scale_shape_manual(values=c(3, 16)) +
     labs(x="", y="") +
     theme(legend.position="bottom")

ggsave(outfile, g, width=6, height=6, units="in")
