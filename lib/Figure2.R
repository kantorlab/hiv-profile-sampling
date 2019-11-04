library(tidyverse)

args <- commandArgs(trailingOnly=TRUE)
mdsfiles <- args[1:(length(args)-1)]
outfile  <- args[length(args)]

Gene <- c()
Consensus <- c()
X <- c()
Y <- c()

for (mdsfile in mdsfiles)
{
    load(mdsfile, verbose=TRUE)

    print(paste0("Percent of variance explained by MDS axes for ", genes[1], ":"))
    print(round(mds$eig*100/sum(mds$eig), 1))

    Gene <- c(Gene, genes)
    Consensus <- c(Consensus, consensus)
    X <- c(X, mds$points[,1])
    Y <- c(Y, mds$points[,2])
}

data <- tibble(Gene=factor(Gene, levels=c("prrt", "int", "env", "wgs")),
               Consensus=Consensus,
               X=X,
               Y=Y) %>%
        arrange(Consensus)

g <- ggplot(data, aes(x=X, y=Y, colour=Gene)) +
     geom_point(shape=3) +
     geom_point(data=filter(data, Consensus==TRUE), shape=16, colour="black") +
     facet_grid(. ~ Gene) +
     xlim(-0.06, 0.06) +
     ylim(-0.06, 0.06) +
     labs(x="MDS Axis 1", y="MDS Axis 2") +
     theme(legend.position="none",
           axis.text=element_text(size=6))

ggsave(outfile, g, width=12, height=3.5, units="in")
