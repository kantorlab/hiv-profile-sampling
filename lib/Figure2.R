library(tidyverse)

args <- commandArgs(trailingOnly=TRUE)
mdsfiles <- args[1:(length(args)-1)]
outfile  <- args[length(args)]

Gene <- c()
Name <- c()
X <- c()
Y <- c()

for (mdsfile in mdsfiles)
{
    load(mdsfile, verbose=TRUE)

    print(paste0("Percent of variance explained by MDS axes for ", genes[1], ":"))
    print(round(mds$eig*100/sum(mds$eig), 1))

    Gene <- c(Gene, genes)
    Name <- c(Name, names)
    X <- c(X, mds$points[,1])
    Y <- c(Y, mds$points[,2])
}

data <- tibble(Gene=factor(Gene, levels=c("prrt", "int", "env", "wgs")),
               Name=Name,
               X=X,
               Y=Y) %>%
        arrange(Name)

g <- ggplot(data, aes(x=X, y=Y)) +
     geom_point(aes(color=Gene), data=filter(data,  is.na(Name)), shape=3, show.legend=FALSE) +
     geom_point(aes(shape=Name), data=filter(data, !is.na(Name)), color="black") +
     facet_grid(. ~ Gene) +
     xlim(-0.06, 0.06) +
     ylim(-0.06, 0.06) +
     labs(x="MDS Axis 1", y="MDS Axis 2") +
     theme(legend.position="right",
           legend.title=element_blank(),
           axis.text=element_text(size=6))

ggsave(outfile, g, width=12, height=3.2, units="in")
