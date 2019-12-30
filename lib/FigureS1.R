library(tidyverse)

args <- commandArgs(trailingOnly=TRUE)
mdsfile <- args[1]
outfile <- args[2]

load(mdsfile, verbose=TRUE)

print("Percent of variance explained by MDS axes:")
print(round(mds$eig*100/sum(mds$eig), 1))

data <- tibble(Gene=factor(genes, levels=c("prrt", "int", "env", "wgs")),
               Name=names,
               X=mds$points[,1],
               Y=mds$points[,2]) %>%
        arrange(Name)

g <- ggplot(data, aes(x=X, y=Y)) +
     geom_point(aes(color=Gene), data=filter(data,  is.na(Name)), shape=3, show.legend=FALSE) +
     geom_point(aes(shape=Name), data=filter(data, !is.na(Name)), color="black") +
     labs(x="MDS Axis 1", y="MDS Axis 2") +
     theme(legend.position="bottom",
           legend.title=element_blank(),
           axis.text=element_text(size=6),
           axis.title=element_text(size=9),
           strip.text=element_text(size=11),
           legend.text=element_text(size=11))

ggsave(outfile, g, width=5, height=5, units="in")
