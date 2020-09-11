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

g <- ggplot(data, aes(x=X, y=Y, color=Gene)) +
     geom_point(data=filter(data,  is.na(Name)), shape=3, alpha=0.2, size=3) +
     geom_point(aes(shape=Name), data=filter(data, !is.na(Name)), size=3) +
     scale_color_manual(values=c("#DF8F44", "#00A1D5", "#374E55", "#B24745")) +
     labs(x="MDS Axis 1", y="MDS Axis 2") +
     theme(legend.position="bottom",
           legend.title=element_blank(),
           axis.text=element_text(size=6),
           axis.title=element_text(size=9),
           strip.text=element_text(size=11),
           legend.text=element_text(size=11))

pdf(outfile, width=6.5, height=6.5)
print(g)
dev.off()
