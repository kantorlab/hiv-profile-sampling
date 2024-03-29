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

g <- ggplot(data, aes(x=X, y=Y, color=Gene)) +
     geom_point(data=filter(data,  is.na(Name)), shape=3, size=3, alpha=0.1) +
     geom_point(aes(shape=Name), data=filter(data, !is.na(Name)), size=3) +
     facet_grid(. ~ Gene) +
     xlim(-0.16, 0.16) +
     ylim(-0.16, 0.16) +
     labs(x="MDS Axis 1", y="MDS Axis 2") +
     scale_color_manual(values=c("#DF8F44", "#00A1D5", "#374E55", "#B24745")) +
     theme_minimal() +
     theme(legend.position="bottom",
           legend.title=element_blank(),
           panel.grid.minor=element_line(size=0.25),
           panel.grid.major=element_line(size=0.25),
           axis.text=element_text(color="black"),
           axis.text.x=element_text(angle=90, size=7, vjust=0.5, hjust=0.95),
           axis.text.y=element_text(size=7))

pdf(outfile, width=6.5, height=3)
print(g)
dev.off()
