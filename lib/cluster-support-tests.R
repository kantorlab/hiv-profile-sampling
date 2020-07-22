library(dunn.test)

# Is cluster support significantly different between genomic regions?

support <- matrix(c(0.28,0.034,0.602,0.704,
                    0,0.894,0.658,0.878,
                    0.28,0,0.602,0.704,
                    0.746,0.59,0.962,0.98,
                    0,0,0,0.49,
                    0.784,0.956,0.18,1,
                    0.992,0.99,0.546,1,
                    0.81,0.98,1,1,
                    0.94,0.238,0.992,1,
                    0.04,0.044,0.896,0.998,
                    0.04,0.04,0.894,0.998,
                    0.994,0.782,0.93,1),
                    nrow=12,
                    byrow=TRUE,
                    dimnames = list(1:12, c("prrt", "int", "env", "wgs")))

apply(support, 2, FUN=median)

friedman.test(support)

support <- as.data.frame(support)
print(support)

dunn.test(as.list(support), method="holm")

# Is cluster support significantly different for clusters detected with Sanger vs. NGS consensus sequences?

sanger <- c(0,0,0.602,0,
            0,0.894,0.658,0,
            0,0,0.602,0,
            0,0.59,0.962,0.98,
            0,0,0,0,
            0,0.956,0.18,1,
            0.992,0.99,0,1,
            0.81,0.98,1,1,
            0.94,0,0.992,1,
            0,0,0.896,0.998,
            0,0,0.894,0.998,
            0.994,0.782,0,1)

ngs    <- c(0,0,0,0.704,
            0,0.894,0,0,
            0.28,0,0.602,0,
            0.746,0.59,0.962,0.98,
            0,0,0,0,
            0.784,0.956,0.18,1,
            0.992,0.99,0.546,1,
            0.81,0.98,1,1,
            0,0.238,0.992,1,
            0,0,0.896,0.998,
            0,0,0.894,0.998,
            0.994,0.782,0.93,1)

print(median(sanger))
print(median(ngs))

wilcox.test(sanger, ngs, paired=TRUE)

# Is the number of clusters detected in each region significantly different for Sanger vs. NGS consensus sequences?

print(sum(sanger != 0))
print(sum(ngs    != 0))

agree <- sanger[sanger != 0 & ngs != 0]
print(agree)
print(median(agree))

disagree <- c(sanger[sanger != 0 & ngs == 0], ngs[ngs != 0 & sanger == 0])
print(disagree)
print(median(disagree))