# PROJET GWAS GENEOM
# Santy
# Bessoul

# Script pour corriger les matrices de distances crées par
# plink --cluster an ajoutant les id de colonnes et de lignes

args <- commandArgs(TRUE)
ind <- read.table(args[1], as.is=T, header=F)
M <- read.table(args[2])

# Ajouter les ID des colonnes et lignes
rownames(M) <- paste(ind$V1)
colnames(M) <- paste(ind$V1)

library(MASS)

min_dist_fact <- 1.75
outliers <- c()

for (iter in 1:10) {
    res <- isoMDS(as.dist(M), k = 2)
    distCentr <- sqrt((res$points[,1] - median(res$points[,1])) ** 2 +
             (res$points[,2] - median(res$points[,2])) ** 2)


    outlier <- distCentr[distCentr > min_dist_fact*mean(distCentr)]
    outlier <- data.frame(outlier)
    outlier <- rownames(outlier)

    # Plot the MDS plot before removing outliers
    if (iter == 1) {
        pdf(file = "MDS_before.pdf")

        first_outliers <- c()
        for(i in 1:length(outlier)) {
            first_index <- grep(outlier[i], colnames(M))
            first_outliers <- c(first_outliers, first_index)
        }

        plot(res$points,
            col = ifelse(distCentr>min_dist_fact*mean(distCentr),"red","black"),
            pch = ifelse(distCentr>min_dist_fact*mean(distCentr), 19, 20),
            cex = 1,
            main = "IBS Multidimensional Scaling, raw data",
            xlab = "axis 1", ylab = "axis 2")

        # Labels outliers
        text(res$points[first_outliers,],
            labels = dimnames(res$points)[[1]][first_outliers],
            pos=2)
        dev.off()
    }

    # If there are no outliers left, break the loop and proceed
    if (length(outlier) == 0) break

    for(i in 1:length(outlier)) {
        index <- grep(outlier[i], colnames(M))
        M <- M[-index, -index]
        outliers <- c(outliers, index)
    }
    iter <- iter + 1
}

# Plot MDS after removing outliers
pdf(file = "MDS_after.pdf")
plot(res$points,
    pch = 20, cex = 1,
    main = "IBS Multidimensional Scaling, after filtering",
    xlab = "axis 1", ylab = "axis 2")
dev.off()


toremove <- cbind(dimnames(res$points)[[1]][outliers], 1)
colnames(toremove) <- c("FID", "IID")


write.table(toremove,"toremove.txt",row.names=FALSE, quote=F)
