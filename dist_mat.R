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

min_dist_fact <- 1.8
outliers <- c()
for (iter in 1:10) {
    res <- isoMDS(as.dist(M), k = 2)
    distCentr <- sqrt((res$points[,1] - median(res$points[,1])) ** 2 +
             (res$points[,2] - median(res$points[,2])) ** 2)


    outlier <- distCentr[distCentr > min_dist_fact*mean(distCentr)]

    if (length(outlier) == 0) break

    outlier <- data.frame(outlier)
    outlier <- rownames(outlier)

    for(i in 1:length(outlier)) {
        index <- grep(outlier[i], colnames(M))
        M <- M[-index, -index]
        outliers <- c(outliers, index)
    }
    iter <- iter + 1
}

toremove <- cbind(dimnames(res$points)[[1]][outliers],dimnames(res$points)[[1]][outliers])
colnames(toremove)<-c("ID","FID")
print("wouhou")


write.table(toremove,"toremove.txt",row.names=FALSE,quote=F)
