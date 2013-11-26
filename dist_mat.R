

args <- commandArgs(TRUE)
ind <- read.table(args[1], as.is=T, header=F)
M <- read.table(args[2])

rownames(M) <- paste(ind$V1)
colnames(M) <- paste(ind$V1)

write.table(M, file=args[2], sep='\t', quote=FALSE)
