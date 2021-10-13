library("data.table")

p <- read.table("BLUP_traits.sorted.txt", header=TRUE)
write.table(p, "BLUP_traits.sorted2.txt", sep="\t", row.names=FALSE, quote=FALSE)

for(i in 1:10){
  infile <- paste0("chr", i, ".SNPs.t.txt")
  c10 <- fread(infile, data.table=FALSE)
#  c10[c10 == 9] <- 0
  nSNP <- ncol(c10)
  c10$sampleID <- p$sampleID
  c10 <- c10[, c(ncol(c10), 1:nSNP)]
  fwrite(c10, paste0("chr", i, ".SNPs.t.sampleid.txt"), sep="\t")
}

