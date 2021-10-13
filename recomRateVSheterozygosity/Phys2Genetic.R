# Usage: Rscript Phys2Genetic.R -i ROH300kb.20snpWindow.1het.hom.simplified -o ROH300kb.20snpWindow.1het.hom.withCM

library(optparse)
library(data.table)

option_list <- list(make_option(c('-i','--in_file'), action='store', type='character', default=NULL, help='Input file'),
make_option(c('-o','--out_file'), action='store', type='character', default=NULL, help='Output file')
                    )
opt <- parse_args(OptionParser(option_list = option_list))



snp <- fread(opt$in_file, header=T)
snp <- as.data.frame(snp)
head(snp)

map<-read.delim("Ogut2015GeneticMap0.2cMresolution.txt") # load map
head(map)
chrLengths <- c(301476924,237917468,232245527,242062272,217959525,169407836,176826311,175377492,157038028,149632204) # these are full lengths for version 3 genome. 

snp$geneticPosition <- NA
for(i in 1:nrow(snp)){
lowerIndex <- which(map$chromosome == snp$CHROM[i] & map$position <= snp$POS[i])
belowPhys <- max(map$position[lowerIndex[length(lowerIndex)]],1,na.rm=T)
belowGen <- max(map$cM[lowerIndex[length(lowerIndex)]],map$cM[which(map$chromosome==snp$CHROM[i])][1]-1,na.rm=T)
higherIndex <- which(map$chromosome == snp$CHROM[i] & map$position >= snp$POS[i])
abovePhys <- min(map$position[higherIndex[1]],chrLengths[snp$CHROM[i]],na.rm=T)
aboveGen <- min(map$cM[higherIndex[1]],map$cM[which(map$chromosome==snp$CHROM[i])][length(which(map$chromosome==snp$CHROM[i]))]+1,na.rm=T)
scale <- {snp$POS[i]-belowPhys}/{abovePhys-belowPhys}
newGen <- {aboveGen-belowGen}*scale + belowGen
snp$geneticPosition[i] <- newGen
}
head(snp)

write.table(snp, file=opt$out_file, sep="\t", quote=F, row.names=F, col.names=F)

