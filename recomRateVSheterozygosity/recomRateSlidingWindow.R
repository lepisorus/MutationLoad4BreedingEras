##making sliding windows
library(plyr)
#library(ggplot2)
library(data.table)
library(stringr)
library(tidyverse)
library(expss)

df <- fread("allChr.120Hybrid.withAncestralAllele.GERP.withCM.sorted.recode.txt", header=T)
df <- as.data.frame(df)
head(df)

#user-defined window size and step size
slidingWindow <- function(df) {

win = 1000000;
step = 1000000;
len=max(df$pos);
breaks=seq(1, len, step);
start <- breaks;
end <- start + win;
chr = rep(unique(df$chr), length(breaks));

#df["bins"] <- cut(df$pos, breaks=c(breaks, len), include.lowest=TRUE)
#df.slim <- df[, c(2:3,128)]
#df.new1 <- df.slim %>% group_by(bins) %>% transmute(genDis=max(cM)-min(cM), phyDis=max(pos)-min(pos), recRate=genDis/phyDis);
#df.new2 <- unique(df.new1)

genDis=rep(NA, length(breaks));
phyDis=rep(NA, length(breaks));
recRate=rep(NA, length(breaks));

hetero=rep(NA, length(breaks)*120);
missing=rep(NA, length(breaks)*120);
percentHetero=rep(NA, length(breaks)*120);


for (i in 1:length(breaks)) {
ipos=NULL;
ipos=which(df$pos>=start[i] & df$pos<=end[i]);

if(length(df[ipos, 1]) > 0) { #if
genDis[i] = max(df[ipos, 3], na.rm=T) - min(df[ipos, 3], na.rm=T);
phyDis[i] = max(df[ipos, 2], na.rm=T) - min(df[ipos, 2], na.rm=T);
recRate[i] = genDis[i]/phyDis[i];
hetero[((i-1)*120+1):(i*120)]=  apply(df[ipos, 8:127], 2, function(x) sum(x, na.rm=T));
missing[((i-1)*120+1):(i*120)]= apply(df[ipos, 8:127], 2, function(x) sum(is.na(x)));
percentHetero[((i-1)*120+1):(i*120)]= hetero[((i-1)*120+1):(i*120)]/(length(df[ipos, 1]) - missing[((i-1)*120+1):(i*120)]);
}#if

else {
genDis[i]=NA;
phyDis[i] =NA;
recRate[i]=NA;
hetero[((i-1)*120+1):(i*120)]=  NA;
missing[((i-1)*120+1):(i*120)]= NA;
percentHetero[((i-1)*120+1):(i*120)]= NA;
}#else

}#for

df.new2 <- matrix(percentHetero, ncol=120, nrow=length(breaks), byrow=TRUE);
df.new2 <- as.data.frame(df.new2);

return(data.frame(chr, start, end, genDis, phyDis, recRate, df.new2));

} #function


df.plyr <- ddply(df, .(chr), slidingWindow)
write.table(df.plyr, file="allChr.120Hybrid.slidingWindow1mb.cM.percHetero.txt", quote=F, sep="\t", row.names=F, col.names=F)


