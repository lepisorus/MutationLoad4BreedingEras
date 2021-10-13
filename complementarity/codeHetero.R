library(data.table)
library(cowplot)
library(dplyr)

#GERP0
df <- fread("120hybridGenotype_allChr_SorghumAllele_GERP0.txt")
head(df)
df2 <- as.data.frame(df[, 1:124])
head(df2)

summary <- apply(df2[, 5:124], 2, table)
head(summary)
summary <- as.data.frame(summary)
apply(summary, 2, length)
subset(summary, is.na(summary))
summary2 <- as.data.frame(t(summary))
head(summary2)
summary2["homo"] <- summary2$AA + summary2$CC + summary2$TT + summary2$GG
summary2["hetero"] <- rowSums(summary2[,1:17]) - summary2$X -summary2$homo
head(summary2)
summary2["ratio"] <- summary2$hetero/(summary2$homo+summary2$hetero)
head(summary2)
write.table(summary2, file="deleteriousAllelesHeteroRatioSummary.txt", quote=F, sep="\t", row.names=F)
deleteriousSummary <- summary2
deleteriousSummary["type"] <- "deleterious"

#GERP0-2
df <- fread("120hybridGenotype_allChr_SorghumAllele_GERP02.txt")
head(df)
df2 <- as.data.frame(df[, 1:124])
head(df2)

summary <- apply(df2[, 5:124], 2, table)
head(summary)
summary <- as.data.frame(summary)
apply(summary, 2, length)
subset(summary, is.na(summary))
summary2 <- as.data.frame(t(summary))
head(summary2)
summary2["homo"] <- summary2$AA + summary2$CC + summary2$TT + summary2$GG
summary2["hetero"] <- rowSums(summary2[,1:17]) - summary2$X -summary2$homo
head(summary2)
summary2["ratio"] <- summary2$hetero/(summary2$homo+summary2$hetero)
head(summary2)
write.table(summary2, file="deleteriousAllelesHeteroRatioSummary_GERP02.txt", quote=F, sep="\t", row.names=F)
deleteriousSummary <- summary2
deleteriousSummary["type"] <- "deleterious"

#GERP2-4

df <- fread("120hybridGenotype_allChr_SorghumAllele_GERP24.txt")
head(df)
df2 <- as.data.frame(df[, 1:124])
head(df2)

summary <- apply(df2[, 5:124], 2, table)
head(summary)
summary <- as.data.frame(summary)
apply(summary, 2, length)
subset(summary, is.na(summary))
summary2 <- as.data.frame(t(summary))
head(summary2)
summary2["homo"] <- summary2$AA + summary2$CC + summary2$TT + summary2$GG
summary2["hetero"] <- rowSums(summary2[,1:17]) - summary2$X -summary2$homo
head(summary2)
summary2["ratio"] <- summary2$hetero/(summary2$homo+summary2$hetero)
head(summary2)
write.table(summary2, file="deleteriousAllelesHeteroRatioSummary_GERP24.txt", quote=F, sep="\t", row.names=F)
deleteriousSummary <- summary2
deleteriousSummary["type"] <- "deleterious"


#nonDeleterious
df <- fread("120hybridGenotype_allChr_SorghumAllele_nonDeleterious.txt")
head(df)
df2 <- as.data.frame(df[, 1:124])
head(df2)

summary <- apply(df2[, 5:124], 2, table)
head(summary)
summary <- as.data.frame(summary)
apply(summary, 2, length)
subset(summary, is.na(summary))
summary2 <- as.data.frame(t(summary))
head(summary2)
summary2["homo"] <- summary2$AA + summary2$CC + summary2$TT + summary2$GG
summary2["hetero"] <- rowSums(summary2[,1:17]) - summary2$X -summary2$homo
head(summary2)
summary2["ratio"] <- summary2$hetero/(summary2$homo+summary2$hetero)
head(summary2)
write.table(summary2, file="nonDeleteriousAllelesHeteroRatioSummary.txt", quote=F, sep="\t", row.names=F)
nonDeleteriousSummary <- summary2
nonDeleteriousSummary["type"] <- "nonDeleterious"

#nonDeleterious_daf_partition

df <- fread("120hybridGenotype_allChr_SorghumAllele_nonDeleterious_daf2.txt")
head(df)
df2 <- as.data.frame(df[, 1:125])
head(df2)
breaks <- seq(0, 1, 0.1)
df2["bin"] <- cut(df2$daf, 10, include.lowest=TRUE, labels=c("0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9"))

df2_0 <- df2 %>% filter(bin==0) 
summary <- apply(df2_0[, 6:125], 2, table)
head(summary)
summary2 <- as.data.frame(t(summary))
head(summary2)
summary2["homo"] <- summary2$AA + summary2$CC + summary2$TT + summary2$GG
summary2["hetero"] <- rowSums(summary2[,1:17]) - summary2$X -summary2$homo
head(summary2)
summary2["ratio"] <- summary2$hetero/(summary2$homo+summary2$hetero)
head(summary2)
nonDeleterious_daf0_Summary <- summary2
nonDeleterious_daf0_Summary["type"] <- "nonDeleterious_daf0"

df2_0.1 <- df2 %>% filter(bin==0.1) 
summary <- apply(df2_0.1[, 6:125], 2, table)
head(summary)
summary2 <- as.data.frame(t(summary))
head(summary2)
summary2["homo"] <- summary2$AA + summary2$CC + summary2$TT + summary2$GG
summary2["hetero"] <- rowSums(summary2[,1:17]) - summary2$X -summary2$homo
head(summary2)
summary2["ratio"] <- summary2$hetero/(summary2$homo+summary2$hetero)
head(summary2)
nonDeleterious_daf0.1_Summary <- summary2
nonDeleterious_daf0.1_Summary["type"] <- "nonDeleterious_daf0.1"

df2_0.2 <- df2 %>% filter(bin==0.2) 
summary <- apply(df2_0.2[, 6:125], 2, table)
head(summary)
summary2 <- as.data.frame(t(summary))
head(summary2)
summary2["homo"] <- summary2$AA + summary2$CC + summary2$TT + summary2$GG
summary2["hetero"] <- rowSums(summary2[,1:17]) - summary2$X -summary2$homo
head(summary2)
summary2["ratio"] <- summary2$hetero/(summary2$homo+summary2$hetero)
head(summary2)
nonDeleterious_daf0.2_Summary <- summary2
nonDeleterious_daf0.2_Summary["type"] <- "nonDeleterious_daf0.2"

df2_0.3 <- df2 %>% filter(bin==0.3) 
summary <- apply(df2_0.3[, 6:125], 2, table)
head(summary)
summary2 <- as.data.frame(t(summary))
head(summary2)
summary2["homo"] <- summary2$AA + summary2$CC + summary2$TT + summary2$GG
summary2["hetero"] <- rowSums(summary2[,1:17]) - summary2$X -summary2$homo
head(summary2)
summary2["ratio"] <- summary2$hetero/(summary2$homo+summary2$hetero)
head(summary2)
nonDeleterious_daf0.3_Summary <- summary2
nonDeleterious_daf0.3_Summary["type"] <- "nonDeleterious_daf0.3"

df2_0.4 <- df2 %>% filter(bin==0.4) 
summary <- apply(df2_0.4[, 6:125], 2, table)
head(summary)
summary2 <- as.data.frame(t(summary))
head(summary2)
summary2["homo"] <- summary2$AA + summary2$CC + summary2$TT + summary2$GG
summary2["hetero"] <- rowSums(summary2[,1:17]) - summary2$X -summary2$homo
head(summary2)
summary2["ratio"] <- summary2$hetero/(summary2$homo+summary2$hetero)
head(summary2)
nonDeleterious_daf0.4_Summary <- summary2
nonDeleterious_daf0.4_Summary["type"] <- "nonDeleterious_daf0.4"

df2_0.5 <- df2 %>% filter(bin==0.5) 
summary <- apply(df2_0.5[, 6:125], 2, table)
head(summary)
summary2 <- as.data.frame(t(summary))
head(summary2)
summary2["homo"] <- summary2$AA + summary2$CC + summary2$TT + summary2$GG
summary2["hetero"] <- rowSums(summary2[,1:17]) - summary2$X -summary2$homo
head(summary2)
summary2["ratio"] <- summary2$hetero/(summary2$homo+summary2$hetero)
head(summary2)
nonDeleterious_daf0.5_Summary <- summary2
nonDeleterious_daf0.5_Summary["type"] <- "nonDeleterious_daf0.5"

df2_0.6 <- df2 %>% filter(bin==0.6) 
summary <- apply(df2_0.6[, 6:125], 2, table)
head(summary)
summary2 <- as.data.frame(t(summary))
head(summary2)
summary2["homo"] <- summary2$AA + summary2$CC + summary2$TT + summary2$GG
summary2["hetero"] <- rowSums(summary2[,1:17]) - summary2$X -summary2$homo
head(summary2)
summary2["ratio"] <- summary2$hetero/(summary2$homo+summary2$hetero)
head(summary2)
nonDeleterious_daf0.6_Summary <- summary2
nonDeleterious_daf0.6_Summary["type"] <- "nonDeleterious_daf0.6"

df2_0.7 <- df2 %>% filter(bin==0.7) 
summary <- apply(df2_0.7[, 6:125], 2, table)
head(summary)
summary2 <- as.data.frame(t(summary))
head(summary2)
summary2["homo"] <- summary2$AA + summary2$CC + summary2$TT + summary2$GG
summary2["hetero"] <- rowSums(summary2[,1:17]) - summary2$X -summary2$homo
head(summary2)
summary2["ratio"] <- summary2$hetero/(summary2$homo+summary2$hetero)
head(summary2)
nonDeleterious_daf0.7_Summary <- summary2
nonDeleterious_daf0.7_Summary["type"] <- "nonDeleterious_daf0.7"

df2_0.8 <- df2 %>% filter(bin==0.8) 
summary <- apply(df2_0.8[, 6:125], 2, table)
head(summary)
summary2 <- as.data.frame(t(summary))
head(summary2)
summary2["homo"] <- summary2$AA + summary2$CC + summary2$TT + summary2$GG
summary2["hetero"] <- rowSums(summary2[,1:17]) - summary2$X -summary2$homo
head(summary2)
summary2["ratio"] <- summary2$hetero/(summary2$homo+summary2$hetero)
head(summary2)
nonDeleterious_daf0.8_Summary <- summary2
nonDeleterious_daf0.8_Summary["type"] <- "nonDeleterious_daf0.8"

df2_0.9 <- df2 %>% filter(bin==0.9) 
summary <- apply(df2_0.9[, 6:125], 2, table)
head(summary)
summary2 <- as.data.frame(t(summary))
head(summary2)
summary2["homo"] <- summary2$AA + summary2$CC + summary2$TT + summary2$GG
summary2["hetero"] <- rowSums(summary2[,1:17]) - summary2$X -summary2$homo
head(summary2)
summary2["ratio"] <- summary2$hetero/(summary2$homo+summary2$hetero)
head(summary2)
nonDeleterious_daf0.9_Summary <- summary2
nonDeleterious_daf0.9_Summary["type"] <- "nonDeleterious_daf0.9"

nonDeleterious <- as.data.frame(rbind(nonDeleterious_daf0_Summary, nonDeleterious_daf0.1_Summary, nonDeleterious_daf0.2_Summary, nonDeleterious_daf0.3_Summary, nonDeleterious_daf0.4_Summary, nonDeleterious_daf0.5_Summary, nonDeleterious_daf0.6_Summary, nonDeleterious_daf0.7_Summary, nonDeleterious_daf0.8_Summary, nonDeleterious_daf0.9_Summary))
write.table(nonDeleterious, file="nonDeleteriousAllelesHeteroRatioSummary_dafPartition.txt", quote=F, sep="\t", row.names=F)



combine <- as.data.frame(rbind(deleteriousSummary, nonDeleteriousSummary))

ggplot(data=combine, aes(x=type, y=ratio, color=type)) +  
geom_boxplot(outlier.shape=NA) +  
geom_jitter(position=position_jitter(width=.1, height=0), aes(color=type)) + 
ylab("Ratio of heterozygous genotypes in hybrid F1s") +
theme(axis.text.x=element_text(angle=90, hjust=1)) 
#coord_cartesian(ylim = ylim1*5)
ggsave("RatioHeteroHybrids.pdf")
