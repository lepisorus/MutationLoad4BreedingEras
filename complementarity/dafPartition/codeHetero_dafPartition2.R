library(data.table)
library(dplyr)
library(tidyverse)


#GERP02

df <- read.delim("120hybridGenotype_allChr_SorghumAllele_GERP02_daf2.txt")
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
summary2["sample"] <- row.names(summary2)
summary2["homo"] <- summary2$AA + summary2$CC + summary2$TT + summary2$GG
summary2["hetero"] <- rowSums(summary2[,1:17]) - summary2$X -summary2$homo
head(summary2)
summary2["ratio"] <- summary2$hetero/(summary2$homo+summary2$hetero)
head(summary2)
summary(summary2$ratio)
GERP02_daf0_Summary <- summary2
GERP02_daf0_Summary["daf"] <- "0"

df2_0.1 <- df2 %>% filter(bin==0.1) 
summary <- apply(df2_0.1[, 6:125], 2, table)
head(summary)
test <- do.call(rbind, lapply(summary, as.data.frame))
test["sample"] <- gsub("[.].*$", "", row.names(test))
head(test)
test_wide <- spread(test, Var1, Freq, fill=0)
row.names(test_wide) <- test_wide$sample
head(test_wide)
summary2 <- test_wide[, 2:18]
head(summary2)
summary2["sample"] <- row.names(summary2)
head(summary2)
summary2["homo"] <- summary2$AA + summary2$CC + summary2$TT + summary2$GG
summary2["hetero"] <- rowSums(summary2[,1:17]) - summary2$X -summary2$homo
head(summary2)
summary2["ratio"] <- summary2$hetero/(summary2$homo+summary2$hetero)
head(summary2)
summary(summary2$ratio)
GERP02_daf0.1_Summary <- summary2
GERP02_daf0.1_Summary["daf"] <- "0.1"

df2_0.2 <- df2 %>% filter(bin==0.2) 
summary <- apply(df2_0.2[, 6:125], 2, table)
head(summary)
summary2 <- as.data.frame(t(summary))
head(summary2)
summary2["sample"] <- row.names(summary2)
summary2["homo"] <- summary2$AA + summary2$CC + summary2$TT + summary2$GG
summary2["hetero"] <- rowSums(summary2[,1:17]) - summary2$X -summary2$homo
head(summary2)
summary2["ratio"] <- summary2$hetero/(summary2$homo+summary2$hetero)
head(summary2)
summary(summary2$ratio)
GERP02_daf0.2_Summary <- summary2
GERP02_daf0.2_Summary["daf"] <- "0.2"

df2_0.3 <- df2 %>% filter(bin==0.3) 
summary <- apply(df2_0.3[, 6:125], 2, table)
head(summary)
summary2 <- as.data.frame(t(summary))
head(summary2)
summary2["sample"] <- row.names(summary2)
summary2["homo"] <- summary2$AA + summary2$CC + summary2$TT + summary2$GG
summary2["hetero"] <- rowSums(summary2[,1:17]) - summary2$X -summary2$homo
head(summary2)
summary2["ratio"] <- summary2$hetero/(summary2$homo+summary2$hetero)
head(summary2)
summary(summary2$ratio)
GERP02_daf0.3_Summary <- summary2
GERP02_daf0.3_Summary["daf"] <- "0.3"


df2_0.4 <- df2 %>% filter(bin==0.4) 
summary <- apply(df2_0.4[, 6:125], 2, table)
head(summary)
summary2 <- as.data.frame(t(summary))
head(summary2)
summary2["sample"] <- row.names(summary2)
summary2["homo"] <- summary2$AA + summary2$CC + summary2$TT + summary2$GG
summary2["hetero"] <- rowSums(summary2[,1:17]) - summary2$X -summary2$homo
head(summary2)
summary2["ratio"] <- summary2$hetero/(summary2$homo+summary2$hetero)
head(summary2)
summary(summary2$ratio)
GERP02_daf0.4_Summary <- summary2
GERP02_daf0.4_Summary["daf"] <- "0.4"

df2_0.5 <- df2 %>% filter(bin==0.5) 
summary <- apply(df2_0.5[, 6:125], 2, table)
head(summary)
summary2 <- as.data.frame(t(summary))
head(summary2)
summary2["sample"] <- row.names(summary2)
summary2["homo"] <- summary2$AA + summary2$CC + summary2$TT + summary2$GG
summary2["hetero"] <- rowSums(summary2[,1:17]) - summary2$X -summary2$homo
head(summary2)
summary2["ratio"] <- summary2$hetero/(summary2$homo+summary2$hetero)
head(summary2)
summary(summary2$ratio)
GERP02_daf0.5_Summary <- summary2
GERP02_daf0.5_Summary["daf"] <- "0.5"

df2_0.6 <- df2 %>% filter(bin==0.6) 
summary <- apply(df2_0.6[, 6:125], 2, table)
head(summary)
summary2 <- as.data.frame(t(summary))
head(summary2)
summary2["sample"] <- row.names(summary2)
summary2["homo"] <- summary2$AA + summary2$CC + summary2$TT + summary2$GG
summary2["hetero"] <- rowSums(summary2[,1:17]) - summary2$X -summary2$homo
head(summary2)
summary2["ratio"] <- summary2$hetero/(summary2$homo+summary2$hetero)
head(summary2)
summary(summary2$ratio)
GERP02_daf0.6_Summary <- summary2
GERP02_daf0.6_Summary["daf"] <- "0.6"

df2_0.7 <- df2 %>% filter(bin==0.7) 
summary <- apply(df2_0.7[, 6:125], 2, table)
head(summary)
summary2 <- as.data.frame(t(summary))
head(summary2)
summary2["sample"] <- row.names(summary2)
summary2["homo"] <- summary2$AA + summary2$CC + summary2$TT + summary2$GG
summary2["hetero"] <- rowSums(summary2[,1:17]) - summary2$X -summary2$homo
head(summary2)
summary2["ratio"] <- summary2$hetero/(summary2$homo+summary2$hetero)
head(summary2)
summary(summary2$ratio)
GERP02_daf0.7_Summary <- summary2
GERP02_daf0.7_Summary["daf"] <- "0.7"

df2_0.8 <- df2 %>% filter(bin==0.8) 
summary <- apply(df2_0.8[, 6:125], 2, table)
head(summary)
summary2 <- as.data.frame(t(summary))
head(summary2)
summary2["sample"] <- row.names(summary2)
head(summary2)
summary2["homo"] <- summary2$AA + summary2$CC + summary2$TT + summary2$GG
summary2["hetero"] <- rowSums(summary2[,1:17]) - summary2$X -summary2$homo
head(summary2)
summary2["ratio"] <- summary2$hetero/(summary2$homo+summary2$hetero)
head(summary2)
summary(summary2$ratio)
GERP02_daf0.8_Summary <- summary2
GERP02_daf0.8_Summary["daf"] <- "0.8"

df2_0.9 <- df2 %>% filter(bin==0.9) 
summary <- apply(df2_0.9[, 6:125], 2, table)
head(summary)
summary2 <- as.data.frame(t(summary))
head(summary2)
summary2["sample"] <- row.names(summary2)
summary2["homo"] <- summary2$AA + summary2$CC + summary2$TT + summary2$GG
summary2["hetero"] <- rowSums(summary2[,1:17]) - summary2$X -summary2$homo
head(summary2)
summary2["ratio"] <- summary2$hetero/(summary2$homo+summary2$hetero)
head(summary2)
summary(summary2$ratio)
GERP02_daf0.9_Summary <- summary2
GERP02_daf0.9_Summary["daf"] <- "0.9"

GERP02 <- as.data.frame(rbind(GERP02_daf0_Summary, GERP02_daf0.1_Summary, GERP02_daf0.2_Summary, GERP02_daf0.3_Summary, GERP02_daf0.4_Summary, GERP02_daf0.5_Summary, GERP02_daf0.6_Summary, GERP02_daf0.7_Summary, GERP02_daf0.8_Summary, GERP02_daf0.9_Summary))
write.table(GERP02, file="GERP02_HeteroRatioSummary_dafPartition.txt", quote=F, sep="\t", row.names=T)

#GERP24

df <- read.delim("120hybridGenotype_allChr_SorghumAllele_GERP24_daf2.txt")
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
summary2["sample"] <- row.names(summary2)
summary2["homo"] <- summary2$AA + summary2$CC + summary2$TT + summary2$GG
summary2["hetero"] <- rowSums(summary2[,1:17]) - summary2$X -summary2$homo
head(summary2)
summary2["ratio"] <- summary2$hetero/(summary2$homo+summary2$hetero)
head(summary2)
GERP24_daf0_Summary <- summary2
GERP24_daf0_Summary["daf"] <- "0"

df2_0.1 <- df2 %>% filter(bin==0.1) 
summary <- apply(df2_0.1[, 6:125], 2, table)
head(summary)
summary2 <- as.data.frame(t(summary))
head(summary2)
summary2["sample"] <- row.names(summary2)
summary2["homo"] <- summary2$AA + summary2$CC + summary2$TT + summary2$GG
summary2["hetero"] <- rowSums(summary2[,1:17]) - summary2$X -summary2$homo
head(summary2)
summary2["ratio"] <- summary2$hetero/(summary2$homo+summary2$hetero)
head(summary2)
GERP24_daf0.1_Summary <- summary2
GERP24_daf0.1_Summary["daf"] <- "0.1"

df2_0.2 <- df2 %>% filter(bin==0.2) 
summary <- apply(df2_0.2[, 6:125], 2, table)
head(summary)
summary2 <- as.data.frame(t(summary))
head(summary2)
summary2["sample"] <- row.names(summary2)
summary2["homo"] <- summary2$AA + summary2$CC + summary2$TT + summary2$GG
summary2["hetero"] <- rowSums(summary2[,1:17]) - summary2$X -summary2$homo
head(summary2)
summary2["ratio"] <- summary2$hetero/(summary2$homo+summary2$hetero)
head(summary2)
GERP24_daf0.2_Summary <- summary2
GERP24_daf0.2_Summary["daf"] <- "0.2"

df2_0.3 <- df2 %>% filter(bin==0.3) 
summary <- apply(df2_0.3[, 6:125], 2, table)
head(summary)
summary2 <- as.data.frame(t(summary))
head(summary2)
summary2["sample"] <- row.names(summary2)
summary2["homo"] <- summary2$AA + summary2$CC + summary2$TT + summary2$GG
summary2["hetero"] <- rowSums(summary2[,1:17]) - summary2$X -summary2$homo
head(summary2)
summary2["ratio"] <- summary2$hetero/(summary2$homo+summary2$hetero)
head(summary2)
GERP24_daf0.3_Summary <- summary2
GERP24_daf0.3_Summary["daf"] <- "0.3"

df2_0.4 <- df2 %>% filter(bin==0.4) 
summary <- apply(df2_0.4[, 6:125], 2, table)
head(summary)
summary2 <- as.data.frame(t(summary))
head(summary2)
summary2["sample"] <- row.names(summary2)
summary2["homo"] <- summary2$AA + summary2$CC + summary2$TT + summary2$GG
summary2["hetero"] <- rowSums(summary2[,1:17]) - summary2$X -summary2$homo
head(summary2)
summary2["ratio"] <- summary2$hetero/(summary2$homo+summary2$hetero)
head(summary2)
GERP24_daf0.4_Summary <- summary2
GERP24_daf0.4_Summary["daf"] <- "0.4"

df2_0.5 <- df2 %>% filter(bin==0.5) 
summary <- apply(df2_0.5[, 6:125], 2, table)
head(summary)
summary2 <- as.data.frame(t(summary))
head(summary2)
summary2["sample"] <- row.names(summary2)
summary2["homo"] <- summary2$AA + summary2$CC + summary2$TT + summary2$GG
summary2["hetero"] <- rowSums(summary2[,1:17]) - summary2$X -summary2$homo
head(summary2)
summary2["ratio"] <- summary2$hetero/(summary2$homo+summary2$hetero)
head(summary2)
GERP24_daf0.5_Summary <- summary2
GERP24_daf0.5_Summary["daf"] <- "0.5"

df2_0.6 <- df2 %>% filter(bin==0.6) 
summary <- apply(df2_0.6[, 6:125], 2, table)
head(summary)
summary2 <- as.data.frame(t(summary))
head(summary2)
summary2["sample"] <- row.names(summary2)
summary2["homo"] <- summary2$AA + summary2$CC + summary2$TT + summary2$GG
summary2["hetero"] <- rowSums(summary2[,1:17]) - summary2$X -summary2$homo
head(summary2)
summary2["ratio"] <- summary2$hetero/(summary2$homo+summary2$hetero)
head(summary2)
GERP24_daf0.6_Summary <- summary2
GERP24_daf0.6_Summary["daf"] <- "0.6"

df2_0.7 <- df2 %>% filter(bin==0.7) 
summary <- apply(df2_0.7[, 6:125], 2, table)
head(summary)
summary2 <- as.data.frame(t(summary))
head(summary2)
summary2["sample"] <- row.names(summary2)
summary2["homo"] <- summary2$AA + summary2$CC + summary2$TT + summary2$GG
summary2["hetero"] <- rowSums(summary2[,1:17]) - summary2$X -summary2$homo
head(summary2)
summary2["ratio"] <- summary2$hetero/(summary2$homo+summary2$hetero)
head(summary2)
GERP24_daf0.7_Summary <- summary2
GERP24_daf0.7_Summary["daf"] <- "0.7"

df2_0.8 <- df2 %>% filter(bin==0.8) 
summary <- apply(df2_0.8[, 6:125], 2, table)
head(summary)
summary2 <- as.data.frame(t(summary))
head(summary2)
summary2["sample"] <- row.names(summary2)
summary2["homo"] <- summary2$AA + summary2$CC + summary2$TT + summary2$GG
summary2["hetero"] <- rowSums(summary2[,1:17]) - summary2$X -summary2$homo
head(summary2)
summary2["ratio"] <- summary2$hetero/(summary2$homo+summary2$hetero)
head(summary2)
GERP24_daf0.8_Summary <- summary2
GERP24_daf0.8_Summary["daf"] <- "0.8"

df2_0.9 <- df2 %>% filter(bin==0.9) 
summary <- apply(df2_0.9[, 6:125], 2, table)
head(summary)
summary2 <- as.data.frame(t(summary))
head(summary2)
summary2["sample"] <- row.names(summary2)
summary2["homo"] <- summary2$AA + summary2$CC + summary2$TT + summary2$GG
summary2["hetero"] <- rowSums(summary2[,1:17]) - summary2$X -summary2$homo
head(summary2)
summary2["ratio"] <- summary2$hetero/(summary2$homo+summary2$hetero)
head(summary2)
GERP24_daf0.9_Summary <- summary2
GERP24_daf0.9_Summary["daf"] <- "0.9"

GERP24 <- as.data.frame(rbind(GERP24_daf0_Summary, GERP24_daf0.1_Summary, GERP24_daf0.2_Summary, GERP24_daf0.3_Summary, GERP24_daf0.4_Summary, GERP24_daf0.5_Summary, GERP24_daf0.6_Summary, GERP24_daf0.7_Summary, GERP24_daf0.8_Summary, GERP24_daf0.9_Summary))
write.table(GERP24, file="GERP24_HeteroRatioSummary_dafPartition.txt", quote=F, sep="\t", row.names=T)

#GERP4

df <- read.delim("120hybridGenotype_allChr_SorghumAllele_GERP4_daf2.txt")
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
summary2["sample"] <- row.names(summary2)
summary2["homo"] <- summary2$AA + summary2$CC + summary2$TT + summary2$GG
summary2["hetero"] <- rowSums(summary2[,1:17]) - summary2$X -summary2$homo
head(summary2)
summary2["ratio"] <- summary2$hetero/(summary2$homo+summary2$hetero)
head(summary2)
GERP4_daf0_Summary <- summary2
GERP4_daf0_Summary["daf"] <- "0"

df2_0.1 <- df2 %>% filter(bin==0.1) 
summary <- apply(df2_0.1[, 6:125], 2, table)
head(summary)
summary2 <- as.data.frame(t(summary))
head(summary2)
summary2["sample"] <- row.names(summary2)
summary2["homo"] <- summary2$AA + summary2$CC + summary2$TT + summary2$GG
summary2["hetero"] <- rowSums(summary2[,1:17]) - summary2$X -summary2$homo
head(summary2)
summary2["ratio"] <- summary2$hetero/(summary2$homo+summary2$hetero)
head(summary2)
GERP4_daf0.1_Summary <- summary2
GERP4_daf0.1_Summary["daf"] <- "0.1"

df2_0.2 <- df2 %>% filter(bin==0.2) 
summary <- apply(df2_0.2[, 6:125], 2, table)
head(summary)
summary2 <- as.data.frame(t(summary))
head(summary2)
summary2["sample"] <- row.names(summary2)
summary2["homo"] <- summary2$AA + summary2$CC + summary2$TT + summary2$GG
summary2["hetero"] <- rowSums(summary2[,1:17]) - summary2$X -summary2$homo
head(summary2)
summary2["ratio"] <- summary2$hetero/(summary2$homo+summary2$hetero)
head(summary2)
GERP4_daf0.2_Summary <- summary2
GERP4_daf0.2_Summary["daf"] <- "0.2"

df2_0.3 <- df2 %>% filter(bin==0.3) 
summary <- apply(df2_0.3[, 6:125], 2, table)
head(summary)
summary2 <- as.data.frame(t(summary))
head(summary2)
summary2["sample"] <- row.names(summary2)
summary2["homo"] <- summary2$AA + summary2$CC + summary2$TT + summary2$GG
summary2["hetero"] <- rowSums(summary2[,1:17]) - summary2$X -summary2$homo
head(summary2)
summary2["ratio"] <- summary2$hetero/(summary2$homo+summary2$hetero)
head(summary2)
GERP4_daf0.3_Summary <- summary2
GERP4_daf0.3_Summary["daf"] <- "0.3"

df2_0.4 <- df2 %>% filter(bin==0.4) 
summary <- apply(df2_0.4[, 6:125], 2, table)
head(summary)
summary2 <- as.data.frame(t(summary))
head(summary2)
summary2["sample"] <- row.names(summary2)
summary2["homo"] <- summary2$AA + summary2$CC + summary2$TT + summary2$GG
summary2["hetero"] <- rowSums(summary2[,1:17]) - summary2$X -summary2$homo
head(summary2)
summary2["ratio"] <- summary2$hetero/(summary2$homo+summary2$hetero)
head(summary2)
GERP4_daf0.4_Summary <- summary2
GERP4_daf0.4_Summary["daf"] <- "0.4"

df2_0.5 <- df2 %>% filter(bin==0.5) 
summary <- apply(df2_0.5[, 6:125], 2, table)
head(summary)
summary2 <- as.data.frame(t(summary))
head(summary2)
summary2["sample"] <- row.names(summary2)
summary2["homo"] <- summary2$AA + summary2$CC + summary2$TT + summary2$GG
summary2["hetero"] <- rowSums(summary2[,1:17]) - summary2$X -summary2$homo
head(summary2)
summary2["ratio"] <- summary2$hetero/(summary2$homo+summary2$hetero)
head(summary2)
GERP4_daf0.5_Summary <- summary2
GERP4_daf0.5_Summary["daf"] <- "0.5"

##missing values for some items
df2_0.6 <- df2 %>% filter(bin==0.6) 
summary <- apply(df2_0.6[, 6:125], 2, table)
head(summary)
test <- do.call(rbind, lapply(summary, as.data.frame))
test["sample"] <- gsub("[.].*$", "", row.names(test))
head(test)
test_wide <- spread(test, Var1, Freq, fill=0)
row.names(test_wide) <- test_wide$sample
head(test_wide)
summary2 <- test_wide[, 2:18]
head(summary2)
summary2["sample"] <- row.names(summary2)
summary2["homo"] <- summary2$AA + summary2$CC + summary2$TT + summary2$GG
summary2["hetero"] <- rowSums(summary2[,1:17]) - summary2$X -summary2$homo
head(summary2)
summary2["ratio"] <- summary2$hetero/(summary2$homo+summary2$hetero)
head(summary2)
GERP4_daf0.6_Summary <- summary2
GERP4_daf0.6_Summary["daf"] <- "0.6"

df2_0.7 <- df2 %>% filter(bin==0.7) 
summary <- apply(df2_0.7[, 6:125], 2, table)
head(summary)
test <- do.call(rbind, lapply(summary, as.data.frame))
test["sample"] <- gsub("[.].*$", "", row.names(test))
head(test)
test_wide <- spread(test, Var1, Freq, fill=0)
row.names(test_wide) <- test_wide$sample
head(test_wide)
summary2 <- test_wide[, 2:18]
head(summary2)
summary2["sample"] <- row.names(summary2)
summary2["homo"] <- summary2$AA + summary2$CC + summary2$TT + summary2$GG
summary2["hetero"] <- rowSums(summary2[,1:17]) - summary2$X -summary2$homo
head(summary2)
summary2["ratio"] <- summary2$hetero/(summary2$homo+summary2$hetero)
head(summary2)
GERP4_daf0.7_Summary <- summary2
GERP4_daf0.7_Summary["daf"] <- "0.7"

df2_0.8 <- df2 %>% filter(bin==0.8) 
summary <- apply(df2_0.8[, 6:125], 2, table)
head(summary)
test <- do.call(rbind, lapply(summary, as.data.frame))
test["sample"] <- gsub("[.].*$", "", row.names(test))
head(test)
test_wide <- spread(test, Var1, Freq, fill=0)
row.names(test_wide) <- test_wide$sample
head(test_wide)
summary2 <- test_wide[, 2:18]
head(summary2)
summary2["sample"] <- row.names(summary2)
summary2["homo"] <- summary2$AA + summary2$CC + summary2$TT + summary2$GG
summary2["hetero"] <- rowSums(summary2[,1:17]) - summary2$X -summary2$homo
head(summary2)
summary2["ratio"] <- summary2$hetero/(summary2$homo+summary2$hetero)
head(summary2)
GERP4_daf0.8_Summary <- summary2
GERP4_daf0.8_Summary["daf"] <- "0.8"

df2_0.9 <- df2 %>% filter(bin==0.9) 
summary <- apply(df2_0.9[, 6:125], 2, table)
head(summary)
test <- do.call(rbind, lapply(summary, as.data.frame))
test["sample"] <- gsub("[.].*$", "", row.names(test))
head(test)
test_wide <- spread(test, Var1, Freq, fill=0)
row.names(test_wide) <- test_wide$sample
head(test_wide)
summary2 <- test_wide[, 2:18]
summary2["sample"] <- row.names(summary2)
summary2["homo"] <- summary2$AA + summary2$CC + summary2$TT + summary2$GG
summary2["hetero"] <- rowSums(summary2[,1:17]) - summary2$X -summary2$homo
head(summary2)
summary2["ratio"] <- summary2$hetero/(summary2$homo+summary2$hetero)
head(summary2)
GERP4_daf0.9_Summary <- summary2
GERP4_daf0.9_Summary["daf"] <- "0.9"

GERP4 <- as.data.frame(rbind(GERP4_daf0_Summary, GERP4_daf0.1_Summary, GERP4_daf0.2_Summary, GERP4_daf0.3_Summary, GERP4_daf0.4_Summary, GERP4_daf0.5_Summary, GERP4_daf0.6_Summary, GERP4_daf0.7_Summary, GERP4_daf0.8_Summary, GERP4_daf0.9_Summary))
write.table(GERP4, file="GERP4_HeteroRatioSummary_dafPartition.txt", quote=F, sep="\t", row.names=T)


#nonDeleterious_daf_partition 重新跑

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
summary2["sample"] <- row.names(summary2)
summary2["homo"] <- summary2$AA + summary2$CC + summary2$TT + summary2$GG
summary2["hetero"] <- rowSums(summary2[,1:17]) - summary2$X -summary2$homo
head(summary2)
summary2["ratio"] <- summary2$hetero/(summary2$homo+summary2$hetero)
head(summary2)
nonDeleterious_daf0_Summary <- summary2
nonDeleterious_daf0_Summary["daf"] <- "0"

df2_0.1 <- df2 %>% filter(bin==0.1) 
summary <- apply(df2_0.1[, 6:125], 2, table)
head(summary)
summary2 <- as.data.frame(t(summary))
head(summary2)
summary2["sample"] <- row.names(summary2)
summary2["homo"] <- summary2$AA + summary2$CC + summary2$TT + summary2$GG
summary2["hetero"] <- rowSums(summary2[,1:17]) - summary2$X -summary2$homo
head(summary2)
summary2["ratio"] <- summary2$hetero/(summary2$homo+summary2$hetero)
head(summary2)
nonDeleterious_daf0.1_Summary <- summary2
nonDeleterious_daf0.1_Summary["daf"] <- "0.1"

df2_0.2 <- df2 %>% filter(bin==0.2) 
summary <- apply(df2_0.2[, 6:125], 2, table)
head(summary)
summary2 <- as.data.frame(t(summary))
head(summary2)
summary2["sample"] <- row.names(summary2)
summary2["homo"] <- summary2$AA + summary2$CC + summary2$TT + summary2$GG
summary2["hetero"] <- rowSums(summary2[,1:17]) - summary2$X -summary2$homo
head(summary2)
summary2["ratio"] <- summary2$hetero/(summary2$homo+summary2$hetero)
head(summary2)
nonDeleterious_daf0.2_Summary <- summary2
nonDeleterious_daf0.2_Summary["daf"] <- "0.2"

df2_0.3 <- df2 %>% filter(bin==0.3) 
summary <- apply(df2_0.3[, 6:125], 2, table)
head(summary)
summary2 <- as.data.frame(t(summary))
head(summary2)
summary2["sample"] <- row.names(summary2)
summary2["homo"] <- summary2$AA + summary2$CC + summary2$TT + summary2$GG
summary2["hetero"] <- rowSums(summary2[,1:17]) - summary2$X -summary2$homo
head(summary2)
summary2["ratio"] <- summary2$hetero/(summary2$homo+summary2$hetero)
head(summary2)
nonDeleterious_daf0.3_Summary <- summary2
nonDeleterious_daf0.3_Summary["daf"] <- "0.3"

df2_0.4 <- df2 %>% filter(bin==0.4) 
summary <- apply(df2_0.4[, 6:125], 2, table)
head(summary)
summary2 <- as.data.frame(t(summary))
head(summary2)
summary2["sample"] <- row.names(summary2)
summary2["homo"] <- summary2$AA + summary2$CC + summary2$TT + summary2$GG
summary2["hetero"] <- rowSums(summary2[,1:17]) - summary2$X -summary2$homo
head(summary2)
summary2["ratio"] <- summary2$hetero/(summary2$homo+summary2$hetero)
head(summary2)
nonDeleterious_daf0.4_Summary <- summary2
nonDeleterious_daf0.4_Summary["daf"] <- "0.4"

df2_0.5 <- df2 %>% filter(bin==0.5) 
summary <- apply(df2_0.5[, 6:125], 2, table)
head(summary)
summary2 <- as.data.frame(t(summary))
head(summary2)
summary2["sample"] <- row.names(summary2)
summary2["homo"] <- summary2$AA + summary2$CC + summary2$TT + summary2$GG
summary2["hetero"] <- rowSums(summary2[,1:17]) - summary2$X -summary2$homo
head(summary2)
summary2["ratio"] <- summary2$hetero/(summary2$homo+summary2$hetero)
head(summary2)
nonDeleterious_daf0.5_Summary <- summary2
nonDeleterious_daf0.5_Summary["daf"] <- "0.5"

df2_0.6 <- df2 %>% filter(bin==0.6) 
summary <- apply(df2_0.6[, 6:125], 2, table)
head(summary)
summary2 <- as.data.frame(t(summary))
head(summary2)
summary2["sample"] <- row.names(summary2)
summary2["homo"] <- summary2$AA + summary2$CC + summary2$TT + summary2$GG
summary2["hetero"] <- rowSums(summary2[,1:17]) - summary2$X -summary2$homo
head(summary2)
summary2["ratio"] <- summary2$hetero/(summary2$homo+summary2$hetero)
head(summary2)
nonDeleterious_daf0.6_Summary <- summary2
nonDeleterious_daf0.6_Summary["daf"] <- "0.6"

df2_0.7 <- df2 %>% filter(bin==0.7) 
summary <- apply(df2_0.7[, 6:125], 2, table)
head(summary)
summary2 <- as.data.frame(t(summary))
head(summary2)
summary2["sample"] <- row.names(summary2)
summary2["homo"] <- summary2$AA + summary2$CC + summary2$TT + summary2$GG
summary2["hetero"] <- rowSums(summary2[,1:17]) - summary2$X -summary2$homo
head(summary2)
summary2["ratio"] <- summary2$hetero/(summary2$homo+summary2$hetero)
head(summary2)
nonDeleterious_daf0.7_Summary <- summary2
nonDeleterious_daf0.7_Summary["daf"] <- "0.7"

df2_0.8 <- df2 %>% filter(bin==0.8) 
summary <- apply(df2_0.8[, 6:125], 2, table)
head(summary)
summary2 <- as.data.frame(t(summary))
head(summary2)
summary2["sample"] <- row.names(summary2)
summary2["homo"] <- summary2$AA + summary2$CC + summary2$TT + summary2$GG
summary2["hetero"] <- rowSums(summary2[,1:17]) - summary2$X -summary2$homo
head(summary2)
summary2["ratio"] <- summary2$hetero/(summary2$homo+summary2$hetero)
head(summary2)
nonDeleterious_daf0.8_Summary <- summary2
nonDeleterious_daf0.8_Summary["daf"] <- "0.8"

df2_0.9 <- df2 %>% filter(bin==0.9) 
summary <- apply(df2_0.9[, 6:125], 2, table)
head(summary)
summary2 <- as.data.frame(t(summary))
head(summary2)
summary2["sample"] <- row.names(summary2)
summary2["homo"] <- summary2$AA + summary2$CC + summary2$TT + summary2$GG
summary2["hetero"] <- rowSums(summary2[,1:17]) - summary2$X -summary2$homo
head(summary2)
summary2["ratio"] <- summary2$hetero/(summary2$homo+summary2$hetero)
head(summary2)
nonDeleterious_daf0.9_Summary <- summary2
nonDeleterious_daf0.9_Summary["daf"] <- "0.9"

nonDeleterious <- as.data.frame(rbind(nonDeleterious_daf0_Summary, nonDeleterious_daf0.1_Summary, nonDeleterious_daf0.2_Summary, nonDeleterious_daf0.3_Summary, nonDeleterious_daf0.4_Summary, nonDeleterious_daf0.5_Summary, nonDeleterious_daf0.6_Summary, nonDeleterious_daf0.7_Summary, nonDeleterious_daf0.8_Summary, nonDeleterious_daf0.9_Summary))
write.table(nonDeleterious, file="nonDeleteriousAllelesHeteroRatioSummary_dafPartition.txt", quote=F, sep="\t", row.names=T)



combine <- as.data.frame(rbind(deleteriousSummary, nonDeleteriousSummary))

ggplot(data=combine, aes(x=type, y=ratio, color=type)) +  
geom_boxplot(outlier.shape=NA) +  
geom_jitter(position=position_jitter(width=.1, height=0), aes(color=type)) + 
ylab("Ratio of heterozygous genotypes in hybrid F1s") +
theme(axis.text.x=element_text(angle=90, hjust=1)) 
#coord_cartesian(ylim = ylim1*5)
ggsave("RatioHeteroHybrids.pdf")
