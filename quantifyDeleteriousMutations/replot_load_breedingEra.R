require(xlsx)
library(data.table)


library(expss)


library(cowplot)
library(dplyr)

df <- read.xlsx("inbred_id_New.xlsx", sheetName="Sheet1", header=T)
head(df)
names(df) <- c("no", "ID", "Inbred", "Era")
df2 <- df[2:length(df$ID), ]
head(df2)
head(df)
era <- df2


SNP <- fread("allChrSNP.GERP.recode2.txt")
names(SNP) <- c(names(SNP)[2:length(names(SNP))], "none")
length(names(SNP))
SNP <- SNP[, 1:353]
names(SNP)[1] <- "ID"
nonSyn <- read.delim("allSNPs.SNPeff.nonSyn.txt", header=F)
head(nonSyn)
names(nonSyn)[1] <- "ID"
DelSNP <- merge(nonSyn, SNP, by="ID")

DelSNP02 <- subset(DelSNP, DelSNP$GERP > 0 & DelSNP$GERP <= 2)



allDelAllele <- apply(DelSNP02[, 5:354], 2, sum, na.rm=TRUE)
homDelAllele <- sum_col_if(2, DelSNP02[, 5:354])
heteroDelAllele <- sum_col_if(1, DelSNP02[, 5:354])




df <- as.data.frame(cbind(allDelAllele, homDelAllele, heteroDelAllele))
df["ID"] <- row.names(df)
combine02 <- merge(df, era, by="ID")
combine02["GERP"] <- "GERP02"

DelSNP24 <- subset(DelSNP, DelSNP$GERP > 2 & DelSNP$GERP <= 4)
length(DelSNP24$Allele)
allDelAllele <- apply(DelSNP24[, 5:354], 2, sum, na.rm=TRUE)
homDelAllele <- sum_col_if(2, DelSNP24[, 5:354])
heteroDelAllele <- sum_col_if(1, DelSNP24[, 5:354])
df <- as.data.frame(cbind(allDelAllele, homDelAllele, heteroDelAllele))
df["ID"] <- row.names(df)
combine24 <- merge(df, era, by="ID")
combine24["GERP"] <- "GERP24"

DelSNP4 <- subset(DelSNP, DelSNP$GERP > 4)
length(DelSNP4$Allele)
allDelAllele <- apply(DelSNP4[, 5:354], 2, sum, na.rm=TRUE)
homDelAllele <- sum_col_if(2, DelSNP4[, 5:354])
heteroDelAllele <- sum_col_if(1, DelSNP4[, 5:354])
df <- as.data.frame(cbind(allDelAllele, homDelAllele, heteroDelAllele))
df["ID"] <- row.names(df)
combine4 <- merge(df, era, by="ID")
combine4["GERP"] <- "GERP4"

combine <- as.data.frame(rbind(combine02, combine24, combine4))
combine["Era2"] <- ifelse(combine$Era=="CN_1960s", "CN_1960&70s", combine$Era)
combine$Era2 <- ifelse(combine$Era=="CN_1970s", "CN_1960&70s", combine$Era2)
combine$Era2 <- ifelse(combine$Era=="CN_1980s", "CN_1980&90s", combine$Era2)
combine$Era2 <- ifelse(combine$Era=="CN_1990s", "CN_1980&90s", combine$Era2)
combine$Era2 <- ifelse(combine$Era=="CN_2000s", "CN_2000&10s", combine$Era2)
combine$Era2 <- ifelse(combine$Era=="CN_2010s", "CN_2000&10s", combine$Era2)
combine$Era2 <- ifelse(combine$Era=="Public_US", "Public_US", combine$Era2)
combine$Era2 <- ifelse(combine$Era=="Ex-PVP", "Ex-PVP", combine$Era2)

write.table(combine, file="Era_DelAllele_nonsynonymous_MetaData.txt", quote=F, sep="\t", row.names=F)



ggplot(data=subset(combine, !is.na(Era2)), aes(x=Era2, y=allDelAllele, color=Era2)) +  
geom_boxplot(outlier.shape=NA) +  
geom_jitter(position=position_jitter(width=.1, height=0), aes(color=Era2)) + 
xlab("") + ylab("count") +
theme(axis.text.x=element_text(angle=90, hjust=1)) +
ggtitle("all deleterious and non-synonymous alleles") +
facet_grid(GERP ~ ., scales="free")+
scale_x_discrete(limits=c("CN_1960&70s","CN_1980&90s","CN_2000&10s", "Public_US", "Ex-PVP"))+
theme(legend.position="none")
ggsave("Era_AllDelAllele_nonsynonymous.pdf")

ggplot(data=subset(combine, !is.na(Era2)), aes(x=Era2, y=homDelAllele, color=Era2)) +  
geom_boxplot(outlier.shape=NA) +  
geom_jitter(position=position_jitter(width=.1, height=0), aes(color=Era2)) + 
xlab("") + ylab("count") +
theme(axis.text.x=element_text(angle=90, hjust=1)) +
ggtitle("homozygous deleterious and non-synonymous alleles") +
facet_grid(GERP ~ ., scales="free")+
scale_x_discrete(limits=c("CN_1960&70s","CN_1980&90s","CN_2000&10s", "Public_US", "Ex-PVP"))+
theme(legend.position="none")
ggsave("Era_HomoDelAllele_nonsynonymous.pdf")

# compute lower and upper whiskers
ylim1 = boxplot.stats(combine$heteroDelAllele)$stats[c(1, 5)]
 
ggplot(data=subset(combine, !is.na(Era2)), aes(x=Era2, y=heteroDelAllele, color=Era2)) +  
geom_boxplot(outlier.shape=NA) +  
geom_jitter(position=position_jitter(width=.1, height=0), aes(color=Era2)) + 
xlab("") + ylab("count") +
theme(axis.text.x=element_text(angle=90, hjust=1)) +
ggtitle("heterozygous deleterious and non-synonymous alleles") +
facet_grid(GERP ~ ., scales="free")+
scale_x_discrete(limits=c("CN_1960&70s","CN_1980&90s","CN_2000&10s", "Public_US", "Ex-PVP"))+
theme(legend.position="none")
ggsave("Era_HeteroDelAllele_nonsynonymous.pdf")

CN1 <- subset(combine, combine$Era2=="CN_1960&70s")
write.table(as.data.frame(CN1$ID), file="CN1_samples.txt", quote=F, sep="\t", row.names=F)

CN2 <- subset(combine, combine$Era2=="CN_1980&90s")
write.table(as.data.frame(CN2$ID), file="CN2_samples.txt", quote=F, sep="\t", row.names=F)

CN3 <- subset(combine, combine$Era2=="CN_2000&10s")
write.table(as.data.frame(CN3$ID), file="CN3_samples.txt", quote=F, sep="\t", row.names=F)

US1 <- subset(combine, combine$Era2=="Public_US")
write.table(as.data.frame(US1$ID), file="US1_samples.txt", quote=F, sep="\t", row.names=F)

US2 <- subset(combine, combine$Era2=="Ex-PVP")
write.table(as.data.frame(US2$ID), file="US2_samples.txt", quote=F, sep="\t", row.names=F)
