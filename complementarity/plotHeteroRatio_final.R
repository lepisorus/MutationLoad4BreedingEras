library(ggplot2)

theme_li <- function () { 
  theme_bw(base_size=15) %+replace% 
    theme(
      axis.text.y=element_text(size=15),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(size=15),
      legend.position=c(0.1, 0.9)
    )
}


deleteriousSummary02 <- read.delim("deleteriousAllelesHeteroRatioSummary_GERP02.txt")
nonDeleteriousSummary <- read.delim("nonDeleteriousAllelesHeteroRatioSummary.txt")

deleteriousSummary02["type"] <- "GERP02"
nonDeleteriousSummary["type"] <- "GERP<=0"

deleteriousSummary24 <- read.delim("deleteriousAllelesHeteroRatioSummary_GERP24.txt")
deleteriousSummary24["type"] <- "GERP24"

deleteriousSummary4 <- read.delim("deleteriousAllelesHeteroRatioSummary_GERP4.txt")
deleteriousSummary4["type"] <- "GERP>4"

##paired t test and new point plot
combine.new <- as.data.frame(cbind(deleteriousSummary02[, 20:21], deleteriousSummary24[, 20:21], deleteriousSummary4[, 20:21], nonDeleteriousSummary[, 20:21]))
names(combine.new)[c(1,3,5,7)] <- c("ratio.GERP02", "ratio.GERP24","ratio.GERP4", "ratio.GERPless0")
t.test(combine.new$ratio.GERP02, combine.new$ratio.GERPless0, alternative="two.sided", mu=0, paired=TRUE)
t.test(combine.new$ratio.GERP24, combine.new$ratio.GERPless0, alternative="two.sided", mu=0, paired=TRUE)
t.test(combine.new$ratio.GERP4, combine.new$ratio.GERPless0, alternative="two.sided", mu=0, paired=TRUE)

names(combine.new)[c(1,3,5)] <- c("ratio", "ratio","ratio")
df <- as.data.frame(rbind(combine.new[, c(1,2,7)], combine.new[, c(3,4,7)], combine.new[, c(5,6,7)]))
df$type <- ifelse(df$type=="GERP>4", "GERP4", df$type)

ggplot(data=df, aes(x=ratio.GERPless0, y=ratio, color=type))+
geom_point(size=2)+
geom_abline(intercept=0, slope=1, linetype="dashed")+
xlab("heterozygous genotypes in non-deleterious genomic regions%")+
ylab ("heterozygous genotypes in deleterious genomic regions%")+
facet_grid(type ~ .)+
theme(legend.position="none")
ggsave("HeteroRatioHybrids_pointPlot.pdf")

ggplot(data=df, aes(x=ratio.GERPless0, y=ratio, color=type, group=type))+
geom_point(size=2)+
geom_abline(intercept=0, slope=1, linetype="dashed")+
xlab("heterozygous genotypes in non-deleterious genomic regions%")+
ylab ("heterozygous genotypes in deleterious genomic regions%")+
geom_smooth(method='lm', formula= y~x)+
theme_li()
ggsave("HeteroRatioHybrids_pointPlot2.pdf")



##partition DAF

theme_li <- function () { 
  theme_bw(base_size=15) %+replace% 
    theme(
      axis.text.y=element_text(size=15),
      axis.text.x=element_text(size=6),
      panel.grid.minor = element_blank(),
      legend.position="none"
    )
}


deleteriousSummary02 <- read.delim("GERP02_HeteroRatioSummary_dafPartition.txt")
deleteriousSummary02["type"] <- "GERP02"
deleteriousSummary02 <- deleteriousSummary02[, c(18,21,22,23)]
deleteriousSummary02["ID"] <- paste(deleteriousSummary02$sample, deleteriousSummary02$daf, sep="_")

deleteriousSummary24 <- read.delim("GERP24_HeteroRatioSummary_dafPartition.txt")
deleteriousSummary24["type"] <- "GERP24"
deleteriousSummary24 <- deleteriousSummary24[, c(18,21,22,23)]
deleteriousSummary24["ID"] <- paste(deleteriousSummary24$sample, deleteriousSummary24$daf, sep="_")

deleteriousSummary4 <- read.delim("GERP4_HeteroRatioSummary_dafPartition.txt")
deleteriousSummary4["type"] <- "GERP4"
deleteriousSummary4 <- deleteriousSummary4[, c(18,21,22,23)]
deleteriousSummary4["ID"] <- paste(deleteriousSummary4$sample, deleteriousSummary4$daf, sep="_")

nonDeleteriousSummary <- read.delim("nonDeleteriousAllelesHeteroRatioSummary_dafPartition.txt")
nonDeleteriousSummary["type"] <- "nonDeleterious"
nonDeleteriousSummary <- nonDeleteriousSummary[, c(18,21,22,23)]
nonDeleteriousSummary["ID"] <- paste(nonDeleteriousSummary$sample, nonDeleteriousSummary$daf, sep="_")

##paired t test and new point plot
combine.new <- merge(deleteriousSummary02, nonDeleteriousSummary, by="ID")
names(combine.new)[c(3,7)] <- c("ratio.GERP02", "ratio.GERPless0")
combine.new2 <- combine.new[, c(3,4,7)]
names(combine.new2)[2] <- "daf"

ggplot(data=combine.new2, aes(x=ratio.GERPless0, y=ratio.GERP02))+
geom_point(size=2)+
geom_abline(intercept=0, slope=1, linetype="dashed")+
xlab("heterozygous genotypes in non-deleterious genomic regions%")+
ylab ("heterozygous genotypes in deleterious genomic regions%")+
facet_grid(daf ~ .)+
scale_y_continuous(breaks=c(0.2, 0.4))+
theme(legend.position="none") + theme_li()
ggsave("HeteroRatioHybrids_dafPartition.pdf")

## we only found significant difference when daf is between 0.6 - 0.8
daf0 <- combine.new2 %>% filter(daf==0)
wilcox.test(daf0$ratio.GERP02, daf0$ratio.GERPless0, alternative="greater")

daf0.1 <- combine.new2 %>% filter(daf==0.1)
wilcox.test(daf0.1$ratio.GERP02, daf0.1$ratio.GERPless0, alternative="greater")

daf0.2 <- combine.new2 %>% filter(daf==0.2)
wilcox.test(daf0.2$ratio.GERP02, daf0.2$ratio.GERPless0, alternative="greater")

daf0.3 <- combine.new2 %>% filter(daf==0.3)
wilcox.test(daf0.3$ratio.GERP02, daf0.3$ratio.GERPless0, alternative="greater")

daf0.4 <- combine.new2 %>% filter(daf==0.4)
wilcox.test(daf0.4$ratio.GERP02, daf0.4$ratio.GERPless0, alternative="greater")

daf0.5 <- combine.new2 %>% filter(daf==0.5)
wilcox.test(daf0.5$ratio.GERP02, daf0.5$ratio.GERPless0, alternative="greater")

daf0.6 <- combine.new2 %>% filter(daf==0.6)
wilcox.test(daf0.6$ratio.GERP02, daf0.6$ratio.GERPless0, alternative="greater")

daf0.7 <- combine.new2 %>% filter(daf==0.7)
wilcox.test(daf0.7$ratio.GERP02, daf0.7$ratio.GERPless0, alternative="greater")

daf0.8 <- combine.new2 %>% filter(daf==0.8)
wilcox.test(daf0.8$ratio.GERP02, daf0.8$ratio.GERPless0, alternative="greater")

daf0.9 <- combine.new2 %>% filter(daf==0.9)
wilcox.test(daf0.9$ratio.GERP02, daf0.9$ratio.GERPless0, alternative="greater")


#boxplot
combine.new <- as.data.frame(rbind(deleteriousSummary02, nonDeleteriousSummary))
ggplot(data=combine.new, aes(x=type, y=ratio, color=type)) +  
geom_boxplot(outlier.shape=NA) +  
geom_jitter(position=position_jitter(width=.1, height=0), aes(color=type)) + 
ylab("Ratio of heterozygous genotypes in hybrid F1s") +
theme(axis.text.x=element_blank()) +
facet_wrap(~daf, nrow=2) + theme_li()
#coord_cartesian(ylim = ylim1*5)
ggsave("RatioHeteroHybrids_GERP02_dafPartition.pdf")

##GERP24  paired t test and new point plot
combine.new <- merge(deleteriousSummary24, nonDeleteriousSummary, by="ID")
names(combine.new)[c(3,7)] <- c("ratio.GERP24", "ratio.GERPless0")
combine.new2 <- combine.new[, c(3,4,7)]
names(combine.new2)[2] <- "daf"

ggplot(data=combine.new2, aes(x=ratio.GERPless0, y=ratio.GERP24))+
geom_point(size=2)+
geom_abline(intercept=0, slope=1, linetype="dashed")+
xlab("heterozygous genotypes in non-deleterious genomic regions%")+
ylab ("heterozygous genotypes in deleterious genomic regions%")+
facet_grid(daf ~ .)+
scale_y_continuous(breaks=c(0.2, 0.4))+
theme(legend.position="none")
ggsave("HeteroRatioHybrids_dafPartition_GERP24.pdf")

daf0 <- combine.new2 %>% filter(daf==0)
wilcox.test(daf0$ratio.GERP24, daf0$ratio.GERPless0, alternative="greater")

daf0.1 <- combine.new2 %>% filter(daf==0.1)
wilcox.test(daf0.1$ratio.GERP24, daf0.1$ratio.GERPless0, alternative="greater")

daf0.2 <- combine.new2 %>% filter(daf==0.2)
wilcox.test(daf0.2$ratio.GERP24, daf0.2$ratio.GERPless0, alternative="greater")

daf0.3 <- combine.new2 %>% filter(daf==0.3)
wilcox.test(daf0.3$ratio.GERP24, daf0.3$ratio.GERPless0, alternative="greater")

daf0.4 <- combine.new2 %>% filter(daf==0.4)
wilcox.test(daf0.4$ratio.GERP24, daf0.4$ratio.GERPless0, alternative="greater")

daf0.5 <- combine.new2 %>% filter(daf==0.5)
wilcox.test(daf0.5$ratio.GERP24, daf0.5$ratio.GERPless0, alternative="greater")

daf0.6 <- combine.new2 %>% filter(daf==0.6)
wilcox.test(daf0.6$ratio.GERP24, daf0.6$ratio.GERPless0, alternative="greater")

daf0.7 <- combine.new2 %>% filter(daf==0.7)
wilcox.test(daf0.7$ratio.GERP24, daf0.7$ratio.GERPless0, alternative="greater")

daf0.8 <- combine.new2 %>% filter(daf==0.8)
wilcox.test(daf0.8$ratio.GERP24, daf0.8$ratio.GERPless0, alternative="greater")

daf0.9 <- combine.new2 %>% filter(daf==0.9)
wilcox.test(daf0.9$ratio.GERP24, daf0.9$ratio.GERPless0, alternative="greater")

#boxplot
combine.new <- as.data.frame(rbind(deleteriousSummary24, nonDeleteriousSummary))
ggplot(data=combine.new, aes(x=type, y=ratio, color=type)) +  
geom_boxplot(outlier.shape=NA) +  
geom_jitter(position=position_jitter(width=.1, height=0), aes(color=type)) + 
ylab("Ratio of heterozygous genotypes in hybrid F1s") +
#theme(axis.text.x=element_blank()) +
facet_wrap(~daf, nrow=2)+theme_li()
#coord_cartesian(ylim = ylim1*5)
ggsave("RatioHeteroHybrids_GERP24_dafPartition.pdf")

##GERP4
combine.new <- merge(deleteriousSummary4, nonDeleteriousSummary, by="ID")
names(combine.new)[c(3,7)] <- c("ratio.GERP4", "ratio.GERPless0")
combine.new2 <- combine.new[, c(3,4,7)]
names(combine.new2)[2] <- "daf"

ggplot(data=combine.new2, aes(x=ratio.GERPless0, y=ratio.GERP4))+
geom_point(size=2)+
geom_abline(intercept=0, slope=1, linetype="dashed")+
xlab("heterozygous genotypes in non-deleterious genomic regions%")+
ylab ("heterozygous genotypes in deleterious genomic regions%")+
facet_grid(daf ~ .)+
scale_y_continuous(breaks=c(0.2, 0.4))+
theme(legend.position="none")
ggsave("HeteroRatioHybrids_dafPartition_GERP4.pdf")

daf0 <- combine.new2 %>% filter(daf==0)
wilcox.test(daf0$ratio.GERP4, daf0$ratio.GERPless0, alternative="greater")

daf0.1 <- combine.new2 %>% filter(daf==0.1)
wilcox.test(daf0.1$ratio.GERP4, daf0.1$ratio.GERPless0, alternative="greater")

daf0.2 <- combine.new2 %>% filter(daf==0.2)
wilcox.test(daf0.2$ratio.GERP4, daf0.2$ratio.GERPless0, alternative="greater")

daf0.3 <- combine.new2 %>% filter(daf==0.3)
wilcox.test(daf0.3$ratio.GERP4, daf0.3$ratio.GERPless0, alternative="greater")

daf0.4 <- combine.new2 %>% filter(daf==0.4)
wilcox.test(daf0.4$ratio.GERP4, daf0.4$ratio.GERPless0, alternative="greater")

daf0.5 <- combine.new2 %>% filter(daf==0.5)
wilcox.test(daf0.5$ratio.GERP4, daf0.5$ratio.GERPless0, alternative="greater")

daf0.6 <- combine.new2 %>% filter(daf==0.6)
wilcox.test(daf0.6$ratio.GERP4, daf0.6$ratio.GERPless0, alternative="greater")

daf0.7 <- combine.new2 %>% filter(daf==0.7)
wilcox.test(daf0.7$ratio.GERP4, daf0.7$ratio.GERPless0, alternative="greater")

daf0.8 <- combine.new2 %>% filter(daf==0.8)
wilcox.test(daf0.8$ratio.GERP4, daf0.8$ratio.GERPless0, alternative="greater")

daf0.9 <- combine.new2 %>% filter(daf==0.9)
wilcox.test(daf0.9$ratio.GERP4, daf0.9$ratio.GERPless0, alternative="greater")

#boxplot
combine.new <- as.data.frame(rbind(deleteriousSummary4, nonDeleteriousSummary))
ggplot(data=combine.new, aes(x=type, y=ratio, color=type)) +  
geom_boxplot(outlier.shape=NA) +  
geom_jitter(position=position_jitter(width=.1, height=0), aes(color=type)) + 
ylab("Ratio of heterozygous genotypes in hybrid F1s") +
#theme(axis.text.x=element_blank()) +
facet_wrap(~daf, nrow=2) + theme_li()
#coord_cartesian(ylim = ylim1*5)
ggsave("RatioHeteroHybrids_GERP4_dafPartition.pdf")

