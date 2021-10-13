library(data.table)
library(cowplot)


#GERPless0
CN1 <- fread("CN1.GERPless0.daf")
head(CN1)
CN1 <- as.data.frame(CN1)
CN1["group"] <- "CN1"

CN2 <- fread("CN2.GERPless0.daf")
CN2 <- as.data.frame(CN2)
CN2["group"] <- "CN2"
CN3 <- fread("CN3.GERPless0.daf")
CN3 <- as.data.frame(CN3)
CN3["group"] <- "CN3"

US1 <- fread("US1.GERPless0.daf")
US1 <- as.data.frame(US1)
US1["group"] <- "US1"
US2 <- fread("US2.GERPless0.daf")
US2 <- as.data.frame(US2)
US2["group"] <- "US2"

CN <- as.data.frame(rbind(CN1, CN2, CN3))
US <- as.data.frame(rbind(US1, US2))
CN["continent"] <- "CN"
US["continent"] <- "US"
CN["daf"] <- CN$Number/60
US["daf"] <- US$Number/148

combine <- as.data.frame(rbind(CN, US))
CN <- CN %>% filter(daf>0 & daf<1)

#cut into 20 maf ranges
breaks <- seq(0, 1, 0.05)
CN["bin"] <- cut(CN$daf, 20, include.lowest=TRUE, labels=c("0", "0.05", "0.1", "0.15", "0.2", "0.25", "0.3", "0.35", "0.4", "0.45", "0.5", "0.55", "0.6", "0.65", "0.7", "0.75", "0.8", "0.85", "0.9", "0.95"))

CN_2 <- CN %>% group_by(group, bin) %>% summarize(count=sum(Observed)) 
CN_2 <- CN_2 %>% group_by(group) %>%  mutate(percent=count/sum())
table(CN$group)
CN_2 <- as.data.frame(CN_2)
CN_2["percent"] <- ifelse(CN_2$group=="CN1", CN_2$count/sum(CN1$Observed), 0)
CN_2$percent <- ifelse(CN_2$group=="CN2", CN_2$count/sum(CN2$Observed), CN_2$percent)
CN_2$percent <- ifelse(CN_2$group=="CN3", CN_2$count/sum(CN3$Observed), CN_2$percent)
head(CN_2)
CN2["GERP"] <- "GERPless0"
CN2_GERPless0 <- CN2

#US
US <- US %>% filter(daf>0 & daf<1)
US["bin"] <- cut(US$daf, 20, include.lowest=TRUE, labels=c("0", "0.05", "0.1", "0.15", "0.2", "0.25", "0.3", "0.35", "0.4", "0.45", "0.5", "0.55", "0.6", "0.65", "0.7", "0.75", "0.8", "0.85", "0.9", "0.95"))

US_2 <- US %>% group_by(group, bin) %>% summarize(count=sum(Observed)) 
US_2 <- US_2 %>% group_by(group) %>%  mutate(percent=count/sum())
table(US$group)
US_2 <- as.data.frame(US_2)
US_2["percent"] <- ifelse(US_2$group=="US1", US_2$count/sum(US1$Observed), 0)
US_2$percent <- ifelse(US_2$group=="US2", US_2$count/sum(US2$Observed), US_2$percent)
US_2$percent <- ifelse(US_2$group=="US3", US_2$count/sum(US3$Observed), US_2$percent)
head(US_2)
US2["GERP"] <- "GERPless0"
US2_GERPless0 <- US2



#GERP02
CN1 <- fread("CN1.GERP02.daf")
CN1 <- as.data.frame(CN1)
CN1["group"] <- "CN1"

CN2 <- fread("CN2.GERP02.daf")
CN2 <- as.data.frame(CN2)
CN2["group"] <- "CN2"
CN3 <- fread("CN3.GERP02.daf")
CN3 <- as.data.frame(CN3)
CN3["group"] <- "CN3"

US1 <- fread("US1.GERP02.daf")
US1 <- as.data.frame(US1)
US1["group"] <- "US1"
US2 <- fread("US2.GERP02.daf")
US2 <- as.data.frame(US2)
US2["group"] <- "US2"

CN <- as.data.frame(rbind(CN1, CN2, CN3))
US <- as.data.frame(rbind(US1, US2))
CN["continent"] <- "CN"
US["continent"] <- "US"
CN["daf"] <- CN$Number/60
US["daf"] <- US$Number/148

combine <- as.data.frame(rbind(CN, US))
CN <- CN %>% filter(daf>0 & daf<1)

#cut into 20 maf ranges
breaks <- seq(0, 1, 0.05)
CN["bin"] <- cut(CN$daf, 20, include.lowest=TRUE, labels=c("0", "0.05", "0.1", "0.15", "0.2", "0.25", "0.3", "0.35", "0.4", "0.45", "0.5", "0.55", "0.6", "0.65", "0.7", "0.75", "0.8", "0.85", "0.9", "0.95"))

CN_2 <- CN %>% group_by(group, bin) %>% summarize(count=sum(Observed)) 
CN_2 <- CN_2 %>% group_by(group) %>%  mutate(percent=count/sum())
table(CN$group)
CN_2 <- as.data.frame(CN_2)
CN_2["percent"] <- ifelse(CN_2$group=="CN1", CN_2$count/sum(CN1$Observed), 0)
CN_2$percent <- ifelse(CN_2$group=="CN2", CN_2$count/sum(CN2$Observed), CN_2$percent)
CN_2$percent <- ifelse(CN_2$group=="CN3", CN_2$count/sum(CN3$Observed), CN_2$percent)
head(CN_2)
CN2["GERP"] <- "GERP02"
CN2_GERP02 <- CN2

#US
US <- US %>% filter(daf>0 & daf<1)
US["bin"] <- cut(US$daf, 20, include.lowest=TRUE, labels=c("0", "0.05", "0.1", "0.15", "0.2", "0.25", "0.3", "0.35", "0.4", "0.45", "0.5", "0.55", "0.6", "0.65", "0.7", "0.75", "0.8", "0.85", "0.9", "0.95"))

US_2 <- US %>% group_by(group, bin) %>% summarize(count=sum(Observed)) 
US_2 <- US_2 %>% group_by(group) %>%  mutate(percent=count/sum())
table(US$group)
US_2 <- as.data.frame(US_2)
US_2["percent"] <- ifelse(US_2$group=="US1", US_2$count/sum(US1$Observed), 0)
US_2$percent <- ifelse(US_2$group=="US2", US_2$count/sum(US2$Observed), US_2$percent)
US_2$percent <- ifelse(US_2$group=="US3", US_2$count/sum(US3$Observed), US_2$percent)
head(US_2)
US2["GERP"] <- "GERP02"
US2_GERP02 <- US2



#GERP24
CN1 <- fread("CN1.GERP24.daf")
CN1 <- as.data.frame(CN1)
CN1["group"] <- "CN1"

CN2 <- fread("CN2.GERP24.daf")
CN2 <- as.data.frame(CN2)
CN2["group"] <- "CN2"
CN3 <- fread("CN3.GERP24.daf")
CN3 <- as.data.frame(CN3)
CN3["group"] <- "CN3"

US1 <- fread("US1.GERP24.daf")
US1 <- as.data.frame(US1)
US1["group"] <- "US1"
US2 <- fread("US2.GERP24.daf")
US2 <- as.data.frame(US2)
US2["group"] <- "US2"

CN <- as.data.frame(rbind(CN1, CN2, CN3))
US <- as.data.frame(rbind(US1, US2))
CN["continent"] <- "CN"
US["continent"] <- "US"
CN["daf"] <- CN$Number/60
US["daf"] <- US$Number/148

combine <- as.data.frame(rbind(CN, US))
CN <- CN %>% filter(daf>0 & daf<1)

#cut into 20 maf ranges
breaks <- seq(0, 1, 0.05)
CN["bin"] <- cut(CN$daf, 20, include.lowest=TRUE, labels=c("0", "0.05", "0.1", "0.15", "0.2", "0.25", "0.3", "0.35", "0.4", "0.45", "0.5", "0.55", "0.6", "0.65", "0.7", "0.75", "0.8", "0.85", "0.9", "0.95"))

CN_2 <- CN %>% group_by(group, bin) %>% summarize(count=sum(Observed)) 
CN_2 <- CN_2 %>% group_by(group) %>%  mutate(percent=count/sum())
table(CN$group)
CN_2 <- as.data.frame(CN_2)
CN_2["percent"] <- ifelse(CN_2$group=="CN1", CN_2$count/sum(CN1$Observed), 0)
CN_2$percent <- ifelse(CN_2$group=="CN2", CN_2$count/sum(CN2$Observed), CN_2$percent)
CN_2$percent <- ifelse(CN_2$group=="CN3", CN_2$count/sum(CN3$Observed), CN_2$percent)
head(CN_2)
CN2["GERP"] <- "GERP24"
CN2_GERP24 <- CN2

#US
US <- US %>% filter(daf>0 & daf<1)
US["bin"] <- cut(US$daf, 20, include.lowest=TRUE, labels=c("0", "0.05", "0.1", "0.15", "0.2", "0.25", "0.3", "0.35", "0.4", "0.45", "0.5", "0.55", "0.6", "0.65", "0.7", "0.75", "0.8", "0.85", "0.9", "0.95"))

US_2 <- US %>% group_by(group, bin) %>% summarize(count=sum(Observed)) 
US_2 <- US_2 %>% group_by(group) %>%  mutate(percent=count/sum())
table(US$group)
US_2 <- as.data.frame(US_2)
US_2["percent"] <- ifelse(US_2$group=="US1", US_2$count/sum(US1$Observed), 0)
US_2$percent <- ifelse(US_2$group=="US2", US_2$count/sum(US2$Observed), US_2$percent)
US_2$percent <- ifelse(US_2$group=="US3", US_2$count/sum(US3$Observed), US_2$percent)
head(US_2)
US2["GERP"] <- "GERP24"
US2_GERP24 <- US2




#GERP4
CN1 <- fread("CN1.GERP4.daf")
CN1 <- as.data.frame(CN1)
CN1["group"] <- "CN1"

CN2 <- fread("CN2.GERP4.daf")
CN2 <- as.data.frame(CN2)
CN2["group"] <- "CN2"
CN3 <- fread("CN3.GERP4.daf")
CN3 <- as.data.frame(CN3)
CN3["group"] <- "CN3"

US1 <- fread("US1.GERP4.daf")
US1 <- as.data.frame(US1)
US1["group"] <- "US1"
US2 <- fread("US2.GERP4.daf")
US2 <- as.data.frame(US2)
US2["group"] <- "US2"

CN <- as.data.frame(rbind(CN1, CN2, CN3))
US <- as.data.frame(rbind(US1, US2))
CN["continent"] <- "CN"
US["continent"] <- "US"
CN["daf"] <- CN$Number/60
US["daf"] <- US$Number/148

combine <- as.data.frame(rbind(CN, US))
CN <- CN %>% filter(daf>0 & daf<1)

#cut into 20 maf ranges
breaks <- seq(0, 1, 0.05)
CN["bin"] <- cut(CN$daf, 20, include.lowest=TRUE, labels=c("0", "0.05", "0.1", "0.15", "0.2", "0.25", "0.3", "0.35", "0.4", "0.45", "0.5", "0.55", "0.6", "0.65", "0.7", "0.75", "0.8", "0.85", "0.9", "0.95"))

CN_2 <- CN %>% group_by(group, bin) %>% summarize(count=sum(Observed)) 
CN_2 <- CN_2 %>% group_by(group) %>%  mutate(percent=count/sum())
table(CN$group)
CN_2 <- as.data.frame(CN_2)
CN_2["percent"] <- ifelse(CN_2$group=="CN1", CN_2$count/sum(CN1$Observed), 0)
CN_2$percent <- ifelse(CN_2$group=="CN2", CN_2$count/sum(CN2$Observed), CN_2$percent)
CN_2$percent <- ifelse(CN_2$group=="CN3", CN_2$count/sum(CN3$Observed), CN_2$percent)
head(CN_2)
CN2["GERP"] <- "GERP4"
CN2_GERP4 <- CN2

#US
US <- US %>% filter(daf>0 & daf<1)
US["bin"] <- cut(US$daf, 20, include.lowest=TRUE, labels=c("0", "0.05", "0.1", "0.15", "0.2", "0.25", "0.3", "0.35", "0.4", "0.45", "0.5", "0.55", "0.6", "0.65", "0.7", "0.75", "0.8", "0.85", "0.9", "0.95"))

US_2 <- US %>% group_by(group, bin) %>% summarize(count=sum(Observed)) 
US_2 <- US_2 %>% group_by(group) %>%  mutate(percent=count/sum())
table(US$group)
US_2 <- as.data.frame(US_2)
US_2["percent"] <- ifelse(US_2$group=="US1", US_2$count/sum(US1$Observed), 0)
US_2$percent <- ifelse(US_2$group=="US2", US_2$count/sum(US2$Observed), US_2$percent)
US_2$percent <- ifelse(US_2$group=="US3", US_2$count/sum(US3$Observed), US_2$percent)
head(US_2)
US2["GERP"] <- "GERP4"
US2_GERP4 <- US2


CN <- as.data.frame(rbind(CN2_GERPless0, CN2_GERP02, CN2_GERP24, CN2_GERP4))
US <- as.data.frame(rbind(US2_GERPless0, US2_GERP02, US2_GERP24, US2_GERP4))

##plot
ggplot(data=CN, aes(x=bin, y=percent, fill=group, group=group)) +   
   geom_bar(stat="identity", width=0.8, position=position_dodge(), alpha=0.8) + 
 #  scale_x_continuous(breaks=c(0, .1,.2,.3,.4,.5,.6,.7,.8, .9)) +
   ylab("Density") + xlab("derived allele frequency") + 
   facet_grid(. ~ GERP) +
   theme_bw() + theme(legend.position=c(0.85, 0.85), plot.title=element_text(size=20),
                      axis.title.y=element_text(size = 16, vjust=+0.2),
                      axis.title.x=element_text(size = 16, vjust=-0.2),
                      axis.text.y=element_text(size = 14),
                      axis.text.x=element_text(size = 14),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank())
ggsave("dafAcrossEra_bar_GERP_CN.pdf")                     


ggplot(data=US, aes(x=bin, y=percent, fill=group, group=group)) +   
   geom_bar(stat="identity", width=0.8, position=position_dodge(), alpha=0.8) + 
 #  scale_x_continuous(breaks=c(0, .1,.2,.3,.4,.5,.6,.7,.8, .9)) +
   ylab("Density") + xlab("derived allele frequency") + 
   facet_grid(. ~ GERP) +
   theme_bw() + theme(legend.position=c(0.85, 0.85), plot.title=element_text(size=20),
                      axis.title.y=element_text(size = 16, vjust=+0.2),
                      axis.title.x=element_text(size = 16, vjust=-0.2),
                      axis.text.y=element_text(size = 14),
                      axis.text.x=element_text(size = 14),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank())
ggsave("dafAcrossEra_bar_GERP_US.pdf")    


