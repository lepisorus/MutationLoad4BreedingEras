library(data.table)
library(cowplot)
library(ggplot2)
library(dplyr)

#GERPless0
CN1 <- fread("CN1.GERPless0.daf")
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

CN["daf"] <- CN$Number/60
US["daf"] <- US$Number/148


#breaks <- seq(0, 1, 0.05)
#CN["bin"] <- cut(CN$daf, 20, include.lowest=TRUE, labels=c("0", "0.05", "0.1", "0.15", "0.2", "0.25", "0.3", "0.35", "0.4", "0.45", "0.5", "0.55", "0.6", "0.65", "0.7", "0.75", "0.8", "0.85", "0.9", "0.95"))

#break daf into 10 bins, each bin is 0.1
breaks <- seq(0, 1, 0.1)
CN["bin"] <- cut(CN$daf, 10, include.lowest=TRUE, labels=c("0", "0.1",  "0.2", "0.3", "0.4", "0.5",  "0.6", "0.7",  "0.8",  "0.9"))


CN_2 <- CN %>% group_by(group, bin) %>% summarize(count=sum(Observed), .groups='drop') 
CN_2 <- CN_2 %>% group_by(group) %>%  mutate(percent=count/sum(count))
CN_3 <- CN_2 %>% group_by(group) %>%  mutate(Sumpercent=sum(percent))

table(CN$group)

head(CN_2,n=40)

ggplot(data=CN_2, aes(x=bin, y=percent, fill=group, group=group)) +   
   geom_bar(stat="identity", width=0.8, position=position_dodge(), alpha=0.8) + 
   scale_x_discrete(breaks=c(0, .1,.2,.3,.4,.5,.6,.7,.8, .9)) +
   ylab("Density") + xlab("derived allele frequency") + 
   theme_bw() + theme(legend.position=c(0.85, 0.85), plot.title=element_text(size=20),
                      axis.title.y=element_text(size = 16, vjust=+0.2),
                      axis.title.x=element_text(size = 16, vjust=-0.2),
                      axis.text.y=element_text(size = 14),
                      axis.text.x=element_text(size = 14),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank())
ggsave("dafAcrossEra_bar_GERPless0_CN.pdf")                     




#US
US <- US %>% filter(daf>0 & daf<1)
US["bin"] <- cut(US$daf, 20, include.lowest=TRUE, labels=c("0", "0.05", "0.1", "0.15", "0.2", "0.25", "0.3", "0.35", "0.4", "0.45", "0.5", "0.55", "0.6", "0.65", "0.7", "0.75", "0.8", "0.85", "0.9", "0.95"))

US_2 <- US %>% group_by(group, bin) %>% summarize(count=sum(Observed)) 
US_2 <- US_2 %>% group_by(group) %>%  mutate(percent=count/sum(count))
table(US$group)

head(US_2)

ggplot(data=US_2, aes(x=bin, y=percent, fill=group, group=group)) +   
   geom_bar(stat="identity", width=0.8, position=position_dodge(), alpha=0.8) + 
   scale_x_discrete(breaks=c(0, .1,.2,.3,.4,.5,.6,.7,.8, .9)) +
   ylab("Density") + xlab("derived allele frequency") + 
   theme_bw() + theme(legend.position=c(0.85, 0.85), plot.title=element_text(size=20),
                      axis.title.y=element_text(size = 16, vjust=+0.2),
                      axis.title.x=element_text(size = 16, vjust=-0.2),
                      axis.text.y=element_text(size = 14),
                      axis.text.x=element_text(size = 14),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank())
ggsave("dafAcrossEra_bar_GERPless0_US.pdf")    


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

CN["daf"] <- CN$Number/60
US["daf"] <- US$Number/148


breaks <- seq(0, 1, 0.05)
CN["bin"] <- cut(CN$daf, 20, include.lowest=TRUE, labels=c("0", "0.05", "0.1", "0.15", "0.2", "0.25", "0.3", "0.35", "0.4", "0.45", "0.5", "0.55", "0.6", "0.65", "0.7", "0.75", "0.8", "0.85", "0.9", "0.95"))


CN_2 <- CN %>% group_by(group, bin) %>% summarize(count=sum(Observed), .groups='drop') 
CN_2 <- CN_2 %>% group_by(group) %>%  mutate(percent=count/sum(count))
CN_3 <- CN_2 %>% group_by(group) %>%  mutate(Sumpercent=sum(percent))

table(CN$group)

head(CN_2,n=40)

ggplot(data=CN_2, aes(x=bin, y=percent, fill=group, group=group)) +   
   geom_bar(stat="identity", width=0.8, position=position_dodge(), alpha=0.8) + 
   scale_x_discrete(breaks=c(0, .1,.2,.3,.4,.5,.6,.7,.8, .9)) +
   ylab("Density") + xlab("derived allele frequency") + 
   theme_bw() + theme(legend.position=c(0.85, 0.85), plot.title=element_text(size=20),
                      axis.title.y=element_text(size = 16, vjust=+0.2),
                      axis.title.x=element_text(size = 16, vjust=-0.2),
                      axis.text.y=element_text(size = 14),
                      axis.text.x=element_text(size = 14),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank())
ggsave("dafAcrossEra_bar_GERP02_CN.pdf")                     




#US
US <- US %>% filter(daf>0 & daf<1)
US["bin"] <- cut(US$daf, 20, include.lowest=TRUE, labels=c("0", "0.05", "0.1", "0.15", "0.2", "0.25", "0.3", "0.35", "0.4", "0.45", "0.5", "0.55", "0.6", "0.65", "0.7", "0.75", "0.8", "0.85", "0.9", "0.95"))

US_2 <- US %>% group_by(group, bin) %>% summarize(count=sum(Observed)) 
US_2 <- US_2 %>% group_by(group) %>%  mutate(percent=count/sum(count))
table(US$group)

head(US_2)

ggplot(data=US_2, aes(x=bin, y=percent, fill=group, group=group)) +   
   geom_bar(stat="identity", width=0.8, position=position_dodge(), alpha=0.8) + 
   scale_x_discrete(breaks=c(0, .1,.2,.3,.4,.5,.6,.7,.8, .9)) +
   ylab("Density") + xlab("derived allele frequency") + 
   theme_bw() + theme(legend.position=c(0.85, 0.85), plot.title=element_text(size=20),
                      axis.title.y=element_text(size = 16, vjust=+0.2),
                      axis.title.x=element_text(size = 16, vjust=-0.2),
                      axis.text.y=element_text(size = 14),
                      axis.text.x=element_text(size = 14),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank())
ggsave("dafAcrossEra_bar_GERP02_US.pdf")    

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

CN["daf"] <- CN$Number/60
US["daf"] <- US$Number/148


breaks <- seq(0, 1, 0.05)
CN["bin"] <- cut(CN$daf, 20, include.lowest=TRUE, labels=c("0", "0.05", "0.1", "0.15", "0.2", "0.25", "0.3", "0.35", "0.4", "0.45", "0.5", "0.55", "0.6", "0.65", "0.7", "0.75", "0.8", "0.85", "0.9", "0.95"))


CN_2 <- CN %>% group_by(group, bin) %>% summarize(count=sum(Observed), .groups='drop') 
CN_2 <- CN_2 %>% group_by(group) %>%  mutate(percent=count/sum(count))
CN_3 <- CN_2 %>% group_by(group) %>%  mutate(Sumpercent=sum(percent))

table(CN$group)

head(CN_2,n=40)

ggplot(data=CN_2, aes(x=bin, y=percent, fill=group, group=group)) +   
   geom_bar(stat="identity", width=0.8, position=position_dodge(), alpha=0.8) + 
   scale_x_discrete(breaks=c(0, .1,.2,.3,.4,.5,.6,.7,.8, .9)) +
   ylab("Density") + xlab("derived allele frequency") + 
   theme_bw() + theme(legend.position=c(0.85, 0.85), plot.title=element_text(size=20),
                      axis.title.y=element_text(size = 16, vjust=+0.2),
                      axis.title.x=element_text(size = 16, vjust=-0.2),
                      axis.text.y=element_text(size = 14),
                      axis.text.x=element_text(size = 14),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank())
ggsave("dafAcrossEra_bar_GERP24_CN.pdf")                     




#US
US <- US %>% filter(daf>0 & daf<1)
US["bin"] <- cut(US$daf, 20, include.lowest=TRUE, labels=c("0", "0.05", "0.1", "0.15", "0.2", "0.25", "0.3", "0.35", "0.4", "0.45", "0.5", "0.55", "0.6", "0.65", "0.7", "0.75", "0.8", "0.85", "0.9", "0.95"))

US_2 <- US %>% group_by(group, bin) %>% summarize(count=sum(Observed)) 
US_2 <- US_2 %>% group_by(group) %>%  mutate(percent=count/sum(count))
table(US$group)

head(US_2)

ggplot(data=US_2, aes(x=bin, y=percent, fill=group, group=group)) +   
   geom_bar(stat="identity", width=0.8, position=position_dodge(), alpha=0.8) + 
   scale_x_discrete(breaks=c(0, .1,.2,.3,.4,.5,.6,.7,.8, .9)) +
   ylab("Density") + xlab("derived allele frequency") + 
   theme_bw() + theme(legend.position=c(0.85, 0.85), plot.title=element_text(size=20),
                      axis.title.y=element_text(size = 16, vjust=+0.2),
                      axis.title.x=element_text(size = 16, vjust=-0.2),
                      axis.text.y=element_text(size = 14),
                      axis.text.x=element_text(size = 14),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank())
ggsave("dafAcrossEra_bar_GERP24_US.pdf")    

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

CN["daf"] <- CN$Number/60
US["daf"] <- US$Number/148


breaks <- seq(0, 1, 0.05)
CN["bin"] <- cut(CN$daf, 20, include.lowest=TRUE, labels=c("0", "0.05", "0.1", "0.15", "0.2", "0.25", "0.3", "0.35", "0.4", "0.45", "0.5", "0.55", "0.6", "0.65", "0.7", "0.75", "0.8", "0.85", "0.9", "0.95"))

CN_2 <- CN %>% group_by(group, bin) %>% summarize(count=sum(Observed), .groups='drop') 
CN_2 <- CN_2 %>% group_by(group) %>%  mutate(percent=count/sum(count))
CN_3 <- CN_2 %>% group_by(group) %>%  mutate(Sumpercent=sum(percent))

table(CN$group)

head(CN_2,n=40)

ggplot(data=CN_2, aes(x=bin, y=percent, fill=group, group=group)) +   
   geom_bar(stat="identity", width=0.8, position=position_dodge(), alpha=0.8) + 
   scale_x_discrete(breaks=c(0, .1,.2,.3,.4,.5,.6,.7,.8, .9)) +
   ylab("Density") + xlab("derived allele frequency") + 
   theme_bw() + theme(legend.position=c(0.85, 0.85), plot.title=element_text(size=20),
                      axis.title.y=element_text(size = 16, vjust=+0.2),
                      axis.title.x=element_text(size = 16, vjust=-0.2),
                      axis.text.y=element_text(size = 14),
                      axis.text.x=element_text(size = 14),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank())
ggsave("dafAcrossEra_bar_GERP4_CN.pdf")                     




#US
US <- US %>% filter(daf>0 & daf<1)
US["bin"] <- cut(US$daf, 20, include.lowest=TRUE, labels=c("0", "0.05", "0.1", "0.15", "0.2", "0.25", "0.3", "0.35", "0.4", "0.45", "0.5", "0.55", "0.6", "0.65", "0.7", "0.75", "0.8", "0.85", "0.9", "0.95"))

US_2 <- US %>% group_by(group, bin) %>% summarize(count=sum(Observed)) 
US_2 <- US_2 %>% group_by(group) %>%  mutate(percent=count/sum(count))
table(US$group)

head(US_2)

ggplot(data=US_2, aes(x=bin, y=percent, fill=group, group=group)) +   
   geom_bar(stat="identity", width=0.8, position=position_dodge(), alpha=0.8) + 
   scale_x_discrete(breaks=c(0, .1,.2,.3,.4,.5,.6,.7,.8, .9)) +
   ylab("Density") + xlab("derived allele frequency") + 
   theme_bw() + theme(legend.position=c(0.85, 0.85), plot.title=element_text(size=20),
                      axis.title.y=element_text(size = 16, vjust=+0.2),
                      axis.title.x=element_text(size = 16, vjust=-0.2),
                      axis.text.y=element_text(size = 14),
                      axis.text.x=element_text(size = 14),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank())
ggsave("dafAcrossEra_bar_GERP4_US.pdf")    


