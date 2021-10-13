library(tidyverse)
library(data.table)

theme_li <- function () { 
  theme_bw(base_size=15) %+replace% 
    theme(
      axis.text.y=element_text(size=15),
      panel.grid.minor = element_blank(),
      legend.position=c(0.15, 0.9)
    )
}


df <- fread("allChr.120Hybrid.slidingWindow1mb.cM.percHetero2.txt")
df <- data.frame(df)


df["bins"] <- cut(df$recRate, quantile(df$recRate, na.rm=T), include.lowest=TRUE, labels=c("0-25%", "25-50%", "50-75%", "75-100%"))
head(df)
names(df)

df["meanHetero"] <- apply(df[, 7:126], 1, mean)

df.new <- df[complete.cases(df$bins), ]

dfa <- df.new %>% filter(bins=="0-25%")
dfb <- df.new %>% filter(bins=="25-50%")
dfc <- df.new %>% filter(bins=="50-75%")
dfd <- df.new %>% filter(bins=="75-100%")

wilcox.test(dfa$meanHetero, dfb$meanHetero, alternative="less") #p=0.0002898
wilcox.test(dfb$meanHetero, dfc$meanHetero, alternative="less") #p=0.1786
wilcox.test(dfc$meanHetero, dfd$meanHetero, alternative="less") #p=0.9529

wilcox.test(dfa$meanHetero, dfd$meanHetero, alternative="less") #p=0.0003

df.new2 <- df.new %>% filter(bins=="0-25%" | bins=="75-100%")

ggplot(data=df.new2, aes(x=bins, y=meanHetero, fill=bins)) + geom_boxplot()+theme_li()+ labs(y="average of percentage of heterozygous sites", x="recombination rate bins")
ggsave("recRateVSheteroPercentPlot2.pdf")

ggplot(data=df.new, aes(x=recRate, y=meanHetero)) + geom_point() + geom_smooth(method="lm")+theme_li()+ylab("average of percentage of heterozygous sites")+xlab("recombination rate(cM/bp)")
ggsave("recRateVSheteroPercentSmoothPlot.pdf")

test <- lm(df.new$meanHetero ~ df.new$recRate)
summary(test)
#p=0.04