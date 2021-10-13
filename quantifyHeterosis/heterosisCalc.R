library(tidyr)
library(ggplot2)
library(tidyverse)

theme_li <- function () { 
  theme_bw(base_size=15) %+replace% 
    theme(
     # axis.text.y=element_text(size=15),
      panel.grid.minor = element_blank(),
      legend.position="none"
    )
}


df <- read.delim("hybrid_traits.txt")
head(df)
parent <- read.csv("BLUPs_Phe_4_Set_no_outliers.csv")
head(parent)
names(parent)[1] <- "Maternal"
head(df)
head(parent)
df1 <- merge(df, parent, by="Maternal")
head(df1)
length(df1$ID)
names(parent)[1] <- "Paternal"
head(parent)
df2 <- merge(df1, parent, by="Paternal")
names(df2)

df2["heterosis.ASI"] <- df2$ASI - (df2$ASI_BLUP.x + df2$ASI_BLUP.y)/2
df2["heterosis.DTA"] <- df2$DTA - (df2$DTA_BLUP.x + df2$DTA_BLUP.y)/2
df2["heterosis.DTS"] <- df2$DTS - (df2$DTS_BLUP.x + df2$DTS_BLUP.y)/2
df2["heterosis.PH"] <- df2$PH - (df2$PH_BLUP.x + df2$PH_BLUP.y)/2
df2["heterosis.EH"] <- df2$EH - (df2$EH_BLUP.x + df2$EH_BLUP.y)/2
df2["heterosis.EP"] <- df2$EP - (df2$EP_BLUP.x + df2$EP_BLUP.y)/2
df2["heterosis.LL"] <- df2$LL - (df2$LL_BLUP.x + df2$LL_BLUP.y)/2
df2["heterosis.LW"] <- df2$LW - (df2$LW_BLUP.x + df2$LW_BLUP.y)/2
df2["heterosis.TL"] <- df2$TL - (df2$TL_BLUP.x + df2$TL_BLUP.y)/2
df2["heterosis.TBN"] <- df2$TBN - (df2$TBN_BLUP.x + df2$TBN_BLUP.y)/2
df2["heterosis.LAL"] <- df2$LAL - (df2$LAL_BLUP.x + df2$LAL_BLUP.y)/2
df2["heterosis.LAU"] <- df2$LAU - (df2$LAU_BLUP.x + df2$LAU_BLUP.y)/2


df2["PercHeterosis.ASI"] <- df2$heterosis.ASI / ((df2$ASI_BLUP.x + df2$ASI_BLUP.y)/2)
df2["PercHeterosis.DTA"] <- df2$heterosis.DTA / ((df2$DTA_BLUP.x + df2$DTA_BLUP.y)/2)
df2["PercHeterosis.DTS"] <- df2$heterosis.DTS / ((df2$DTS_BLUP.x + df2$DTS_BLUP.y)/2)
df2["PercHeterosis.PH"] <- df2$heterosis.PH/((df2$PH_BLUP.x + df2$PH_BLUP.y)/2)
df2["PercHeterosis.EH"] <- df2$heterosis.EH / ((df2$EH_BLUP.x + df2$EH_BLUP.y)/2)
df2["PercHeterosis.EP"] <- df2$heterosis.EP /((df2$EP_BLUP.x + df2$EP_BLUP.y)/2)
df2["PercHeterosis.LL"] <- df2$heterosis.LL / ((df2$LL_BLUP.x + df2$LL_BLUP.y)/2)
df2["PercHeterosis.LW"] <- df2$heterosis.LW / ((df2$LW_BLUP.x + df2$LW_BLUP.y)/2)
df2["PercHeterosis.TL"] <- df2$heterosis.TL /((df2$TL_BLUP.x + df2$TL_BLUP.y)/2)
df2["PercHeterosis.TBN"] <- df2$heterosis.TBN / ((df2$TBN_BLUP.x + df2$TBN_BLUP.y)/2)
df2["PercHeterosis.LAL"] <- df2$heterosis.LAL / ((df2$LAL_BLUP.x + df2$LAL_BLUP.y)/2)
df2["PercHeterosis.LAU"] <- df2$heterosis.LAU /((df2$LAU_BLUP.x + df2$LAU_BLUP.y)/2)

df.new <- df2[, 1:20]

ASI <- lm(df2$ASI ~ df2$Time)
summary(ASI)

y= -0.03* X + 63.4
R2 = 0.07， p = 0.002

DTA <- lm(df2$DTA ~ df2$Time)
summary(DTA)

y= 0.06* X - 62.49
R2 = 0.03， p = 0.029

DTS <- lm(df2$DTS ~ df2$Time)
summary(DTS)

y= 0.03 * X - 5.22
R2 = 0.002， p = 0.26

PH <- lm(df2$PH ~ df2$Time)
summary(PH)

y= 0.16 * X - 11.18
R2 = 0.004， p = 0.49

EH <- lm(df2$EH ~ df2$Time)
summary(EH)
y= -0.42 * X - 976.09
R2 = 0.038， p = 0.02

EP <- lm(df2$EP ~ df2$Time)
summary(EP)

y= -0.17 * X + 382.95
R2 = 0.155， p = 7.383e-06

LL <- lm(df2$LL ~ df2$Time)
summary(LL)

y= 0.0005 * X + 89.96
R2 = 5.756e-05， p = 0.994

LW <- lm(df2$LW ~ df2$Time)
summary(LW)

y= 0.004 * X + 1.53
R2 = 0.3705， p = 0.5439

TL <- lm(df2$TL ~ df2$Time)
summary(TL)

y= -0.012 * X + 63.51
R2 = 0.0013， p = 0.7054

TBN <- lm(df2$TBN ~ df2$Time)
summary(TBN)

y= -0.26 * X + 530.09
R2 = 0.2609， p = 2.803e-09

LAL <- lm(df2$LAL ~ df2$Time)
summary(LAL)

y= 0.35944 * X - 556.8
R2 = 0.226， p = 4.087e-08


LAU <- lm(df2$LAU ~ df2$Time)
summary(LAU)

y= -0.35 * X + 728.88
R2 = 0.2249， p = 4.511e-08


df.new.long <- gather(df.new, trait, measurement, ASI:LAU, factor_key=TRUE)

df.new.long["Era"] <- ifelse(df.new.long$Time>=1960 & df.new.long$Time<=1979, "CN_1960&70s", df.new.long$Time)
df.new.long$Era <- ifelse(df.new.long$Time>=1980 & df.new.long$Time<=1999, "CN_1980&90s", df.new.long$Era)
df.new.long$Era <- ifelse(df.new.long$Time>=2000 & df.new.long$Time<=2020, "CN_2000&10s", df.new.long$Era)

#subset the traits which showed the signficant correlation with Time
df.new.long2 <- df.new.long %>% filter (trait == "ASI" | trait== "EP" | trait=="LAU" | trait== "TBN" )

##
df <- read.delim("finalTraitValuePlot.txt")
df.new.long2 <- df %>% mutate(trait = fct_relevel(trait, c("ASI", "EP", "LAU", "TBN")))

 
ggplot(data=df.new.long2, aes(x=Time, y=measurement, color=trait, group=trait)) + geom_point() + facet_grid(trait ~., scales="free")+
scale_x_continuous(breaks=c(1970, 1980, 1990, 2000, 2010)) +  
geom_smooth(method='lm',formula=y~x)+
theme_li() +
#geom_text(color="black", aes(label = lm_eqn(lm(measurement ~ Time, df.new.long2))), parse = TRUE) + 
theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ylab("trait values")
ggsave("hybridTraitChange.pdf")

write.table(df.new.long2, file="finalTraitValuePlot.txt", quote=F, row.names=F, sep="\t")


##df2 percentage of heterosis

df.new <- df2[, c(1:6, 51:74)]
df.new.long <- gather(df.new, trait, measurement,PercHeterosis.ASI:PercHeterosis.LAU, factor_key=TRUE)

df.new.long["Era"] <- ifelse(df.new.long$Time>=1960 & df.new.long$Time<=1979, "CN_1960&70s", df.new.long$Time)
df.new.long$Era <- ifelse(df.new.long$Time>=1980 & df.new.long$Time<=1999, "CN_1980&90s", df.new.long$Era)
df.new.long$Era <- ifelse(df.new.long$Time>=2000 & df.new.long$Time<=2020, "CN_2000&10s", df.new.long$Era)

head(df.new.long)
table(df.new.long$trait)

 df.new.long["trait2"] <- str_split(df.new.long$trait, "\\.", simplify=TRUE)[,2]
head(df.new.long)
table(df.new.long$trait2)

df.new.long2 <- df.new.long %>% filter (trait2 == "ASI" |  trait2== "EP" | trait2=="LAU" | trait2== "TBN" )

ggplot(data=df.new.long2, aes(x=Time, y=measurement, color=trait2, group=trait2)) + geom_point() + facet_grid(trait2 ~., scales="free")+
scale_x_continuous(breaks=c(1970, 1980, 1990, 2000, 2010)) +  
geom_smooth(method='lm',formula=y~x)+
theme_li() +
#geom_text(color="black", aes(label = lm_eqn(lm(measurement ~ Time, df.new.long2))), parse = TRUE) + 
theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ylab("percentage of heterosis")
ggsave("percentageHeterosis4traitsWithSignificantChangeOverTime.pdf")

write.table(df.new.long2, file="finalPercHeterosisChangPlot.txt", quote=F, row.names=F, sep="\t")

ASI <- df.new.long2 %>% filter (trait2 == "ASI")
summary(lm(ASI$measurement ~ ASI$Time))
#R2 = 0.00297; p=0.5883

EH <- df.new.long2 %>% filter (trait2 == "EH")
summary(lm(EH$measurement ~ EH$Time))
#R2 = 0.006391; p=0.3936

EP <- df.new.long2 %>% filter (trait2 == "EP")
summary(lm(EP$measurement ~ EP$Time))
#R2 = 0.007558; p=0.3534

TBN <- df.new.long2 %>% filter (trait2 == "TBN")
summary(lm(TBN$measurement ~ TBN$Time))
#R2 = 0.02958; p=0.06488

LAL <- df.new.long2 %>% filter (trait2 == "LAL")
summary(lm(LAL$measurement ~ LAL$Time))
#R2 = 0.01728; p=0.1596

LAU <- df.new.long2 %>% filter (trait2 == "LAU")
summary(lm(LAU$measurement ~ LAU$Time))
#R2 = 0.05233; p=0.01352

