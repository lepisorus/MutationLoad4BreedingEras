library(data.table)
library(tidyverse)

ed <- fread("cache/gerp_snps_effects.csv", data.table=FALSE)
ed$type <- "GERP0"
ed[ed$GERP > 4, ]$type <- "GERP4"
ed[ed$GERP > 2 & ed$GERP <= 4, ]$type <- "GERP24"
ed[ed$GERP > 0 & ed$GERP <= 2, ]$type <- "GERP02"

library(plyr)
stat <- ddply(ed, .(trait, type), summarise,
              meff = mean(A1Effect),
              mfrq = mean(A1Frq))

ed$ancfrq <- -9
#if A1 == sorghum
ed[ed$A1 == ed$sorghum, ]$ancfrq <- ed[ed$A1 == ed$sorghum, ]$A1Frq
ed[ed$A2 == ed$sorghum, ]$ancfrq <- 1- ed[ed$A2 == ed$sorghum, ]$A1Frq

ed$dtype <- "GERP>0"
ed[ed$GERP <=0, ]$dtype <- "GERP<=0"       
              
lau <- ed %>% filter(trait=="LAU")

p0 <- ggplot(lau, aes(x=dtype, y=log10(A1Effect), fill=dtype)) +
    geom_boxplot() +
    #geom_jitter() +
    facet_wrap(.~trait, ncol=2, scales = "free_y") +
    xlab("") +
    ylab("Effect Size (log10)") +
    ggtitle("") +
    theme_classic() +
    #labs(color = "") +
    #scale_y_continuous(limits = c(-0.01, 1)) +
    #geom_vline(xintercept=0, linetype="dashed", color = "red") +
    theme(plot.title = element_text(size=20, face = "bold"), 
          axis.text=element_text(size=12, face="bold"),
          strip.text.x = element_text(size = 16, face = "bold"),
          axis.title=element_text(size=18, face="bold"),
          #axis.text.x = element_text(angle = 20, hjust=0.6),
          legend.position = "none", 
          #axis.text.x = element_text(angle = 20, hjust=0.6),
          legend.title = element_text(size=18, face="bold"),
          legend.text = element_text(size=18))

p0
ggsave("SNPeffectSizeLAU.pdf")

lau.nondel <- lau %>% filter(dtype == "GERP<=0")
lau.del <- lau %>% filter(dtype == "GERP>0")
wilcox.test(lau.nondel$A1Effect, lau.del$A1Effect, alternative="less")  #p=0.037

##dta
DTA <- ed %>% filter(trait=="DTA" & log10(A1Effect) < -1.8)
p0 <- ggplot(DTA, aes(x=dtype, y=log10(A1Effect), fill=dtype)) +
    geom_boxplot() +
    #geom_jitter() +
    facet_wrap(.~trait, ncol=2, scales = "free_y") +
    xlab("") +
    ylab("Effect Size (log10)") +
    ggtitle("") +
    theme_classic() +
    #labs(color = "") +
    scale_y_continuous(limits = c(-2.595, -1.8)) +
    #geom_vline(xintercept=0, linetype="dashed", color = "red") +
    theme(plot.title = element_text(size=20, face = "bold"), 
          axis.text=element_text(size=12, face="bold"),
          strip.text.x = element_text(size = 16, face = "bold"),
          axis.title=element_text(size=18, face="bold"),
          #axis.text.x = element_text(angle = 20, hjust=0.6),
          legend.position = "none", 
          #axis.text.x = element_text(angle = 20, hjust=0.6),
          legend.title = element_text(size=18, face="bold"),
          legend.text = element_text(size=18))
p0
ggsave("SNPeffectSizeDTA.pdf")

dta.nondel <- DTA %>% filter(dtype == "GERP<=0")
dta.del <- DTA %>% filter(dtype == "GERP>0")
wilcox.test(dta.nondel$A1Effect, dta.del$A1Effect, alternative="less")  #p=0.037

##dts
DTS <- ed %>% filter(trait=="DTS" & log10(A1Effect) < -1.8)
p0 <- ggplot(DTS, aes(x=dtype, y=log10(A1Effect), fill=dtype)) +
    geom_boxplot() +
    #geom_jitter() +
    facet_wrap(.~trait, ncol=2, scales = "free_y") +
    xlab("") +
    ylab("Effect Size (log10)") +
    ggtitle("") +
    theme_classic() +
    #labs(color = "") +
    #scale_y_continuous(limits = c(-2.595, -1.8)) +
    #geom_vline(xintercept=0, linetype="dashed", color = "red") +
    theme(plot.title = element_text(size=20, face = "bold"), 
          axis.text=element_text(size=12, face="bold"),
          strip.text.x = element_text(size = 16, face = "bold"),
          axis.title=element_text(size=18, face="bold"),
          #axis.text.x = element_text(angle = 20, hjust=0.6),
          legend.position = "none", 
          #axis.text.x = element_text(angle = 20, hjust=0.6),
          legend.title = element_text(size=18, face="bold"),
          legend.text = element_text(size=18))
p0
ggsave("SNPeffectSizeDTS.pdf")

DTS.nondel <- DTS %>% filter(dtype == "GERP<=0")
DTS.del <- DTS %>% filter(dtype == "GERP>0")
wilcox.test(DTS.nondel$A1Effect, DTS.del$A1Effect, alternative="less")  #p=0.03

##asi
ASI <- ed %>% filter(trait=="ASI" & log10(A1Effect) < -1.8)
p0 <- ggplot(ASI, aes(x=dtype, y=log10(A1Effect), fill=dtype)) +
    geom_boxplot() +
    #geom_jitter() +
    facet_wrap(.~trait, ncol=2, scales = "free_y") +
    xlab("") +
    ylab("Effect Size (log10)") +
    ggtitle("") +
    theme_classic() +
    #labs(color = "") +
    scale_y_continuous(limits = c(-2.595, -1.8)) +
    #geom_vline(xintercept=0, linetype="dashed", color = "red") +
    theme(plot.title = element_text(size=20, face = "bold"), 
          axis.text=element_text(size=12, face="bold"),
          strip.text.x = element_text(size = 16, face = "bold"),
          axis.title=element_text(size=18, face="bold"),
          #axis.text.x = element_text(angle = 20, hjust=0.6),
          legend.position = "none", 
          #axis.text.x = element_text(angle = 20, hjust=0.6),
          legend.title = element_text(size=18, face="bold"),
          legend.text = element_text(size=18))
p0
ggsave("SNPeffectSizeASI.pdf")

ASI.nondel <- ASI %>% filter(dtype == "GERP<=0")
ASI.del <- ASI %>% filter(dtype == "GERP>0")
wilcox.test(ASI.nondel$A1Effect, ASI.del$A1Effect, alternative="less")  #p=0.037





p0 <- ggplot(ed, aes(x=type, y=log10(A1Effect), fill=type)) +
    geom_boxplot() +
    #geom_jitter() +
    facet_wrap(.~trait, ncol=4, scales = "free_y") +
    xlab("") +
    ylab("Effect Size (log10)") +
    ggtitle("") +
    theme_classic() +
    #labs(color = "") +
    #scale_y_continuous(limits = c(-0.01, 1)) +
    #geom_vline(xintercept=0, linetype="dashed", color = "red") +
    theme(plot.title = element_text(size=20, face = "bold"), 
          axis.text=element_text(size=16, face="bold"),
          strip.text.x = element_text(size = 16, face = "bold"),
          axis.title=element_text(size=18, face="bold"),
          axis.text.x=element_blank(),
          legend.position = c(0.85, 0.15), 
          #axis.text.x = element_text(angle = 20, hjust=0.6),
          legend.title = element_text(size=18, face="bold"),
          legend.text = element_text(size=18))

p0
ggsave("SNPeffect4AllTraits.pdf")