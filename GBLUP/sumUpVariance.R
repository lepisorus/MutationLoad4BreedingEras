library(cowplot)
library(plyr)
library(ggplot2)
library(dplyr)

library(tidyr)
library(data.table)


theme_li <- function () { 
  theme_bw(base_size=15) %+replace% 
    theme(
      axis.text.y=element_text(size=15),
      panel.grid.minor = element_blank(),
      legend.position=c(0.15, 0.9)
    )
}

#df <- fread("allSNPs.frq.GERP.region.uniq.txt")
#df <- as.data.frame(df)
#names(df) <- c("pos", "frq", "GERP", "annotation")
#head(df)

for(g in c(2:15)){

SNP <- fread(paste0("output_snpeff_ce_", g, ".snpe"))
head(SNP)
names(SNP)

SNP <- as.data.frame(SNP)
names(SNP)[1] <- "pos"
SNP_slim <- SNP[, c(1,6,7,10)]
length(SNP_slim$pos)
SNP_slim <- SNP_slim %>% filter(Effect_A2 != 0) 
length(SNP_slim$pos)
SNP_combine <- merge(df, SNP_slim, by="pos")

SNP_combine["SNPfreq"] <- ifelse(SNP_combine$frq <= 0.1, 0.1, SNP_combine$frq)
SNP_combine["SNPfreq"] <- ifelse(SNP_combine$frq <= 0.2 & SNP_combine$frq > 0.1, 0.2, SNP_combine$SNPfreq)
SNP_combine["SNPfreq"] <- ifelse(SNP_combine$frq <= 0.3 & SNP_combine$frq > 0.2, 0.3, SNP_combine$SNPfreq)
SNP_combine["SNPfreq"] <- ifelse(SNP_combine$frq <= 0.4 & SNP_combine$frq > 0.3, 0.4, SNP_combine$SNPfreq)
SNP_combine["SNPfreq"] <- ifelse(SNP_combine$frq <= 0.5 & SNP_combine$frq > 0.4, 0.5, SNP_combine$SNPfreq)
table(SNP_combine $ SNPfreq)
names(SNP_combine)[7] <- "varExp"
 

SNP_combine2 <- SNP_combine %>% filter(GERP>0) %>% group_by(SNPfreq) %>% summarize(nSNP=length(SNPfreq), varExp2=sum(varExp))
SNP_combine2["group"] <- "GERP>0"
write.table(SNP_combine2, file=paste0("GBLUP_trait", g, ".sum.txt"), quote=F, sep="\t", row.names=F)

#############################################################################
##### generate same number of random SNPs

getRandomSNP <- function(SNP_combine, verbose=FALSE, outfile="./randomSNP/rsnp1.txt"){
  
  ###
  sub <- SNP_combine %>% filter(GERP>0)
  myr <- subset(SNP_combine, SNP_combine$GERP <= 0)
  out <- data.frame()
  a <- c(0.10, 0.20, 0.30, 0.40, 0.50)
    for (f in a){
      num <- nrow(subset(sub,  SNPfreq == f))
      subr <- subset(myr, SNPfreq == f)
      
      idx <- sample(1:nrow(subr), num, replace=FALSE)
      out <- rbind(out, subr[idx, ])
    }
out2 <- out %>% group_by(SNPfreq) %>% summarize(nSNP=length(SNPfreq), varExp2=sum(varExp))
out2["group"] <- "random"
   write.table(out2, outfile, sep="\t", row.names=FALSE, quote=FALSE)
} #end of function

for(i in 1:10){
  getRandomSNP(SNP_combine, verbose=FALSE, 
               outfile= paste0("./randomSNP", g, "/rsnp", i, ".sum.txt"))
}

#############################################################################
##### generate same number of random SNPs from the genic regions

getRandomSNP <- function(SNP_combine, verbose=FALSE, outfile="./randomSNP2_genic/rsnp1.txt"){
  
  ###
  sub <- SNP_combine %>% filter(GERP>0)
  myr <- subset(SNP_combine, SNP_combine$GERP <= 0 & SNP_combine$annotation != "repeat_region")
  out <- data.frame()
  a <- c(0.10, 0.20, 0.30, 0.40, 0.50)
    for (f in a){
      num <- nrow(subset(sub,  SNPfreq == f))
      subr <- subset(myr, SNPfreq == f)
      
      idx <- sample(1:nrow(subr), num, replace=FALSE)
      out <- rbind(out, subr[idx, ])
    }
out2 <- out %>% group_by(SNPfreq) %>% summarize(nSNP=length(SNPfreq), varExp2=sum(varExp))
out2["group"] <- "random_genic"
   write.table(out2, outfile, sep="\t", row.names=FALSE, quote=FALSE)
} #end of function

for(i in 1:10){
  getRandomSNP(SNP_combine, verbose=FALSE, 
               outfile= paste0("./randomSNP", g, "_genic/", "rsnp", i, ".sum.txt"))
}

}#for


#trait2 to trait16

for(g in c(2:15)){

#random1 total random from the genome

path <- paste0("./randomSNP", g, "/")
filenames <- paste0(path, list.files(path, pattern = "sum.txt"))

random = ldply(filenames, function(filename) {
  dum = read.table(filename, header=T)
  return(dum)
})
random1 <- random  

#random2 total random from the genic regions

path <- paste0("./randomSNP", g, "_genic/")
filenames <- paste0(path, list.files(path, pattern = "sum.txt"))

random = ldply(filenames, function(filename) {
  dum = read.table(filename, header=T)
  return(dum)
})
random2 <- random 
random2$group <- "random_genic"

df <- read.delim(paste0("GBLUP_trait", g, ".sum.txt"))
df$group <- "GERP>0"
combine <- as.data.frame(rbind(df, random1, random2))

ggplot(data=combine, aes(x=SNPfreq, y=varExp2, color=group)) + geom_point() +xlab("maf") + ylab("variance explained")  + theme_li()
ggsave(paste0("trait", g, "TotalVarExp_allInOne_sum.pdf"))
}

#trait 2,3,5,6,11,14 see the difference
##trait 12,13,15  difference in random from the whole genome but no difference with random from genic regions


