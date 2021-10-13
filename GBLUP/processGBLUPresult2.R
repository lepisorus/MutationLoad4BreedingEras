library(dplyr)
library(data.table)
library(tidyr)
library(ggplot2)
#library(optparse)
 
#option_list <- list(make_option(c('-i','--in_file'), action='store', type='character', default=NULL, help='Input file containing per site gerp scores'),
#                    make_option(c('-o1','--out_file'), action='store', type='character', default=NULL, help='Output file'),
#                    make_option(c('-o2','--out_dir'), action='store', type='character', default=NULL, help='Output directory')
#                    )
#opt <- parse_args(OptionParser(option_list = option_list))
 


df <- fread("allSNPs.frq.GERP.txt")
df <- as.data.frame(df)
names(df) <- c("pos", "frq", "GERP")
head(df)

#SNP
SNP <- fread("output_snpeff_ce_2.snpe")
head(SNP)
names(SNP)

SNP <- as.data.frame(SNP)
names(SNP)[1] <- "pos"
SNP_slim <- SNP[, c(1,6,7,10)]
length(SNP_slim$pos)
SNP_slim <- SNP_slim %>% filter(Effect_A2 != 0) %>% mutate(k=Effect_D2/Effect_A2) 

length(SNP_slim$pos)
SNP_combine1 <- merge(df, SNP_slim, by="pos")

SNP_combine <- SNP_combine1 %>% filter(k > -2 & k < 2)
SNP_combine["SNPfreq"] <- ifelse(SNP_combine$frq <= 0.1, 0.1, SNP_combine$frq)
SNP_combine["SNPfreq"] <- ifelse(SNP_combine$frq <= 0.2 & SNP_combine$frq > 0.1, 0.2, SNP_combine$SNPfreq)
SNP_combine["SNPfreq"] <- ifelse(SNP_combine$frq <= 0.3 & SNP_combine$frq > 0.2, 0.3, SNP_combine$SNPfreq)
SNP_combine["SNPfreq"] <- ifelse(SNP_combine$frq <= 0.4 & SNP_combine$frq > 0.3, 0.4, SNP_combine$SNPfreq)
SNP_combine["SNPfreq"] <- ifelse(SNP_combine$frq <= 0.5 & SNP_combine$frq > 0.4, 0.5, SNP_combine$SNPfreq)
table(SNP_combine $ SNPfreq)


names(SNP_combine)[6] <- "varExp"

#deleterious SNPs
SNP_combine2 <- SNP_combine %>% filter(GERP>0) %>% group_by(SNPfreq) %>% summarize(nSNP=length(SNPfreq), varExp2=mean(varExp))
SNP_combine2["group"] <- "GERP0"

write.table(SNP_combine1, file="GBLUP_trait2_full.txt", quote=F, sep="\t", row.names=F)
write.table(SNP_combine2, file="GBLUP_trait2.txt", quote=F, sep="\t", row.names=F)

#SNP_combine24 <- SNP_combine %>% filter(GERP>=2 &  GERP<4) %>% group_by(SNPfreq) %>% summarize(nSNP=length(SNPfreq), varExp2=mean(varExp))
#SNP_combine24["group"] <- "GERP24"

#SNP_combine4 <- SNP_combine %>% filter(GERP>=4 ) %>% group_by(SNPfreq) %>% summarize(varExp2=mean(varExp))
#SNP_combine4["group"] <- "GERP4"



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
out2 <- out %>% group_by(SNPfreq) %>% summarize(nSNP=length(SNPfreq), varExp2=mean(varExp))
out2["group"] <- "random"
   write.table(out2, outfile, sep="\t", row.names=FALSE, quote=FALSE)
}

for(i in 1:10){
  getRandomSNP(SNP_combine, verbose=FALSE, 
               outfile= paste0("./randomSNP2/", "rsnp", i, ".txt"))
}

