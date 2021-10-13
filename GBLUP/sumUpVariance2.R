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
      legend.position=c(0.85, 0.9)
    )
}


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


