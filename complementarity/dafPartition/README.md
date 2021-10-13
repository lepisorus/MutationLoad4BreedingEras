#non-deleterious
> `join -1 1 -2 1 <(cut -f1,5 allSNPs.daf.txt | sed 's/ /\t/g' | sed '1d' | sort -k1,1) <(sed '1d' 120hybridGenotype_allChr_SorghumAllele_nonDeleterious.txt | sort -k1,1 ) | sed 's/ /\t/g' > 120hybridGenotype_allChr_SorghumAllele_nonDeleterious_daf.txt`
> 
> `cat header.hybrid.txt 120hybridGenotype_allChr_SorghumAllele_nonDeleterious_daf.txt > 120hybridGenotype_allChr_SorghumAllele_nonDeleterious_daf2.txt`

#GERP02
> `join -1 1 -2 1 <(cut -f1,5 allSNPs.daf.txt | sed 's/ /\t/g' | sed '1d' | sort -k1,1) <(sed '1d' 120hybridGenotype_allChr_SorghumAllele_GERP02.txt | sort -k1,1 ) | sed 's/ /\t/g' > 120hybridGenotype_allChr_SorghumAllele_GERP02_daf.txt`
> 
> `cat header.hybrid.txt 120hybridGenotype_allChr_SorghumAllele_GERP02_daf.txt > 120hybridGenotype_allChr_SorghumAllele_GERP02_daf2.txt`

#GERP24
> `join -1 1 -2 1 <(cut -f1,5 allSNPs.daf.txt | sed 's/ /\t/g' | sed '1d' | sort -k1,1) <(sed '1d' 120hybridGenotype_allChr_SorghumAllele_GERP24.txt | sort -k1,1 ) | sed 's/ /\t/g' > 120hybridGenotype_allChr_SorghumAllele_GERP24_daf.txt`
> 
> `cat header.hybrid.txt 120hybridGenotype_allChr_SorghumAllele_GERP24_daf.txt > 120hybridGenotype_allChr_SorghumAllele_GERP24_daf2.txt`

#GERP4
> `join -1 1 -2 1 <(cut -f1,5 allSNPs.daf.txt | sed 's/ /\t/g' | sed '1d' | sort -k1,1) <(sed '1d' 120hybridGenotype_allChr_SorghumAllele_GERP4.txt | sort -k1,1 ) | sed 's/ /\t/g' > 120hybridGenotype_allChr_SorghumAllele_GERP4_daf.txt`
> 
> `cat header.hybrid.txt 120hybridGenotype_allChr_SorghumAllele_GERP4_daf.txt > 120hybridGenotype_allChr_SorghumAllele_GERP4_daf2.txt`

#plot the results
> `R CMD BATCH codeHetero_dafPartition2.R`

