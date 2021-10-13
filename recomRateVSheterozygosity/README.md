#recombination rate versus heterozygosity

## if we run the R script to turn the physical position to genetic position for the whole genome altogether, it will take forever. so we split it into 10 chromosomes, and for each chromosome, we split into two parts, and run the R script in parallel for each file.
> Rscript Phys2Genetic.R -i 120hybridGenotype_chr10_pos.txt -o 120hybridGenotype_chr10_pos_cm.txt
> Rscript Phys2Genetic.R -i 120hybridGenotype_chr1_posaa -o 120hybridGenotype_chr1_posaa_cm.txt
> Rscript Phys2Genetic.R -i 120hybridGenotype_chr1_posab -o 120hybridGenotype_chr1_posab_cm.txt
> 
> Rscript Phys2Genetic.R -i 120hybridGenotype_chr2_posaa -o 120hybridGenotype_chr2_posaa_cm.txt
> Rscript Phys2Genetic.R -i 120hybridGenotype_chr2_posab -o 120hybridGenotype_chr2_posab_cm.txt
> 
> Rscript Phys2Genetic.R -i 120hybridGenotype_chr3_posaa -o 120hybridGenotype_chr3_posaa_cm.txt
> Rscript Phys2Genetic.R -i 120hybridGenotype_chr3_posab -o 120hybridGenotype_chr3_posab_cm.txt
> 
> Rscript Phys2Genetic.R -i 120hybridGenotype_chr4_posaa -o 120hybridGenotype_chr4_posaa_cm.txt
> Rscript Phys2Genetic.R -i 120hybridGenotype_chr4_posab -o 120hybridGenotype_chr4_posab_cm.txt
> 
> Rscript Phys2Genetic.R -i 120hybridGenotype_chr5_posaa -o 120hybridGenotype_chr5_posaa_cm.txt
> Rscript Phys2Genetic.R -i 120hybridGenotype_chr5_posab -o 120hybridGenotype_chr5_posab_cm.txt
> 
> Rscript Phys2Genetic.R -i 120hybridGenotype_chr6_posaa -o 120hybridGenotype_chr6_posaa_cm.txt
> Rscript Phys2Genetic.R -i 120hybridGenotype_chr6_posab -o 120hybridGenotype_chr6_posab_cm.txt
> 
> Rscript Phys2Genetic.R -i 120hybridGenotype_chr7_posaa -o 120hybridGenotype_chr7_posaa_cm.txt
> Rscript Phys2Genetic.R -i 120hybridGenotype_chr7_posab -o 120hybridGenotype_chr7_posab_cm.txt
> 
> Rscript Phys2Genetic.R -i 120hybridGenotype_chr8_posaa -o 120hybridGenotype_chr8_posaa_cm.txt
> Rscript Phys2Genetic.R -i 120hybridGenotype_chr8_posab -o 120hybridGenotype_chr8_posab_cm.txt
> 
> Rscript Phys2Genetic.R -i 120hybridGenotype_chr9_posaa -o 120hybridGenotype_chr9_posaa_cm.txt
> Rscript Phys2Genetic.R -i 120hybridGenotype_chr9_posab -o 120hybridGenotype_chr9_posab_cm.txt
> 
> for f in 120hybridGenotype_chr1*_pos*_cm.txt 
> do cat $f >> allSNPs.cm.txt
> done

##prepare the file
> join -1 1 -2 1 <(awk -v OFS="\t" '{print $1"_"$2, $3}' allSNPs.cm.txt | sort -k1,1) <(paste <(awk -v OFS="\t" '{print $1"_"$2}' allChr.120Hybrid.withAncestralAllele.GERP.txt) <(cut -f3- allChr.120Hybrid.withAncestralAllele.GERP.txt) | sort -k1,1) | sed 's/ /\t/g' | sed 's/\_/\t/g'  | sort -n -k1,1 -n -k2,2 > allChr.120Hybrid.withAncestralAllele.GERP.withCM.txt
> cat header.txt allChr.120Hybrid.withAncestralAllele.GERP.withCM.txt > allChr.120Hybrid.withAncestralAllele.GERP.withCM.sorted.txt
> 
> perl recodeHybridGenotype.pl  allChr.120Hybrid.withAncestralAllele.GERP.withCM.sorted.txt allChr.120Hybrid.withAncestralAllele.GERP.withCM.sorted.recode.txt

##do the sliding window for the recombination rate and the mean heterozygosity 
> R CMD BATCH recomRateSlidingWindow.R

##plot the results
> R CMD BATCH recRateVSheteroPercentPlot.R


