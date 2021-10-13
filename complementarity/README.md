#It lists the steps for figure 3 in the manuscript

*  produce hybrid genotypes and calculate the percentage of heterozygous sites


> `perl produceHybridGenotype.pl allChr.jvcf 120hybridDetails.txt 120hybridGenotype.txt`

> `cat <(cut -f1-126 header.txt) <(zcat 120hybridGenotype_allChr_SorghumAllele_GERP.txt.gz | sed 's/\//\t/g' | sed 's/\_/\t/g') > allChr.120Hybrid.withAncestralAllele.GERP.txt`
> 
> `join -1 1 -2 1 120hybridGenotype_allChr_SorghumAllele_GERP.txt <(awk -v OFS="\t"
 '{$2=""; print $0}' allChr_sorted.jvcf) | sed 's/ /\t/g' > allChr.parentAndHybr
id.txt`

> `awk -v OFS="\t" '$2 != "N" && $3 > 0'  allChr.parentAndHybrid.txt | sed 's/\//\t
> /g' > allChr.parentAndHybrid.withAncestralAllele.GERPO.txt`

> `cat header.txt allChr.parentAndHybrid.withAncestralAllele.GERPO.txt > allChr.parentAndHybrid.withAncestralAllele.GERPO.2.txt`
> 
> sed -i 's/AA/A/g' allChr.parentAndHybrid.withAncestralAllele.GERPO.2.txt 
> 
> sed -i 's/CC/C/g' allChr.parentAndHybrid.withAncestralAllele.GERPO.2.txt 
> 
> sed -i 's/GG/G/g' allChr.parentAndHybrid.withAncestralAllele.GERPO.2.txt 
> 
> sed -i 's/TT/T/g' allChr.parentAndHybrid.withAncestralAllele.GERPO.2.txt 
> 
> sed -i 's/X/x/g' allChr.parentAndHybrid.withAncestralAllele.GERPO.2.txt
> 
> sed -i 's/AC/X/g' allChr.parentAndHybrid.withAncestralAllele.GERPO.2.txt 
> 
> sed -i 's/AG/X/g' allChr.parentAndHybrid.withAncestralAllele.GERPO.2.txt 
> 
> sed -i 's/AT/X/g' allChr.parentAndHybrid.withAncestralAllele.GERPO.2.txt 
> 
> sed -i 's/CA/X/g' allChr.parentAndHybrid.withAncestralAllele.GERPO.2.txt 
> 
> sed -i 's/CG/X/g' allChr.parentAndHybrid.withAncestralAllele.GERPO.2.txt 
> 
> sed -i 's/CT/X/g' allChr.parentAndHybrid.withAncestralAllele.GERPO.2.txt 
> 
> sed -i 's/GA/X/g' allChr.parentAndHybrid.withAncestralAllele.GERPO.2.txt 
> 
> sed -i 's/GC/X/g' allChr.parentAndHybrid.withAncestralAllele.GERPO.2.txt 
> 
> sed -i 's/GT/X/g' allChr.parentAndHybrid.withAncestralAllele.GERPO.2.txt 
> 
> sed -i 's/TA/X/g' allChr.parentAndHybrid.withAncestralAllele.GERPO.2.txt 
> 
> sed -i 's/TC/X/g' allChr.parentAndHybrid.withAncestralAllele.GERPO.2.txt 
> 
> sed -i 's/TG/X/g' allChr.parentAndHybrid.withAncestralAllele.GERPO.2.txt 
> 
> `perl recodeHybridParentGenotype.pl allChr.parentAndHybrid.withAncestralAllele.GERPO.2.txt allChr.parentAndHybrid.withAncestralAllele.GERPO.recode.txt`
> 
> module load r
> 
> `R CMD BATCH codeHetero.R`

* plot the results

> `R CMD BATCH plotHeteroRatio.R`


