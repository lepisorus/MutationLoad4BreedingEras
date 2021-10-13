join -1 1 -2 1 120hybridGenotype_allChr_SorghumAllele_GERP.txt <(awk -v OFS="\t" '{$2=""; print $0}' allChr_sorted.jvcf) | sed 's/ /\t/g' > allChr.parentAndHybrid.txt
awk -v OFS="\t" '$2 != "N" && $3 > 0'  allChr.parentAndHybrid.txt | sed 's/\//\t/g' > allChr.parentAndHybrid.withAncestralAllele.GERPO.txt

cat header.txt allChr.parentAndHybrid.withAncestralAllele.GERPO.txt > allChr.parentAndHybrid.withAncestralAllele.GERPO.2.txt

sed -i 's/AA/A/g' allChr.parentAndHybrid.withAncestralAllele.GERPO.2.txt 
sed -i 's/CC/C/g' allChr.parentAndHybrid.withAncestralAllele.GERPO.2.txt 
sed -i 's/GG/G/g' allChr.parentAndHybrid.withAncestralAllele.GERPO.2.txt 
sed -i 's/TT/T/g' allChr.parentAndHybrid.withAncestralAllele.GERPO.2.txt 

#sed -i 's/X/x/g' allChr.parentAndHybrid.withAncestralAllele.GERPO.2.txt
sed -i 's/AC/X/g' allChr.parentAndHybrid.withAncestralAllele.GERPO.2.txt 
sed -i 's/AG/X/g' allChr.parentAndHybrid.withAncestralAllele.GERPO.2.txt 
sed -i 's/AT/X/g' allChr.parentAndHybrid.withAncestralAllele.GERPO.2.txt 
sed -i 's/CA/X/g' allChr.parentAndHybrid.withAncestralAllele.GERPO.2.txt 
sed -i 's/CG/X/g' allChr.parentAndHybrid.withAncestralAllele.GERPO.2.txt 
sed -i 's/CT/X/g' allChr.parentAndHybrid.withAncestralAllele.GERPO.2.txt 
sed -i 's/GA/X/g' allChr.parentAndHybrid.withAncestralAllele.GERPO.2.txt 
sed -i 's/GC/X/g' allChr.parentAndHybrid.withAncestralAllele.GERPO.2.txt 
sed -i 's/GT/X/g' allChr.parentAndHybrid.withAncestralAllele.GERPO.2.txt 
sed -i 's/TA/X/g' allChr.parentAndHybrid.withAncestralAllele.GERPO.2.txt 
sed -i 's/TC/X/g' allChr.parentAndHybrid.withAncestralAllele.GERPO.2.txt 
sed -i 's/TG/X/g' allChr.parentAndHybrid.withAncestralAllele.GERPO.2.txt 

perl recodeHybridParentGenotype.pl allChr.parentAndHybrid.withAncestralAllele.GERPO.2.txt allChr.parentAndHybrid.withAncestralAllele.GERPO.recode.txt
join -1 1 -2 1 <(sort -k1,1 hm32SNP.correctedGERP.txt | cut -f1,2,7 | sed 's/ /\t/g') <(sort -k1,1 allChr.parentAndHybrid.withAncestralAllele.GERPO.recode.txt) | sed 's/ /\t/g' | sed '1h;1d;$!H;$!d;G' > allChr.parentAndHybrid.withAncestralAllele.GERPO.recode3.txt

module load r
R CMD BATCH parentHybridDeleteriousAllele.R



