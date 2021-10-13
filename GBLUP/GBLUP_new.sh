# convert all data into traditional vcf format
perl jvcf2vcf.pl allChr_seperateAllele.jvcf allSNPs.vcf

#calculate the allele frequency of all SNPs in the vcf file. 
#It is used to compare the explained variance of SNPs at categories of maf (allele frequency).
# Bias: SNPs with higher frequency usually explained more variation.
vcftools --vcf allSNPs.vcf --freq --out allSNPs

#remove SNPs with maf < 0.05
sed 's/\:/\t/g' allSNPs.frq | awk -v OFS="\t" '$6>0.05 && $6 < 0.95' > allSNPs.maf0.05.frq

#get the GERP value for all SNP sites
join -1 1 -2 1 <(sed 's/ /\t/g' ../AllChr.SorghumAllele.GERP | cut -f1,3)  <(cut -f3,10- allSNPs.vcf | sort -k1,1)  > allSNPs.withGERP.txt

#convert 0/0 0/1 1/1 with 0,1,2 code
perl recodeVCF4GBLUP.pl allSNPs.withGERP.txt  allSNPs.withGERP.recode.txt

#join the allele frequency file and the genotype file together; and processing file format etc.
join -1 1 -2 1 <(awk '{print $1"_"$2, $5, $6, $7, $8}' allSNPs.maf0.05.frq | sort -k1,1) <(sort -k1,1 allSNPs.withGERP.recode.txt) > allSNPs.withGERP.recode2.txt
sed -i -e 's/ /\t/g' -e 's/\_/\t/g' allSNPs.withGERP.recode2.txt

#determine which is the minor allele and write down the minor allele frequency (maf)
perl mafFrq.pl allSNPs.withGERP.recode2.txt allSNPs.withGERP.frq.txt

#the list of deleterious allele with GERP>0 form the map.txt file
awk -v OFS="\t" '$4 >0 {print $1, $2, $1"_"$2}' allSNPs.withGERP.frq.txt | sort -n -k1,1 -n -k2,2 | sed 's/ /\t/g' > map.txt

## we will calculate the variance explained, dominance effect and additive effect for all SNPs together (delterious + nonDelterious)
##and then compare the difference between the two categories


##get the snp for each chromosome; and move to each sub-directory
paste <(awk -v OFS="\t" '$0 ~ /^1\t/ {print $1"_"$2}' allSNPs.withGERP.frq.txt)  <(awk -v OFS="\t" '$0 ~ /^1\t/ ' allSNPs.withGERP.frq.txt | cut -f5-) | sed 's/ /\t/g' > chr1.SNPs.txt
paste <(awk -v OFS="\t" '$0 ~ /^2\t/ {print $1"_"$2}' allSNPs.withGERP.frq.txt)  <(awk -v OFS="\t" '$0 ~ /^2\t/ ' allSNPs.withGERP.frq.txt | cut -f5-) | sed 's/ /\t/g' > chr2.SNPs.txt
paste <(awk -v OFS="\t" '$0 ~ /^3\t/ {print $1"_"$2}' allSNPs.withGERP.frq.txt)  <(awk -v OFS="\t" '$0 ~ /^3\t/ ' allSNPs.withGERP.frq.txt | cut -f5-) | sed 's/ /\t/g' > chr3.SNPs.txt
paste <(awk -v OFS="\t" '$0 ~ /^4\t/ {print $1"_"$2}' allSNPs.withGERP.frq.txt)  <(awk -v OFS="\t" '$0 ~ /^4\t/ ' allSNPs.withGERP.frq.txt | cut -f5-) | sed 's/ /\t/g' > chr4.SNPs.txt
paste <(awk -v OFS="\t" '$0 ~ /^5\t/ {print $1"_"$2}' allSNPs.withGERP.frq.txt)  <(awk -v OFS="\t" '$0 ~ /^5\t/ ' allSNPs.withGERP.frq.txt | cut -f5-) | sed 's/ /\t/g' > chr5.SNPs.txt
paste <(awk -v OFS="\t" '$0 ~ /^6\t/ {print $1"_"$2}' allSNPs.withGERP.frq.txt)  <(awk -v OFS="\t" '$0 ~ /^6\t/ ' allSNPs.withGERP.frq.txt | cut -f5-) | sed 's/ /\t/g' > chr6.SNPs.txt
paste <(awk -v OFS="\t" '$0 ~ /^7\t/ {print $1"_"$2}' allSNPs.withGERP.frq.txt)  <(awk -v OFS="\t" '$0 ~ /^7\t/ ' allSNPs.withGERP.frq.txt | cut -f5-) | sed 's/ /\t/g' > chr7.SNPs.txt
paste <(awk -v OFS="\t" '$0 ~ /^8\t/ {print $1"_"$2}' allSNPs.withGERP.frq.txt)  <(awk -v OFS="\t" '$0 ~ /^8\t/ ' allSNPs.withGERP.frq.txt | cut -f5-) | sed 's/ /\t/g' > chr8.SNPs.txt
paste <(awk -v OFS="\t" '$0 ~ /^9\t/ {print $1"_"$2}' allSNPs.withGERP.frq.txt)  <(awk -v OFS="\t" '$0 ~ /^9\t/ ' allSNPs.withGERP.frq.txt | cut -f5-) | sed 's/ /\t/g' > chr9.SNPs.txt
paste <(awk -v OFS="\t" '$0 ~ /^10\t/ {print $1"_"$2}' allSNPs.withGERP.frq.txt)  <(awk -v OFS="\t" '$0 ~ /^10\t/ ' allSNPs.withGERP.frq.txt | cut -f5-) | sed 's/ /\t/g' > chr10.SNPs.txt

mv chr1.SNPs.txt ./chr1
mv chr2.SNPs.txt ./chr2
mv chr3.SNPs.txt ./chr3
mv chr4.SNPs.txt ./chr4
mv chr5.SNPs.txt ./chr5
mv chr6.SNPs.txt ./chr6
mv chr7.SNPs.txt ./chr7
mv chr8.SNPs.txt ./chr8
mv chr9.SNPs.txt ./chr9
mv chr10.SNPs.txt ./chr10

##we will do the transposition of each snp file (required format for GVC_BLUP). As the snp file for each chromosome is still very big, we will split each file into multiple files with 300000 lines.
split -l 300000 chr1.SNPs.txt
for f in x*; do awk '{for (i=1; i<=NF; i++) a[i, NR]=$i; max=(max < NF?NF:max)} END {for (i=1;i<=max; i++) {for (j=1; j<= NR; j++) printf "%s%s", a[i, j], (j==NR?RS:FS)}}' $f | sed 's/ /\t/g' > $(basename $f | sed 's/x/trans\./g'); done 
for f in trans.*; do cat chr1.SNPs.t.txt | paste - $f >temp; cp temp chr1.SNPs.t.txt ; done; rm temp
##determine the number of deleterious SNPs in each chromosome, which is an element of parameter settings in the gp3.dat file
awk -v FS=" " '{print NF}' chr1.SNPs.t.txt | head -3

split -l 300000 chr2.SNPs.txt
for f in x*; do awk '{for (i=1; i<=NF; i++) a[i, NR]=$i; max=(max < NF?NF:max)} END {for (i=1;i<=max; i++) {for (j=1; j<= NR; j++) printf "%s%s", a[i, j], (j==NR?RS:FS)}}' $f | sed 's/ /\t/g' > $(basename $f | sed 's/x/trans\./g'); done 
for f in trans.*; do cat chr2.SNPs.t.txt | paste - $f >temp; cp temp chr2.SNPs.t.txt ; done; rm temp
awk -v FS=" " '{print NF}' chr2.SNPs.t.txt | head -3

split -l 300000 chr3.SNPs.txt
for f in x*; do awk '{for (i=1; i<=NF; i++) a[i, NR]=$i; max=(max < NF?NF:max)} END {for (i=1;i<=max; i++) {for (j=1; j<= NR; j++) printf "%s%s", a[i, j], (j==NR?RS:FS)}}' $f | sed 's/ /\t/g' > $(basename $f | sed 's/x/trans\./g'); done 
for f in trans.*; do cat chr3.SNPs.t.txt | paste - $f >temp; cp temp chr3.SNPs.t.txt ; done; rm temp
awk -v FS=" " '{print NF}' chr3.SNPs.t.txt | head -3

split -l 300000 chr4.SNPs.txt
for f in x*; do awk '{for (i=1; i<=NF; i++) a[i, NR]=$i; max=(max < NF?NF:max)} END {for (i=1;i<=max; i++) {for (j=1; j<= NR; j++) printf "%s%s", a[i, j], (j==NR?RS:FS)}}' $f | sed 's/ /\t/g' > $(basename $f | sed 's/x/trans\./g'); done 
for f in trans.*; do cat chr4.SNPs.t.txt | paste - $f >temp; cp temp chr4.SNPs.t.txt ; done; rm temp
awk -v FS=" " '{print NF}' chr4.SNPs.t.txt | head -3

split -l 300000 chr5.SNPs.txt
for f in x*; do awk '{for (i=1; i<=NF; i++) a[i, NR]=$i; max=(max < NF?NF:max)} END {for (i=1;i<=max; i++) {for (j=1; j<= NR; j++) printf "%s%s", a[i, j], (j==NR?RS:FS)}}' $f | sed 's/ /\t/g' > $(basename $f | sed 's/x/trans\./g'); done 
for f in trans.*; do cat chr5.SNPs.t.txt | paste - $f >temp; cp temp chr5.SNPs.t.txt ; done; rm temp
awk -v FS=" " '{print NF}' chr5.SNPs.t.txt | head -3

split -l 300000 chr6.SNPs.txt
for f in x*; do awk '{for (i=1; i<=NF; i++) a[i, NR]=$i; max=(max < NF?NF:max)} END {for (i=1;i<=max; i++) {for (j=1; j<= NR; j++) printf "%s%s", a[i, j], (j==NR?RS:FS)}}' $f | sed 's/ /\t/g' > $(basename $f | sed 's/x/trans\./g'); done 
for f in trans.*; do cat chr6.SNPs.t.txt | paste - $f >temp; cp temp chr6.SNPs.t.txt ; done; rm temp
awk -v FS=" " '{print NF}' chr6.SNPs.t.txt | head -3


split -l 300000 chr7.SNPs.txt
for f in x*; do awk '{for (i=1; i<=NF; i++) a[i, NR]=$i; max=(max < NF?NF:max)} END {for (i=1;i<=max; i++) {for (j=1; j<= NR; j++) printf "%s%s", a[i, j], (j==NR?RS:FS)}}' $f | sed 's/ /\t/g' > $(basename $f | sed 's/x/trans\./g'); done 
for f in trans.*; do cat chr7.SNPs.t.txt | paste - $f >temp; cp temp chr7.SNPs.t.txt ; done; rm temp
awk -v FS=" " '{print NF}' chr7.SNPs.t.txt | head -3


split -l 300000 chr8.SNPs.txt
for f in x*; do awk '{for (i=1; i<=NF; i++) a[i, NR]=$i; max=(max < NF?NF:max)} END {for (i=1;i<=max; i++) {for (j=1; j<= NR; j++) printf "%s%s", a[i, j], (j==NR?RS:FS)}}' $f | sed 's/ /\t/g' > $(basename $f | sed 's/x/trans\./g'); done 
for f in trans.*; do cat chr8.SNPs.t.txt | paste - $f >temp; cp temp chr8.SNPs.t.txt ; done; rm temp
awk -v FS=" " '{print NF}' chr8.SNPs.t.txt | head -3


split -l 300000 chr9.SNPs.txt
for f in x*; do awk '{for (i=1; i<=NF; i++) a[i, NR]=$i; max=(max < NF?NF:max)} END {for (i=1;i<=max; i++) {for (j=1; j<= NR; j++) printf "%s%s", a[i, j], (j==NR?RS:FS)}}' $f | sed 's/ /\t/g' > $(basename $f | sed 's/x/trans\./g'); done 
for f in trans.*; do cat chr9.SNPs.t.txt | paste - $f >temp; cp temp chr9.SNPs.t.txt ; done; rm temp
awk -v FS=" " '{print NF}' chr9.SNPs.t.txt | head -3


split -l 300000 chr10.SNPs.txt
for f in x*; do awk '{for (i=1; i<=NF; i++) a[i, NR]=$i; max=(max < NF?NF:max)} END {for (i=1;i<=max; i++) {for (j=1; j<= NR; j++) printf "%s%s", a[i, j], (j==NR?RS:FS)}}' $f | sed 's/ /\t/g' > $(basename $f | sed 's/x/trans\./g'); done 
for f in trans.*; do cat chr10.SNPs.t.txt | paste - $f >temp; cp temp chr10.SNPs.t.txt ; done; rm temp
awk -v FS=" " '{print NF}' chr10.SNPs.t.txt | head -3

##utilize script.R to convert the above files into the format GVC_BLUP requires
R CMD BATCH script.R
##from the above steps, we get the genotype files as the input in the GVC_BLUP, those were indicated in gp3.dat.


##then we ran the gvc_blup for each trait consecutively; greml_ce was utilized; 
##note: here we only gave one parameter setting file "gp3_2.dat", it takes into account of the trait on the second column of the input trait file 
## if you have multiple traits, you have to prepare one parameter setting file for each trait. 
OMP_NUM_THREADS=16 
for i in {2..16}; do ./greml_ce gp3_$i.dat; done

##the above command gave us the result file "output_snpeff_ce_2.snpe", we utilize the R script to generate what the deleterious mutation and random SNPs tells us.

Rscript processGBLUPresult2.R -i output_snpeff_ce_2.snpe -o1 GBLUP_trait2.txt -o2 ./randomSNP2/
Rscript processGBLUPresult2.R -i output_snpeff_ce_3.snpe -o1 GBLUP_trait3.txt -o2 ./randomSNP3/
Rscript processGBLUPresult2.R -i output_snpeff_ce_4.snpe -o1 GBLUP_trait4.txt -o2 ./randomSNP4/
Rscript processGBLUPresult2.R -i output_snpeff_ce_5.snpe -o1 GBLUP_trait5.txt -o2 ./randomSNP5/
Rscript processGBLUPresult2.R -i output_snpeff_ce_6.snpe -o1 GBLUP_trait6.txt -o2 ./randomSNP6/
Rscript processGBLUPresult2.R -i output_snpeff_ce_7.snpe -o1 GBLUP_trait7.txt -o2 ./randomSNP7/
Rscript processGBLUPresult2.R -i output_snpeff_ce_8.snpe -o1 GBLUP_trait8.txt -o2 ./randomSNP8/
Rscript processGBLUPresult2.R -i output_snpeff_ce_9.snpe -o1 GBLUP_trait9.txt -o2 ./randomSNP9/
Rscript processGBLUPresult2.R -i output_snpeff_ce_10.snpe -o1 GBLUP_trait10.txt -o2 ./randomSNP10/
Rscript processGBLUPresult2.R -i output_snpeff_ce_11.snpe -o1 GBLUP_trait11.txt -o2 ./randomSNP11/
Rscript processGBLUPresult2.R -i output_snpeff_ce_12.snpe -o1 GBLUP_trait12.txt -o2 ./randomSNP12/
Rscript processGBLUPresult2.R -i output_snpeff_ce_13.snpe -o1 GBLUP_trait13.txt -o2 ./randomSNP13/
Rscript processGBLUPresult2.R -i output_snpeff_ce_14.snpe -o1 GBLUP_trait14.txt -o2 ./randomSNP14/
Rscript processGBLUPresult2.R -i output_snpeff_ce_15.snpe -o1 GBLUP_trait15.txt -o2 ./randomSNP15/
Rscript processGBLUPresult2.R -i output_snpeff_ce_16.snpe -o1 GBLUP_trait16.txt -o2 ./randomSNP16/

##generate the random snps from genic regions
##get the annotation of the snps
bedtools intersect -a <(sed -e 's/\_/\t/g' allSNPs.frq.GERP.txt | awk -v OFS="\t" '{print $1, $2, $2, $3, $4}')  -b /work/LAS/mhufford-lab/lwang/AGPv3/Zea_mays.AGPv3.21.bed -wo | cut -f1,2,4,5,9 | awk -v OFS="\t" '{print $1"_"$2, $3, $4, $5}' >  allSNPs.frq.GERP.region.txt
#get the uniq 
sort -k1,1 allSNPs.frq.GERP.region.txt | awk -F"\t" '!seen[$1, $2, $3]++' > allSNPs.frq.GERP.region.uniq.txt

Rscript processGBLUPresult2_genic.R -i output_snpeff_ce_2.snpe -o1 GBLUP_trait2.txt -o2 ./randomSNP2_genic/
Rscript processGBLUPresult2_genic.R -i output_snpeff_ce_3.snpe -o1 GBLUP_trait3.txt -o2 ./randomSNP3_genic/
Rscript processGBLUPresult2_genic.R -i output_snpeff_ce_4.snpe -o1 GBLUP_trait4.txt -o2 ./randomSNP4_genic/
Rscript processGBLUPresult2_genic.R -i output_snpeff_ce_5.snpe -o1 GBLUP_trait5.txt -o2 ./randomSNP5_genic/
Rscript processGBLUPresult2_genic.R -i output_snpeff_ce_6.snpe -o1 GBLUP_trait6.txt -o2 ./randomSNP6_genic/
Rscript processGBLUPresult2_genic.R -i output_snpeff_ce_7.snpe -o1 GBLUP_trait7.txt -o2 ./randomSNP7_genic/
Rscript processGBLUPresult2_genic.R -i output_snpeff_ce_8.snpe -o1 GBLUP_trait8.txt -o2 ./randomSNP8_genic/
Rscript processGBLUPresult2_genic.R -i output_snpeff_ce_9.snpe -o1 GBLUP_trait9.txt -o2 ./randomSNP9_genic/
Rscript processGBLUPresult2_genic.R -i output_snpeff_ce_10.snpe -o1 GBLUP_trait10.txt -o2 ./randomSNP10_genic/
Rscript processGBLUPresult2_genic.R -i output_snpeff_ce_11.snpe -o1 GBLUP_trait11.txt -o2 ./randomSNP11_genic/
Rscript processGBLUPresult2_genic.R -i output_snpeff_ce_12.snpe -o1 GBLUP_trait12.txt -o2 ./randomSNP12_genic/
Rscript processGBLUPresult2_genic.R -i output_snpeff_ce_13.snpe -o1 GBLUP_trait13.txt -o2 ./randomSNP13_genic/
Rscript processGBLUPresult2_genic.R -i output_snpeff_ce_14.snpe -o1 GBLUP_trait14.txt -o2 ./randomSNP14_genic/
Rscript processGBLUPresult2_genic.R -i output_snpeff_ce_15.snpe -o1 GBLUP_trait15.txt -o2 ./randomSNP15_genic/
Rscript processGBLUPresult2_genic.R -i output_snpeff_ce_16.snpe -o1 GBLUP_trait16.txt -o2 ./randomSNP16_genic/

#the above files are utilized to generate the mean variance explained  by each snp

#in order to get the sum of variance explained:
#please run the following:
R CMD BATCH sumUpVariance2.R
R CMD BATCH sumUpVariance2.R

##plus: plotK.R is utilized to plot the dominance effect of the snps






