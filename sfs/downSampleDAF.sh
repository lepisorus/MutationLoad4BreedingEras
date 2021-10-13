#
uniq -c CN2_samples.txt | sed 's/^[ \t]*//g' | cut -d ' ' -f2 > CN2.txt
uniq -c CN3_samples.txt | sed 's/^[ \t]*//g' | cut -d ' ' -f2 > CN3.txt
uniq -c US1_samples.txt | sed 's/^[ \t]*//g' | cut -d ' ' -f2 > US1.txt
uniq -c US2_samples.txt | sed 's/^[ \t]*//g' | cut -d ' ' -f2 > US2.txt

#as we are calculating the derived allele frequency, we got to know the ancestral allele (AA) state. the AA flag was indicated in the ANNO tab of vcf file. 
perl addAA_info_in_vcf.pl allSNPs.SNPeff.GERP.vcf allSNPs.SNPeff.GERP.addAA.vcf


cat header.txt <(awk -v OFS="\t" '$1>2 && $1 < 4 {$1=""; print $0}' allSNPs.SNPeff.GERP.addAA.vcf) > allSNPs.SNPeff.GERP24.addAA.vcf
cat header.txt <(awk -v OFS="\t" '$1>0 && $1 <= 2{$1=""; print $0}' allSNPs.SNPeff.GERP.addAA.vcf) > allSNPs.SNPeff.GERP02.addAA.vcf
cat header.txt <(awk -v OFS="\t" '$1<=0 {$1=""; print $0}' allSNPs.SNPeff.GERP.addAA.vcf) > allSNPs.SNPeff.GERPless0.addAA.vcf

sed -i -e 's/^[ \t]*//g' allSNPs.SNPeff.GERPless0.addAA.vcf
sed -i -e 's/^[ \t]*//g' allSNPs.SNPeff.GERP02.addAA.vcf
sed -i -e 's/^[ \t]*//g' allSNPs.SNPeff.GERP24.addAA.vcf

#get the snps, whose ancestral state is known, not unknown as indicated by "N". And prepare the file as the input of sofos
sed 's/\=/\t/g' allSNPs.SNPeff.GERP4.addAA.vcf | sed 's/\;/\t/g' | sed '1d' | awk -v OFS="\t" '$9 == $4 || $9==$5 {print $3}' > allSNPs.SNPeff.GERP4.knownAncestral.txt
join -1 1 -2 3 <(sort -k1,1 allSNPs.SNPeff.GERP4.knownAncestral.txt)  <(sort -k3,3 allSNPs.SNPeff.GERP4.addAA.vcf ) | sed 's/ /\t/g' | cut -f2,3,1,4- | sed 's/ /\t/g' | sort -n -k1,1 -n -k2,2 > allSNPs.SNPeff.GERP4.addAA2.vcf
cat header.txt <(paste <(cut -f2,3 allSNPs.SNPeff.GERP4.addAA2.vcf) <(cut -f1 allSNPs.SNPeff.GERP4.addAA2.vcf) <(cut -f4- allSNPs.SNPeff.GERP4.addAA2.vcf) | sed 's/ /\t/g')  > allSNPs.SNPeff.GERP4.addAA3.vcf

sed 's/\=/\t/g' allSNPs.SNPeff.GERPless0.addAA.vcf | sed 's/\;/\t/g' | sed '1d' | awk -v OFS="\t" '$9 == $4 || $9==$5 {print $3}' > allSNPs.SNPeff.GERPless0.knownAncestral.txt
join -1 1 -2 3 <(sort -k1,1 allSNPs.SNPeff.GERPless0.knownAncestral.txt)  <(sort -k3,3 allSNPs.SNPeff.GERPless0.addAA.vcf ) | sed 's/ /\t/g' | cut -f2,3,1,4- | sed 's/ /\t/g' | sort -n -k1,1 -n -k2,2 > allSNPs.SNPeff.GERPless0.addAA2.vcf
cat header.txt <(paste <(cut -f2,3 allSNPs.SNPeff.GERPless0.addAA2.vcf) <(cut -f1 allSNPs.SNPeff.GERPless0.addAA2.vcf) <(cut -f4- allSNPs.SNPeff.GERPless0.addAA2.vcf) | sed 's/ /\t/g')  > allSNPs.SNPeff.GERPless0.addAA3.vcf

sed 's/\=/\t/g' allSNPs.SNPeff.GERP02.addAA.vcf | sed 's/\;/\t/g' | sed '1d' | awk -v OFS="\t" '$9 == $4 || $9==$5 {print $3}' > allSNPs.SNPeff.GERP02.knownAncestral.txt
join -1 1 -2 3 <(sort -k1,1 allSNPs.SNPeff.GERP02.knownAncestral.txt)  <(sort -k3,3 allSNPs.SNPeff.GERP02.addAA.vcf ) | sed 's/ /\t/g' | cut -f2,3,1,4- | sed 's/ /\t/g' | sort -n -k1,1 -n -k2,2 > allSNPs.SNPeff.GERP02.addAA2.vcf
cat header.txt <(paste <(cut -f2,3 allSNPs.SNPeff.GERP02.addAA2.vcf) <(cut -f1 allSNPs.SNPeff.GERP02.addAA2.vcf) <(cut -f4- allSNPs.SNPeff.GERP02.addAA2.vcf) | sed 's/ /\t/g')  > allSNPs.SNPeff.GERP02.addAA3.vcf


sed 's/\=/\t/g' allSNPs.SNPeff.GERP24.addAA.vcf | sed 's/\;/\t/g' | sed '1d' | awk -v OFS="\t" '$9 == $4 || $9==$5 {print $3}' > allSNPs.SNPeff.GERP24.knownAncestral.txt
join -1 1 -2 3 <(sort -k1,1 allSNPs.SNPeff.GERP24.knownAncestral.txt)  <(sort -k3,3 allSNPs.SNPeff.GERP24.addAA.vcf ) | sed 's/ /\t/g' | cut -f2,3,1,4- | sed 's/ /\t/g' | sort -n -k1,1 -n -k2,2 > allSNPs.SNPeff.GERP24.addAA2.vcf
cat header.txt <(paste <(cut -f2,3 allSNPs.SNPeff.GERP24.addAA2.vcf) <(cut -f1 allSNPs.SNPeff.GERP24.addAA2.vcf) <(cut -f4- allSNPs.SNPeff.GERP24.addAA2.vcf) | sed 's/ /\t/g')  > allSNPs.SNPeff.GERP24.addAA3.vcf



#GERPless0
vcftools --vcf allSNPs.SNPeff.GERPless0.addAA3.vcf --keep CN1.txt --recode --out CN1.GERPless0
vcftools --vcf allSNPs.SNPeff.GERPless0.addAA3.vcf --keep CN2.txt --recode --out CN2.GERPless0
vcftools --vcf allSNPs.SNPeff.GERPless0.addAA3.vcf --keep CN3.txt --recode --out CN3.GERPless0
vcftools --vcf allSNPs.SNPeff.GERPless0.addAA3.vcf --keep US1.txt --recode --out US1.GERPless0
vcftools --vcf allSNPs.SNPeff.GERPless0.addAA3.vcf --keep US2.txt --recode --out US2.GERPless0

##get the annotation column
cut -f8 allSNPs.SNPeff.GERPless0.addAA3.vcf > anno.GERPless0.txt

##paste the three parts of the vcf
paste <(cut -f1-7 CN1.GERPless0.recode.vcf) anno.GERPless0.txt <(cut -f9- CN1.GERPless0.recode.vcf) | sed 's/ /\t/g' > CN1.GERPless0.recode2.vcf
paste <(cut -f1-7 CN2.GERPless0.recode.vcf) anno.GERPless0.txt <(cut -f9- CN2.GERPless0.recode.vcf) | sed 's/ /\t/g' > CN2.GERPless0.recode2.vcf
paste <(cut -f1-7 CN3.GERPless0.recode.vcf) anno.GERPless0.txt <(cut -f9- CN3.GERPless0.recode.vcf) | sed 's/ /\t/g' > CN3.GERPless0.recode2.vcf

paste <(cut -f1-7 US1.GERPless0.recode.vcf) anno.GERPless0.txt <(cut -f9- US1.GERPless0.recode.vcf) | sed 's/ /\t/g' > US1.GERPless0.recode2.vcf
paste <(cut -f1-7 US2.GERPless0.recode.vcf) anno.GERPless0.txt <(cut -f9- US2.GERPless0.recode.vcf) | sed 's/ /\t/g' > US2.GERPless0.recode2.vcf

cat VCFheader.txt <(sed -e 's/^[ \t]*//g' CN1.GERPless0.recode2.vcf) > CN1.GERPless0.recode3.vcf
cat VCFheader.txt <(sed -e 's/^[ \t]*//g' CN2.GERPless0.recode2.vcf) > CN2.GERPless0.recode3.vcf
cat VCFheader.txt <(sed -e 's/^[ \t]*//g' CN3.GERPless0.recode2.vcf) > CN3.GERPless0.recode3.vcf

cat VCFheader.txt <(sed -e 's/^[ \t]*//g' US1.GERPless0.recode2.vcf) > US1.GERPless0.recode3.vcf
cat VCFheader.txt <(sed -e 's/^[ \t]*//g' US2.GERPless0.recode2.vcf) > US2.GERPless0.recode3.vcf


/work/LAS/mhufford-lab/lwang/bin/SoFoS/sofos -n 60 -u -t -a 1.0 -b 1.0 CN1.GERPless0.recode3.vcf >  CN1.GERPless0.daf
/work/LAS/mhufford-lab/lwang/bin/SoFoS/sofos -n 60 -u -t -a 1.0 -b 1.0 CN2.GERPless0.recode3.vcf >  CN2.GERPless0.daf
/work/LAS/mhufford-lab/lwang/bin/SoFoS/sofos -n 60 -u -t -a 1.0 -b 1.0 CN3.GERPless0.recode3.vcf >  CN3.GERPless0.daf

/work/LAS/mhufford-lab/lwang/bin/SoFoS/sofos -n 148 -u -t -a 1.0 -b 1.0 US1.GERPless0.recode3.vcf >  US1.GERPless0.daf
/work/LAS/mhufford-lab/lwang/bin/SoFoS/sofos -n 148 -u -t -a 1.0 -b 1.0 US2.GERPless0.recode3.vcf >  US2.GERPless0.daf

#GERP02
vcftools --vcf allSNPs.SNPeff.GERP02.addAA3.vcf --keep CN1.txt --recode --out CN1.GERP02
vcftools --vcf allSNPs.SNPeff.GERP02.addAA3.vcf --keep CN2.txt --recode --out CN2.GERP02
vcftools --vcf allSNPs.SNPeff.GERP02.addAA3.vcf --keep CN3.txt --recode --out CN3.GERP02
vcftools --vcf allSNPs.SNPeff.GERP02.addAA3.vcf --keep US1.txt --recode --out US1.GERP02
vcftools --vcf allSNPs.SNPeff.GERP02.addAA3.vcf --keep US2.txt --recode --out US2.GERP02

cut -f8 allSNPs.SNPeff.GERP02.addAA3.vcf > anno.GERP02.txt

paste <(cut -f1-7 CN1.GERP02.recode.vcf) anno.GERP02.txt <(cut -f9- CN1.GERP02.recode.vcf) | sed 's/ /\t/g' > CN1.GERP02.recode2.vcf
paste <(cut -f1-7 CN2.GERP02.recode.vcf) anno.GERP02.txt <(cut -f9- CN2.GERP02.recode.vcf) | sed 's/ /\t/g' > CN2.GERP02.recode2.vcf
paste <(cut -f1-7 CN3.GERP02.recode.vcf) anno.GERP02.txt <(cut -f9- CN3.GERP02.recode.vcf) | sed 's/ /\t/g' > CN3.GERP02.recode2.vcf

paste <(cut -f1-7 US1.GERP02.recode.vcf) anno.GERP02.txt <(cut -f9- US1.GERP02.recode.vcf) | sed 's/ /\t/g' > US1.GERP02.recode2.vcf
paste <(cut -f1-7 US2.GERP02.recode.vcf) anno.GERP02.txt <(cut -f9- US2.GERP02.recode.vcf) | sed 's/ /\t/g' > US2.GERP02.recode2.vcf



cat VCFheader.txt <(sed -e 's/^[ \t]*//g' CN1.GERP02.recode2.vcf) > CN1.GERP02.recode3.vcf
cat VCFheader.txt <(sed -e 's/^[ \t]*//g' CN2.GERP02.recode2.vcf) > CN2.GERP02.recode3.vcf
cat VCFheader.txt <(sed -e 's/^[ \t]*//g' CN3.GERP02.recode2.vcf) > CN3.GERP02.recode3.vcf

cat VCFheader.txt <(sed -e 's/^[ \t]*//g' US1.GERP02.recode2.vcf) > US1.GERP02.recode3.vcf
cat VCFheader.txt <(sed -e 's/^[ \t]*//g' US2.GERP02.recode2.vcf) > US2.GERP02.recode3.vcf


/work/LAS/mhufford-lab/lwang/bin/SoFoS/sofos -n 60 -u -t -a 1.0 -b 1.0 CN1.GERP02.recode3.vcf >  CN1.GERP02.daf
/work/LAS/mhufford-lab/lwang/bin/SoFoS/sofos -n 60 -u -t -a 1.0 -b 1.0 CN2.GERP02.recode3.vcf >  CN2.GERP02.daf
/work/LAS/mhufford-lab/lwang/bin/SoFoS/sofos -n 60 -u -t -a 1.0 -b 1.0 CN3.GERP02.recode3.vcf >  CN3.GERP02.daf

/work/LAS/mhufford-lab/lwang/bin/SoFoS/sofos -n 148 -u -t -a 1.0 -b 1.0 US1.GERP02.recode3.vcf >  US1.GERP02.daf
/work/LAS/mhufford-lab/lwang/bin/SoFoS/sofos -n 148 -u -t -a 1.0 -b 1.0 US2.GERP02.recode3.vcf >  US2.GERP02.daf

#GERP24
vcftools --vcf allSNPs.SNPeff.GERP24.addAA3.vcf --keep CN1.txt --recode --out CN1.GERP24
vcftools --vcf allSNPs.SNPeff.GERP24.addAA3.vcf --keep CN2.txt --recode --out CN2.GERP24
vcftools --vcf allSNPs.SNPeff.GERP24.addAA3.vcf --keep CN3.txt --recode --out CN3.GERP24
vcftools --vcf allSNPs.SNPeff.GERP24.addAA3.vcf --keep US1.txt --recode --out US1.GERP24
vcftools --vcf allSNPs.SNPeff.GERP24.addAA3.vcf --keep US2.txt --recode --out US2.GERP24

cut -f8 allSNPs.SNPeff.GERP24.addAA3.vcf > anno.GERP24.txt

paste <(cut -f1-7 CN1.GERP24.recode.vcf) anno.GERP24.txt <(cut -f9- CN1.GERP24.recode.vcf) | sed 's/ /\t/g' > CN1.GERP24.recode2.vcf
paste <(cut -f1-7 CN2.GERP24.recode.vcf) anno.GERP24.txt <(cut -f9- CN2.GERP24.recode.vcf) | sed 's/ /\t/g' > CN2.GERP24.recode2.vcf
paste <(cut -f1-7 CN3.GERP24.recode.vcf) anno.GERP24.txt <(cut -f9- CN3.GERP24.recode.vcf) | sed 's/ /\t/g' > CN3.GERP24.recode2.vcf

paste <(cut -f1-7 US1.GERP24.recode.vcf) anno.GERP24.txt <(cut -f9- US1.GERP24.recode.vcf) | sed 's/ /\t/g' > US1.GERP24.recode2.vcf
paste <(cut -f1-7 US2.GERP24.recode.vcf) anno.GERP24.txt <(cut -f9- US2.GERP24.recode.vcf) | sed 's/ /\t/g' > US2.GERP24.recode2.vcf

cat VCFheader.txt <(sed -e 's/^[ \t]*//g' CN1.GERP24.recode2.vcf) > CN1.GERP24.recode3.vcf
cat VCFheader.txt <(sed -e 's/^[ \t]*//g' CN2.GERP24.recode2.vcf) > CN2.GERP24.recode3.vcf
cat VCFheader.txt <(sed -e 's/^[ \t]*//g' CN3.GERP24.recode2.vcf) > CN3.GERP24.recode3.vcf

cat VCFheader.txt <(sed -e 's/^[ \t]*//g' US1.GERP24.recode2.vcf) > US1.GERP24.recode3.vcf
cat VCFheader.txt <(sed -e 's/^[ \t]*//g' US2.GERP24.recode2.vcf) > US2.GERP24.recode3.vcf

/work/LAS/mhufford-lab/lwang/bin/SoFoS/sofos -n 60 -u -t -a 1.0 -b 1.0 CN1.GERP24.recode3.vcf >  CN1.GERP24.daf
/work/LAS/mhufford-lab/lwang/bin/SoFoS/sofos -n 60 -u -t -a 1.0 -b 1.0 CN2.GERP24.recode3.vcf >  CN2.GERP24.daf
/work/LAS/mhufford-lab/lwang/bin/SoFoS/sofos -n 60 -u -t -a 1.0 -b 1.0 CN3.GERP24.recode3.vcf >  CN3.GERP24.daf

/work/LAS/mhufford-lab/lwang/bin/SoFoS/sofos -n 148 -u -t -a 1.0 -b 1.0 US1.GERP24.recode3.vcf >  US1.GERP24.daf
/work/LAS/mhufford-lab/lwang/bin/SoFoS/sofos -n 148 -u -t -a 1.0 -b 1.0 US2.GERP24.recode3.vcf >  US2.GERP24.daf


#GERP4

vcftools --vcf allSNPs.SNPeff.GERP4.addAA3.vcf --keep CN1.txt --recode --out CN1.GERP4
vcftools --vcf allSNPs.SNPeff.GERP4.addAA3.vcf --keep CN2.txt --recode --out CN2.GERP4
vcftools --vcf allSNPs.SNPeff.GERP4.addAA3.vcf --keep CN3.txt --recode --out CN3.GERP4
vcftools --vcf allSNPs.SNPeff.GERP4.addAA3.vcf --keep US1.txt --recode --out US1.GERP4
vcftools --vcf allSNPs.SNPeff.GERP4.addAA3.vcf --keep US2.txt --recode --out US2.GERP4

cut -f8 allSNPs.SNPeff.GERP4.addAA3.vcf > anno.GERP4.txt

paste <(cut -f1-7 CN1.GERP4.recode.vcf) anno.GERP4.txt <(cut -f9- CN1.GERP4.recode.vcf) | sed 's/ /\t/g' > CN1.GERP4.recode2.vcf
paste <(cut -f1-7 CN2.GERP4.recode.vcf) anno.GERP4.txt <(cut -f9- CN2.GERP4.recode.vcf) | sed 's/ /\t/g' > CN2.GERP4.recode2.vcf
paste <(cut -f1-7 CN3.GERP4.recode.vcf) anno.GERP4.txt <(cut -f9- CN3.GERP4.recode.vcf) | sed 's/ /\t/g' > CN3.GERP4.recode2.vcf

paste <(cut -f1-7 US1.GERP4.recode.vcf) anno.GERP4.txt <(cut -f9- US1.GERP4.recode.vcf) | sed 's/ /\t/g' > US1.GERP4.recode2.vcf
paste <(cut -f1-7 US2.GERP4.recode.vcf) anno.GERP4.txt <(cut -f9- US2.GERP4.recode.vcf) | sed 's/ /\t/g' > US2.GERP4.recode2.vcf


cat VCFheader.txt CN1.GERP4.recode2.vcf > CN1.GERP4.recode3.vcf
cat VCFheader.txt CN2.GERP4.recode2.vcf > CN2.GERP4.recode3.vcf
cat VCFheader.txt CN3.GERP4.recode2.vcf > CN3.GERP4.recode3.vcf

cat VCFheader.txt US1.GERP4.recode2.vcf > US1.GERP4.recode3.vcf
cat VCFheader.txt US2.GERP4.recode2.vcf > US2.GERP4.recode3.vcf


/work/LAS/mhufford-lab/lwang/bin/SoFoS/sofos -n 60 -u -t -a 1.0 -b 1.0 CN1.GERP4.recode3.vcf >  CN1.GERP4.daf
/work/LAS/mhufford-lab/lwang/bin/SoFoS/sofos -n 60 -u -t -a 1.0 -b 1.0 CN2.GERP4.recode3.vcf >  CN2.GERP4.daf
/work/LAS/mhufford-lab/lwang/bin/SoFoS/sofos -n 60 -u -t -a 1.0 -b 1.0 CN3.GERP4.recode3.vcf >  CN3.GERP4.daf

/work/LAS/mhufford-lab/lwang/bin/SoFoS/sofos -n 148 -u -t -a 1.0 -b 1.0 US1.GERP4.recode3.vcf >  US1.GERP4.daf
/work/LAS/mhufford-lab/lwang/bin/SoFoS/sofos -n 148 -u -t -a 1.0 -b 1.0 US2.GERP4.recode3.vcf >  US2.GERP4.daf

#we utilize the "DAFacrossEraPlot.R" to plot the result out.
#"CN1.GERP02.daf" and "dafAcrossEra_bar_GERP02_CN.pdf" gave the examples of result file.

