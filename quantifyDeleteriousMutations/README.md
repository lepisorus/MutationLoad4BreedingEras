#processing the data
> `zcat 25M_chr1.jvcf.gz | head -3 | sed '1d' | sed '1d' > header.txt`

> `cat header.txt <(zcat 25M_chr1.jvcf.gz | sed '1d' | sed '1d' | sed '1d') <(zcat 25M_chr2.jvcf.gz | sed '1d' | sed '1d' | sed '1d') <(zcat 25M_chr3.jvcf.gz | sed '1d' | sed '1d' | sed '1d') <(zcat 25M_chr4.jvcf.gz | sed '1d' | sed '1d' | sed '1d') <(zcat 25M_chr5.jvcf.gz | sed '1d' | sed '1d' | sed '1d') <(zcat 25M_chr6.jvcf.gz | sed '1d' | sed '1d' | sed '1d') <(zcat 25M_chr7.jvcf.gz | sed '1d' | sed '1d' | sed '1d') <(zcat 25M_chr8.jvcf.gz | sed '1d' | sed '1d' | sed '1d') <(zcat 25M_chr9.jvcf.gz | sed '1d' | sed '1d' | sed '1d') <(zcat 25M_chr10.jvcf.gz | sed '1d' | sed '1d' | sed '1d') > allChr.jvcf`


> `paste <(awk -v OFS="\t" '{print $1"_"$2}' allChr.jvcf) <(awk -v OFS="\t" '{$1=$2=""; print $0}' allChr.jvcf) | sort -k1,1 > allChr_sorted.jvcf`

> `sed -i 's/chr//g' allChr_sorted.jvcf`

> `cd ../breedingEra`

> `ln -s ../recalGERP/allChr.GERP.noB73.txt`

> `join -1 1 -2 1 <(sort -k1,1 allChr.GERP.noB73.txt) allChr_sorted.jvcf | sed 's/ /\t/g' > allChrSNP.GERP.txt`

> `cat <(tail -1 allChrSNP.GERP.txt) <(sed '$ d' allChrSNP.GERP.txt)  > allChrSNP.GERP2.txt`
> 
> `perl recodeGenotype.pl allChrSNP.GERP2.txt allChrSNP.GERP.recode.txt`

> `head -1 allChrSNP.GERP2.txt > longHeader.txt`

> `cat longHeader.txt allChrSNP.GERP.recode.txt > allChrSNP.GERP.recode2.txt`
> 

#plot the results
> `R CMD BATCH replot_load_breedingEra.R`





